"""Synchronous adapter for the dynamically discovered PDBe MCP tools.

The upstream API server generates tool names from its OpenAPI specification, so names are
intentionally never encoded here. Each discovery/call starts a short-lived stdio server;
this is slower than a persistent session but gives the synchronous CLI deterministic
cleanup and avoids leaking subprocesses in Passenger workers.
"""

from __future__ import annotations

import asyncio
import json
from collections.abc import Callable
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from typing import Any

from .config import Settings

API_SERVER = "pdbe_api_server"
SEARCH_SERVER = "pdbe_search_server"
_SERVER_TYPES = (API_SERVER, SEARCH_SERVER)


class PdbeMcpError(RuntimeError):
    """PDBe MCP discovery or invocation failed."""


@dataclass(frozen=True, slots=True)
class DiscoveredTool:
    name: str
    server_type: str
    description: str
    input_schema: dict[str, Any]

    def openai_schema(self) -> dict[str, Any]:
        return {
            "type": "function",
            "function": {
                "name": self.name,
                "description": self.description,
                "parameters": self.input_schema,
            },
        }


Runner = Callable[[str, str, str | None, dict[str, Any] | None], Any]


class PdbeMcpClient:
    """Discover and execute the official PDBe API and Search MCP tools over stdio."""

    def __init__(self, settings: Settings, *, runner: Runner | None = None) -> None:
        self.settings = settings
        self._runner = runner or self._run_stdio
        self._tools: dict[str, DiscoveredTool] | None = None
        self._search_schema_tool: str | None = None
        self._search_query_tool: str | None = None
        self._search_schema_seen = False

    def discover(self, *, refresh: bool = False) -> list[DiscoveredTool]:
        if self._tools is not None and not refresh:
            return list(self._tools.values())

        discovered: dict[str, DiscoveredTool] = {}
        for server_type in _SERVER_TYPES:
            payload = self._runner(server_type, "list", None, None)
            for item in payload:
                name = str(item.get("name") or "").strip()
                schema = item.get("inputSchema") or item.get("input_schema") or {}
                if not name or not isinstance(schema, dict):
                    raise PdbeMcpError(f"{server_type} advertised an invalid tool schema")
                if name in discovered:
                    raise PdbeMcpError(
                        f"PDBe MCP tool name {name!r} is duplicated across servers"
                    )
                discovered[name] = DiscoveredTool(
                    name=name,
                    server_type=server_type,
                    description=str(item.get("description") or "PDBe retrieval tool"),
                    input_schema=schema,
                )
        if not discovered:
            raise PdbeMcpError("PDBe MCP servers advertised no tools")
        search_tools = [
            item for item in discovered.values() if item.server_type == SEARCH_SERVER
        ]
        schema_tools = [
            item for item in search_tools
            if not (item.input_schema.get("properties") or {})
            and not (item.input_schema.get("required") or [])
        ]
        query_tools = [
            item for item in search_tools
            if "query" in (item.input_schema.get("properties") or {})
            and "query" in (item.input_schema.get("required") or [])
        ]
        if len(schema_tools) != 1 or len(query_tools) != 1:
            raise PdbeMcpError(
                "PDBe Search MCP tools could not be classified unambiguously from their schemas"
            )
        self._search_schema_tool = schema_tools[0].name
        self._search_query_tool = query_tools[0].name
        self._tools = discovered
        return list(discovered.values())

    def openai_tools(self) -> list[dict[str, Any]]:
        return [item.openai_schema() for item in self.discover()]

    def supports(self, name: str) -> bool:
        return name in {item.name for item in self.discover()}

    def call(self, name: str, arguments: dict[str, Any]) -> dict[str, Any]:
        tools = {item.name: item for item in self.discover()}
        tool = tools.get(name)
        if tool is None:
            raise PdbeMcpError(f"PDBe MCP did not advertise tool {name!r}")
        if name == self._search_query_tool and not self._search_schema_seen:
            field_options = {"fl", "filters", "fq", "sort", "facet", "facet_fields",
                             "facet_queries", "group", "group_field", "params"}
            if field_options.intersection(arguments):
                raise PdbeMcpError(
                    "Call the discovered PDBe Search schema tool before using field, "
                    "filter, sort, facet, group, or advanced query options"
                )
        bounded = _bounded_arguments(name, arguments, tool)
        result = self._runner(tool.server_type, "call", name, bounded)
        if not isinstance(result, dict):
            result = {"content": result}
        encoded_size = len(json.dumps(result, ensure_ascii=False).encode("utf-8"))
        if encoded_size > self.settings.pdbe_mcp_max_result_bytes:
            raise PdbeMcpError(
                f"PDBe MCP result was {encoded_size} bytes, above the configured "
                f"{self.settings.pdbe_mcp_max_result_bytes}-byte limit; request fewer fields"
            )
        if name == self._search_schema_tool and result.get("isError") is not True:
            self._search_schema_seen = True
        result["_provenance"] = {
            "source": "PDBe MCP",
            "server": tool.server_type,
            "tool": name,
            "arguments": bounded,
            "retrieved_at": datetime.now(UTC).isoformat(timespec="seconds"),
        }
        return result

    def close(self) -> None:
        """Short-lived MCP sessions leave no persistent resources to close."""

    def _run_stdio(
        self,
        server_type: str,
        operation: str,
        tool_name: str | None,
        arguments: dict[str, Any] | None,
    ) -> Any:
        try:
            return asyncio.run(
                _stdio_request(
                    command=self.settings.pdbe_mcp_command,
                    server_type=server_type,
                    operation=operation,
                    tool_name=tool_name,
                    arguments=arguments,
                    timeout=self.settings.pdbe_mcp_timeout,
                )
            )
        except PdbeMcpError:
            raise
        except (OSError, RuntimeError, TimeoutError) as exc:
            raise PdbeMcpError(
                f"PDBe MCP {server_type} {operation} failed: {type(exc).__name__}: {exc}"
            ) from exc


async def _stdio_request(
    *,
    command: str,
    server_type: str,
    operation: str,
    tool_name: str | None,
    arguments: dict[str, Any] | None,
    timeout: float,
) -> Any:
    try:
        from mcp import ClientSession, StdioServerParameters
        from mcp.client.stdio import stdio_client
    except ModuleNotFoundError as exc:  # pragma: no cover - packaging failure
        raise PdbeMcpError("The 'mcp' package is required for PDBe retrieval") from exc

    params = StdioServerParameters(
        command=command,
        args=[
            "--transport",
            "stdio",
            "--server-type",
            server_type,
        ],
        # The MCP SDK merges this with its safe environment allowlist. In particular,
        # RESTECH_API_KEY is not inherited by the third-party subprocess.
        env={"UV_NO_PROGRESS": "1"},
    )
    try:
        async with (
            asyncio.timeout(timeout),
            stdio_client(params) as (read_stream, write_stream),
            ClientSession(
                read_stream,
                write_stream,
                read_timeout_seconds=timedelta(seconds=timeout),
            ) as session,
        ):
            await session.initialize()
            if operation == "list":
                tools: list[dict[str, Any]] = []
                cursor: str | None = None
                while True:
                    response = await session.list_tools(cursor=cursor)
                    tools.extend(tool.model_dump(mode="json") for tool in response.tools)
                    cursor = response.nextCursor
                    if not cursor:
                        return tools
            if operation != "call" or not tool_name:
                raise PdbeMcpError(f"Unsupported MCP operation {operation!r}")
            response = await session.call_tool(tool_name, arguments or {})
            return response.model_dump(mode="json", exclude_none=True)
    except TimeoutError as exc:
        raise PdbeMcpError(
            f"PDBe MCP {server_type} exceeded the {timeout:g}-second timeout"
        ) from exc


def _bounded_arguments(
    name: str, arguments: dict[str, Any], tool: DiscoveredTool
) -> dict[str, Any]:
    """Bound search fan-out without mutating arbitrary API-tool arguments."""
    out = dict(arguments)
    properties = (tool.input_schema.get("properties") or {}) if tool.input_schema else {}
    if "rows" in properties:
        try:
            out["rows"] = max(1, min(int(out.get("rows", 10)), 25))
        except (TypeError, ValueError) as exc:
            raise PdbeMcpError(f"PDBe tool {name!r} received a non-integer rows value") from exc
    if "start" in out:
        try:
            out["start"] = max(0, int(out["start"]))
        except (TypeError, ValueError) as exc:
            raise PdbeMcpError(f"PDBe tool {name!r} received a non-integer start value") from exc
    return out


__all__ = [
    "API_SERVER",
    "SEARCH_SERVER",
    "DiscoveredTool",
    "PdbeMcpClient",
    "PdbeMcpError",
]
