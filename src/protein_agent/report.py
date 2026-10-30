"""Offline, deterministic single-sequence report for the local web app.

No UniProt identity, structure, expression, or experimental claims are inferred from a
pasted sequence. This entry point never invokes the inference gateway or requires a key.
"""

# HTML template is intentionally readable rather than wrapped to Python's line limit.
# ruff: noqa: E501
from __future__ import annotations

import argparse
import html
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

from .disorder import predict_disorder
from .pdbe_report import fetch_pdbe_structures
from .plots import plot_hydropathy
from .properties import STANDARD_AA, SequenceError, chromophore_candidates, compute_properties
from .uniprot_report import fetch_annotations, find_sequence_annotations

MAX_INPUT_BYTES = 100_000
MAX_RESIDUES = 10_000


def parse_sequence(text: str) -> tuple[str, str]:
    """Accept one FASTA record or bare sequence; reject ambiguous/lossy inputs."""
    if len(text.encode("utf-8")) > MAX_INPUT_BYTES:
        raise SequenceError("input exceeds the 100 KB limit")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise SequenceError("no sequence provided")
    if any(line.startswith(">") for line in lines[1:]):
        raise SequenceError("provide exactly one FASTA record")
    header = lines.pop(0)[1:].strip() if lines[0].startswith(">") else "User sequence"
    if header != "User sequence" and not header:
        raise SequenceError("FASTA header is empty")
    if not lines:
        raise SequenceError("no residues found")
    sequence = "".join("".join(line.split()) for line in lines).upper()
    if len(sequence) > MAX_RESIDUES:
        raise SequenceError("sequence exceeds the 10,000-residue limit")
    invalid = sorted(set(sequence) - STANDARD_AA)
    if invalid:
        raise SequenceError(
            "unsupported characters: " + ", ".join(repr(char) for char in invalid)
            + "; use the 20 standard amino-acid letters (no gaps or ambiguity codes)"
        )
    return header[:200], sequence


def _row(name: str, value: object) -> str:
    return f"<tr><th scope='row'>{html.escape(name)}</th><td>{html.escape(str(value))}</td></tr>"


def _named_value(node: object) -> str | None:
    if isinstance(node, dict):
        value = node.get("value")
        if isinstance(value, str) and value:
            return value
    return None


def _protein_name(description: object) -> str | None:
    if not isinstance(description, dict):
        return None
    recommended = description.get("recommendedName") or {}
    name = _named_value(recommended.get("fullName"))
    if name:
        return name
    for key in ("submissionNames", "alternativeNames"):
        for item in description.get(key) or []:
            name = _named_value((item or {}).get("fullName"))
            if name:
                return name
    return None


def _evidence_label(items: list[dict] | None) -> str:
    return ", ".join(
        " ".join(str(item[key]) for key in ("code", "source", "id") if item.get(key))
        for item in items or []
        if isinstance(item, dict)
    )


def _feature_position(feature: dict) -> str:
    start, end = feature.get("start"), feature.get("end")
    if start is None and end is None:
        return "position not returned"
    start_text = str(start if start is not None else end)
    end_text = str(end if end is not None else start)
    if feature.get("start_modifier") not in {None, "EXACT"}:
        start_text = str(feature["start_modifier"]) + " " + start_text
    if feature.get("end_modifier") not in {None, "EXACT"}:
        end_text = str(feature["end_modifier"]) + " " + end_text
    return start_text if start_text == end_text else f"{start_text}-{end_text}"


def _feature_detail(feature: dict) -> str:
    parts = [f"Residue(s) {_feature_position(feature)}"]
    if feature.get("description"):
        parts.append(str(feature["description"]))
    ligand = feature.get("ligand") or {}
    ligand_name = ligand.get("name") if isinstance(ligand, dict) else None
    if ligand_name:
        parts.append("ligand: " + str(ligand_name))
    alternative = feature.get("alternative_sequence") or {}
    if isinstance(alternative, dict):
        original = alternative.get("originalSequence")
        replacements = alternative.get("alternativeSequences") or []
        if original and replacements:
            parts.append(f"{original} → {', '.join(map(str, replacements))}")
    if feature.get("feature_id"):
        parts.append("feature: " + str(feature["feature_id"]))
    cross_references = feature.get("cross_references") or []
    if cross_references:
        parts.append(
            "cross-references: "
            + ", ".join(
                " ".join(
                    str(value)
                    for value in (item.get("database"), item.get("id"))
                    if value
                )
                for item in cross_references
                if isinstance(item, dict)
            )
        )
    evidence = _evidence_label(feature.get("evidence"))
    if evidence:
        parts.append("evidence: " + evidence)
    return "; ".join(parts)


def _feature_section(
    title: str, features: list[dict], types: set[str], *, empty: str
) -> str:
    selected = [item for item in features if item.get("type") in types]
    if not selected:
        return f"<h3>{html.escape(title)}</h3><p>{html.escape(empty)}</p>"
    rows = "".join(
        _row(str(item.get("type") or "Feature"), _feature_detail(item))
        for item in selected
    )
    return f"<h3>{html.escape(title)}</h3><table>{rows}</table>"


def _comment_texts(comment: dict) -> list[str]:
    values: list[str] = []
    for item in comment.get("texts") or []:
        if isinstance(item, dict) and item.get("value"):
            text = str(item["value"])
            evidence = _evidence_label([
                {
                    "code": evidence_item.get("evidenceCode"),
                    "source": evidence_item.get("source"),
                    "id": evidence_item.get("id"),
                }
                for evidence_item in item.get("evidences") or []
                if isinstance(evidence_item, dict)
            ])
            values.append(text + (f" (evidence: {evidence})" if evidence else ""))
    return values


def _text_comment_section(
    title: str, comments: list[dict], comment_types: set[str], *, empty: str
) -> str:
    values = [
        value
        for comment in comments
        if comment.get("commentType") in comment_types
        for value in _comment_texts(comment)
    ]
    if not values:
        return f"<h3>{html.escape(title)}</h3><p>{html.escape(empty)}</p>"
    return (
        f"<h3>{html.escape(title)}</h3><ul>"
        + "".join(f"<li>{html.escape(value)}</li>" for value in values)
        + "</ul>"
    )


def _pdb_section(cross_references: list[dict], comments: list[dict]) -> str:
    pdb_refs = [item for item in cross_references if item.get("database") == "PDB"]
    shown_refs = pdb_refs[:20]
    rows = []
    for item in shown_refs:
        properties = {
            str(prop.get("key")): prop.get("value")
            for prop in item.get("properties") or []
            if isinstance(prop, dict) and prop.get("key")
        }
        pdb_id = str(item.get("id") or "")
        detail_values = []
        for key in ("Method", "Resolution", "Chains"):
            value = properties.get(key)
            if value and value != "-":
                detail_values.append(str(value))
        details = " · ".join(detail_values) or "No method, resolution, or chain coverage returned"
        url = html.escape(
            "https://www.ebi.ac.uk/pdbe/entry/pdb/" + pdb_id.lower(), quote=True
        )
        rows.append(
            f"<tr><th scope='row'><a href=\"{url}\">{html.escape(pdb_id)}</a></th>"
            f"<td>{html.escape(details)}</td></tr>"
        )
    output = "<h3>Experimental structures</h3>"
    if rows:
        output += (
            "<p>UniProt PDB cross-references; chain ranges may cover only a construct "
            "or domain.</p><table>" + "".join(rows) + "</table>"
        )
        if len(pdb_refs) > len(shown_refs):
            output += (
                f"<p>Showing {len(shown_refs)} of {len(pdb_refs)} PDB cross-references; "
                "all are retained in report.json.</p>"
            )
    else:
        output += "<p>No experimental structure is cross-referenced in this UniProt entry.</p>"

    mass_spec: list[str] = []
    for comment in comments:
        if comment.get("commentType") != "MASS SPECTROMETRY":
            continue
        parts = []
        for key, label in (
            ("method", "method"), ("mass", "mass"), ("molWeight", "mass"),
            ("error", "error"), ("moleculeType", "molecule type"),
        ):
            if comment.get(key) is not None:
                parts.append(f"{label}: {comment[key]}")
        raw_positions = comment.get("positions") or comment.get("range") or []
        positions = [raw_positions] if isinstance(raw_positions, dict) else raw_positions
        position_text = [
            _feature_position({
                "start": (position.get("start") or {}).get("value"),
                "end": (position.get("end") or {}).get("value"),
            })
            for position in positions
            if isinstance(position, dict)
        ]
        if position_text:
            parts.append("positions: " + ", ".join(position_text))
        evidence = _evidence_label([
            {
                "code": item.get("evidenceCode"),
                "source": item.get("source"),
                "id": item.get("id"),
            }
            for item in comment.get("evidences") or []
            if isinstance(item, dict)
        ])
        if evidence:
            parts.append("evidence: " + evidence)
        mass_spec.append("; ".join(parts) or "Mass-spectrometry annotation returned")
    output += "<h3>Mass-spectrometry protocol</h3>"
    if mass_spec:
        output += "<ul>" + "".join(
            f"<li>{html.escape(value)}</li>" for value in mass_spec
        ) + "</ul>"
    else:
        output += "<p>Mass-spectrometry protocol details are not annotated in this entry.</p>"
    return output


def _localisation_section(comments: list[dict], features: list[dict]) -> str:
    locations: list[str] = []
    for comment in comments:
        if comment.get("commentType") != "SUBCELLULAR LOCATION":
            continue
        for item in comment.get("subcellularLocations") or []:
            parts = []
            evidence_items = []
            for key in ("location", "topology", "orientation"):
                node = item.get(key) or {} if isinstance(item, dict) else {}
                if node.get("value"):
                    parts.append(str(node["value"]))
                evidence_items.extend(node.get("evidences") or [])
            evidence = _evidence_label([
                {
                    "code": evidence_item.get("evidenceCode"),
                    "source": evidence_item.get("source"),
                    "id": evidence_item.get("id"),
                }
                for evidence_item in evidence_items
                if isinstance(evidence_item, dict)
            ])
            if evidence:
                parts.append("evidence: " + evidence)
            if parts:
                locations.append("; ".join(parts))
        locations.extend(_comment_texts(comment.get("note") or {}))
    topology_types = {"Transmembrane", "Intramembrane", "Topological domain", "Signal",
                      "Transit peptide"}
    topology = [item for item in features if item.get("type") in topology_types]
    if not locations and not topology:
        return (
            "<h3>Localisation</h3>"
            "<p>Subcellular localisation and membrane topology are not annotated in this entry.</p>"
        )
    output = "<h3>Localisation</h3>"
    if locations:
        output += "<ul>" + "".join(
            f"<li>{html.escape(value)}</li>" for value in locations
        ) + "</ul>"
    if topology:
        output += "<table>" + "".join(
            _row(str(item["type"]), _feature_detail(item)) for item in topology
        ) + "</table>"
    return output


def _catalytic_section(comments: list[dict]) -> str:
    reactions: list[str] = []
    for comment in comments:
        if comment.get("commentType") != "CATALYTIC ACTIVITY":
            continue
        reaction = comment.get("reaction") or {}
        parts = []
        if reaction.get("name"):
            parts.append(str(reaction["name"]))
        if reaction.get("ecNumber"):
            parts.append("EC " + str(reaction["ecNumber"]))
        reaction_refs = reaction.get("reactionCrossReferences") or []
        if reaction_refs:
            parts.append(
                "reaction cross-references: "
                + ", ".join(
                    " ".join(
                        str(value)
                        for value in (item.get("database"), item.get("id"))
                        if value
                    )
                    for item in reaction_refs
                    if isinstance(item, dict)
                )
            )
        physiological = comment.get("physiologicalReactions") or []
        for item in physiological:
            direction = item.get("directionType") if isinstance(item, dict) else None
            if direction:
                parts.append("physiological direction: " + str(direction))
        evidence = _evidence_label([
            {
                "code": item.get("evidenceCode"),
                "source": item.get("source"),
                "id": item.get("id"),
            }
            for item in reaction.get("evidences") or []
            if isinstance(item, dict)
        ])
        if evidence:
            parts.append("evidence: " + evidence)
        if parts:
            reactions.append("; ".join(parts))
    if not reactions:
        return (
            "<h3>Catalytic activity</h3>"
            "<p>Catalytic activity is not annotated in this UniProt entry.</p>"
        )
    return "<h3>Catalytic activity</h3><ul>" + "".join(
        f"<li>{html.escape(value)}</li>" for value in reactions
    ) + "</ul>"


def _ppi_section(interactions: list[dict], cross_references: list[dict]) -> str:
    shown_interactions = interactions[:50]
    rows = []
    for interaction in shown_interactions:
        partner = interaction.get("interactantTwo") or {}
        partner_label = partner.get("geneName") or partner.get("uniProtKBAccession")
        if not partner_label:
            partner_label = "interaction partner not named"
        details = []
        if partner.get("uniProtKBAccession"):
            details.append(str(partner["uniProtKBAccession"]))
        if partner.get("intActId"):
            details.append("IntAct " + str(partner["intActId"]))
        if interaction.get("numberOfExperiments") is not None:
            details.append(f"{interaction['numberOfExperiments']} experiment(s)")
        rows.append(_row(str(partner_label), "; ".join(details)))
    ppi_databases = {"IntAct", "BioGRID", "MINT", "DIP", "STRING", "ComplexPortal", "CORUM"}
    database_refs = [item for item in cross_references if item.get("database") in ppi_databases]
    if not rows and not database_refs:
        return (
            "<h3>Protein-protein interactions</h3>"
            "<p>No protein-interaction record or supported PPI database cross-reference is "
            "annotated in this entry.</p>"
        )
    output = "<h3>Protein-protein interactions</h3>"
    if rows:
        output += f"<table>{''.join(rows)}</table>"
        if len(interactions) > len(shown_interactions):
            output += (
                f"<p>Showing {len(shown_interactions)} of {len(interactions)} interaction "
                "records; all are retained in report.json.</p>"
            )
    if database_refs:
        output += "<p><strong>Database cross-references:</strong> " + html.escape(
            ", ".join(
                f"{item.get('database')} {item.get('id')}" for item in database_refs
            )
        ) + "</p>"
    return output


def _domain_section(features: list[dict], cross_references: list[dict], comments: list[dict]) -> str:
    boundary_types = {"Domain", "Repeat", "Zinc finger", "Coiled coil", "Region", "Motif"}
    boundaries = [item for item in features if item.get("type") in boundary_types]
    domain_refs = [
        item for item in cross_references if item.get("database") in {"InterPro", "Pfam"}
    ]
    domain_comments = [
        value
        for comment in comments
        if comment.get("commentType") == "DOMAIN"
        for value in _comment_texts(comment)
    ]
    if not boundaries and not domain_refs and not domain_comments:
        return (
            "<h3>Domain boundaries</h3>"
            "<p>Domain boundaries are not annotated in this UniProt entry.</p>"
        )
    output = "<h3>Domain boundaries</h3>"
    if boundaries:
        output += "<table>" + "".join(
            _row(str(item["type"]), _feature_detail(item)) for item in boundaries
        ) + "</table>"
    else:
        output += "<p>No residue-level domain boundaries are annotated in this entry.</p>"
    if domain_comments:
        output += "<ul>" + "".join(
            f"<li>{html.escape(value)}</li>" for value in domain_comments
        ) + "</ul>"
    if domain_refs:
        output += "<p><strong>Domain database cross-references:</strong> " + html.escape(
            ", ".join(f"{item.get('database')} {item.get('id')}" for item in domain_refs)
        ) + "</p>"
    return output


def _variant_section(record: dict, features: list[dict]) -> str:
    variants = [
        item for item in features
        if item.get("type") in {
            "Natural variant", "Mutagenesis", "Sequence conflict", "Alternative sequence"
        }
    ]
    output = "<h3>Variants</h3>"
    if variants:
        shown = variants[:50]
        output += "<table>" + "".join(
            _row(str(item["type"]), _feature_detail(item)) for item in shown
        ) + "</table>"
        if len(variants) > len(shown):
            output += (
                f"<p>Showing {len(shown)} of {len(variants)} UniProt variant features; "
                "all are retained in report.json.</p>"
            )
    else:
        output += "<p>Variants are not annotated in this UniProt entry.</p>"
    alphamissense = [
        item for item in record.get("cross_references") or []
        if item.get("database") == "AlphaMissense"
    ]
    if alphamissense:
        output += "<p><strong>AlphaMissense cross-references:</strong> " + html.escape(
            ", ".join(str(item.get("id")) for item in alphamissense if item.get("id"))
        ) + "</p>"
    elif (record.get("organism_details") or {}).get("taxonId") == 9606:
        output += (
            "<p><strong>AlphaMissense:</strong> AlphaMissense classifications were not "
            "returned by this UniProtKB REST entry, so they are not shown or inferred here.</p>"
        )
    return output


def _uniprot_section(record: dict) -> str:
    """Display fetched fields only after an exact canonical-sequence match."""
    if record["status"] == "not_requested":
        return (
            "<h2>UniProt annotations</h2><p>No accession was supplied. "
            "To add annotations, enter a UniProt accession with your sequence.</p>"
        )
    if record["status"] != "verified":
        candidates = record.get("candidates") or []
        suggestion = (
            " Candidate accessions: " + ", ".join(html.escape(acc) for acc in candidates) + "."
            if candidates else ""
        )
        return (
            "<h2>UniProt annotations</h2><p>Not included: "
            + html.escape(record["reason"]) + suggestion + "</p>"
        )
    features = record.get("features") or []
    comments = record.get("comments") or []
    cross_references = record.get("cross_references") or []
    annotation_sections = "".join([
        _pdb_section(cross_references, comments),
        _localisation_section(comments, features),
        _feature_section(
            "Glycosylation sites",
            features,
            {"Glycosylation"},
            empty="Glycosylation sites are not annotated in this UniProt entry.",
        ),
        _feature_section(
            "Post-translational modifications (PTMs)",
            features,
            {
                "Modified residue", "Lipidation", "Disulfide bond", "Cross-link",
                "Non-standard residue",
            },
            empty="Residue-level PTMs are not annotated in this UniProt entry.",
        ),
        _text_comment_section(
            "PTM annotation notes",
            comments,
            {"PTM"},
            empty="No additional PTM comment is annotated in this entry.",
        ),
        _text_comment_section(
            "Stoichiometry / subunit organisation",
            comments,
            {"SUBUNIT"},
            empty="Subunit organisation or stoichiometry is not annotated in this entry.",
        ),
        _domain_section(features, cross_references, comments),
        _feature_section(
            "Functional / active sites",
            features,
            {
                "Active site", "Binding site", "Site", "DNA binding", "Nucleotide binding",
                "Metal binding", "Calcium binding",
            },
            empty="Functional or active-site residues are not annotated in this entry.",
        ),
        _text_comment_section(
            "Function annotation",
            comments,
            {"FUNCTION"},
            empty="A free-text function comment is not annotated in this entry.",
        ),
        _catalytic_section(comments),
        _ppi_section(record.get("protein_interactions") or [], cross_references),
        _variant_section(record, features),
    ])
    absorption = record.get("absorption")
    if absorption:
        maximum = absorption.get("max_nm")
        notes = "; ".join(absorption.get("note") or [])
        absorption_text = (
            f"Annotated absorption maximum: {maximum} nm" if maximum is not None
            else "Absorption annotated (no maximum returned)"
        )
        if notes:
            absorption_text += ". " + notes
        evidence = _evidence_label(absorption.get("evidence"))
        if evidence:
            absorption_text += " (evidence: " + evidence + ")"
        annotation_sections += (
            "<h3>Other biophysical annotations</h3><p>"
            + html.escape(absorption_text)
            + "</p>"
        )
    else:
        annotation_sections += (
            "<h3>Other biophysical annotations</h3>"
            "<p>Absorption is not annotated in this UniProt entry.</p>"
        )
    annotation_sections += (
        "<p><small>Residue-level features are displayed selectively by category; "
        "the full normalized feature list is retained in report.json.</small></p>"
    )
    url = html.escape(record["source_url"], quote=True)
    identity = html.escape(
        f"{record['accession']} · {record['organism'] or 'organism not provided'} "
        f"· {record['entry_type'] or 'review status not provided'}"
    )
    date = html.escape(record["retrieved_at"] or "retrieval date not available")
    protein_name = _protein_name(record.get("protein_description"))
    gene_names = [
        _named_value((gene or {}).get("geneName")) for gene in record.get("genes") or []
    ]
    gene_names = [name for name in gene_names if name]
    identity_details = ""
    if protein_name:
        identity_details += "<p><strong>Protein:</strong> " + html.escape(protein_name) + "</p>"
    if gene_names:
        identity_details += (
            "<p><strong>Gene name(s):</strong> "
            + html.escape(", ".join(gene_names))
            + "</p>"
        )
    comment_types = sorted({
        str(item.get("commentType")) for item in record.get("comments") or []
        if isinstance(item, dict) and item.get("commentType")
    })
    if comment_types:
        identity_details += (
            "<p><strong>Annotation sections retrieved:</strong> "
            + html.escape(", ".join(comment_types))
            + ". Full fetched fields are included in report.json.</p>"
        )
    if record.get("identity_basis") == "supplied_accession":
        identity_statement = (
            f"<p>Supplied accession resolved to: <a href=\"{url}\">{identity}</a> "
            f"(retrieved {date}). Calculations and positions use its canonical sequence.</p>"
        )
    else:
        identity_statement = (
            f"<p>Canonical sequence matched exactly: <a href=\"{url}\">{identity}</a> "
            f"(retrieved {date}). Positions refer to that canonical sequence.</p>"
        )
    if record.get("input_note"):
        identity_statement += (
            "<p><strong>Input note:</strong> "
            + html.escape(str(record["input_note"]))
            + "</p>"
        )
    return (
        "<h2>UniProt annotations</h2>"
        + identity_statement
        + identity_details
        + annotation_sections
    )


def _pdbe_section(record: dict) -> str:
    """Render only structures mapped by PDBe to the sequence-verified accession."""
    if record["status"] != "verified":
        return (
            "<h2>Experimental structures (PDBe)</h2><p>Not included: "
            + html.escape(record.get("reason") or "UniProt identity was not verified.")
            + "</p>"
        )
    rows = []
    for item in record["structures"]:
        pdb_id = html.escape(item["pdb_id"])
        url = html.escape(item["source_url"], quote=True)
        method = item.get("experimental_method") or "method not returned"
        resolution = item.get("resolution_angstrom")
        detail = html.escape(str(method))
        if resolution is not None:
            detail += f" · {html.escape(str(resolution))} Å"
        chain = item.get("chain_id")
        if chain:
            detail += " · chain " + html.escape(str(chain))
        start, end = item.get("uniprot_start"), item.get("uniprot_end")
        if start is not None and end is not None:
            detail += f" · UniProt residues {html.escape(str(start))}-{html.escape(str(end))}"
        coverage = item.get("coverage_fraction")
        if isinstance(coverage, (int, float)):
            detail += f" · {coverage:.0%} coverage"
        rows.append(
            f"<tr><th scope='row'><a href=\"{url}\">{pdb_id}</a></th><td>{detail}</td></tr>"
        )
    return (
        "<h2>Experimental structures (PDBe)</h2>"
        "<p>Ranked PDBe mappings for the exact sequence-verified UniProt accession. "
        "Coverage can represent a construct or domain rather than the full protein.</p>"
        f"<table>{''.join(rows)}</table>"
    )


def render_html(data: dict, *, generated_at: str) -> str:
    """Render only escaped input/data; never execute arbitrary HTML from a FASTA header."""
    ext = data["molar_extinction_coefficient_M-1cm-1"]
    rows = [
        _row("Length", f"{data['length_residues']} residues"),
        _row("Molecular weight", f"{data['molecular_weight_da']} Da"),
        _row("Isoelectric point (pI)", data["isoelectric_point"]),
        _row("Net charge at pH 7", data["net_charge_pH_7"]),
        _row("Average hydropathy (GRAVY)", data["gravy_kyte_doolittle"]),
        _row("Aromatic amino-acid fraction", f"{data['aromaticity_iron_dos_santos']:.1%}"),
        _row("Instability index (sequence estimate)", data["instability_index_guruprasad"]),
        _row("A280 extinction coefficient (reduced)", f"{ext['reduced']} M⁻¹ cm⁻¹"),
        _row("A280 extinction coefficient (oxidised)", f"{ext['oxidised']} M⁻¹ cm⁻¹"),
    ]
    composition = "".join(
        _row(aa, f"{item['number']} residues · {item['mole_percent']}% · "
             f"{item['dayhoff_stat']}x reference frequency")
        for aa, item in sorted(data["amino_acid_dayhoff_stats"].items())
        if aa in STANDARD_AA
    )
    charge = "".join(
        _row(f"pH {ph}", value)
        for ph, value in data["charge_vs_ph"]["ph_to_net_charge"].items()
    )
    classes = "".join(
        _row(name, f"{values['mole_percent']}%")
        for name, values in data["physico_chemical_classes"].items()
    )
    disorder = data["disorder_prediction"]
    if disorder["status"] == "predicted":
        regions = disorder["disordered_regions"]
        range_text = ", ".join(f"{item['start']}-{item['end']}" for item in regions)
        region_text = f"Predicted regions: {range_text}." if regions else "No regions called by this predictor."
        disorder_section = (
            "<h2>Predicted disorder</h2>"
            f"<p>{html.escape(region_text)} This is a sequence-based prediction, not a "
            "measured structure.</p>"
            '<img src="disorder.png" alt="Predicted disorder scores along the sequence">'
            '<p><a href="disorder.png" download>Download the disorder plot</a></p>'
        )
    else:
        disorder_section = (
            "<h2>Predicted disorder</h2>"
            f"<p>Not available: {html.escape(disorder['reason'])}</p>"
        )
    motifs = data["chromophore_candidates"]
    motif_section = ""
    if motifs["status"] == "candidates_found":
        positions = sorted({pos for group in motifs["candidates"].values() for pos in group})
        motif_section = (
            "<details><summary>Possible chromophore motif matches</summary>"
            f"<p>Pattern matches start at residue(s) {', '.join(map(str, positions))}. "
            "A match alone does not establish that the protein absorbs light or forms "
            "a chromophore. Verify against an identified, annotated protein entry.</p></details>"
        )
    ext_mg = data["a280_extinction_per_mg_ml_M-1cm-1"]
    estimates = "".join([
        _row("Average residue weight", f"{data['average_residue_weight_da']} Da"),
        _row("A280 per 1 mg/mL (reduced / oxidised)",
             f"{ext_mg['reduced']} / {ext_mg['oxidised']}"),
    ])
    structure = data["secondary_structure_composition_estimate"]
    structure_rows = "".join(
        _row(label.title(), f"{value:.1%}") for label, value in structure.items()
    )
    expression = data["inclusion_body_propensity"]
    expression_flag = "higher" if expression["class"] == "likely_insoluble" else "lower"
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Sequence analysis report</title><style>
body{{font:1rem/1.5 system-ui,sans-serif;max-width:850px;margin:2rem auto;padding:0 1rem;color:#17243a}}
h1,h2{{line-height:1.2}}table{{border-collapse:collapse;width:100%;margin-bottom:2rem}}
th,td{{padding:.45rem;border-bottom:1px solid #ccd;text-align:left;vertical-align:top}}
th{{width:48%}}img{{max-width:100%;height:auto}}.note{{background:#eef2f7;padding:1rem}}
details{{margin:1rem 0;padding:.5rem;border:1px solid #ccd}}summary{{cursor:pointer;font-weight:600}}
a{{color:#0057b8}}
</style></head><body><h1>Sequence analysis report</h1>
<p><strong>Your sequence:</strong> {html.escape(data['header'])} · {data['length_residues']} residues</p>
<p>Report created {html.escape(generated_at)} UTC. <a href="../../">Analyse another sequence</a></p>
<h2>At a glance</h2><table>{''.join(rows)}</table>
<h2>Hydropathy along the sequence</h2>
<p>Positive values are more hydrophobic; negative values are more hydrophilic.</p>
<img src="hydropathy.png" alt="Hydropathy profile across the submitted sequence">
<p><a href="hydropathy.png" download>Download the hydropathy image</a></p>
{disorder_section}
<details><summary>Amino-acid composition</summary>
<p>Counts, share of the sequence, and enrichment relative to a reference composition.</p>
<table>{composition}</table></details>
<details><summary>Charge at different pH values</summary><table>{charge}</table></details>
<details><summary>Residue groups</summary><table>{classes}</table>
<p>Groups overlap, so percentages do not add to 100%.</p></details>
<details><summary>Other sequence estimates</summary><table>{estimates}</table>
<p>Secondary-structure composition estimate (not a structure prediction):</p>
<table>{structure_rows}</table>
<p>Exploratory E. coli inclusion-body composition flag: {expression_flag}. This simple
sequence-only flag does not establish actual expression or solubility.</p></details>
{motif_section}
{_uniprot_section(data['uniprot'])}
{_pdbe_section(data['pdbe'])}
<h2>Download all calculated results</h2>
<p><a href="report.json" download>Download full results (JSON)</a> — includes additional
composition statistics, disorder scores (if available), and calculation details.</p>
<div class="note"><strong>What this report cannot tell you:</strong> These are sequence-derived
estimates, not measurements. A UniProt entry is shown after resolving a supplied accession
or verifying an exact sequence match; neither establishes an expression outcome.
A missing annotation or prediction is not evidence that a property is absent.</div></body></html>
"""


def generate_report(
    input_path: Path | None, output_dir: Path, *, accession: str = ""
) -> Path:
    raw_input = ""
    if input_path is not None:
        if input_path.stat().st_size > MAX_INPUT_BYTES:
            raise SequenceError("input exceeds the 100 KB limit")
        raw_input = input_path.read_text(encoding="utf-8")

    accession = accession.strip().upper()
    if accession:
        uniprot = fetch_annotations("", accession, accession_precedence=True)
        if uniprot.get("status") != "verified":
            raise SequenceError(
                "the supplied UniProt accession could not provide a canonical sequence: "
                + str(uniprot.get("reason") or uniprot.get("status"))
            )
        sequence = str(uniprot.pop("canonical_sequence"))
        header = f"UniProt {uniprot['accession']} canonical sequence"
        input_source = "canonical sequence retrieved from the supplied UniProt accession"
        if raw_input.strip():
            try:
                submitted_header, submitted_sequence = parse_sequence(raw_input)
            except SequenceError:
                uniprot["submitted_sequence_match"] = None
                uniprot["input_note"] = (
                    "The submitted sequence could not be parsed; accession precedence was "
                    "applied and calculations use the UniProt canonical sequence."
                )
                input_source += "; unparseable submitted sequence was not analysed"
            else:
                matches = submitted_sequence == sequence
                uniprot["submitted_sequence_match"] = matches
                if matches:
                    header = submitted_header
                    input_source += "; submitted sequence matched exactly"
                else:
                    uniprot["input_note"] = (
                        "The submitted sequence differed from the supplied accession; accession "
                        "precedence was applied and calculations use the UniProt canonical sequence."
                    )
                    input_source += "; differing submitted sequence was not analysed"
    else:
        if not raw_input.strip():
            raise SequenceError("provide either a protein sequence or a UniProt accession")
        header, sequence = parse_sequence(raw_input)
        uniprot = find_sequence_annotations(sequence)
        input_source = "user-supplied sequence; exact UniProt identity assessed separately"

    output_dir.mkdir(parents=True, exist_ok=True)
    data = compute_properties(sequence, header=header).as_dict()
    data["input_source"] = input_source
    data["generated_at_utc"] = datetime.now(UTC).isoformat(timespec="seconds")
    plot_hydropathy(sequence, label=header, out_path=output_dir / "hydropathy.png")
    data["disorder_prediction"] = predict_disorder(sequence, output_dir / "disorder.png")
    data["chromophore_candidates"] = chromophore_candidates(sequence)
    data["uniprot"] = uniprot
    data["pdbe"] = (
        fetch_pdbe_structures(data["uniprot"]["accession"])
        if data["uniprot"]["status"] == "verified"
        else {
            "status": "not_requested",
            "reason": "A verified UniProt accession is required before PDBe lookup.",
        }
    )
    (output_dir / "report.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    report = output_dir / "report.html"
    report.write_text(render_html(data, generated_at=data["generated_at_utc"]), encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a report from a protein sequence or UniProt accession"
    )
    parser.add_argument("--input", type=Path, help="optional single FASTA record or bare sequence")
    parser.add_argument("--output", type=Path, required=True, help="report directory")
    parser.add_argument(
        "--accession",
        default="",
        help="UniProt accession; takes precedence over --input when both are supplied",
    )
    args = parser.parse_args()
    try:
        report = generate_report(args.input, args.output, accession=args.accession)
    except (SequenceError, OSError, UnicodeError, ValueError) as exc:
        print(f"Report not generated: {exc}", file=sys.stderr)
        return 2
    print(f"Report generated: {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
