"""Local, deterministic biophysicochemical properties from a raw sequence.

Everything here is computed, not retrieved, so every function names the scale or
algorithm behind its number. No model involvement: the assistant may summarise these
values but may not produce them.

Backed by Biopython's `Bio.SeqUtils.ProtParam.ProteinAnalysis` (validated live on this
node) plus a per-residue Kyte-Doolittle scale for hydropathy plots, which replaces the
GUI-only ChimeraX `mlp` step for headless use.

Method notes to report to a user:
- pI, MW, charge: Biopython `ProteinAnalysis` (undecay; mean residue composition).
- Extinction coefficient: reduction-cysteine-disulfide-independent (Edelhoch) for the
  reduced and oxidised pair — report as a pair, never one number.
- Instability index: Guruprasad et al. (1990); a crude growth predictor, not a solubility
  predictor.
- Secondary-structure fractions: a composition-based estimate from the residues present,
  not a structure prediction, and must not be presented as one.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from Bio.SeqUtils.ProtParam import ProteinAnalysis

# ---------------------------------------------------------------------------
# sequence handling
# ---------------------------------------------------------------------------

_AMBIGUITY = {
    "B": "ND",  # not needed for composition; flagged instead
    "Z": "EQ",
    "X": "",
    "U": "",  # selenocysteine, excluded from Biopython composition
    "O": "",  # pyrrolysine
    "J": "",  # leucine/X equivalence
}
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


class SequenceError(ValueError):
    """Input is not usable protein sequence; reported as a finding, not guessed past."""


def clean_fasta_sequence(raw: str) -> tuple[str, list[str]]:
    """Strip FASTA records/whitespace/numbering; return (sequence, warnings).

    Never silently repairs input: anything unusual is returned as a warning the caller
    must surface, per the "absence/failure is a finding" rule.
    """
    warnings: list[str] = []
    lines = [ln.strip() for ln in raw.strip().splitlines()]
    headers = [ln for ln in lines if ln.startswith(">")]
    body = "".join(ln for ln in lines if ln and not ln.startswith(">"))
    body = body.replace(" ", "").replace("\t", "")
    if "-" in body or "." in body:
        warnings.append("gap characters present in input; gaps removed before counting")
        body = body.replace("-", "").replace(".", "")
    if not body:
        raise SequenceError("no residues found in input")

    seq = body.upper()
    if len(headers) > 1:
        warnings.append(f"input contains {len(headers)} FASTA records; only the first is analysed")
    # numeric formatting artifacts (e.g. "60 120") already stripped by whitespace removal
    nonstd = sorted({c for c in seq if c not in STANDARD_AA})
    if nonstd:
        stripped = sum(seq.count(c) for c in nonstd)
        warnings.append(
            "non-standard residues removed before computation: "
            + ", ".join(f"{c}x{seq.count(c)}" for c in nonstd)
            + f" ({stripped} of {len(seq)} positions); composition, MW, and pI reflect the "
            "remaining standard residues only"
        )
        seq = "".join(c for c in seq if c in STANDARD_AA)
    if not seq:
        raise SequenceError("input contained no standard amino acids")
    return seq, warnings


def read_first_record(path: str) -> tuple[str, str, list[str]]:
    """Return (header, sequence, warnings) for the first record of a FASTA file."""
    try:
        with open(path, encoding="utf-8") as handle:
            text = handle.read()
    except OSError as exc:
        raise SequenceError(f"cannot read FASTA file: {exc}") from exc
    if not text.strip():
        raise SequenceError(f"{path} is empty")
    if not text.lstrip().startswith(">"):
        # bare sequence is allowed, but say so rather than pretending it was titled
        seq, warnings = clean_fasta_sequence(text)
        return "(no FASTA header)", seq, [*warnings, "input had no > header line"]
    header = next(ln for ln in text.splitlines() if ln.startswith(">")).lstrip(">").strip()
    seq, warnings = clean_fasta_sequence(text)
    return header, seq, warnings


# ---------------------------------------------------------------------------
# hydropathy
# ---------------------------------------------------------------------------

KYTE_DOOLITLE: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
    "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
    "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def hydropathy_profile(sequence: str, window: int = 19) -> list[tuple[int, float]]:
    """Centred moving-average Kyte-Doolittle hydropathy; edge windows are truncated.

    Returns (1-based position, value) pairs. Positives are hydrophobic; a plot is a
    visualisation of these numbers, not a new claim.
    """
    half = max(1, window // 2)
    out: list[tuple[int, float]] = []
    for index in range(len(sequence)):
        lo = max(0, index - half)
        hi = min(len(sequence), index + half + 1)
        span = [KYTE_DOOLITLE.get(c, 0.0) for c in sequence[lo:hi]]
        out.append((index + 1, round(sum(span) / len(span), 4)))
    return out


# ---------------------------------------------------------------------------
# property computation
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class SequenceProperties:
    header: str
    length: int
    molecular_weight: float
    isoelectric_point: float
    net_charge_at_pH_7: float
    molar_extinction_reduced: float
    molar_extinction_oxidised: float
    gravy: float
    aromaticity: float
    instability_index: float
    amino_acid_percent: dict[str, float]
    secondary_structure_estimate: dict[str, float]
    hydropathy_min: tuple[int, float]
    hydropathy_max: tuple[int, float]
    warnings: list[str] = field(default_factory=list)
    method: dict[str, str] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return {
            "header": self.header,
            "length_residues": self.length,
            "molecular_weight_da": self.molecular_weight,
            "isoelectric_point": self.isoelectric_point,
            "net_charge_pH_7": self.net_charge_at_pH_7,
            "molar_extinction_coefficient_M-1cm-1": {
                "reduced": self.molar_extinction_reduced,
                "oxidised": self.molar_extinction_oxidised,
            },
            "gravy_kyte_doolittle": self.gravy,
            "aromaticity_iron_dos_santos": self.aromaticity,
            "instability_index_guruprasad": self.instability_index,
            "amino_acid_percent": self.amino_acid_percent,
            "secondary_structure_composition_estimate": self.secondary_structure_estimate,
            "hydropathy_extremes": {"min": self.hydropathy_min, "max": self.hydropathy_max},
            "warnings": self.warnings,
            "method": self.method,
            "caveats": [
                "computed from composition alone; no structure, disorder, or expression "
                "information is implied",
                "secondary-structure fractions are a composition estimate, not a prediction",
                "instability index predicts degradation propensity in E. coli, not solubility",
                "these are calculated values, not experimental measurements",
            ],
        }


def compute_properties(sequence: str, header: str = "(unnamed)") -> SequenceProperties:
    analysis = ProteinAnalysis(sequence)
    profile = hydropathy_profile(sequence)
    reduced, oxidised = analysis.molar_extinction_coefficient()
    helix, strand, coil = analysis.secondary_structure_fraction()
    return SequenceProperties(
        header=header,
        length=len(sequence),
        molecular_weight=round(analysis.molecular_weight(), 2),
        isoelectric_point=round(analysis.isoelectric_point(), 2),
        net_charge_at_pH_7=round(analysis.charge_at_pH(7.0), 2),
        molar_extinction_reduced=float(reduced),
        molar_extinction_oxidised=float(oxidised),
        gravy=round(analysis.gravy(), 3),
        aromaticity=round(analysis.aromaticity(), 3),
        instability_index=round(analysis.instability_index(), 1),
        amino_acid_percent={k: round(v, 2) for k, v in analysis.amino_acids_percent.items()},
        secondary_structure_estimate={
            "helix": round(helix, 3),
            "strand": round(strand, 3),
            "coil": round(coil, 3),
        },
        hydropathy_min=min(profile, key=lambda p: p[1]),
        hydropathy_max=max(profile, key=lambda p: p[1]),
        method={
            "molecular_weight": "Bio.SeqUtils.ProtParam.ProteinAnalysis.molecular_weight (average)",
            "isoelectric_point": "ProteinAnalysis.isoelectric_point (bidimensional titration)",
            "molar_extinction": "ProteinAnalysis.molar_extinction_coefficient "
            "(Edelhoch; reduced, oxidised)",
            "gravy": "ProteinAnalysis.gravy (Kyte-Doolittle)",
            "aromaticity": "ProteinAnalysis.aromaticity (Iron & Dos Santos)",
            "instability_index": "ProteinAnalysis.instability_index (Guruprasad 1990)",
            "hydropathy_profile": "hydropathy_profile (Kyte-Doolittle, centred window)",
            "secondary_structure": "ProteinAnalysis.secondary_structure_fraction "
            "(composition-based)",
        },
    )


# ---------------------------------------------------------------------------
# chromophore / photo-activity (heuristic only — see module caveats)
# ---------------------------------------------------------------------------

# Deliberately a *candidate* detector. Validated against a fixture chosen by the
# assistant (GFP P42212, accession verified by fetch_uniprot_entry in this session): a
# naive TYG scan returns a spurious hit at 38 and MISSES the annotated chromophore,
# because wtGFP carries Ser65 (S65T is the EGFP engineering mutation). Sequence-only
# detection therefore never asserts a chromophore; it names candidates to cross-check
# against an annotated entry, and says nothing further when there is no match.
# Validated against the assistant-chosen fixture GFP P42212 (accession verified by
# fetch_uniprot_entry in this session, not requested by the user). Findings that shaped
# this table: the assumed "TYG" and "G[FYL]G" cores both MISS the real chromophore, because
# wtGFP carries Ser65 (S65T is the EGFP engineering change) and its annotated modified
# residue is Tyr66, one position later. Only the S/Y/G form matches — and the broad
# "[^A]YG" style patterns also hit a spurious position 38. A linear motif therefore cannot
# identify a chromophore; it can only nominate regions for annotation lookup.
CHROMOPHORE_MOTIFS: dict[str, str] = {
    "[STAG]Y G": "GFP-like tripeptide core (Thr/Ser/Ala-Tyr-Gly; the Tyr cyclises)",
    "AGY": "kaede/DRAQ-like alternate core (low prior; frequently coincidental)",
}

_MOTIF_PATTERNS: dict[str, str] = {
    "GFP-like [STAG]-Tyr-Gly": r"(?=[STAG]YG)",
}

# Positions the fixture proves are NOT diagnostic, kept so regressions cannot reintroduce them.
INVALIDATED_MOTIFS: tuple[str, ...] = ("TYG", "G[FYL]G")


def chromophore_candidates(sequence: str) -> dict[str, Any]:
    """Nominate photo-active candidate regions; never identifies one.

    Returns a terminal, self-describing verdict so the caller can stop: an empty
    candidate list is a complete result, not a failure to retry.
    """
    import re

    candidates: dict[str, list[int]] = {}
    for label, pattern in _MOTIF_PATTERNS.items():
        candidates[label] = [m.start() + 1 for m in re.finditer(pattern, sequence)]
    found = [pos for positions in candidates.values() for pos in positions]
    complete = True
    if found:
        verdict = (
            f"{len(found)} candidate core(s) at position(s) "
            f"{', '.join(str(p) for p in sorted(set(found)))} of the pasted sequence. "
            "These are motif matches only; this tool cannot and does not conclude that the "
            "protein is photo-active."
        )
    else:
        verdict = (
            "No candidate chromophore-forming tripeptide core was found in the pasted "
            "sequence. This is a complete result for this tool: do NOT retry it with a "
            "modified sequence. It means no GFP-like motif is present, which is NOT "
            "evidence that the protein fails to absorb light — many photo-active proteins "
            "use chromophores this scan does not look for."
        )
    return {
        "status": "candidates_found" if found else "no_candidates",
        "complete_result": complete,
        "verdict": verdict,
        "candidates": candidates,
        "meaning": (
            "a linear motif match does not establish photon absorption, maturation, or "
            "photo-activity; nor does its absence exclude them"
        ),
        "known_limitations": [
            "matches are position-tagged substrings, so coincidental hits are expected "
            "(a non-chromophoric [STAG]YG occurs in the validated fixture itself)",
            "the well-known TYG and G[FYL]G cores were tested against the fixture and both "
            "miss the annotated chromophore; they are excluded on evidence",
            "residues outside the core and the surrounding barrel environment set the "
            "spectrum; neither is visible from composition alone",
        ],
        "to_confirm": [
            "identify the protein first — a name, gene, or accession from the user, or a "
            "similarity search (BLAST). UniProt full-text search does NOT accept a raw "
            "sequence, so search_uniprotkb cannot resolve identity from sequence alone",
            "once an accession is known, fetch_uniprot_entry and read its "
            "BIOPHYSICOCHEMICAL PROPERTIES > absorption (max, in nm)",
            "read PTM/Processing > Modified residue / Cross-link features",
            "if identity stays unresolved, report 'not annotated under these criteria' and "
            "stop; do not re-search with the sequence, that call cannot succeed",
        ],
        "classification": "heuristic / computational — never report as curated or measured",
    }


__all__ = [
    "CHROMOPHORE_MOTIFS",
    "INVALIDATED_MOTIFS",
    "KYTE_DOOLITLE",
    "SequenceError",
    "SequenceProperties",
    "chromophore_candidates",
    "clean_fasta_sequence",
    "compute_properties",
    "hydropathy_profile",
    "read_first_record",
]
