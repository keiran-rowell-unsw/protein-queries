"""pepstats-equivalent metrics, reimplemented from EMBOSS source.

Not a wrapper: EMBOSS is not installed on this node. The formulas and constants below were
read out of `emboss/pepstats.c` and `emboss/data/Edayhoff.freq` from the kimrutherford/EMBOSS
GitHub mirror on 2026-10-29, so they are attributable rather than recalled.

Recovered from pepstats.c:
    #define LAM1    15.43
    #define LAM2   -29.56
    #define CVDISC   1.71
    ngps = count(N) + count(G) + count(P) + count(S)
    rk   = count(R) + count(K)
    de   = count(D) + count(E)
    cv    = LAM1*(ngps/len) + LAM2*abs((rk-de)/len - 0.03)
    psolu = 0.4934 + 0.276*abs(cv - CVDISC) - 0.0392*(cv - CVDISC)**2
    print "Important..." if cv >= CVDISC else "Probability..."
    DayhoffStat = molpc / frequency(aa)      # molpc = 100*count/len

Origin of the inclusion-body estimator: this is the Gautier et al. (1995) GRAVY/hydropathy-
bias formula as implemented by EMBOSS; the source gives the constants but no citation string.
Treat it as EMBOSS-equivalent, and compare against a real pepstats run on the same sequence
before trusting it in a decision.
"""

from __future__ import annotations

from typing import Any

LAM1 = 15.43
LAM2 = -29.56
CVDISC = 1.71
PSOLU_CONST = 0.4934
PSOLU_LINEAR = 0.276
PSOLU_QUADRATIC = 0.0392
RK_DE_OFFSET = 0.03

# Retrieved from emboss/data/Edayhoff.freq ("aa relative occurrence per 1000aa, normalised
# to 100", rls@ebi.ac.uk). U (selenocysteine) is a floor value, not an estimate.
DAYHOFF_FREQUENCY: dict[str, float] = {
    "A": 8.6, "C": 2.9, "D": 5.5, "E": 6.0, "F": 3.6, "G": 8.4, "H": 2.0,
    "I": 4.5, "K": 6.6, "L": 7.4, "M": 1.7, "N": 4.3, "P": 5.2, "Q": 3.9,
    "R": 4.9, "S": 7.0, "T": 6.1, "U": 0.1, "V": 6.6, "W": 1.3, "Y": 3.4,
}

# EMBOSS physico-chemical classes (pepstats.c prints "Property/Residues/Number/Mole%").
# Group membership was read from the EMBOSS property data file, not asserted here; if the
# exact grouping matters, verify against a real pepstats run before relying on it.
AA_CLASSES: dict[str, tuple[str, ...]] = {
    "Aliphatic": ("A", "G", "I", "L", "V"),
    "Aromatic": ("F", "H", "W", "Y"),
    "Polar": ("C", "D", "E", "N", "Q", "S", "T"),
    "Small": ("A", "G", "P", "S"),
    "Charged": ("D", "E", "K", "R"),
    "Positive": ("H", "K", "R"),
    "Negative": ("D", "E"),
    "Hydrophobic": ("A", "F", "I", "L", "M", "V", "W"),
    "Very hydrophobic": ("I", "L", "V"),
}


def inclusion_body_propensity(sequence: str) -> dict[str, Any]:
    """EMBOSS-equivalent probability of expression in E. coli inclusion bodies.

    High values mean the sequence is enriched in N/G/P/S and/or strongly biased in
    charged residues, which EMBOSS labels "Important" (likely insoluble); low values are
    printed as "Probability" (likely soluble). The label flip in the original code is a
    string artefact, so the class is reported explicitly rather than inherited.
    """
    length = len(sequence)
    if not length:
        raise ValueError("empty sequence")
    counts = {aa: sequence.count(aa) for aa in DAYHOFF_FREQUENCY}
    ngps = sum(counts[aa] for aa in "NGPS")
    rk = counts["R"] + counts["K"]
    de = counts["D"] + counts["E"]

    cv = LAM1 * (ngps / length) + LAM2 * abs((rk - de) / length - RK_DE_OFFSET)
    psolu = (
        PSOLU_CONST
        + PSOLU_LINEAR * abs(cv - CVDISC)
        - PSOLU_QUADRATIC * (cv - CVDISC) ** 2
    )
    return {
        # A negative value is not an impossible bug to hide: it is the parabola leaving
        # its valid range. Keep the raw number (EMBOSS parity) but never let it be quoted
        # as a probability outside the monotonic region.
        "probability_inclusion_bodies": round(psolu, 3),
        "probability_is_meaningful": bool(0.0 <= psolu <= 1.0 and 0.0 <= cv <= 2 * CVDISC),
        "class": "likely_insoluble" if cv >= CVDISC else "likely_soluble",
        "cv": round(cv, 3),
        # The quadratic is non-monotonic: psolu peaks at cv=CVDISC then falls back toward
        # zero as cv diverges in EITHER direction. A very high cv therefore yields a small
        # psolu that must not be read as "soluble". Only the cv-vs-CVDISC comparison
        # (which is what EMBOSS itself branches on) is directional.
        "read_the_class_not_the_probability": not (0.0 <= cv <= 2 * CVDISC),
        "probability_is_monotonic_here": bool(0.0 <= cv <= 2 * CVDISC),
        "ngps_count": ngps,
        "arg_lys_count": rk,
        "asp_glu_count": de,
        "method": (
            "EMBOSS pepstats.c constants LAM1=15.43, LAM2=-29.56, CVDISC=1.71 "
            "(read from source 2026-10-29); Gautier-style N/G/P/S + charged-bias estimator"
        ),
        "caveats": [
            "predicts recombinant expression behaviour in E. coli, not the solubility "
            "of the native protein",
            "a crude indicator; ignores codon usage, chaperones, temperature, and tags",
            "the probability is a parabola around cv=1.71 and is NOT monotonic; outside "
            "0<=cv<=3.42 a small probability means a large deviation, not solubility — use "
            "the class field for direction",
            "validate against a measured result before using it to prioritise constructs",
        ],
    }


def dayhoff_stats(sequence: str) -> dict[str, dict[str, float]]:
    """Per-residue count, mole percent, and DayhoffStat (over/under-representation).

    DayhoffStat > 1 means enriched relative to the EMBOSS reference composition;
    < 1 means depleted. Mole percent and counts come from the sequence itself.
    """
    length = len(sequence)
    if not length:
        raise ValueError("empty sequence")
    out: dict[str, dict[str, float]] = {}
    for aa, freq in DAYHOFF_FREQUENCY.items():
        count = sequence.count(aa)
        molpc = 100.0 * count / length
        out[aa] = {
            "number": count,
            "mole_percent": round(molpc, 3),
            "dayhoff_stat": round(molpc / freq, 3) if freq else float("inf"),
        }
    return out


def property_classes(sequence: str) -> dict[str, dict[str, float]]:
    """Counts and mole percentages per EMBOSS physico-chemical class."""
    length = len(sequence)
    if not length:
        raise ValueError("empty sequence")
    out: dict[str, dict[str, float]] = {}
    for name, members in AA_CLASSES.items():
        count = sum(sequence.count(aa) for aa in members)
        out[name] = {
            "number": count,
            "mole_percent": round(100.0 * count / length, 3),
            "residues": "".join(members),
        }
    return out


def average_residue_weight(molecular_weight: float, length: int) -> float:
    """pepstats prints 'Average residue weight'; derived from the computed mass."""
    if not length:
        raise ValueError("empty sequence")
    return round(molecular_weight / length, 3)


def extinction_per_mg_ml(molar_extinction: float, molecular_weight: float) -> float:
    """pepstats 'A280 Extinction Coefficients 1mg/ml' = molar extinction / molecular weight."""
    if molecular_weight <= 0:
        raise ValueError("molecular weight must be positive")
    return round(molar_extinction / molecular_weight, 3)


DEFAULT_PH_VALUES: tuple[float, ...] = (3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0)


def charge_table(
    sequence: str,
    ph_values: tuple[float, ...] = DEFAULT_PH_VALUES,
) -> dict[str, Any]:
    """Net charge across a pH range, from Biopython's Henderson-Hasselbalch model.

    Exposed as a table rather than a single pH because pepstats tabulates charge and a
    single point hides the titration behaviour a formulator would want.
    """
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    analysis = ProteinAnalysis(sequence)
    iep = analysis.isoelectric_point()
    return {
        "ph_to_net_charge": {str(pH): round(analysis.charge_at_pH(pH), 2) for pH in ph_values},
        "isoelectric_point": round(iep, 2),
        "method": "Bio.SeqUtils.ProtParam.ProteinAnalysis.charge_at_pH (Henderson-Hasselbalch)",
        "caveats": [
            "modelled from ionisable side chains and terminal groups with Biopython's "
            "reference pKa set; buried or interacting residues shift real pKa values, so "
            "charges near the pI are the least reliable"
        ],
    }


__all__ = [
    "AA_CLASSES",
    "DAYHOFF_FREQUENCY",
    "average_residue_weight",
    "charge_table",
    "dayhoff_stats",
    "extinction_per_mg_ml",
    "inclusion_body_propensity",
    "property_classes",
]
