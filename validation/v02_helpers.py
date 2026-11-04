"""Independent SBML reading for cross-checks, using regular expressions rather than the engine."""
from __future__ import annotations

import gzip
import re


def reaction_table(path) -> dict[str, dict[str, float]]:
    """Reaction identifier to metabolite coefficient."""
    path = str(path)
    text = gzip.open(path, "rt").read() if path.endswith(".gz") else open(path).read()
    table = {}
    for block in re.findall(r"<reaction .*?</reaction>", text, flags=re.S):
        rid = re.search(r'id="R_([^"]+)"', block).group(1)
        stoich: dict[str, float] = {}
        for tag, sign in (("listOfReactants", -1.0), ("listOfProducts", 1.0)):
            section = re.search(rf"<{tag}>(.*?)</{tag}>", block, flags=re.S)
            if not section:
                continue
            for ref in re.findall(r"<speciesReference [^>]*>", section.group(1)):
                met = re.search(r'species="M_([^"]+)"', ref).group(1)
                coef = re.search(r'stoichiometry="([^"]+)"', ref)
                stoich[met] = stoich.get(met, 0.0) + sign * (float(coef.group(1)) if coef else 1.0)
        table[rid] = stoich
    return table
