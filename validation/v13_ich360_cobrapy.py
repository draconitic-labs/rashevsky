#!/usr/bin/env python3
"""iCH360, the bundled model: every comparison of v01 against COBRApy on the same SBML file.

The authors publish no growth rate to match, so COBRApy is the reference here. Tolerances as in v01.
"""
from __future__ import annotations

import sys

from common import ICH360
from v01_core_cobrapy import compare, record

if __name__ == "__main__":
    sys.exit(record("v13_ich360_cobrapy", compare(ICH360)))
