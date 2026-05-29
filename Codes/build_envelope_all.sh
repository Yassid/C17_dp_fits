#!/usr/bin/env bash
# Convenience: regenerate every envelope artifact in one shot.
# Assumes results/<binning>[_<variant>]/ are all populated.
set -euo pipefail
BINNING=${1:-fine}
HERE="$(cd "$(dirname "$0")" && pwd)"

python3 "$HERE/plot_yield_envelope.py" --binning "$BINNING"
python3 "$HERE/plot_dsdo_envelope.py"  --binning "$BINNING"
python3 "$HERE/build_systematics_summary.py" --binning "$BINNING"
echo "envelope artifacts in plots/ and results/"
