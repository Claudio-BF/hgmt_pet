SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
REPO_ROOT=$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)
python3 $REPO_ROOT/simulation_materials/make_hgmt_detector_derenzo.py $REPO_ROOT/simulation_materials/hgmt_detector_derenzo.topas
topas $REPO_ROOT/simulation_materials/hgmt_detector_derenzo.topas
mv $REPO_ROOT/HGMTDerenzo $REPO_ROOT/data/
rm $REPO_ROOT/HGMTDerenzo.phsp
rm $REPO_ROOT/HGMTDerenzo.header
