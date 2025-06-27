SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
REPO_ROOT=$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)
$REPO_ROOT/hgmt_lor_creator $REPO_ROOT/data/HGMTDerenzo $REPO_ROOT/simulation_materials/kapton_effs.csv $REPO_ROOT/data/
