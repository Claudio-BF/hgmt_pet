SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
REPO_ROOT=$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)
$REPO_ROOT/neural_network/training_data $REPO_ROOT/data/HGMTDerenzo $REPO_ROOT/simulation_materials/kapton_effs.csv
python3 trainer.py
