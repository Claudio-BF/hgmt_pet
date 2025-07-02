SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
REPO_ROOT=$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)
$REPO_ROOT/rendering/imager $REPO_ROOT/data/HGMTDerenzo.lor $REPO_ROOT/data/ 150
python3 $REPO_ROOT/rendering/render.py $REPO_ROOT/data/image.voxels 0
