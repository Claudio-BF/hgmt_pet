SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
REPO_ROOT=$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)
echo "Running Diagnostics"
$REPO_ROOT/hgmt_lor_creator $REPO_ROOT/data/HGMTDerenzo $REPO_ROOT/simulation_materials/kapton_effs.csv $REPO_ROOT/data/ -e0 -e1 -e2 -e3 -e4 | tee $REPO_ROOT/plots/detector_diagnostics.txt
$REPO_ROOT/scripts/hgmt_debug $REPO_ROOT/data/debug0.data 12 12 -hi | tee $REPO_ROOT/plots/scatter_detector.txt
python3 $REPO_ROOT/scripts/plot_bars.py $REPO_ROOT/plots/scatter_detector.txt Detector\ ID $REPO_ROOT/plots/scatters_by_detector.png
$REPO_ROOT/scripts/hgmt_debug $REPO_ROOT/data/debug1.data 12 12 -hi | tee $REPO_ROOT/plots/first_scatter_detector.txt
python3 $REPO_ROOT/scripts/plot_bars.py $REPO_ROOT/plots/first_scatter_detector.txt Detector\ ID $REPO_ROOT/plots/first_scatters_by_detector.png
python3 $REPO_ROOT/scripts/plot_histogram.py $REPO_ROOT/data/debug2.data Lor\ Center\ Error\ Transverse\ \(cm\) $REPO_ROOT/plots/lor_center_error_transverse.png 10 0.4
$REPO_ROOT/scripts/hgmt_debug $REPO_ROOT/data/debug2.data 5 20 -hi | tee $REPO_ROOT/plots/lor_center_error_transverse.txt
python3 $REPO_ROOT/scripts/plot_histogram.py $REPO_ROOT/data/debug3.data Lor\ Center\ Error\ Longitudinal\ \(cm\) $REPO_ROOT/plots/lor_center_error_longitudinal.png 10 0.2
$REPO_ROOT/scripts/hgmt_debug $REPO_ROOT/data/debug3.data 5 20 -hi | tee $REPO_ROOT/plots/lor_center_error_longitudinal.txt
python3 $REPO_ROOT/scripts/plot_multiple.py $REPO_ROOT/data/debug4.data Impact\ Parameter\ \(cm\) $REPO_ROOT/plots/lor_center_error_transverse_multiple.png 1 1
echo "All Tasks Complete"
