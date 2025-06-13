echo "Running Diagnostics"
../hgmt_lor_creator ../data/HGMTDerenzo ../simulation_materials/kapton_effs.csv ../data/ -e0 -e1 -e2 -e3 -e4 -d | tee ../plots/detector_diagnostics.txt
./hgmt_debug ../data/debug0.data 12 12 -hi | tee ../plots/scatter_detector.txt
python3 plot_bars.py ../plots/scatter_detector.txt Detector\ ID Scatters\ Occured
./hgmt_debug ../data/debug1.data 12 12 -hi | tee ../plots/first_scatter_detector.txt
python3 plot_bars.py ../plots/first_scatter_detector.txt Detector\ ID First\ Scatters\ Occured
python3 plot_histogram.py ../data/debug2.data Lor\ Center\ Error\ Transverse\ \(cm\) lor_center_error_transverse.png 10 0.4
./hgmt_debug ../data/debug2.data 5 20 -hi | tee ../plots/lor_center_error_transverse.txt
python3 plot_histogram.py ../data/debug3.data Lor\ Center\ Error\ Longitudinal\ \(cm\) lor_center_error_longitudinal.png 10 0.2
./hgmt_debug ../data/debug3.data 5 20 -hi | tee ../plots/lor_center_error_longitudinal.txt
python3 plot_multiple.py ../data/debug4.data Impact\ Parameter\ \(cm\) lor_center_error_transverse_mutliple.png 1 1
echo "All Tasks Complete"
