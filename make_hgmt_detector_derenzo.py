# Helper script to generate a hgmt_detector_derenzo file with a custom detector
import numpy as np
import os
import shutil


####################################
###          PARAMETERS          ###
####################################
file_name = "hgmt_detector_derenzo.topas"
dead_material = "Air"
active_material = "B33"
detector_thickness = 2.54  # cm
detector_inner_radii = np.array([45] * 12) + 5 * np.array(range(12))  # MUST BE SORTED
detector_volume_inner_rad = 45  # cm
detector_volume_outer_rad = 105  # cm
####################################
###             CODE             ###
####################################
with open(file_name, "w") as f:
    f.write(
        f"""
b:Ts/UseQt              = "False"
b:Gr/ViewA/Active	= "False"
b:Gr/ViewB/Active	= "False"
b:Gr/ViewC/Active	= "False"
b:Gr/ViewD/Active	= "False"

#sv:Ph/Default/Modules = 1 "g4em-standard_opt1"
#d:Ph/Default/CutForElectron = 0.5 um # overrides CutForAllParticles for Electron
sv:Ph/Default/Modules = 2 "g4em-standard_opt4" "g4em-penelope"

Ge/CheckForOverlaps = "True"
i:Ge/CheckForOverlapsResolution = 100000

i:Ts/ShowHistoryCountAtInterval = 100000
i:Ts/NumberOfThreads = 1


b:Ts/PauseBeforeQuit = "false"

b:Ma/B33/BuildFromMaterials = "True"
sv:Ma/B33/Components = 4 "G4_SILICON_DIOXIDE" "G4_BORON_OXIDE" "G4_SODIUM_MONOXIDE" "G4_ALUMINUM_OXIDE"
uv:Ma/B33/Fractions = 4 0.81 0.13 0.04 0.02
d:Ma/B33/Density = 2.23 g/cm3#1.67 g/cm3
s:Ge/DetectorVolume/Type	= "TsCylinder"
s:Ge/DetectorVolume/Parent = "AirBox"
d:Ge/DetectorVolume/HL		= 1 m
d:Ge/DetectorVolume/RMin	= {detector_volume_inner_rad} cm
d:Ge/DetectorVolume/RMax	= {detector_volume_outer_rad} cm
s:Ge/DetectorVolume/Material	= "{dead_material}"
d:Ge/DetectorVolume/MinStepSize 	= 0.01 mm
"""
    )
    for i in range(len(detector_inner_radii)):
        f.write(
            f"""
s:Ge/Detector_{i}/Type = "TsCylinder"
s:Ge/Detector_{i}/Material = "{active_material}"
s:Ge/Detector_{i}/Parent = "DetectorVolume"
d:Ge/Detector_{i}/Rmin = {detector_inner_radii[i]} cm 
d:Ge/Detector_{i}/RMax	= {detector_inner_radii[i] + detector_thickness} cm
d:Ge/Detector_{i}/HL = 1 m
d:Ge/Detector_{i}/MinStepSize 	= 0.01 mm
"""
        )
    f.write(
        """
includeFile = ./derenzo_300keV_100.topas

s:Ge/Derenzo/Parent = "AirBox"
background_e = 30
rod_e = 90



s:Ge/World/Material  = "Vacuum"
d:Ge/World/HLX       = 10.0 m
d:Ge/World/HLY       = 10.0 m
d:Ge/World/HLZ       = 10.0 m
b:Ge/World/Invisible = "TRUE"

s:Ge/AirBox/Material	= "Air"
s:Ge/AirBox/Type    	= "TsBox"
d:Ge/AirBox/HLX	= 10.0 m
d:Ge/AirBox/HLY	= 10.0 m
d:Ge/AirBox/HLZ	= 10.0 m
s:Ge/AirBox/Parent	= "world"
b:Ge/AirBox/Invisible	= "true"


i:Ts/Seed = 26

s:Sc/HGMTDerenzo/Quantity			= "AnnihilScorer"
s:Sc/HGMTDerenzo/Component			= "AirBox"
b:Sc/HGMTDerenzo/PropagateToChildren	= "True"
s:Sc/HGMTDerenzo/OutputType			= "binary"
s:Sc/HGMTDerenzo/IfOutputFileAlreadyExists	= "Overwrite"

s:Gr/ViewA/Type           = "OpenGL"
sv:Gr/ViewA/VisibleWorlds = 1 "All"
i:Gr/ViewA/WindowSizeX    = 900
i:Gr/ViewA/WindowSizeY    = 900
d:Gr/ViewA/Theta          = 45 deg
d:Gr/ViewA/Phi            = 0 deg
u:Gr/ViewA/Zoom	   = 5

s:Gr/ViewB/Type           = "OpenGL"
sv:Gr/ViewB/VisibleWorlds = 1 "All"
i:Gr/ViewB/WindowSizeX    = 900
i:Gr/ViewB/WindowSizeY    = 900
d:Gr/ViewB/Theta          = 45 deg
d:Gr/ViewB/Phi            = 45 deg
u:Gr/ViewB/Zoom	   = 1

s:Gr/ViewC/Type           = "OpenGL"
sv:Gr/ViewC/VisibleWorlds = 1 "All"
i:Gr/ViewC/WindowSizeX    = 900
i:Gr/ViewC/WindowSizeY    = 900
d:Gr/ViewC/Theta          = 90 deg
d:Gr/ViewC/Phi            = 0 deg
u:Gr/ViewC/Zoom	   = 40

s:Gr/ViewD/Type           = "OpenGL"
sv:Gr/ViewD/VisibleWorlds = 1 "All"
i:Gr/ViewD/WindowSizeX    = 900
i:Gr/ViewD/WindowSizeY    = 900
d:Gr/ViewD/Theta          = 0 deg
d:Gr/ViewD/Phi            = 0 deg
u:Gr/ViewD/Zoom	   = 2
"""
    )
