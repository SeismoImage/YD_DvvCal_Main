## YD DvvCal Main
# Phase1_ScorrCal
These scripts are used for SCF calculation and stacking

step1_scorr_cal.jl: this script is used for SCF calculation in julia
step2_scorr_stack.py: this script is used for SCF stacking

# Phase2_DvvCal
These scripts are used for $\Delta v/v$ measurement, you need to install the pycwt, follow the cross-wavelet-transform (https://github.com/Qhig/cross-wavelet-transform).

step1_WCS_measure.py: this script is used for $\Delta v/v$ measurements
step2_dvov_obatin.py: this script is used for $\Delta v/v$ stacking in different component-pairs
visual_dvov.ipynb: this script is used for plotting the $\Delta v/v$ measurements result
xwt.py: core function for the cross-wavelet-transform

# Phase3_DvvAnalysis
These scripts are used for the further analysis of $\Delta v/v$ measurement

step1_pore_pressure.ipynb: this script is used for pore pressure simulation