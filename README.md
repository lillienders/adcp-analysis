# adcp-analysis
Analysis (wake law, logarithmic law, spectra) for 5-Beam Nortek ADCPs (MatLab version)

**Create Index of Tide Phases (Divide_Phases.m)**

Adds a field to data structure in QC file containing a matrix of indices with the tidal phase of each ensembles (of size [1 x NEns]. Phases are currently set as: (0) Decelerating Ebb (1) Accelerating Flood (2) Peak Flood (3) Decelerating Flood (4) Accelerating Ebb (5) Peak Ebb  

**Fit Vertical Profiles using Logarithmic Law  (LogLaw_Fits.m)**

Fits vertical profiles to velocity bin averaged along stream velocities using the Logarithmic Law of the Wall. Currently set up to fit one tide phase and wind condition (ex accelerating flood, windy ensembles) at a time. 

(1) Performs a linear fit at each depth bin, such that the first iteration included the bottom three depth bins and each subsequent iteration included one additional depth bin, until all depth bins were used for the fit. Calculates parameters 𝑢∗ and 𝑧𝑜 using the best (maximized with respect to R2) fit for each averaged profile. 

(2) Plots the fit in logarithmic scale along with the observed profiles 

**Fit Vertical Profiles using Wake Law (Wake_Fits.m)**

Fits vertical profiles to velocity bin averaged along stream velocities using an adapted Logarithmic Law of the Wall, which includes a cubic wake term to capture near surface behaviour. Currently set up to fit one tide phase and wind condition (ex accelerating flood, windy ensembles) at a time. Also calculates fit error (RMSE, normalized by flow speed). Draws figures showing observational and fit data. 

Dependency: wake_function.m - contains the wake function which is fit

**Calculate Turbulence Metrics (Turbulence_Metrics.m)**

Calculates and creates colourplots of the following turbulence metrics: Turbulent kinetic energy (TKE) (using variance method), horizontal variance, turbulence intensity (TI), vertical (beam 5) variance. Generates plots that show all data, and also plots that remove data above the sigma = -0.15 threshold.

**Perform Spectral Analysis (Spectral_Analysis.m)**

Computes spectra of vertical velocities in frequency domain for all ensembles individually at one depth level, option to then average all ensembles or a subset of ensembles by tide stage. 

Dependency: get_ensembles.m - Grabs concatenated file of vertical velocities and sorts it into 5-minute ensembles (retains individual measurements) 

            get_specs.m - Uses pwelch.m to calculate spectra of vertical velocities 

