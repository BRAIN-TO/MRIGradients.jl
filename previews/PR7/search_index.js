var documenterSearchIndex = {"docs":
[{"location":"Utilities/","page":"Utilities","title":"Utilities","text":"Modules = [MRIGradients]\nPages   = [\"Util.jl\"]","category":"page"},{"location":"Utilities/#MRIGradients.time2freq-Tuple{Any}","page":"Utilities","title":"MRIGradients.time2freq","text":"time2freq(t)\n\nconverts time domain vector into a vector of frequencies according to the Nyquist criterion Port of Johanna Vannesjo's time2freq method: https://github.com/MRI-gradient/GIRF/blob/main/utils/time2freq.m\n\nArguments\n\nt - input time vector\n\n\n\n\n\n","category":"method"},{"location":"GIRFApplier/","page":"GIRFApplier","title":"GIRFApplier","text":"Modules = [MRIGradients]\nPages   = [\"GIRFApplier.jl\"]","category":"page"},{"location":"GIRFApplier/#MRIGradients.apply_girf-Tuple{GirfApplier, Any, Any, Any, Any}","page":"GIRFApplier","title":"MRIGradients.apply_girf","text":"apply_girf(g::GirfApplier, gradient_in, time_in, time_out, direction)\n\nFunction for predicting the gradient by applying the GIRF\nThis is done in one dimension only. Repeat this function for multiple dimensions\n\nArguments\n\ng::GirfApplier - GirfApplier containing GIRF data\ngradient_in - input gradient waveform corresponding to time_in vector\ntime_in - time vector corresponding to input gradient waveform\ntime_out - output time vector desired (decimation)\ndirection - direction of the GIRF to use (x,y,z etc...)\n\n\n\n\n\n","category":"method"},{"location":"#Welcome-to-the-MRIGradients.jl-Documentation!","page":"Home","title":"Welcome to the MRIGradients.jl Documentation!","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Follow the side pages for information about the provided functionality of the package :) ","category":"page"},{"location":"GIRFEssential/","page":"GIRFEssential","title":"GIRFEssential","text":"Modules = [MRIGradients]\nPages   = [\"GIRFEssential.jl\"]","category":"page"},{"location":"GIRFEssential/#MRIGradients.displayGirf-Tuple{GirfEssential}","page":"GIRFEssential","title":"MRIGradients.displayGirf","text":"displayGirf(g::GirfEssential)\n\nFunction to print out the girf data in a human-readable format as a summary to make sure data is loaded properly\n\nArguments\n\ng::GirfEssential - GirfEssential structure containing GIRF data\n\n\n\n\n\n","category":"method"},{"location":"GIRFEssential/#MRIGradients.loadGIRFMatlabTim","page":"GIRFEssential","title":"MRIGradients.loadGIRFMatlabTim","text":"GIRFfreq, GIRFdata, GIRF_length = loadGIRFTimMatlab() loads different measured GIRFs from Siemens Prisma, as computed in Matlab (Zhe \"Tim\" Wu's code)\n\nArguments\n\nidGirf              - identifier for selected GIRF                        0 = example GIRF in this repo                        1 = Nov 2020 pos blips                        2 = June 2021 pos blips                        3 = June 2021 pos/neg blips (4 averages)                         4 = November 2021, pos/neg blips (4 averages)\n\nOutputs\n\nGIRF_freq           - \nGIRF_data           - x 3 data\nGIRF_length           -                        \n\n\n\n\n\n","category":"function"},{"location":"GIRFEssential/#MRIGradients.loadGirf","page":"GIRFEssential","title":"MRIGradients.loadGirf","text":"loadGirf(degree, id)\n\nFunction to load girf with specified id according to hardcoded file structure (internal use only)\n\nArguments\n\ndegree - specifies degree of the girf (1st or 0th order correction)\nid - specifies girf measurement id since there have been multiple measurements\n\n\n\n\n\n","category":"function"},{"location":"GIRFEssential/#MRIGradients.readGIRFFile-Tuple{String, String, String, String, Bool}","page":"GIRFEssential","title":"MRIGradients.readGIRFFile","text":"readGIRFFile(g::GirfEssential, pathX::String, pathY::String, pathZ::String, varName::String, doFilter::Bool)\n\nFunction to read GIRF files with specific variable name and return a corresponding GIRFEssential object.\n\nArguments\n\npathX::String - Full path for GIRF MAT-file on X axis\npathY::String - Full path for GIRF MAT-file on Y axis\npathZ::String - Full path for GIRF MAT-file on Z axis\nvarName::String - Variable name that needs to be read as GIRF data from the MAT-files.\ndoFilter::Bool - Whether we filter the GIRF data using Tukey filter\n\n\n\n\n\n","category":"method"}]
}
