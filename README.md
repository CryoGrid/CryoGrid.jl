# CryoGridJulia
Custom Julia implementation of the CryoGrid permafrost land surface model.

Author: Moritz Langer (moritz.langer@awi.de)

Alfred Wegener Institute, Helmholtz Center for Polar and Marine Research (AWI)
Telegrafenberg A45 | 14473 Potsdam

requires Julia Version 0.6.4 which can be downloaded at: https://julialang.org/downloads/oldreleases/
Additional packages are required which can be installed with Pkg.add("GLOB")
succefully tested on a Linux Server running Ubuntu 18.04

Required climate forcing data (JSON file format) are located in the folder: FORCING_JSONfiles
Associated model parameter (JSON file format) are located in the folder: PARA_JSONfiles

For ruining the model for different locations and climate forcing change the following lines in
CryoGrid_shell_landonly.jl accordingly including brackets []:

#### Set site coordinates as upper left corner (grid: 0-360E, 55-90N)
ULC_lat = [insert latitude as integer] 
ULC_lon = [insert longitude as integer]

#### load forcing data -------------------------------------------------------------
forcing_path = "FORCING_JSONfiles/[change file prefix accordingly]_ULC_" * string(ULC_lon) * "_" * string(ULC_lat) * ".json"

Either run the shell script: CryoGrid_shell_landonly.jl directly from within Julia or within a terminal

