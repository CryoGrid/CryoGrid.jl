module Implicit
#Start script for CryoGrid (adapted to standalone module)
#Created on Fri Jul 17 19:06:49 2020
#author: Moritz Langer (moritz.langer@awi.de) Brian Groenke (brian.groenke@awi.de)
#runs with Julia Version 0.6.2 (2017-12-13 18:08 UTC)
#https://julialang.org/downloads/oldreleases/
#note that the script might require to install additional packages e.g. Pkg.add("Interpolations")
# ========================================================================
#ArcticLakes shell script to run the CryoGridImplicit Model
using MAT
using Glob
using JLD
using JSON
using Dates
using Statistics
using ProgressMeter
include("json2type.jl")
include("DataLoader.jl")
include("matlab.jl")
include("ArcticLakesSetup.jl")
include("CryoGridTyps.jl")

export load_forcing, load_para, run_landonly, save
export HeatCapacity, ThermalConductivity, Enthalpy, EnthalpyInv, Model, MakeGrid

function load_forcing(ULC_lat=72, ULC_lon=126)
    # ==============================================================================
    #Set site coordinates as upper left corner (grid: 0-360E, 55-90N)
    # ==============================================================================
    #load forcing data -------------------------------------------------------------
    forcing_path = "input/FORCING_JSONfiles/FORCING_ULC_" * string(ULC_lon) * "_" * string(ULC_lat) * ".json"
    jsontxt = open(forcing_path,"r") do file
            dicttxt = read(file, String)
    end
    FORCING = json2type.typenarrow!(JSON.parse(jsontxt))
    #in the case nan values exisit in Ptot
    for i = 1:length(FORCING["data"]["Ptot"])
            if FORCING["data"]["Ptot"][i]==nothing
                    FORCING["data"]["Ptot"][i] = 0.0
            end
    end
    FORCING["data"]["Ptot"] = convert(Array{Float64,1},FORCING["data"]["Ptot"])
    #snow fall in winter (snow fall is in m/day - snow water equivalent - SWE!!!)
    FORCING["data"]["snowfall"] = (FORCING["data"]["Tair"].<=0.0).*FORCING["data"]["Ptot"]/1000.0;
    FORCING
end

function load_para(ULC_lat=72,ULC_lon=126)
    #load parameter ----------------------------------------------------------------
    para_path = "input/PARA_JSONfiles/PARA_ULC_" * string(ULC_lon) * "_" * string(ULC_lat) * ".json"
    jsontxt = open(para_path,"r") do file
        dicttxt = read(file, String)
    end
    PARA = json2type.typenarrow!(JSON.parse(jsontxt))
end

function run_landonly(FORCING,PARA,t0)
    GRID = PARA["GRID"] #soil grid with depth
    PARA = PARA["PARA"] #soil and snow parameters

    #Set to land only (no lakes)
    tile_number = 1 #use only one tile for reference run
    LAKESTAT = Dict()
    LAKESTAT["LandMaskArea"]=[1.0]
    LAKESTAT["LakeDistance"]=[1.0]
    LAKESTAT["LakeNumber"]=[0.0]
    LAKESTAT["LakeDepth"]=[0.0]
    LAKESTAT["SnowDepth"]=[NaN]
    LAKESTAT["VorArea"]=[1.0]
    LAKESTAT["LandMaskArea"]=[1.0]
    LAKESTAT["LakeArea"]=[0.0]
    LAKESTAT["LakePerimeter"]=[NaN]

    #BUILD input structures --------------------------------------------------------
    FOR, PAR, TEM, GRI, STRA, STAT, OUT = ArcticLakesSetup.SetUpInputStructs(FORCING, GRID, PARA, LAKESTAT, tile_number)

    #RUN model ---------------------------------------------------------------------
    @time Model(FOR, GRI, STRA, PAR, STAT, TEM, OUT; start=t0)
end

function save(OUT)
    #REDUCE output precision -------------------------------------------------------
    SAVE = CryoGridTyps.save(OUT.Date, OUT.H_av, OUT.T_av, OUT.T_min, OUT.T_max, OUT.W_av, OUT.W_min, OUT.W_max, OUT.Q_lat, OUT.FDD, OUT.TDD, OUT.FrostDays, OUT.SnowDepth_av, OUT.SnowDepth_max, OUT.SnowDays)
    #SAVE result as JSON for each year ---------------------------------------------
    subfolder = "RESULT_JSONfiles"
    mkpath(subfolder);
    subsubfolder = "ULC_" * string(ULC_lon) * "_" * string(ULC_lat)
    mkpath(subfolder * "/" * subsubfolder);
    for idx = 1:length(SAVE.Date)
            out_dict = Dict()
            year = SAVE.Date[idx]
            FieldsInStruct = fieldnames(typeof(SAVE))
            for i = 1:length(FieldsInStruct)
                    #Check fields
                    Value = getfield(SAVE, FieldsInStruct[i])
                    if size(Value)[2]>=idx
                            value = Value[:,idx,:]
                    else
                            value = Value[:,end,:]
                    end
                    name = string(FieldsInStruct[i])
                    out_dict[name] = value
            end
            #save as json uncompressed
            saveloc = subfolder * "/" * subsubfolder * "/" * "RESULT_" * string(year) * ".zjson"
            open(saveloc,"w") do file
                    write(file,JSON.json(out_dict))
            end
    end
end

#implcit heat transfere with phase change following an approch by:
#C.R. SWAMINATHAN and V.R. VOLLER (1992) in METALLURGICAL TRANSACTION
#this version uses a free water freeze curve
#single side lateral heat flux is included
# ========================================================================
function HeatCapacity(WaterIce, Water, Mineral, Organic, HC)
    cs = 2000; #[J/kgK]  heat capacity solid
    cl = 2500; #[J/kgK]  heat capacity liquid
    cm = 1000; #[J/kgK]  heat capacity matrix
    cp = 1.0;  #[J/kgK]  heat capacity pore space

    cw = 4.2*10^6; #[J/m^3K] heat capacity water
    co = 2.5*10^6; #[J/m^3K]  heat capacity organic
    cm = 2*10^6; #[J/m^3K]  heat capacity mineral
    ca = 0.00125*10^6;#[J/m^3K]  heat capacity pore space
    ci = 1.9*10^6;#[J/m^3K]  heat capacity ice

    @inbounds @fastmath for i=1:length(HC)
        air = 1.0 - WaterIce[i] - Mineral[i] - Organic[i];
        ice = WaterIce[i] - Water[i];
        HC[i] = Water[i].*cw + ice.*ci + Mineral[i].*cm + Organic[i].*co + air.*ca;
    end

    return HC
end
#--------------------------------------------------------------------------
function ThermalConductivity(WaterIce, Water, Mineral, Organic, TC)
    ka = 0.025;       #air [Hillel(1982)]
    kw = 0.57;        #water [Hillel(1982)]
    ko = 0.25;        #organic [Hillel(1982)]
    km = 3.8;         #mineral [Hillel(1982)]
    ki = 2.2;         #ice [Hillel(1982)]

    @inbounds @fastmath for i=1:length(TC)
        ice = WaterIce[i] - Water[i];
        air = 1.0 - WaterIce[i] - Mineral[i] - Organic[i];
        TC[i] = (Water[i].* kw.^0.5 + ice.* ki.^0.5 + Mineral[i].* km.^0.5 + Organic[i].* ko.^0.5 + air.* ka.^0.5).^2.0;
    end

    return TC
end
#--------------------------------------------------------------------------
function Enthalpy(T, Water, HC, H, dHdT)
    rho = 1000.; #[kg/m^3]
    Lsl = 334000.; #[J/kg]
    L = rho*Lsl; #[J/m^3]
    theta = Water; #[Vol. fraction]

    @inbounds @fastmath for i = 1:length(H)
        H[i] = T[i]*HC[i] + theta[i]*L
        dHdT[i] = HC[i]
        if T[i]==0.0
            dHdT[i] = dHdT[i] + 1.0e9
        end
    end

    return H, dHdT
end
#--------------------------------------------------------------------------
function EnthalpyInv(H, WaterIce, HC, T, fl)

    rho = 1000.; #[kg/m^3]
    Lsl = 334000.; #[J/kg]
    L = rho*Lsl;#[J/m^3]

    theta = WaterIce;
    @inbounds @fastmath for i = 1:length(T)
        if theta[i]<=0.0
            theta[i] = 1e-8
        end
        if H[i]>0.0 && H[i]<=L*theta[i]
            fl[i] = H[i]./(L*theta[i])
            T[i] = 0.0
        elseif H[i]>L*theta[i]
            fl[i] = 1.0
            T[i] = (H[i]-L*theta[i])./HC[i]
        elseif H[i]<0.0
            fl[i] = 0.0
            T[i] = H[i]./HC[i]
        end
    end

    return T, fl
end
#--------------------------------------------------------------------------
function TDMAsolver(a, b, c, d, x)

    #a, b, c are the column vectors for the compressed tridiagonal matrix, d is the right vector
    n = length(b); # n is the number of rows
    #x = zeros(n,1);
    # Modify the first-row coefficients
    c[1] = c[1] / b[1];    # Division by zero risk.
    d[1] = d[1] / b[1];    # Division by zero would imply a singular matrix.

    @inbounds @fastmath for i = 2:n-1
        temp = b[i] - a[i] * c[i-1];
        c[i] = c[i] / temp;
        d[i] = (d[i] - a[i] * d[i-1]) / temp;
    end

    d[n] = (d[n] - a[n] * d[n-1])/( b[n] - a[n] * c[n-1]);

    # Now back substitute.
    x[n] = d[n];
    @inbounds @fastmath for i = n-1:-1:1
        x[i] = d[i] - c[i] * x[i + 1];
    end

    return x
end
#--------------------------------------------------------------------------
function MakeGrid(Zp, dxp, dxn, dxs, kp, kn, ks)
    N = length(Zp);

    Zn = Zp-dxp/2.;
    Zs = Zp+dxp/2.;

    #dxn=ones(N,1);
    #dxs=ones(N,1);

    dxs[1,1] = Zp[2] - Zp[1];
    dxs[N,1] = dxp[N]/2.;

    dxn[1,1] = dxp[1,1]/2.;
    dxn[N,1] = Zp[N] - Zp[N-1];

    #kn=ones(N,1);
    #ks=ones(N,1);

    kn[1,1] = kp[1];
    kn[N,1] = kp[N];
    ks[1,1] = kp[1];
    ks[N,1] = kp[N];

    @inbounds @fastmath for i=2:N-1
        dxn[i,1] = Zp[i] - Zp[i-1];
        dxs[i,1] = Zp[i+1] - Zp[i];

        kn[i,1] = (dxp[i,1]/(2*dxn[i])*kp[i,1].^-1 + dxp[i-1,1]/(2*dxn[i])*kp[i-1].^-1).^-1
        ks[i,1] = (dxp[i,1]/(2*dxs[i])*kp[i,1].^-1 + dxp[i+1,1]/(2*dxs[i])*kp[i+1].^-1).^-1
    end

    return Zn, Zs, dxn, dxs, kn, ks
end

#--------------------------------------------------------------------------
function SnowCover2(T0, Water, WaterIce, Zs, snow_depth, snow_max, snow_fall, SnowDensity, SnowDensityMax, SnowDensityK1, SnowDensityK2, WaterDensity, current_date)

    #get current frozen snow cover dept
    snow_dw = sum(Zs.<=0.0); #lower most snow cell
    snow_up = max(snow_dw-sum((WaterIce.>0.0) .& (Zs.<=0.0)),1); #cell above last snow cell
    old_snow_up = snow_up;

    delta_snow_depth = snow_depth - max(-Zs[snow_up],0.0) #the difference between real and grided snow depth

    #snow depletion
    if snow_up < snow_dw #snow cover exists
        #route melt water downward
        #water_holding_cap = 0.1;
        meltwater = 0.0; #initlize melt trem
        snow_rm = 0.65; #together with snow grid cell size this value scales the positive degree day factor of snow melt (see e.g. VICTOR C. TSAI & XIAOZHOU RUAN 2018)
        water_hold = 0.025;
        for i=snow_up+1:snow_dw-1
            water_holding_cap = water_hold*(WaterIce[i]-Water[i]); #2.5 percent of the water equivalent of the snowpack existing in the ice phase (U.S. Army, 1956; Leaf, 1966)
            #route melt water downward
            if Water[i] > water_holding_cap
                meltwater = max(0.0, Water[i]-water_holding_cap);
                WaterIce[i] = WaterIce[i] - meltwater;
                Water[i] = Water[i] - meltwater;
                WaterIce[i+1] = WaterIce[i+1] + meltwater;
                Water[i+1] = Water[i+1] + meltwater;

                #remove snow cell if it contains less ice than 50% of inituial SWE
                if (WaterIce[i]-Water[i])<=(snow_rm*SnowDensity[1]/WaterDensity[1])
                    #treat rest as additional meltwater
                    meltwater = WaterIce[i];
                    WaterIce[i] = 0.0;
                    Water[i] = 0.0;
                    WaterIce[i+1] = WaterIce[i+1] + meltwater;
                    Water[i+1] = Water[i+1] + meltwater;

                    #update upper snow indx and snow depth
                    snow_up = snow_up+1;
                    snow_depth = max(-Zs[snow_up],0.0) + delta_snow_depth;
                end
            end
        end

        #last snow cell
        i = snow_dw;
        #remove last snow cell and/or determine routetable water
        if (WaterIce[i]-Water[i])<=(snow_rm*SnowDensity[1]/WaterDensity[1])
            meltwater = WaterIce[i];
            WaterIce[i] = 0.0;
            Water[i] = 0.0;
            #update upper snow indx and snow depth
            snow_up = snow_up+1;
            snow_depth = max(-Zs[snow_up],0.0) + delta_snow_depth;
            if snow_up == snow_dw
                #if the snow cover is completely melted set snow_depth to zero
                snow_depth = 0.0
            end
        else
            meltwater = Water[i];
            #consider this water as free routetable water for upward routing
            Water[i] = Water[i]-meltwater;
            WaterIce[i] = WaterIce[i]-meltwater;
        end

        #route metlwater upward
        if (meltwater>0.0) & (snow_up<snow_dw)
            #copy water and ice arrays
            WaterIce_cp = WaterIce;
            Water_cp = Water;
            j = snow_dw;
            for k = snow_dw:-1:snow_up

                #step over previously removed snow cells
                while (WaterIce_cp[j]==0.0) & (j>=old_snow_up+1)
                    j = j-1;
                end

                #melt water routeting
                if j>0 #only for safety reason
                    WaterIce[k] = WaterIce_cp[j] + meltwater;
                    Water[k] = Water_cp[j] + meltwater;
                    meltwater = 0.0;

                    #further route rest if cell is filled until an uptake
                    #maximum
                    UTmax = 0.9; #1.0 would be completely filled
                    if WaterIce[k]>UTmax
                        meltwater = WaterIce[k]-UTmax;
                        WaterIce[k] = UTmax;
                        Water[k] = Water[k]-meltwater;
                    end
                else #only for safety reason
                    WaterIce[k] = 0.0;
                    Water[k] = 0.0;
                end
                j = j-1;
            end
        end
    end
    #remove all water and ice rests above actual snow cover
    WaterIce[1:snow_up] .= 0.0;
    Water[1:snow_up] .= 0.0;

    #snow compactionmodule according to Sturm et al.(2010)
    if snow_up < snow_dw #snow cover exists
        #current_date = DateTime(2019,1,30)
        current_year = year(current_date)
        if current_date < DateTime(current_year,10,1)
            delta_doy = 0
        else
            delta_doy = dayofyear(DateTime(current_year,10,1))+92
        end
        ##=
        DOY = dayofyear(current_date)
        DOY_snowseason = DOY-delta_doy
        DOY_snowseason = max(-92,min(181,DOY_snowseason)) #limints according to Sturm et al. 2010
        k1 = SnowDensityK1[1]
        k2 = SnowDensityK2[1]
        SnowDensity_0 = SnowDensity[1]
        SnowDensity_max = SnowDensityMax[1]
        BulkSnowDensity_new = (SnowDensity_max - SnowDensity_0).*(1.0 - exp(-k1*snow_depth - k2*DOY_snowseason)) + SnowDensity_0
        #current snow bulk density according to gridded ice content (only takes Ice into account)
        BulkSnowDensity_current = max(SnowDensity_0, WaterDensity[1]*mean(WaterIce[snow_up+1:snow_dw]))

        if (current_date >= DateTime(current_year,8,1)) .& (current_date < DateTime(current_year,9,1))
            snow_depth = 0.0 #set snow depth to zero in the period from August 1 to August 31 to prevent glaciers
        else
            snow_depth = max(min(snow_depth, snow_max), 0.0); #avoid values smaller 0.0 and larger max snow depth
        end
        #new snow depth according to new snow bulk snow density this emulates compaction only
        snow_depth_compact = snow_depth*BulkSnowDensity_current/BulkSnowDensity_new
        snow_depth_compact = max(min(snow_depth_compact, snow_max), 0.0); #avoid values smaller 0.0 and larger max snow depth
        new_snow_up = argmin(abs.(snow_depth_compact.-(-Zs)));
        #remove snow cover if snow depth is lower than gridded

        if new_snow_up>snow_up
            #println("compaction from: " * string(snow_depth) * " to: " * string(snow_depth_compact))
            #tot_WaterIce_current = sum(WaterIce[snow_up+1:snow_dw])
            #tot_WaterIce_new = sum(WaterIce[new_snow_up+1:snow_dw])
            #delta_WaterIce = tot_WaterIce_current-tot_WaterIce_new;
            #println("delta total Ice Water is: " * string(delta_WaterIce))

            WaterIce[snow_up:new_snow_up] .= 0.0;
            Water[snow_up:new_snow_up] .= 0.0;
            snow_up = new_snow_up;
            snow_depth = max(-Zs[snow_up],0.0)
            #adjust ice contend to new Bulk Snow Density
            #WaterIceSnow = WaterIce[snow_up+1:snow_dw]
            #low_density_idx = (WaterIceSnow.*WaterDensity[1] .< BulkSnowDensity_new)
            #BulkSnowDensity_low_density = mean(WaterIceSnow[low_density_idx])*WaterDensity[1]
            #WaterIceSnow_new = WaterIceSnow[low_density_idx]*(BulkSnowDensity_new/BulkSnowDensity_low_density)
            #WaterIceSnow[low_density_idx] = WaterIceSnow_new
            #WaterIce[snow_up+1:snow_dw] = WaterIceSnow

            WaterIce_current = WaterIce[snow_up+1:snow_dw]
            WaterIce_new = WaterIce[snow_up+1:snow_dw]*(BulkSnowDensity_new/BulkSnowDensity_current)
            idx_smaller_new_BD = WaterIce_current .< BulkSnowDensity_new/WaterDensity[1]
            #WaterIce_current[idx_smaller_new_BD] = WaterIce_new[idx_smaller_new_BD]
            WaterIce_current = WaterIce_new
            WaterIce_current[WaterIce_current.>0.9] .= 0.9
            WaterIce[snow_up+1:snow_dw] = WaterIce_current

        end
    end

    #snow accumulation
    snow_depth = snow_depth + snow_fall*WaterDensity[1]/SnowDensity[1]; #increses snow depth by snow fall
    snow_depth = max(min(snow_depth, snow_max), 0.0); #avoid values smaller 0.0 and larger max snow depth
    new_snow_up = argmin(abs.(snow_depth.-(-Zs)));
    #build a new snow cell on top if snow depth fills a complet new cell
    if new_snow_up<snow_up
        WaterIce[new_snow_up+1:snow_up] .= SnowDensity[1]/WaterDensity[1];
        Water[new_snow_up+1:snow_up] .= 0.0;
        snow_up = new_snow_up;
    end

    #get new air-ground interface index
    idx = snow_up+1;

    #snow cover flag indicating if a fully developed snow cover exists
    if snow_up < snow_dw
        snow_flag = true
    else
        snow_flag = false
    end

    return Water, WaterIce, snow_depth, idx, snow_flag
end
# ========================================================================
function Model(FORCING, GRID, STRAT, PARA, STATVAR, TEMP, OUT; start::DateTime=nothing)
    dt = 60.0*60.0*24.0;
    t_span = FORCING.t_span;
    Tair = FORCING.Tair;
    snowfall = FORCING.snowfall;

    dxp = GRID.dxp;
    Zp = GRID.Zp;
    Zn = GRID.Zn;
    Zs = GRID.Zs;
    dxn = GRID.dxn;
    dxs = GRID.dxs;
    dxo = GRID.dxo;
    An = GRID.An;
    As = GRID.As;
    Ao = GRID.Ao;
    Vp = GRID.Vp;

    cp = TEMP.cp;
    kp = TEMP.kp;
    kn = TEMP.kn;
    ks = TEMP.ks;
    lat_flux = TEMP.lat_flux;
    SnowDepth = TEMP.SnowDepth;

    WaterIce = STRAT.WaterIce;
    Mineral = STRAT.Mineral;
    Organic = STRAT.Organic;
    Water = STRAT.Water;

    WaterDepth = PARA.WaterDepth;
    SnowDensity = PARA.SnowDensityZero;
    SnowDensityMax = PARA.SnowDensityMax;
    SnowDensityK1 = PARA.SnowDensityK1;
    SnowDensityK2 = PARA.SnowDensityK2;
    WaterDensity = PARA.WaterDensity;
    Qgeo = PARA.Qgeo;
    SnowDepthMax = PARA.SnowCoverMax;

    T = STATVAR.T;
    H = STATVAR.H;

    #preallocate temp variables
    dHdT = zeros(size(T,1))
    fl = zeros(size(T,1))

    N = size(T,1);
    N_tiles = size(T,2);
    bp_lat = zeros(N,size(T,2));
    bp = zeros(N,1);

    #preallocate output arrays
    out_years = unique(year.(matlab.datestr(t_span)))
    out_length = length(out_years)
    out_intervall = Year(1)
    out_date = DateTime(DateTime(out_years[1])+out_intervall)
    out_N = 1.0
    out_idx = 1

    #local output variables
    H_av = OUT.H_av
    T_av = OUT.T_av
    T_min = OUT.T_min
    T_max = OUT.T_max
    T_out = OUT.T

    W_av = OUT.W_av
    W_min = OUT.W_min
    W_max = OUT.W_max
    Water_out = OUT.Water
    WaterIce_out = OUT.WaterIce

    Qlat_out = OUT.Q_lat

    TDD = OUT.TDD;
    FDD = OUT.FDD;
    FrostDays = OUT.FrostDays;

    SnowDepth_av = OUT.SnowDepth_av;
    SnowDepth_max = OUT.SnowDepth_max;
    SnowDays = OUT.SnowDays;

    @showprogress for t = 1:length(t_span)
        current_date = matlab.datestr([t_span[t]])[1]
        if start != nothing && current_date < start
            continue
        end
        #upper boundary temperature [°C] and snow fall rate [m/day]
        T0 = Tair[t];
        snow_fall = snowfall[t];

        dp = 1.0/dt;
        T_old = copy(T);
        H_old = copy(H);

        for j = N_tiles:-1:1 #loop backward through tiles, tile1 is matrix
            H_old[:,j] = copy(H[:,j]);
            T_old[:,j] = copy(T[:,j]);

            #Prognostic Snow Scheme
            Water[:,j], WaterIce[:,j], SnowDepth[1,j], idx, snow_flag = SnowCover2(T0, Water[:,j], WaterIce[:,j], Zs, SnowDepth[1,j], SnowDepthMax[1,j], snow_fall, SnowDensity, SnowDensityMax, SnowDensityK1, SnowDensityK2, WaterDensity[1], current_date);
            #update properties and state variables
            idx_snowmelt = (T[:,j].<0.0) .& (Water[:,j] .>0.0) .& (Zs.<=0.0)
            T[idx_snowmelt,j] .= 0.0;
            cp[:,j] = HeatCapacity(WaterIce[:,j], Water[:,j], Mineral[:,j], Organic[:,j], cp[:,j]);
            kp[:,j] = ThermalConductivity(WaterIce[:,j], Water[:,j], Mineral[:,j], Organic[:,j], kp[:,j]);
            Zn, Zs, dxn, dxs, kn, ks = MakeGrid(Zp, dxp, dxn, dxs, kp[:,j], kn, ks);
            Hⱼ, dHdT = Enthalpy(T[:,j], Water[:,j], cp[:,j], H[:,j], dHdT);
            H[:,j] .= Hⱼ
            Tⱼ, fl = EnthalpyInv(H[:,j], WaterIce[:,j], cp[:,j], T[:,j], fl);
            T[:,j] .= Tⱼ

            #simple snow cover module
            #Water[:,j], SnowDepth[1,j], idx = SnowCover(Water[:,j], WaterIce[:,j], Zp, T0, SnowDepth[1,j], SnowDepthMax[1,j],snow_fall);

            #simple lake scheme
            if T0>=0.0 #Thawing, in the case of freezing upper boundary will be at lake top
                #get the last cell that consists of 100% liquid water
                SolidMatrix = Mineral[:,j] + Organic[:,j]
                lic_lake = Water[idx,j]/(WaterIce[idx,j]+SolidMatrix[idx]);
                while lic_lake == 1.0
                    lic_lake = Water[idx,j]/(WaterIce[idx,j]+SolidMatrix[idx]);
                    if lic_lake == 1.0
                        idx = idx+1;
                    end
                end
            end
            ubc_idx = idx;

            temp_offset_max = 1.;
            iter_count = 1;
            while ((temp_offset_max>=1e-3) || (iter_count<10)) && (iter_count<=500)
                if iter_count>499
                    println("warning: did not reach convergence!")
                    println(temp_offset_max)
                end
                iter_count=iter_count+1;

                #implicit update thermal conductivity and heat capacity
                kp[:,j] = ThermalConductivity(WaterIce[:,j], Water[:,j], Mineral[:,j], Organic[:,j], kp[:,j])
                cp[:,j] = HeatCapacity(WaterIce[:,j], Water[:,j], Mineral[:,j], Organic[:,j], cp[:,j])

                #update grid
                Zn, Zs, dxn, dxs, kn, ks = MakeGrid(Zp, dxp, dxn, dxs, kp[:,j], kn, ks);
                ko = kp[:,1];#ks; #for testing lat flux

                #pre-factores acording to grid cell sizes and thermal conductivities
                kn[ubc_idx] = kp[ubc_idx];#ensures full thermal conductivity at upper boundary
                anpn = An[:,j]./Vp[:,j].*kn./dxn;
                anps = As[:,j]./Vp[:,j].*ks./dxs;

                #Additional heat fluxes----------------------------------------
                bp = bp.*0.0;
                #flux from upper boundary
                bp[1] = (An[1,j]./Vp[1,j].*kn[1]./dxn[1]).*T0; #[W/m³]
                #flux from lower boundary
                bp[end] = (An[1,j]./Vp[end,j]) * Qgeo[1]; #[W/m³]

                #if upper boundary is somewhere within the grid
                if ubc_idx > 1
                    bp[1:ubc_idx] = (An[:,j]./Vp[1:ubc_idx,j].*kn[1:ubc_idx]./dxn[1:ubc_idx]) * T0;
                    anps[1:ubc_idx] .= 0.0;
                end

                if j!=1
                    ap = anpn + anps + Ao[:,j]./Vp[:,j].*ko[:]/dxo[j];
                    ap[end] = anpn[end] + Ao[end,j]./Vp[end,j].*ko[end]/dxo[j];
                    anpn[1:ubc_idx] = 0.0;
                    #lateral heat flux (relative [W/m³])
                    bp_lat[:,j] = (Ao[:,j]./Vp[:,j].*ko[:]./dxo[j]) .* T_old[:,1];
                    #update total lateral heat fluxes [W]
                    lat_flux[:,j] = (Ao[:,j].*ko[:]./dxo[j]).* (T[:,j]-T_old[:,1]);
                    #println(sum(lat_flux[:,j]))
                    #println(sum(bp_lat[:,j]))
                else
                    #for tile 1 the lateral flux is assumed to be a static
                    #external flux reulting from the sum of the lateral heat
                    #fluxes from the other tiles.
                    ap = anpn + anps;
                    ap[end] = anpn[end];
                    anpn[1:ubc_idx] .= 0.0;
                    bp_lat[:,j] = sum(lat_flux,dims=2)./Vp[:,j]; #[W/m³];
                    #println(sum(bp_lat[:,j]))
                    #println(' ')
                end
                bp = bp + bp_lat[:,j];

                #--------------------------------------------------------------
                Hinv, fl = EnthalpyInv(H[:,j], WaterIce[:,j], cp[:,j], T[:,j], fl);
                #Hinv = temperature from inverse T-H function
                #Water = total vol fraction of liquid water
                dummy, dHdT = Enthalpy(Hinv, Water[:,j], cp[:,j], H[:,j], dHdT);
                #dHdT = approximated deriviative of H at T
                #fl = fraction of liquid water
                Water[:,j] = fl.*WaterIce[:,j];
                Sp = -dp*dHdT;
                Sc = dp*(H_old[:,j] - H[:,j]) - Sp.*Hinv;

                #TDMA solver --------------------------------------------------
                alpha = -anpn;
                beta = (ap - Sp);
                gamma = -anps;
                delta = Sc + bp;
                T[:,j] = TDMAsolver(alpha, beta, gamma, delta, T[:,j]);

                #---------------------------------------------------------------
                #update current state of H
                H[:,j] = H[:,j] + dHdT.*(T[:,j] - Hinv);

                Hinv_check, dummy = EnthalpyInv(H[:,j], WaterIce[:,j], cp[:,j], T[:,j], fl);
                R = Hinv_check;
                B = T[:,j];
                temp_offset=R-B;
                temp_offset_max = maximum(temp_offset[2:end-1,:]);

            end

            #-------------------------------------------------------------------
            #Model OUTPUT
            if (current_date == out_date) .& (j == N_tiles)
                #update to next output date
                out_date = DateTime(current_date+out_intervall)
                #increase save index
                out_idx += 1
                #reset counter
                out_N = 1.0
                println(out_date)
            end
            if out_N == 1.0 #first values of the averaging period
                H_av[:,out_idx,j] = copy(H[:,j])
                T_av[:,out_idx,j] = copy(T[:,j])
                T_min[:,out_idx,j] = copy(T[:,j])
                T_max[:,out_idx,j] = copy(T[:,j])

                W_av[:,out_idx,j] = copy(Water[:,j])
                W_min[:,out_idx,j] = copy(Water[:,j])
                W_max[:,out_idx,j] = copy(Water[:,j])

                Qlat_out[:,out_idx,j] = -copy(lat_flux[:,j]);

                TDD[:,out_idx,j] = T[:,j].*(T[:,j].>=0.0);
                FDD[:,out_idx,j] = T[:,j].*(T[:,j].<0.0);
                FrostDays[:,out_idx,j] = 1.0*(T[:,j].<0.0);

                SnowDays[1,out_idx,j] = 1 * (SnowDepth[1,j]>0.0);
                SnowDepth_av[1,out_idx,j] = copy(SnowDepth[1,j]);
                SnowDepth_max[1,out_idx,j] = copy(SnowDepth[1,j]);


            else #cumulative averaging and min max
                H_av[:,out_idx,j] = H_av[:,out_idx,j] + (H[:,j]-H_av[:,out_idx,j])/out_N
                T_av[:,out_idx,j] = T_av[:,out_idx,j] + (T[:,j]-T_av[:,out_idx,j])/out_N
                T_min[:,out_idx,j], T_max[:,out_idx,j] = matlab.fastminmax(T_min[:,out_idx,j], T_max[:,out_idx,j], T[:,j])

                W_av[:,out_idx,j] = W_av[:,out_idx,j] + (Water[:,j]-W_av[:,out_idx,j])/out_N
                W_min[:,out_idx,j], W_max[:,out_idx,j] = matlab.fastminmax(W_min[:,out_idx,j], W_max[:,out_idx,j], Water[:,j])

                if snow_flag
                    SnowDays[1,out_idx,j] = SnowDays[1,out_idx,j] + 1
                    SnowDepth_av[1,out_idx,j] = SnowDepth_av[1,out_idx,j] + (SnowDepth[1,j]-SnowDepth_av[1,out_idx,j])/SnowDays[1,out_idx,j];
                    SnowDepth_max[1,out_idx,j] = max(SnowDepth[1,j], SnowDepth_max[1,out_idx,j]);
                end

                if j!=1
                    Qlat_out[:,out_idx,j] = Qlat_out[:,out_idx,j] - lat_flux[:,j];
                else
                    #different Q_lat of first tile which is the sum of all other lateral fluxes
                    Qlat_out[:,out_idx,j] = Qlat_out[:,out_idx,j] + sum(lat_flux,dims=2);
                end

                TDD[:,out_idx,j] = TDD[:,out_idx,j] + T[:,j].*(T[:,j].>=0.0);
                FDD[:,out_idx,j] = FDD[:,out_idx,j] + T[:,j].*(T[:,j].<0.0);
                FrostDays[:,out_idx,j] = FrostDays[:,out_idx,j] + 1.0*(T[:,j].<0.0);
            end
            #increase counter (number)
            if j == 1
                out_N = out_N + 1.0
            end
            if size(T_out,2)>=length(t_span)
                T_out[:,t,j] = copy(T[:,j]);
                Water_out[:,t,j] = copy(Water[:,j]);
                WaterIce_out[:,t,j] = copy(WaterIce[:,j]);
            end


            #-------------------------------------------------------------------
        end
    end

    return OUT
end
end
