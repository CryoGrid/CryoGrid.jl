module ArcticLakesSetup
using MAT
include("CryoGridTyps.jl")
include("CryoGridImplicit.jl")
include("matlab.jl")

function SetUpInputStructs(FORCING, GRID, PARA, lakestat, tile_number)
    #-------------------------------------------------------------------------------
    STRAT=Dict()
    STATVAR=Dict()
    TEMP=Dict()
    lc_idx = (lakestat["LakeNumber"].>0)

    Z = [collect(-2.10:0.01:-0.01); PARA["technical"]["subsurfaceGrid"]] #grid with soil depth (upper edge of grid cell!)
    GRID["dxp"] = diff(Z); #grid cell sizes
    GRID["Zp"] = Z[1:end-1]+GRID["dxp"]/2.; #grid cell midpoints
    N = length(GRID["Zp"]); #number of grid cells

    PARA["WaterDepth"] = [0.0 lakestat["LakeDepth"][lc_idx]']
    PARA["WaterDensity"] = 1000.0; #[kg/m³]
    PARA["SnowDensity"] = PARA["snow"]["rho_snow"]; #[kg/m³]
    PARA["Qgeo"] = PARA["soil"]["Qgeo"]; #[W/m²] geothermal heat flux
    PARA["SnowCoverMax"] = [2.0 lakestat["SnowDepth"][lc_idx]']; #maximum snow cover thickness on landscape unit (if empty snow cover is limited by snowfall and Z domain only)

    #Initialize soil stratigraphie, status variables, temporary variables
    STRAT["WaterIce"] = zeros(N,tile_number); #volume fraction Water+Ice
    STRAT["Water"] = zeros(N,tile_number); #volume fraction liquid Water
    STRAT["Mineral"] = zeros(N,tile_number); #volume fraction other solid soil components
    STRAT["Organic"] = zeros(N,tile_number); #volume fraction other solid soil components
    STATVAR["T"] = zeros(N,tile_number); #[°C] status variable soil temperature
    STATVAR["H"] = zeros(N,tile_number); #[J/m^3] status variable Enthalpy
    TEMP["cp"] = zeros(N,tile_number);
    TEMP["kp"] = zeros(N,tile_number);

    for j = 1:tile_number
        STRAT["WaterIce"][GRID["Zp"].>0.0,j] = GRID["soil"]["cT_water"]; #1e-8; #at the moment WaterIce=0 not possible, use a very small number e.g. 1e-8!!!!
        STRAT["Mineral"][GRID["Zp"].>0.0,j] = GRID["soil"]["cT_mineral"];# [GRID["Zp"].>PARA["WaterDepth"][j],j] = 0.45; #1.0-1e-8;
        STRAT["Organic"][GRID["Zp"].>0.0,j] = GRID["soil"]["cT_organic"];# [GRID["Zp"].>PARA["WaterDepth"][j],j] = 0.05; #1.0-1e-8;

        STRAT["WaterIce"][GRID["Zp"].<=PARA["WaterDepth"][j],j] = 1.0;
        STRAT["Mineral"][GRID["Zp"].<=PARA["WaterDepth"][j],j] = 0.0;
        STRAT["Organic"][GRID["Zp"].<=PARA["WaterDepth"][j],j] = 0.0;

        if PARA["WaterDepth"][j]>0.0
            STRAT["WaterIce"][GRID["Zp"].>PARA["WaterDepth"][j],j] = 1.0 - STRAT["Mineral"][GRID["Zp"].>PARA["WaterDepth"][j],j] -STRAT["Organic"][GRID["Zp"].>PARA["WaterDepth"][j],j]
        end

        STRAT["Water"][:,j] = 0.0;
        STRAT["WaterIce"][GRID["Zp"].<=0.0,j] = 0.0 #no snow #PARA["SnowDensity"]/PARA["WaterDensity"];
        STRAT["Mineral"][GRID["Zp"].<=0.0,j] = 0.0;
        STRAT["Organic"][GRID["Zp"].<=0.0,j] = 0.0;

        TEMP["cp"][:,j] = CryoGridImplicit.HeatCapacity(STRAT["WaterIce"][:,j], STRAT["Water"][:,j], STRAT["Mineral"][:,j],STRAT["Organic"][:,j],TEMP["cp"][:,j]);
        TEMP["kp"][:,j] = CryoGridImplicit.ThermalConductivity(STRAT["WaterIce"][:,j], STRAT["Water"][:,j], STRAT["Mineral"][:,j],STRAT["Organic"][:,j],TEMP["kp"][:,j]);
    end
    TEMP["lat_flux"] = zeros(N,tile_number);
    TEMP["SnowDepth"] = zeros(1,tile_number);

    #Grid
    #MakeGrid(Zp, dxp, dxn, dxs, kp, kn, ks)
    GRID["dxn"] = copy(GRID["dxp"])
    GRID["dxs"] = copy(GRID["dxp"])
    TEMP["kn"] = copy(TEMP["kp"][:,1])
    TEMP["ks"] = copy(TEMP["kp"][:,1])

    GRID["Zn"], GRID["Zs"], GRID["dxn"], GRID["dxs"], TEMP["kn"], TEMP["ks"] = CryoGridImplicit.MakeGrid(GRID["Zp"],GRID["dxp"],GRID["dxn"],GRID["dxs"],TEMP["kp"][:,1],TEMP["kn"][:,1],TEMP["ks"][:,1]);
    GRID["dxo"] = [NaN lakestat["LakeDistance"][lc_idx]'/2.0]#half the average distance between lake class #[NaN 300.0 300.0 300.0 300.0 300.0 300.0 300.0 300.0 300.0 300.0 300.0 300.0 300.0 300.0 300.0 300.0]; #lateral distance between landscape units[m]
    LandArea = sum(lakestat["VorArea"][lc_idx])-sum(lakestat["LakeArea"][lc_idx])

    #if no lakes exisit
    if LandArea == 0.0 || tile_number <= 1
        LandArea = lakestat["LandMaskArea"]
        tile_number = 1
    end

    GRID["An"] = [LandArea lakestat["LakeArea"][lc_idx]'] # [10000.0 500.0 500.0 500.0 500.0 200.0 200.0 200.0 200.0 100.0 100.0 100.0 100.0 50.0 50.0 50.0 50.0]; #area of landscape unit [m²]
    GRID["As"] = GRID["An"];
    GRID["Ao"] = GRID["dxp"].*[NaN lakestat["LakePerimeter"][lc_idx]']#GRID["dxp"].*(GRID["An"].^0.5); #lateral interface area of each layer [m²]
    GRID["Vp"] = GRID["An"].*GRID["dxp"]; #volume of each layer [m³]

    #State variables
    for j = 1:tile_number
        forcing_days = length(FORCING["data"]["Tair"])
        if forcing_days <= 100 *365
            a = forcing_days
        else
            a = 100*365
        end
        #STATVAR["T"][:,j] = mean(FORCING["data"]["Tair"][1:a])+1.5 + PARA["Qgeo"]/2.0*GRID["Zp"];

        idx = indmin(abs(GRID["Zp"]-0.0))
        STATVAR["T"][1:idx,j] = mean(FORCING["data"]["Tair"][1:a])+1.5
        for i = idx:length(GRID["Zp"])-1
            STATVAR["T"][i+1,j] = PARA["Qgeo"]*(GRID["Zp"][i+1]-GRID["Zp"][i])/TEMP["kn"][i] + STATVAR["T"][i,j]
        end

        STATVAR["H"][:,j], dummy = CryoGridImplicit.Enthalpy(STATVAR["T"][:,j], STRAT["Water"][:,j], TEMP["cp"][:,j],STATVAR["H"][:,j],STATVAR["H"][:,j]);
    end

    TEM = CryoGridTyps.temporary(TEMP["lat_flux"],TEMP["SnowDepth"],TEMP["kp"],TEMP["cp"],TEMP["kn"],TEMP["ks"])
    FOR = CryoGridTyps.forcing(FORCING["data"]["t_span"],FORCING["data"]["Tair"],FORCING["data"]["snowfall"])
    GRI = CryoGridTyps.grid(GRID["Zp"],GRID["Zn"],GRID["Zs"],GRID["dxo"],GRID["dxp"],GRID["dxn"],GRID["dxs"],GRID["An"],GRID["As"],GRID["Ao"],GRID["Vp"])
    STRA = CryoGridTyps.stratigraphy(STRAT["Water"],STRAT["Mineral"],STRAT["Organic"],STRAT["WaterIce"])
    STAT = CryoGridTyps.statvar(STATVAR["T"],STATVAR["H"])
    PAR = CryoGridTyps.para([PARA["snow"]["rho_snow"]],[PARA["snow"]["rho_snow_max"]],[PARA["snow"]["k1_snow"]],[PARA["snow"]["k2_snow"]],PARA["SnowCoverMax"],[PARA["WaterDensity"]],PARA["WaterDepth"],[PARA["Qgeo"]])

    #Preallocate output arrays
    out_years = unique(Dates.year.(matlab.datestr(FOR.t_span)))'
    out_length = length(out_years)
    out_intervall = Dates.Year(1)
    out_date = Dates.DateTime(Dates.DateTime(out_years[1])+out_intervall)
    out_length = length(out_years)

    T = ones(1,1,tile_number).*NaN #dummy for reduced output
    Water = ones(1,1,tile_number).*NaN #dummy for reduced output
    WaterIce = ones(1,1,tile_number).*NaN #dummy for reduced output
    #T = ones(N,length(FOR.t_span),size(STAT.T,2)).*NaN; #output for full temperature field
    #Water = ones(N,length(FOR.t_span),size(STAT.T,2)).*NaN; #output for full temperature field
    #WaterIce = ones(N,length(FOR.t_span),size(STAT.T,2)).*NaN; #output for full temperature field

    H_av = zeros(N,out_length+1,tile_number);
    T_av = zeros(N,out_length+1,tile_number);
    T_min = zeros(N,out_length+1,tile_number);
    T_max = zeros(N,out_length+1,tile_number);
    W_av = zeros(N,out_length+1,tile_number);
    W_min = zeros(N,out_length+1,tile_number);
    W_max = zeros(N,out_length+1,tile_number);

    Q_lat = zeros(N,out_length+1,tile_number);

    FDD = zeros(N,out_length+1,tile_number);
    TDD = zeros(N,out_length+1,tile_number);
    FrostDays = round.(zeros(N,out_length+1,tile_number));

    SnowDepth_av = zeros(1,out_length+1,tile_number);
    SnowDepth_max = zeros(1,out_length+1,tile_number);
    SnowDays = zeros(1,out_length+1,tile_number);

    OUT = CryoGridTyps.out(out_years, T, Water, WaterIce, H_av, T_av, T_min, T_max, W_av, W_min, W_max, Q_lat, FDD, TDD, FrostDays, SnowDepth_av, SnowDepth_max, SnowDays)

    return FOR, PAR, TEM, GRI, STRA, STAT, OUT
end

end
