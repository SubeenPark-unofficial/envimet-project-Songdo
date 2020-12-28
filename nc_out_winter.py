from matplotlib import pyplot as plt
import pandas as pd
import netCDF4 as nc
from scipy.io import netcdf
import numpy as np
import matplotlib as mpl


################### SIMPLE STATISTICS #####################
mins, maxes, meanss = dict(), dict(), dict()

for i in range(7,18):
    
    ################################ READ NC FILES ############################
    f_dir= '/home/subeen/ENVI/output/WINTER/'

    f_PMV1 = '/PMV/exp1_BIO_PMV_2018-12-21_' # exp1_summer_AT_2018-12-21_07.00.01
    f_PMV2 = '.00.01.nc'

    f_AT1 = '/ATMO/exp1_AT_2018-12-21_' # exp1_summer_AT_2018-12-21_07.00.01
    f_AT2 = '.00.01.nc'



    f_PMV = f_dir + f_PMV1 + f"{i:02}" + f_PMV2 
    f_AT = f_dir + f_AT1 + f"{i:02}" + f_AT2 

    ################################ READ VARIABLE VALUES ############################
    data_PMV = nc.Dataset(f_PMV)
    data_AT = nc.Dataset(f_AT)

    text_PMV = 'float64 nx(nx), float64 ny(ny), float64 nz(nz), float64 Objects(nz,ny,nx), float64 PMV(nz,ny,nx), float64 PPD(nz,ny,nx), float64 TCloths(nz,ny,nx), float64 Windspeed(nz,ny,nx), float64 Airtemperature(nz,ny,nx), float64 MeanRadiantTemperature(nz,ny,nx), float64 Specifichumidity(nz,ny,nx)'
    vars_PMV = [word.split('(')[0] for word in text_PMV.split('float64 ')[1:]]

    text_AT= 'float64 nx(nx), float64 ny(ny), float64 nz(nz), float64 Objects(nz,ny,nx), float64 Flowu(nz,ny,nx), float64 Flowv(nz,ny,nx), float64 Floww(nz,ny,nx), float64 WindSpeed(nz,ny,nx), float64 WindSpeedChange(nz,ny,nx), float64 WindDirection(nz,ny,nx), float64 PressurePerturbation(nz,ny,nx), float64 PotentialAirTemperature(nz,ny,nx), float64 AirTemperatureDelta(nz,ny,nx), float64 AirTemperatureChange(nz,ny,nx), float64 SpecHumidity(nz,ny,nx), float64 RelativeHumidity(nz,ny,nx), float64 TKE(nz,ny,nx), float64 Dissipation(nz,ny,nx), float64 VerticalExchangeCoefImpuls(nz,ny,nx), float64 HorizontalExchangeCoefImpuls(nz,ny,nx), float64 VegetationLAD(nz,ny,nx), float64 DirectSwRadiation(nz,ny,nx), float64 DiffuseSwRadiation(nz,ny,nx), float64 ReflectedSwRadiation(nz,ny,nx), float64 TemperatureFlux(nz,ny,nx), float64 VapourFlux(nz,ny,nx), float64 WateronLeafes(nz,ny,nx), float64 LeafTemperature(nz,ny,nx), float64 LocalMixingLength(nz,ny,nx), float64 MeanRadiantTemp(nz,ny,nx), float64 TKEnormalised1D(nz,ny,nx), float64 Dissipationnormalised1D(nz,ny,nx), float64 Kmnormalised1D(nz,ny,nx), float64 TKEMechanicalTurbulenceProd(nz,ny,nx), float64 StomataResistance(nz,ny,nx), float64 CO2.(nz,ny,nx), float64 CO2(nz,ny,nx), float64 PlantCO2Flux(nz,ny,nx), float64 DivRlwTempchange(nz,ny,nx), float64 BuildingNumber(nz,ny,nx)'
    vars_AT = [word.split('(')[0] for word in text_AT.split('float64 ')[1:]]


    ################################ SEPERATE VALUES BY DIMENSIONS ############################
    data1d , data3d = dict(), dict()
    for var in vars_PMV:
        nd = len(np.shape(data_PMV[var]))
        if nd == 1:
            data1d[var] = data_PMV[var]
        elif nd == 3:
            data3d[var] = data_PMV[var]

    for var in vars_AT:
        nd = len(np.shape(data_AT[var]))
        if nd == 1:
            data1d[var] = data_AT[var]
        elif nd == 3:
            data3d[var] = data_AT[var]

    ################################ MODEL DIMENSIONS/VARIABLE LISTS ############################
    nx, ny, nz = len(data1d['nx']), len(data1d['ny']), len(data1d['nz'])
    key3d = data3d.keys()

    ################################ BUILDING LOCATIONS ############################
    BLDG = np.empty((nz, ny, nx))
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                if data3d['PMV'][z, y, x] <= -800:
                    BLDG[z, y, x] = 1
                else:
                    BLDG[z, y, x] = 0

    ############################### add_nan FUNCTION ############################
    def add_nan(data3d, bldg = BLDG):
        res = np.empty((nz, ny, nx))
        for z in range(nz):
            for y in range(ny):
                for x in range(nx):
                    if bldg[z, y, x] == 1:
                        res[z,y,x] = np.NaN
                    else:
                        res[z,y,x] = 1.0*data3d[z,y,x]     
        return res

    ############################### DATA WITH BUILDING LOC NANNED ############################
    data3d_nan = dict()
    for var in key3d:
        data3d_nan[var] = add_nan(data3d[var])
        
        
    ############################# EXPORT DATA AS NC FILE #####################################
    f_dir_out= '/home/subeen/ENVI/output/WINTER/NaN/'
    f_out1 = 'NaN_exp1_winter_2018-12-21_' # exp1_summer_AT_2018-12-21_07.00.01
    f_out2 = '.00.01.nc'
    f_out = f_dir_out + f_out1 + f"{i:02}" + f_out2 
    nc_out =  nc.Dataset(f_out, "w", format="NETCDF4")
    
    # Dimensions
    nx = nc_out.createDimension('nx', nx)
    ny = nc_out.createDimension('ny', ny)
    nz = nc_out.createDimension('nz', nz)
    
    # Variables
    for var in key3d:
        data3d[var] = nc_out.createVariable(var, "f8", ("nz", "ny", "nz"))        

        
