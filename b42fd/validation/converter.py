import numpy as np
from scipy.io import loadmat
import json
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
import cartopy

# change this when we have our flight data
PATH = "FTISxprt-20200310_flight2.mat"
indexs = ["vane_AOA","elevator_dte","column_fe","lh_engine_FMF","rh_engine_FMF","lh_engine_itt","rh_engine_itt","lh_engine_OP","rh_engine_OP", "dcoc", "lh_engine_fan_N1","lh_engine_turbine_N2","rh_engine_fan_N1","rh_engine_turbine_N2","lh_engine_FU","rh_engine_FU","delta_a","delta_e","delta_r","Gps_date","Gps_utcSec","Ahrs1_Roll","Ahrs1_Pitch","Fms1_trueHeading","Gps_lat","Gps_long","Ahrs1_bRollRate","Ahrs1_bPitchRate","Ahrs1_bYawRate","Ahrs1_bLongAcc","Ahrs1_bLatAcc","Ahrs1_bNormAcc","Ahrs1_aHdgAcc","Ahrs1_xHdgAcc","Ahrs1_VertAcc","Dadc1_sat","Dadc1_tat","Dadc1_alt","Dadc1_bcAlt","Dadc1_bcAltMb","Dadc1_mach","Dadc1_cas","Dadc1_tas","Dadc1_altRate","measurement_running","measurement_n_rdy","display_graph_state","display_active_screen","time"]

# comment both path and indexs bc they have different indexs
# PATH = "matlab.mat"
# indexs = ["vane_AOA","elevator_dte","column_fe","lh_engine_FMF","rh_engine_FMF","lh_engine_itt","rh_engine_itt","lh_engine_OP","rh_engine_OP","lh_engine_fan_N1","lh_engine_turbine_N2","rh_engine_fan_N1","rh_engine_turbine_N2","lh_engine_FU","rh_engine_FU","delta_a","delta_e","delta_r","Gps_date","Gps_utcSec","Ahrs1_Roll","Ahrs1_Pitch","Fms1_trueHeading","Gps_lat","Gps_long","Ahrs1_bRollRate","Ahrs1_bPitchRate","Ahrs1_bYawRate","Ahrs1_bLongAcc","Ahrs1_bLatAcc","Ahrs1_bNormAcc","Ahrs1_aHdgAcc","Ahrs1_xHdgAcc","Ahrs1_VertAcc","Dadc1_sat","Dadc1_tat","Dadc1_alt","Dadc1_bcAlt","Dadc1_bcAltMb","Dadc1_mach","Dadc1_cas","Dadc1_tas","Dadc1_altRate","measurement_running","measurement_n_rdy","display_graph_state","display_active_screen","time"]


conversions = { "deg/s" : { "unit" : "rad/s",
                          "factor" : lambda x:x*0.0174533},
                "deg" : { "unit" : "rad",
                         "factor" : lambda x:np.radians(x)},
                "lbs/hr" : { "unit" : "kg/s",
                            "factor" : lambda x:x/7936.64},
                "ft" : { "unit" : "m",
                            "factor" : lambda x:x*0.3048},
                "knots" : { "unit" : "m/s",
                            "factor" : lambda x:x*0.514444},
                "ft/min" : { "unit" : "m/s",
                            "factor" : lambda x:x*0.00508},
                "deg C" : { "unit" : "K",
                            "factor" : lambda x:x+273.15}}

mat = loadmat(PATH)
mat = mat["flightdata"][0][0]

dataDict = {}

for index, dat in zip(indexs, mat):
    dat = dat[0][0]
    try:
       lDict = { "data" : dat[0].flatten().tolist(),
                            "unit" : dat[1][0][0][0],
                            "name" : dat[2][0][0][0]}
    except IndexError:
        lDict = { "data" : dat[0].flatten().tolist(),
                            "unit" : "nan",
                            "name" : dat[2][0][0][0]}
    
    if lDict["unit"] in conversions.keys(): # check for conversion

        func = conversions[lDict["unit"]]["factor"] # get the conversion function
        lDict["data"] = func(np.array(lDict["data"])).tolist() # apply it and put the array back in the dict
        lDict["unit"] = conversions[lDict["unit"]]["unit"] # change the unit

    dataDict[index] = lDict

with open("flight_data.json", "w") as f:
    f.write(json.dumps(dataDict,indent=2))