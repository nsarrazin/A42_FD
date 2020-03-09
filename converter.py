import numpy as np
from scipy.io import loadmat
import json
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
import cartopy

# change this when we have our flight data
PATH = "matlab.mat"


mat = loadmat(PATH)
indexs = ["vane_AOA","elevator_dte","column_fe","lh_engine_FMF","rh_engine_FMF","lh_engine_itt","rh_engine_itt","lh_engine_OP","rh_engine_OP","lh_engine_fan_N1","lh_engine_turbine_N2","rh_engine_fan_N1","rh_engine_turbine_N2","lh_engine_FU","rh_engine_FU","delta_a","delta_e","delta_r","Gps_date","Gps_utcSec","Ahrs1_Roll","Ahrs1_Pitch","Fms1_trueHeading","Gps_lat","Gps_long","Ahrs1_bRollRate","Ahrs1_bPitchRate","Ahrs1_bYawRate","Ahrs1_bLongAcc","Ahrs1_bLatAcc","Ahrs1_bNormAcc","Ahrs1_aHdgAcc","Ahrs1_xHdgAcc","Ahrs1_VertAcc","Dadc1_sat","Dadc1_tat","Dadc1_alt","Dadc1_bcAlt","Dadc1_bcAltMb","Dadc1_mach","Dadc1_cas","Dadc1_tas","Dadc1_altRate","measurement_running","measurement_n_rdy","display_graph_state","display_active_screen","time"]
mat = mat["flightdata"][0][0]

dataDict = {}

for index, dat in zip(indexs, mat):
    dat = dat[0][0]
    try:
        dataDict[index] = { "data" : dat[0].flatten().tolist(),
                            "unit" : dat[1][0][0][0],
                            "name" : dat[2][0][0][0]}
    except IndexError:
        dataDict[index] = { "data" : dat[0].flatten().tolist(),
                            "unit" : "nan",
                            "name" : dat[2][0][0][0]}
    
    
with open("data.json", "w") as f:
    f.write(json.dumps(dataDict,indent=2))