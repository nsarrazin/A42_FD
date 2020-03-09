import json
import numpy as np
import matplotlib.pyplot as plt

# load the formatted data into a dictionary
with open("data.json", "r") as f:
    raw = f.read()
    dataDict = json.loads(raw)

# structure is as follow
# you have a dictionary whose keys represent a dataset each
#
# for each dataset you have another dictionary which uses the following keys :
# - data (all the data in a list)
# - unit (the unit used to represent the dat)
# - name ( a more verbose description of what the data is)
#
# for example 
#  dataDict["vane_AOA"] = { "data" : [0.0, 0.1, ... etc],
#                           "unit" : "deg",
#                           "name" : "Angle of attack"}

# if you don't know what a dictionary is in python, you might want
# to read up on that
# you can get a list of all the datasets using dataDict.keys()

for key in dataDict.keys():
    print("[\""+key+"\"] = "+dataDict[key]["name"]+ " - " + dataDict[key]["unit"])

from scipy.integrate import cumtrapz

plt.figure(dpi=150)
plt.plot(dataDict["time"]["data"], dataDict["lh_engine_FMF"]["data"], label="Left Engine")
plt.plot(dataDict["time"]["data"], dataDict["rh_engine_FMF"]["data"], label="Right Engine")
plt.xlabel("Time [s]")
plt.ylabel("Fuel mass flow (lbs/hr)")
plt.legend()
plt.show()


plt.figure(dpi=150)
plt.plot(dataDict["time"]["data"], dataDict["lh_engine_itt"]["data"], label="Left")
plt.plot(dataDict["time"]["data"], dataDict["rh_engine_itt"]["data"], label="Right")
plt.xlabel("Time [s]")
plt.ylabel("Engine internal temperature (deg C)")
plt.legend()
plt.show()

lh_fueltot = cumtrapz(np.array(dataDict["lh_engine_FMF"]["data"])/3600, dataDict["time"]["data"])
# i think the unit value for recorded fuel use says L but it should be in lbs
# check with the TA maybe idk
plt.figure(dpi=150)
plt.plot(lh_fueltot, label="Integrated fuel use [lbs]")
plt.plot(np.array(dataDict["lh_engine_FU"]["data"]), label="Recorded fuel use [L???]")
plt.legend()
plt.show()


# sick map plot 
# to get it you have to run 
# conda install -c conda-forge cartopy
# into the anaconda prompt
has_cartopy=True
try:
    import cartopy.crs as ccrs
    import cartopy
except ModuleNotFoundError:
    has_cartopy=False

if has_cartopy:
    alt = dataDict["Dadc1_alt"]["data"]
    plt.figure(dpi=140)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
    ax.add_feature(cartopy.feature.LAKES, alpha=0.8)
    ax.add_feature(cartopy.feature.RIVERS)
    ax.set_xlim([4., 6.])
    ax.set_ylim([51., 53.])
    ax.scatter(dataDict["Gps_long"]["data"],dataDict["Gps_lat"]["data"], c=alt, s=0.05, cmap='coolwarm')
    plt.show()


# sick 3d trajectory plot
# conda install -c conda-forge mayavi

from helpers import wgs84_to_ecef, TEC

lon = np.radians(np.array(dataDict["Gps_long"]["data"]))
lat = np.radians(np.array(dataDict["Gps_lat"]["data"]))
x,y,z = [], [], []
for lo,la,al in zip(lon,lat,alt):
    x_ecef = wgs84_to_ecef(np.array([lo, la, al]))
    x_ned = x_ecef @ TEC(lo, la)
    x.append(x_ecef[0])
    y.append(x_ecef[1])
    z.append(x_ecef[2])

v = np.array(dataDict["Dadc1_tas"]["data"])
has_mayavi = True
try:
    from mayavi import mlab
except ModuleNotFoundError:
    has_mayavi=False

if has_mayavi:
    mlab.show(mlab.plot3d(x,y,z,v, tube_radius=100,colormap='Spectral'))
    input()