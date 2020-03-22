import numpy as np
import matplotlib.pyplot as plt
import json
from b42fd.helpers import load_data, find_nearest

class PostProcessing:
    def __init__(self, path, events):
        self.data = load_data(path)
        self.events = events
        self.t = self.data["time"]["data"]

        self.times = (120, 30, 60, 60, 60)
        self.keys = ("phugoid", "spm", "dutchroll", "ape_roll", "ape_spiral")

    def plotAngles(self, key_event, t_plot = 60, angles = ["roll", "pitch"]):
        """
        key_event [str] : the key used to get to the flight event like phugoid, dutchroll, etc.
        t_plot [float] : how long to plot for
        angles [list] : contains the keys of which angles you want to plot
        """
        t_0 = self.events[key_event] - 10
        t_max = t_0 + t_plot + 10 


        plt.figure()

        if "roll" in angles:
            self._plotData("Ahrs1_Roll", t_0, t_max)
        if "pitch" in angles:
            self._plotData("Ahrs1_Pitch", t_0, t_max)
        # if "yaw" in angles:
        #     self._plotData("Fms1_trueHeading", t_0, t_max)

        plt.axvline(x=self.events[key_event], label=f"Time of {key_event}", c="black", linestyle="dashed")

        plt.xlabel("Time [s]")
        plt.ylabel("Angle [rad]")
        plt.legend()
        plt.show(block=False)

    def plotRates(self, key_event, t_plot = 60, angles = ["roll", "pitch", "yaw"]):
        """
        key_event [str] : the key used to get to the flight event like phugoid, dutchroll, etc.
        t_plot [float] : how long to plot for
        angles [list] : contains the keys of which angles you want to plot
        """
        t_0 = self.events[key_event] - 10
        t_max = t_0 + t_plot + 10 


        plt.figure()

        if "roll" in angles:
            self._plotData("Ahrs1_bRollRate", t_0, t_max)
        if "pitch" in angles:
            self._plotData("Ahrs1_bPitchRate", t_0, t_max)
        if "yaw" in angles:
            self._plotData("Ahrs1_bYawRate", t_0, t_max)

        plt.axvline(x=self.events[key_event], label=f"Time of {key_event}", c="black", linestyle="dashed")

        plt.xlabel("Time [s]")
        plt.ylabel("Angular rate [rad/s]")
        plt.legend()
        plt.show(block=False)

    def plotInputs(self, key_event, t_plot = 60, inputs = ["delta_a", "delta_e", "delta_r"]):
        t_0 = self.events[key_event] - 10
        t_max = t_0 + t_plot + 10

        plt.figure()

        if "delta_a" in inputs:
            self._plotData("delta_a", t_0, t_max)

        if "delta_e" in inputs:
            self._plotData("delta_e", t_0, t_max)

        if "delta_r" in inputs:
            self._plotData("delta_r", t_0, t_max)

        if "elevator_dte" in inputs:
            self._plotData("elevator_dte", t_0, t_max)
        plt.axvline(x=self.events[key_event], label=f"Time of {key_event}", c="black", linestyle="dashed")

        plt.xlabel("Time [s]")
        plt.ylabel("Angle [rad]")
        plt.legend()
        plt.show(block=False)

    def _plotData(self, key, t_0, t_max):
        t_0_i = np.searchsorted(self.t, t_0)
        t_max_i = np.searchsorted(self.t, t_max)
        plt.plot(self.t[t_0_i:t_max_i], self.data[key]["data"][t_0_i:t_max_i], label=self.data[key]["name"])

    def plotAll(self):
        for key,time in zip(self.keys, self.times):
            pp.plotAngles(key, t_plot =time)
            pp.plotRates(key, t_plot=time)
            pp.plotInputs(key, t_plot = time)
            input()
    
    def _getDeflections(self, t_0, t_max):
        t_0_i = np.searchsorted(self.t, t_0)
        t_max_i = np.searchsorted(self.t, t_max)
        
        delta_a = self.data["delta_a"]["data"][t_0_i:t_max_i]
        delta_e = self.data["delta_e"]["data"][t_0_i:t_max_i]
        delta_r = self.data["delta_r"]["data"][t_0_i:t_max_i]
        t = self.data["time"]["data"][t_0_i:t_max_i]

        return {"t" : t, 
                "delta_a" : delta_a, 
                "delta_e" : delta_e, 
                "delta_r" : delta_r
                }

    def controlJSON(self, name):
        fullDict = {}
        for key_event, time in zip(self.keys, self.times):
            t_0 = self.events[key_event] - 5
            t_max = t_0 + time + 5

            dictinputs = self._getDeflections(t_0, t_max)
            fullDict[key_event] = dictinputs

        json_raw = json.dumps(fullDict, indent=2)
        with open(name, "w") as f:
            f.write(json_raw)

    def computeParams(self, key):
        Thalf, P = 0,0
        return Thalf, P

if __name__ == "__main__":
    events_ref = {"phugoid" : 53*60+57,
                     "spm" : 60*60+35,
                     "dutchroll" : 60*60,
                     "ape_roll" : 59*60+10,
                     "ape_spiral" : 62*60+20}
    
    events_flight = {"phugoid" : 53*60+38,
                     "spm" : 58*60+39,
                     "dutchroll" : 60*60+10,
                     "ape_roll" : 57*60+12.5,
                     "ape_spiral" : 62*60+30}
    
    pp = PostProcessing("data/flight_data/flight_data.json", events_flight)
    # pp.plotAll()
    pp.controlJSON("input.json")


    # for flight data :
    # PHUGOID
    # P =
    # T_1/2 = 

    # SPM
    # P =
    # T_1/2 = 

    # DUTCHROLL
    # P =
    # T_1/2 = 

    # APERIODIC ROLL
    # P =
    # T_1/2 = 

    # SPIRAL
    # P =
    # T_1/2 = 