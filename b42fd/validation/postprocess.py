import numpy as np
import matplotlib.pyplot as plt
import json, control

from b42fd.helpers import load_data, find_nearest
from b42fd.numerical_model.model import A_s_h, B_s_h, C_s_h, D_s_h, A_a_h, B_a_h, C_a_h, D_a_h

class PostProcessing:
    def __init__(self, path, events):
        self.data = load_data(path)
        self.events = events
        self.t = self.data["time"]["data"]

        # self.times = (120, 30, 60, 60, 60)
        self.times = {"phugoid" : 120, 
                      "spm" : 30,
                      "dutchroll" : 60, 
                      "ape_roll" : 60,
                      "ape_spiral" : 60}


    def _plotData(self, key, t_0, t_max):
        t_0_i = np.searchsorted(self.t, t_0)
        t_max_i = np.searchsorted(self.t, t_max)
        plt.plot(self.t[t_0_i:t_max_i], self.data[key]["data"][t_0_i:t_max_i], label=self.data[key]["name"])

    def _getDeflections(self, t_0, t_max):
        t_0_i = np.searchsorted(self.t, t_0)
        t_max_i = np.searchsorted(self.t, t_max)
        
        delta_a = self.data["delta_a"]["data"][t_0_i:t_max_i]
        delta_e = np.array(self.data["delta_e"]["data"][t_0_i:t_max_i])
        delta_r = self.data["delta_r"]["data"][t_0_i:t_max_i]
        t = self.data["time"]["data"][t_0_i:t_max_i]

        return {"t" : t, 
                "delta_a" : delta_a, 
                "delta_e" : delta_e, 
                "delta_r" : delta_r
                }

    def _getInitialCondition(self, key):
        t_0 = self.events[key] - 5
        t_0_i = np.searchsorted(self.t, t_0)

        u_0 = self.data["Dadc1_tas"]["data"][t_0_i]*0
        alpha_0 = self.data["Ahrs1_Pitch"]["data"][t_0_i]
        theta_0 = self.data["Ahrs1_Roll"]["data"][t_0_i]
        q_0 = self.data["Ahrs1_bPitchRate"]["data"][t_0_i]
        return [u_0, alpha_0, theta_0, q_0]

    def plotAngles(self, key_event, t_plot = 60, angles = ["roll", "pitch"], show=True):
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
        if show:
            plt.show(block=False)

    def plotRates(self, key_event, t_plot = 60, angles = ["roll", "pitch", "yaw"], show=True):
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
        
        if show:
            plt.show(block=False)

    def plotInputs(self, key_event, t_plot = 60, inputs = ["delta_a", "delta_e", "delta_r"], show=True):
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
        if show:
            plt.show(block=False)

    def plotAll(self):
        for key,time in zip(self.times.keys(), self.times.values()):
            pp.plotAngles(key, t_plot =time)
            pp.plotRates(key, t_plot=time)
            pp.plotInputs(key, t_plot = time)
            input()
    
    def controlJSON(self, name):
        fullDict = {}
        for key_event, time in zip(self.times.keys(), self.times.values()):
            t_0 = self.events[key_event] - 5
            t_max = t_0 + time + 5

            dictinputs = self._getDeflections(t_0, t_max)
            fullDict[key_event] = dictinputs

        json_raw = json.dumps(fullDict, indent=2)
        with open(name, "w") as f:
            f.write(json_raw)

    def compareNumerical(self, key):
        t_0 = self.events[key] -5
        t_max = t_0 + self.times[key] + 5

        defl = self._getDeflections(t_0, t_max)
        x0 = self._getInitialCondition(key)

        sys_s = control.ss(A_s_h, B_s_h, C_s_h, D_s_h)
        sys_a = control.ss(A_a_h, B_a_h, C_a_h, D_a_h)

        T, yout, bla = control.forced_response(sys_s, defl["t"], \
            defl["delta_e"], X0 = x0)

        self.plotAngles(key, t_plot=self.times[key], show=False)
        plt.plot(T, yout[1, :], label="pitch numerical model")
        plt.plot(T, yout[2, :], label="roll numerical model")
        plt.legend()
        plt.show()

        self.plotRates(key, t_plot=self.times[key], show=False, angles=["roll"])
        plt.plot(T, yout[3, :], label="pitch rate numerical model")
        plt.legend()
        plt.show()
        # # for i in range(yout.shape[0]):

    # def computeParams(self, key):
    #     Thalf, P = 0,0
    #     return Thalf, P

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
    pp.compareNumerical("phugoid")
    # pp.controlJSON("input.json")


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