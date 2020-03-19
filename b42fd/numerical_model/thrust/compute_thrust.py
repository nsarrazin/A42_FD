import numpy as np
import os 
from sys import platform as _platform

def compute_thrust(inputs):
    input_str = " ".join([str(i) for i in inputs.tolist()]) # build the input string

    # go to the folder where the script is located
    old_dir = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # write the input string to the file
    with open("matlab.dat", "w") as f:
        f.write(input_str)

    # run the program
    if _platform == "linux" or _platform == "linux2" or _platform == "darwin":
        os.system("wine thrust.exe")
    elif _platform == "win32":
        os.system("thrust.exe")
    else:
        raise RuntimeError("Operating system not supported.")

    # read the output
    with open("thrust.dat", "r") as f:
        thrusts = f.read().split()
    
    # get the summed thrust
    thrust = sum([float(i) for i in thrusts])

    # go back to the old directory
    os.chdir(old_dir)
    
    return thrust
