from math import pi, sin, cos
from pathlib import Path
import numpy as np


def calculate_eigenvalues(A, B, C, V, c):
    return complex((-B / (2*A)) * V / c, np.sqrt(4 * A * C - B ** 2) / (2*A) * V / c), \
           complex(-B / (2*A) * V / c, -np.sqrt(4 * A * C - B ** 2) / (2*A) * V / c)


def short_period_eigenvalues(data):
    "data given as dict"
    A = 2*data['muc']*data['KY2'] * (2*data['muc'] - data['CZadot'])
    B = -2*data['muc']*data['KY2']*data['CZa'] - (2*data['muc'] + data['CZq'])
    C = -2*data['CZa']*data['Cmq'] - (2*data['muc'] + data['CZq']) * data['Cma']
    return calculate_eigenvalues(A, B, C, data['V0'], data['c'])


def phugoid_eigenvalues(data):
    "data given as dict"
    A = 2*data['muc']*(data['CZa']*data['Cmq'] - 2*data['muc']*data['Cma'])
    B = 2*data['muc']*(data['CXu']*data['Cma'] - data['Cmu']*data['CXa']) + data['Cmq']*(data['CZu']*data['CXa']-data['CXu']*data['CZa'])
    C = data['CZ0']*(data['Cmu']*data['CZa']-data['Cma']*data['CZu'])
    return calculate_eigenvalues(A, B, C, data['V0'], data['c'])

def dutch_roll_eigenvalues(data):
    A = -2*data['mub']*data['KX2']
    B = 0.5*data['Cnr']
    C = -1*data['Cnb']
    return calculate_eigenvalues(A, B, C, data['V0'], data['c'])

