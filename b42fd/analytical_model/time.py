#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 16:45:39 2020

@author: wenweidong
"""

from b42fd.validation.fuelmass import data
from b42fd.helpers import load_data
from pathlib import Path

class TimeTool:
    
    def __init__(self):
        self.data=data    #flight data
        #sefl.data=load_data("data/ref_data/ref_data.json")
        self.time=self.data["time"]["data"]
        self.altitude=self.data["Dadc1_alt"]["data"]
        self.rh_fu=self.data["rh_engine_FU"]["data"]
        self.rh_fu=self.data["rh_engine_FU"]["data"]
        self.lf_fu=self.data["lh_engine_FU"]["data"]
        self.V_TAS=self.data['Dadc_tas']["data"]
        self.alpha=self.data['']["data"]

    
