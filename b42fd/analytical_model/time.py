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
        self.mtime_flat = self.time.flatten()
        self.mtime_list = self.mtime_flat.tolist()
        self.mtime_list_rounded = [round(t, 1) for t in self.mtime_list]
        
    def get_mdat_tstep_list_idx_for_matching_pdat_tstep(self, pdat_t):
        idx = self.mtime_list_rounded.index(round(pdat_t,1))
        return idx

    def get_t_specific_mdat_values(self, pdat_t):
        pdat_t_idx = self.get_mdat_tstep_list_idx_for_matching_pdat_tstep(pdat_t)
        vars_keys_list = list(self.mdat.keys())
        t_specific_mdat = {}
        for var_key in vars_keys_list:
            t_specific_mdat[var_key] = self.mdat[var_key][pdat_t_idx]
        print("At t= {0} the corresponding recorded 'black-box' data is:\n {1}".format(pdat_t, t_specific_mdat))
        return t_specific_mdat
    
if __name__ == "__main__":
      ts_tool = TimeTool()
      t = 5171
      specific_t_mdat_vals = ts_tool.get_t_specific_mdat_values(t)
      print("At t= {0} the corresponding recorded 'black-box' data is:\n {1}".format(t, specific_t_mdat_vals))
    # print(ts_tool.get_t_specific_mdat_values(1665))
