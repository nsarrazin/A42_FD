#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import numpy as np
from b42fd.helpers import load_data
from pathlib import Path
import os

path = "data/flight_data/flight_data.json"

data = load_data(path)
FFl=data["lh_engine_FMF"]["data"]
FFr=data["lh_engine_FMF"]["data"]



