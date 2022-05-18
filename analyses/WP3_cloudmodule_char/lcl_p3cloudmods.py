# -*- coding: utf-8 -*-
"""
Created on Tue May 17 13:59:41 2022

@author: Dean

Get LCL for each P-3 cloud module.

LCL calcs will use near-surface data from the dropsondes released during each 
cloud module.
"""



# Built in
import os

# Local code
import thermo



path_drpsnd = "../WP3_cloudmodule_char/cldmod_datafiles/" # Path to cloud mod. dropsondes folder.
fnames = [f for f in os.listdir(path_drpsnd) 
          if f.startswith("p3cld_dropsondes") and f.endswith(".nc")
          ]





















