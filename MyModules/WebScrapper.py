#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 12:03:31 2019

@author: cmosbeux
"""

import sys
import os
import requests

url = "https://daacdata.apps.nsidc.org/pub/DATASETS/ICEBRIDGE/IAMET1B_NSERCmetXyMet_v01/"
url = 'https://www.nature.com/articles/s41561-019-0510-82'
source = requests.get(url)
source.text

#%%
mycommand = 'wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies  --cut-dirs=4 --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r --reject "index.html*" -np -e robots=off'+address
os.system(mycommand)