# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 11:40:05 2022

@author: cgomezco
"""
# Assign new roughness from hazen williams friction factor to the Darcy Weisbach friction factor equation
import wntr
import pandas as pnd
import numpy as np
import random as rnd
import math

nomb_red = "MOD"
net_name = "C:/Users/cgomezco/CurrentResearch/Code/BilzenKortessem.inp"
real_name = "C:/Users/cgomezco/CurrentResearch/Code/BilzenKortessem2.inp"
format_cat = "C:/Users/cgomezco/CurrentResearch/Code/Formatos/PipeInfo_BilzenKortessem.xlsx"

# roughness ranges
coef_rug = {"DCI":{"C1":10.533,"C2":1.523,"C3":1.084},
            "GCI":{"C1":22.760,"C2":-4.180,"C3":1.2},
            "PEH":{"C1":0.302,"C2":0.274,"C3":0.993},
            "PVC":{"C1":0.302,"C2":0.274,"C3":0.993},
            "STE":{"C1":74.681,"C2":-5.436,"C3":1.311},
            "FCE":{"C1":7.627,"C2":-0.171,"C3":1.211},
            "CON":{"C1":7.627,"C2":-0.171,"C3":1.211}}

rug_ini = {"DCI":0.26,
           "GCI":0.8,
           "PEH":0.007,
           "PVC":0.0015,
           "STE":0.045,
           "FCE":0.3,
           "CON":0.6}
           
pre_rug = {"DCI":2,
           "GCI":1,
           "PEH":3,
           "PVC":4,
           "STE":3,
           "FCE":1,
           "CON":1}
rang_rug = {}
materiales = coef_rug.keys()
for material in materiales:
    rang_rug[material] = {}
    for i in range(0,50,1):
        rang_rug[material][str(i)] =round(rug_ini[material]*((coef_rug[material]["C1"]*((i/50)**2))+(coef_rug[material]["C2"]*(i/50))+(coef_rug[material]["C3"])),pre_rug[material])

wn = wntr.network.WaterNetworkModel(net_name)
# Load pipe info
sheet = pnd.ExcelFile(format_cat).sheet_names
xl = pnd.read_excel(format_cat,sheet_name=None)
for n in sheet:
    if xl[n].columns[0] != 'Pipes' or xl[n].columns[2] != 'Age':
        format_cat = ""
        print("Error pipe info #1")
        break
    if len(xl[n].iloc[:,0]) != len(wn.pipe_name_list):
        format_cat = ""
        print("Error pipe info #2")
        break
    if pnd.isnull(xl[n].iloc[:,1]).any() == True or np.isnan(xl[n].iloc[:,2]).any() == True:
        format_cat = ""
        print("Error pipe info #3")
        break
    t = 0
    for pipe_name, pipe in wn.pipes():
        if xl[n].iloc[t,1] in materiales and xl[n].iloc[t,2] <= 50:
            pipe.tag = xl[n].iloc[t,1]+"_"+str(int(xl[n].iloc[t,2]))
            t += 1
        else:
            format_cat = ""
            print("Error pipe info #4")
            break
    if format_cat != "":
        print("Load pipe info successfully")
        wn.write_inpfile(wn.name)

wn.options.hydraulic.headloss='D-W'
for pipe_name, pipe in wn.pipes():
    eq = pipe.roughness*(pipe.diameter**0.01)
    if eq > 120:
        pipe.roughness = 0.0015
    else:
        pipe.roughness = round((pipe.diameter)*(math.exp(-0.04125*eq))*((3.32-(0.021*eq))**(2.173)),4)
    if pipe.roughness < 0.0015:
        pipe.roughness = 0.0015
    
    expex_rough = round(rang_rug[pipe.tag[0:3]][pipe.tag[4:6]],pre_rug[pipe.tag[0:3]])
    if (expex_rough/pipe.roughness)>10:
        print("Warning in roughness link"+str(pipe_name))
wn.write_inpfile(real_name)




