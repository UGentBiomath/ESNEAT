# -*- coding: utf-8 -*-
"""
Created on Sun Aug 14 15:22:25 2022

@author: Cristian Camilo Gómez Cortés

### Description:
Create a shp file with the warning issues raised before the calibration process

### The input parameters are the following:
    - network: (str) - File path to the .inp file of the current network
    - format_measurement: (str) - File path to the .xlsx file which contains all the measurement related information.
### The output parameters are the following:
    - Print list of current warning of the network
    - Link and pipe warning layer in a .html and .shp file
"""
# Load required packages
import wntr 
import pandas as pnd
import numpy as np
from statistics import mean
import shapefile
from spotpy.objectivefunctions import nashsutcliffe as nse

class warning:
    def __init__(self, network, format_measurement):
        #Load current network and measured information
        elem = pnd.DataFrame(columns=["Type","ID"])
        Current_measure = []
        wn = wntr.network.WaterNetworkModel(network)
        wnp = wntr.network.WaterNetworkModel(network)
        wnf = wntr.network.WaterNetworkModel(network)
        xl = pnd.read_excel(format_measurement,sheet_name=None,header=None)
        inter = int(wn.options.time.hydraulic_timestep/60)
        leng = int(wn.options.time.duration/3600+inter/60)
        sheet = pnd.ExcelFile(format_measurement).sheet_names
        # Verify the measured format and validate the information
        x =[]
        y = []
        if xl[sheet[0]].iloc[0,0] != 'Type' or xl[sheet[0]].iloc[1,0] != 'ID':
            format_measurement = ""
            print("Error measurements #1: Incorrect format")
            
        if len(xl[sheet[0]].iloc[:,0]) != leng/(inter/60)+3:
            format_measurement = ""
            print("Error measurements #2: Incorrect timestep or time lenght")
        
        if all(p == "T" or p == "N" for p in xl[sheet[0]].iloc[0,1:]):
            x.extend(xl[sheet[0]].iloc[0,1:])
        else:
            format_measurement = ""
            print("Error measurements #3: Incorrect element type. Must be T or N")
        
        if all(str(p) in wn.pipe_name_list or str(p) in wn.junction_name_list for p in xl[sheet[0]].iloc[1,1:]):
            y.extend([str(i) for i in xl[sheet[0]].iloc[1,1:]])
        else:
            format_measurement = ""
            print("Error measurements #4: One or more of the current elements is not included in the network")
        
        for a in xl[sheet[0]].columns[1:]:
            Current_measure.extend(xl[sheet[0]][a][2:])
        
        if np.isnan(Current_measure).any() == True:
            format_measurement = ""
            print("Error measurements #5: There are missing data in the current measurement file")
        if format_measurement != "":
            elem["Type"] = x
            elem["ID"] = y
            print("Load measurements successfully")
        
        #Simulate starting conditions
        sim_test = wntr.sim.EpanetSimulator(wn)
        results = sim_test.run_sim()
        
        #Create the conditions to test the flow case
        for valve_name, valve in wnf.valves():
            if valve.valve_type == "TCV":                                  
                valve.initial_status = wntr.network.base.LinkStatus["Active"]
        for node_name, node in wnf.junctions():
            node.demand_timeseries_list[0].base_value = node.demand_timeseries_list[0].base_value*2
        
        sim_testf = wntr.sim.EpanetSimulator(wnf)
        resultsf = sim_testf.run_sim()
        
        #Create the conditions to test the pressure case
        for node_name, node in wnp.junctions():
            node.emitter_coefficient = 0
            node.demand_timeseries_list[0].base_value = 0
        for pipe_name, pipe in wnp.pipes():
            pipe.roughness = 0
            pipe.minor_loss = 0
        sim_testp = wntr.sim.EpanetSimulator(wnp)
        resultsp = sim_testp.run_sim()
        
        #List of errors per category (high flow, high mean flow, pattern flow, high pressure, high mean pressure, pattern pressure)
        hf = []
        hmf = []
        pf = []
        hp = []
        hmp = []
        pp = []
        
        #Pressure review
        d = 0
        div = len(Current_measure)/len(elem["Type"])
        for i, j in elem.iterrows():
            if j["Type"] == "T":
                if max(Current_measure[int(d*div):int(((d+1)*div))]) > max(resultsf.link['flowrate'].loc[:,j["ID"]].tolist()):
                    print("Warning high max flow link "+str(j["ID"]))
                    hf.append(str(j["ID"]))
                if mean(Current_measure[int(d*div):int(((d+1)*div))]) > mean(resultsf.link['flowrate'].loc[:,j["ID"]].tolist()):
                    print("Warning high mean flow link "+str(j["ID"]))
                    hmf.append(str(j["ID"]))
                if (nse([float(i)/max(Current_measure[int(d*div):int(((d+1)*div))]) for i in Current_measure[int(d*div):int(((d+1)*div))]],[float(i)/max(results.link['flowrate'].loc[:,j["ID"]].tolist()) for i in results.link['flowrate'].loc[:,j["ID"]].tolist()])) < -1:
                    print("Warning pattern flow link "+str(j["ID"]))
                    pf.append(str(j["ID"]))    
                d += 1
            if j["Type"] == "N":
                if max(Current_measure[int(d*div):int(((d+1)*div))]) > max(resultsp.node['pressure'].loc[:,j["ID"]].tolist()):
                    print("Warning high max pressure node "+str(j["ID"]))
                    hp.append(str(j["ID"]))
                if mean(Current_measure[int(d*div):int(((d+1)*div))]) > mean(resultsp.node['pressure'].loc[:,j["ID"]].tolist()):
                    print("Warning high mean pressure node "+str(j["ID"]))
                    hmp.append(str(j["ID"]))
                if (nse([float(i)/max(Current_measure[int(d*div):int(((d+1)*div))]) for i in Current_measure[int(d*div):int(((d+1)*div))]],[float(i)/max(results.node['pressure'].loc[:,j["ID"]].tolist()) for i in results.node['pressure'].loc[:,j["ID"]].tolist()])) < -1:
                    print("Warning pattern pressure link "+str(j["ID"]))
                    pp.append(str(j["ID"])) 
                d += 1
                
        #Raise pressure and flow errors
        attf = pnd.Series([0]*len(wnf.pipe_name_list), index = wnf.pipe_name_list,name="Warning #")
        popupf = pnd.DataFrame(columns=['Warning high max flow', 'Warning high mean flow', 'Warning pattern flow'],index = wnf.pipe_name_list)
        attp = pnd.Series([0]*len(wnp.junction_name_list), index = wnp.junction_name_list,name="Warning #")
        popupp = pnd.DataFrame(columns=['Warning high max pressure','Warning high mean pressure', 'Warning pattern pressure'],index = wnp.junction_name_list)
        for a in hf:
            attf[a] += 1
            popupf['Warning high max flow'][a] = 1
        for a in hmf:
            attf[a] += 1
            popupf['Warning high mean flow'][a] = 1
        for a in pf:
            attf[a] += 1
            popupf['Warning pattern flow'][a] = 1
            
        for a in hmp:
            attp[a] += 1
            popupp['Warning high max pressure'][a] = 1
        for a in hp:
            attp[a] += 1
            popupp['Warning high mean pressure'][a] = 1
        for a in pp:
            attp[a] += 1
            popupp['Warning pattern pressure'][a] = 1
        #Plot warnings
        #Use coordinates of 2 nodes
        longlat_map_f = {'774942':(5.35906, 50.8789), 'WT_Blz': (5.5020, 50.8749)}
        wn3 = wntr.morph.convert_node_coordinates_to_longlat(wnf, longlat_map_f)
        wntr.graphics.plot_leaflet_network(wn3,link_attribute_name="Warnings #",link_attribute=attf,filename='Warnings_links.html',link_range=[0,2],node_size=1,link_cmap=['grey','gold','firebrick'],link_width=4,add_to_link_popup=popupf,add_legend=True)
        
        #Use coordinates of 2 nodes
        longlat_map_f = {'774942':(5.35906, 50.8789), 'WT_Blz': (5.5020, 50.8749)}
        wn3 = wntr.morph.convert_node_coordinates_to_longlat(wnp, longlat_map_f)
        wntr.graphics.plot_leaflet_network(wn3,node_attribute_name="Warnings #",node_attribute=attp,filename='Warnings_nodes.html',node_size=4,link_cmap=['grey','gold','firebrick'],link_width=1,add_to_node_popup=popupp,add_legend=True)
        
        #Create shp file
        p = shapefile.Writer('Warning_links.shp')
        p.field('name', 'C')
        p.field('Warn max flow', 'N')
        p.field('Warn mean flow', 'N')
        p.field('Warn pattern flow', 'N')
        for dot_name, dot in wn.links():
            p.line([[[dot.start_node.coordinates[0],dot.start_node.coordinates[1]],[dot.end_node.coordinates[0],dot.end_node.coordinates[1]]]])
            p.record(str(dot_name),1 if dot_name in hf else 0,1 if dot_name in hmf else 0,1 if dot_name in pf else 0)
        p.close

        w = shapefile.Writer('Warning_nodes.shp')
        w.field('name', 'C')
        w.field('Warn max pressure', 'N')
        w.field('Warn mean pressure', 'N')
        w.field('Warn pattern pressure', 'N')
        for dot_name, dot in wn.junctions():
            w.point(dot.coordinates[0],dot.coordinates[1])
            w.record(str(dot_name),1 if dot_name in hp else 0,1 if dot_name in hmp else 0,1 if dot_name in pp else 0)
        w.close
        
        print("Validation analisys successfully")