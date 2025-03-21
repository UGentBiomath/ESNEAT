# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 01:35:19 2022

@author: Cristian Camilo Gómez Cortés

### Description:
Resample the data file, clean the noise and adjust the timestep and lengthen list of elements.

### The input parameters are the following:
    - network: (str) - File path to the .inp file of the current network.
    - data_file: (str) - File path to the .csv raw data that you would like to resample.
    - flow_unit_factor: (int) - Conversion factor between actual flow data file and network units
    - pressure_unit_factor: (int) - Conversion factor between actual pressure data file and network units
    - Time_config: (str) Time step in which the information will be resampled. ['5T':5 minutes,'10T':10 minutes,'15T':15 minutes,'30T':30 minutes,'H':1 hour,'2H':2 hour]
### The output parameters are the following:
    - xlsx file with the resampled data
    - image file comparing the resampled data and the simulated by the current inp file
"""
# Load required packages
import wntr
import pandas as pnd
from matplotlib import pyplot
import matplotlib.pyplot as plt
import os

class data_resample:
    #Initialize the elements of the function
    def __init__(self, network, data_file, flow_unit_factor = 3600, pressure_unit_factor = 10.199773339984, time_config= '5T'):
        #Equivalent of each time configuration in seconds
        equivalence = {'5T':300,'10T':600,'15T':600,'30T':1800,'H':3600,'2H':7200}
        # Create results folder
        u = network.split("/")
        u[-1] = "Results"
        results = "/".join(u)
        if not os.path.exists(results):
            os.makedirs(results)
        #Load inp file
        wn = wntr.network.WaterNetworkModel(network)
        #Load data to resample
        df = pnd.read_csv(data_file)
        num = df._get_numeric_data()
        num[num < 1] = None
        # Clean and resample data based on time configuration
        df.index = pnd.to_datetime(df[df.keys()[0]])
        df = df.resample(time_config).mean()
        c = df.groupby(by = pnd.to_datetime(pnd.to_datetime(df.index).hour,unit='h')+pnd.to_timedelta(pnd.to_datetime(df.index).minute,unit='m')).agg('mean')
        #Plot all the data resampled for a 24 hours period
        c.plot()
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        pyplot.show()
        # Create a list that will contain the type of element information
        y=[]
        #Identify the type of element
        for p in list(c):    
            if str(p) in wn.pipe_name_list:
                y.append('T')
            if str(p) in wn.junction_name_list:
                y.append('N')
            if str(p) not in wn.junction_name_list and str(p) not in wn.pipe_name_list:
                print('Error in '+str(p))
                y.append('')
        #Create xlsx file and load the real data provided
        node = pnd.DataFrame(columns=y)
        with pnd.ExcelWriter('Results/Measurements_format_resampled.xlsx') as writer:  
            node.to_excel(writer,startrow=0)
            c.to_excel(writer,startrow=1)
        # Simulate the current flow and pressure patterns of the inp network to compare it later    
        wn.options.time.report_timestep = equivalence[time_config]
        wn.options.time.hydraulic_timestep = equivalence[time_config]
        sim2 = wntr.sim.EpanetSimulator(wn)
        resex = sim2.run_sim()
        #Plot the data and compare it with the current pressure and flow in each element
        o=0
        for j in list(c):
            #Plot measured node
            if str(j) in wn.junction_name_list:    
                pinr = pnd.DataFrame(columns = ["Real", "Calibrated"])
                pinr["Calibrated"]= resex.node['pressure'].loc[:(len(list(c[j]))-1)*equivalence[time_config],j]
                pinr["Real"]= list(c[j])
                #Pressure in bar to m
                pinr["Real"]=pinr["Real"]*pressure_unit_factor
                pinr.index /= 3600 # convert time to hours
                lab = pinr.plot(kind="line",figsize=(8,4),grid = True)
                lab.set_xlabel("Time [hr]",fontsize=10)
                lab.set_ylabel("Pressure [m]",fontsize=10)
                lab.set_title("Pressure_node_"+str(j))
                plt.savefig('Results/'+'Resampled_'+'Compare'+"_Element"+str(o)+"_Pressure.png",dpi=800)
                plt.close()
                o+=1
            #Plot measured link
            if str(j) in wn.pipe_name_list:
                
                qinr = pnd.DataFrame(columns = ["Real", "Calibrated"])
                qinr["Calibrated"]= resex.link['flowrate'].loc[:(len(list(c[j]))-1)*equivalence[time_config],j]
                qinr["Real"] = list(c[j])
                qinr["Calibrated"]= qinr["Calibrated"]*flow_unit_factor
                qinr.index /= 3600 # convert time to hours
                lab = qinr.plot(kind="line",figsize=(8,4),grid = True)
                lab.set_xlabel("Time [hr]",fontsize=10)
                lab.set_ylabel("Flow [CMH]",fontsize=10)
                lab.set_title("Flow_pipe_"+str(j))
                plt.savefig('Results/'+'Resampled_'+'Compare'+"_Element"+str(o)+"_Flow.png",dpi=800)
                plt.close()
                o+=1
        print("data resample successfully")
         