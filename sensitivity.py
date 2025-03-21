# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 09:11:15 2020

@author: Cristian Camilo Gómez Cortés

### Description:
Perform a sensitivity analysis of the elements to detect the most sensitive measurement points.
For the current analysis the following parameters were considered: demand, leaks, roughness

### The input parameters are the following:
    - network: (str) - File path to the .inp file of the current network.
    - net_name: (str) - Name of the current network just for saving purposes
### The output parameters are the following:
    - xlsx file with the ranked most sensitive elements by parameter
    - image file comparing the sensitivity across the network by parameter
    - xlsx file with the ranked most sensitive nodes and links
    - image file comparing the sensitivity across the network for each node and link
"""
# Load required packages
import wntr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pnd
import os

class sensitivity:
    def __init__(self, network, net_name):
        # Create sensitivity folder
        u = network.split("/")
        u[-1] = "Sensitivity"
        sensitiv = "/".join(u)
        if not os.path.exists(sensitiv):
            os.makedirs(sensitiv)
        #Load network
        wn = wntr.network.WaterNetworkModel(network)
        #Set simulation as static
        wn.options.time.duration = 0
        #define delta to change each parameter
        dt_dem = 0.5
        dt_roug = 0.5
        dt_leak = 1
        #Read the number of pipes and nodes
        numnod = len(wn.junction_name_list)
        numtub = len(wn.pipe_name_list)
        #Create a matrix for each parameter
        #Pressure and flow for roughness matrices
        pre_rough = np.zeros((numnod,numtub))
        flow_rough = np.zeros((numtub,numtub))
        #Pressure and flow for demand matrices
        pre_dem = np.zeros((numnod,numnod))
        flow_dem = np.zeros((numtub,numnod))
        #Pressure and flow for leak matrices
        pre_leak = np.zeros((numnod,numnod))
        flow_leak = np.zeros((numtub,numnod))
              
        #Remove demand patterns in order to set the simulation as static
        for node_name, node in wn.junctions():
            if wn.get_node(node_name).demand_timeseries_list[0].pattern_name != "None":
                wn.get_node(node_name).demand_timeseries_list[0].pattern_name = "'None'"
        # Perform sensitivity analysis for the roughness
        g = 0
        for pipe_name, pipe in wn.pipes():
            if pipe.roughness > 0:
                #Simulation without changes
                t = pipe.roughness
                sim = wntr.sim.EpanetSimulator(wn)
                results = sim.run_sim()
                p_temp = results.node['pressure'].values.tolist()
                q_temp = results.link['flowrate'].values.tolist()
                #Simulation with changes
                pipe.roughness = pipe.roughness*(1+dt_roug)
                sim2 = wntr.sim.EpanetSimulator(wn)
                results2 = sim2.run_sim()
                p_temp2 = results2.node['pressure'].values.tolist()
                q_temp2 = results2.link['flowrate'].values.tolist()
                #restore previous value
                pipe.roughness = t
                #Comparation between the 2 simulations
                for x in range(numnod):
                    pre_rough[x][g] = abs(p_temp[0][x]-p_temp2[0][x])/(t*dt_roug)
                for y in range(numtub):
                    flow_rough[y][g] = abs(q_temp[0][y]-q_temp2[0][y])/(t*dt_roug)
            g += 1
            
        
        #Perform sensitivity analysis for the demand
        g = 0
        for node_name, node in wn.junctions():
            if node.demand_timeseries_list[0].base_value > 0:
                #Simulation without changes
                t = node.demand_timeseries_list[0].base_value
                sim = wntr.sim.EpanetSimulator(wn)
                results = sim.run_sim()
                p_temp = results.node['pressure'].values.tolist()
                q_temp = results.link['flowrate'].values.tolist()
                #Simulation with changes
                node.demand_timeseries_list[0].base_value = node.demand_timeseries_list[0].base_value*(1+dt_dem)
                sim2 = wntr.sim.EpanetSimulator(wn)
                results2 = sim2.run_sim()
                p_temp2 = results2.node['pressure'].values.tolist()
                q_temp2 = results2.link['flowrate'].values.tolist()
                #restore previous value
                node.demand_timeseries_list[0].base_value = t
                #Comparation between the 2 simulations
                for x in range(numnod):
                        pre_dem[x][g] = abs(p_temp[0][x]-p_temp2[0][x])/(t*dt_dem)
                for y in range(numtub):
                        flow_dem[y][g] = abs(q_temp[0][y]-q_temp2[0][y])/(t*dt_dem)
            g+=1
            
        #Perform sensitivity analysis for the leaks
        g = 0
        for node_name, node in wn.junctions():
            #Simulation without changes
            t = 0
            node.emitter_coefficient = 0
            sim = wntr.sim.EpanetSimulator(wn)
            results = sim.run_sim()
            p_temp = results.node['pressure'].values.tolist()
            q_temp = results.link['flowrate'].values.tolist()
            #Simulation with changes
            node.emitter_coefficient = dt_leak
            sim2 = wntr.sim.EpanetSimulator(wn)
            results2 = sim2.run_sim()
            p_temp2 = results2.node['pressure'].values.tolist()
            q_temp2 = results2.link['flowrate'].values.tolist()
            #restore previous value
            node.demand_timeseries_list[0].base_value = t
            #Comparation between the 2 simulations
            for x in range(numnod):
                pre_leak[x][g] = abs(p_temp[0][x]-p_temp2[0][x])/(dt_leak)
            for y in range(numtub):
                flow_leak[y][g] = abs(q_temp[0][y]-q_temp2[0][y])/(dt_leak)
            g+=1
            
            
        #Create a list to save al the sensitivity information by parameter
        S_pre_rough = np.zeros((numnod,2))
        S_flow_rough = np.zeros((numtub,2))
        S_pre_dem = np.zeros((numnod,2))
        S_flow_dem = np.zeros((numtub,2))
        S_pre_leak = np.zeros((numnod,2))
        S_flow_leak = np.zeros((numtub,2))
        
        # Create 3 networks to store the sensitivity information just with plot purposes
        wn1 = wntr.network.WaterNetworkModel(network)
        wn2 = wntr.network.WaterNetworkModel(network)
        wn3 = wntr.network.WaterNetworkModel(network)
        # Save information with plot purposes
        # Save roughness information
        for z in range(numnod):
            S_pre_rough[z][0] = z
            S_pre_rough[z][1] = round(np.sum(pre_rough[z]),2)
            wn1.get_node(wn.node_name_list[z]).elevation = np.sum(pre_rough[z])
        for w in range(numtub):
            S_flow_rough[w][0] = w
            S_flow_rough[w][1] = round(np.sum(flow_rough[w]),2)
            wn1.get_link(wn.pipe_name_list[w]).length = np.sum(flow_rough[w])
        #Plot roughness sensitivity
        wntr.graphics.plot_network(wn1,title="Sensitivity nodes by roughness", node_attribute='elevation',node_colorbar_label='Σ(Δp/Δks)', node_size=10)
        #plt.savefig('Sensitivity/'+net_name+"_Sensitivity nodes byroughness",dpi=300,orientation="landscape")
        wntr.graphics.plot_network(wn1,title="Sensitivity pipes by roughness", link_attribute='length',link_colorbar_label='Σ(ΔQ/Δks)', node_size=0.1)
        #plt.savefig('Sensitivity/'+net_name+"_Sensitivity pipes by roughness",dpi=300,orientation="landscape")
        
        # Save demand information
        for m in range(numnod):
            S_pre_dem[m][0] = m
            S_pre_dem[m][1] = round(np.sum(pre_dem[m]),2)
            wn2.get_node(wn.node_name_list[m]).elevation = np.sum(pre_dem[m])
        for n in range(numtub):
            S_flow_dem[n][0] = n
            S_flow_dem[n][1] = round(np.sum(flow_dem[n]),2)
            wn2.get_link(wn.pipe_name_list[n]).length = np.sum(flow_dem[n])
        #Plot demand sensitivity
        wntr.graphics.plot_network(wn2,title="Sensitivity nodes by demand", node_attribute='elevation',node_colorbar_label='Σ(Δp/Δdem)', node_size=10)
        # plt.savefig('Sensitivity/'+net_name+"_Sensitivity nodes by demand",dpi=300,orientation="landscape")
        wntr.graphics.plot_network(wn2,title="Sensitivity pipes by demand", link_attribute="length",link_colorbar_label='Σ(ΔQ/Δdem)', node_size=0.1)
        # plt.savefig('Sensitivity/'+net_name+"_Sensitivity pipes by demand",dpi=300,orientation="landscape")
        
        # Save leak information
        for m in range(numnod):
            S_pre_leak[m][0] = m
            S_pre_leak[m][1] = round(np.sum(pre_leak[m]),2)
            wn2.get_node(wn.node_name_list[m]).elevation = np.sum(pre_leak[m])
        for n in range(numtub):
            S_flow_leak[n][0] = n
            S_flow_leak[n][1] = round(np.sum(flow_leak[n]),2)
            wn2.get_link(wn.pipe_name_list[n]).length = np.sum(flow_leak[n])
        #Plot leak sensitivity
        wntr.graphics.plot_network(wn2,title="Sensitivity nodes by leak", node_attribute='elevation',node_colorbar_label='Σ(Δp/Δleak)', node_size=10)
        # plt.savefig('Sensitivity/'+net_name+"_Sensitivity nodes by leak",dpi=300,orientation="landscape")
        wntr.graphics.plot_network(wn2,title="Sensitivity pipes by leak", link_attribute="length",link_colorbar_label='Σ(ΔQ/Δleak)', node_size=0.1)
        # plt.savefig('Sensitivity/'+net_name+"_Sensitivity pipes by leak",dpi=300,orientation="landscape")
        
        #Save results in xlsx format
        #Save roughness info in xlsx
        EnR = pnd.DataFrame(data=np.flip(S_pre_rough[np.argsort(S_pre_rough[:, 1])],0), columns=["Node","S_pR"])
        EtR = pnd.DataFrame(np.flip(S_flow_rough[np.argsort(S_flow_rough[:, 1])],0), columns=["Pipe","S_qR"])
        #Save demand info in xlsx
        EnD = pnd.DataFrame(np.flip(S_pre_dem[np.argsort(S_pre_dem[:, 1])],0), columns=["Node","S_pD"])
        EtD = pnd.DataFrame(np.flip(S_flow_dem[np.argsort(S_flow_dem[:, 1])],0), columns=["Pipe","S_qD"])
        #Save leak info in xlsx
        EnL = pnd.DataFrame(np.flip(S_pre_leak[np.argsort(S_pre_leak[:, 1])],0), columns=["Node","S_pL"])
        EtL = pnd.DataFrame(np.flip(S_flow_leak[np.argsort(S_flow_leak[:, 1])],0), columns=["Pipe","S_qL"])
        
        #Include index of element in each list of results and strucuture the printing format
        v = 0
        for z in EnR["Node"]:
            EnR.loc[v, "Node"] = wn.node_name_list[int(z)]
            v += 1
        v = 0
        for w in EtR["Pipe"]:
            EtR.loc[v, "Pipe"] = wn.pipe_name_list[int(w)]
            v += 1
        v = 0
        for m in EnD["Node"]:
            EnD.loc[v,"Node"] = wn.node_name_list[int(m)]
            v += 1
        v = 0
        for n in EtD["Pipe"]:
            EtD.loc[v, "Pipe"] = wn.pipe_name_list[int(n)]
            v += 1
        v = 0
        for m in EnL["Node"]:
            EnL.loc[v,"Node"] = wn.node_name_list[int(m)]
            v += 1
        v = 0
        for n in EtL["Pipe"]:
            EtL.loc[v, "Pipe"] = wn.pipe_name_list[int(n)]
            v += 1

        #Create the sensitivity file including each parameter
        with pnd.ExcelWriter('Sensitivity/'+'Sensitivity_parameters'+net_name+'.xlsx') as writer:  
            EnR.to_excel(writer)
            EtR.to_excel(writer,startcol=3)
            EnD.to_excel(writer,startcol=6)
            EtD.to_excel(writer,startcol=9)
            EnL.to_excel(writer,startcol=12)
            EtL.to_excel(writer,startcol=15)
        
        # Normalize the results to compare all parameters as equal
        EnR["S"] = round(EnR["S_pR"].div(EnR["S_pR"].max()),2)
        EtR["S"] = round(EtR["S_qR"].div(EtR["S_qR"].max()),2)
        EnD["S"] = round(EnD["S_pD"].div(EnD["S_pD"].max()),2)
        EtD["S"] = round(EtD["S_qD"].div(EtD["S_qD"].max()),2)
        EnL["S"] = round(EnL["S_pL"].div(EnL["S_pL"].max()),2)
        EtL["S"] = round(EtL["S_qL"].div(EtL["S_qL"].max()),2)
        
        #Delete the extra column previously created
        del EnR["S_pR"]
        del EtR["S_qR"]
        del EnD["S_pD"]
        del EtD["S_qD"]
        del EnL["S_pL"]
        del EtL["S_qL"]
        
        #Create a list for nodes and links that compile the results of all the parameters
        #Create nodes list
        N_results = EnR
        for n in range(numnod):
            N_results.loc[len(N_results.index)] = EnD.loc[n]
        for n in range(numnod):
            N_results.loc[len(N_results.index)] = EnL.loc[n]
        #Create pipe list
        P_results = EtR
        for t in range(numtub):
            P_results.loc[len(P_results.index)] = EtD.loc[t]
        for t in range(numtub):
            P_results.loc[len(P_results.index)] = EtL.loc[t]
        
        #Sum all the results by element to obtain the sensitivity of each element
        Node_sens = N_results.groupby(["Node"]).sum()
        Pipe_sens = P_results.groupby(["Pipe"]).sum()
        
        #Create and the results in a network just for plot purposes
        wn4 = wntr.network.WaterNetworkModel(network)
        for m in Node_sens.index:
            wn4.get_node(m).elevation = Node_sens.loc[m,"S"]
        for n in Pipe_sens.index:
            wn4.get_link(n).length = Pipe_sens.loc[n,"S"]
        #Plot sensitivity by element
        wntr.graphics.plot_network(wn4,title="Node Sensitivity", node_attribute='elevation',node_colorbar_label='Sensitivity\nNormalized\nNodes',node_size=10)
        # plt.savefig('Sensitivity/'+"Node_Sensitivity_"+net_name,dpi=300,orientation="landscape")
        wntr.graphics.plot_network(wn4,title="Pipe Sensitivity", link_attribute="length",link_colorbar_label='Sensitivity\nNormalized\nPipes',node_size=0.1)
        # plt.savefig('Sensitivity/'+"Pipe_Sensitivity_"+net_name,dpi=300,orientation="landscape")
        #Save results in a xlsx file
        #Create format and sort
        Node_sens = Node_sens.sort_values("S",ascending=False)
        Pipe_sens = Pipe_sens.sort_values("S",ascending=False)
        #
        with pnd.ExcelWriter('Sensitivity/'+'Sensitivity_elements'+net_name+'.xlsx') as writer:
            Node_sens.to_excel(writer)
            Pipe_sens.to_excel(writer, startcol=2)
        
        print("Sensitivity analysis successfully")
red = "C:/Users/cgomezco/OneDrive - UGent/Desktop/Articles/AutomaticCalibration/WaterResearch/fossolo.inp"
sens = sensitivity(network=red, net_name="Fossolo")