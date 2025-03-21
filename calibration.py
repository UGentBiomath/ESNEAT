# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 10:07:07 2020

@author: Cristian Camilo Gómez Cortés

### Description:
Perform the automatic calibration process of a EPANET network gathering all the expert knowledge and implementing the NEAT algorithm
For the current analysis the following parameters were considered: demand, leaks, roughness, minor loss, valve status, valve coefficient, demand patterns

### The input parameters are the following:
    - network: (str) - File path to the .inp file of the current network.
    - net_name: (str) - Name of the current network just for saving purposes
    - Real_data: (list) - List with all the measured information saved loaded from the measurement file
    - Elements: (DataFrame) - Dataframe with the id and type of each measured element
    - calib_roughness: (bool) - Boolean parameter to enable the calibration of the roughness
    - calib_minorloss: (bool) - Boolean parameter to enable the calibration of the minorloss
    - calib_demand: (bool) - Boolean parameter to enable the calibration of the demand
    - calib_leak: (bool) - Boolean parameter to enable the calibration of the leak
    - calib_valves_coeff: (bool) - Boolean parameter to enable the calibration of the valves coefficient
    - calib_valve_status: (bool) - Boolean parameter to enable the calibration of the valve status
    - calib_demand_pattern: (bool) - Boolean parameter to enable the calibration of the demand pattern
    - generations: (int) - Number of generations executed in each mass and energy calibration process
    - population: (int) - Number of individual per generation executed in each mass and energy calibration process
    - loops: (int) - Number of mass and energy calibration process (Recomended in case of convergence problems). Recomended: 1
    - hidden_layers: (int) - Number of initial hidden layers of each ANN created for the algorithm. Recomended: 0
    - disc_mass: (int) - Number of discretization steps on the mass ranges parameters. Recomended: 10
    - disc_energy: (int) - Number of discretization steps on the energy ranges parameters. Recomended: 4
    - percent_closure_valves: (int) - Percent of closed valves at the same time between 0 and 100. Recomended: 20
    - rang_valve_coeff: (list) - Search range of the valve coefficient (valve setting). Recomended: [0.2,10]
    - rang_minorloss: (list) - Search range of the minor loss coefficient of each pipe. Recomended:[0.1,2]
    - rang_leak: (list) - Search range of the leak coefficient of each node (emitter coefficient). Recomended: [0,1]
    - rang_demand: (list) - Search range of the base demand of each node. Recomended: [0.5,1.5]
### The output parameters are the following:
    - inp network with the final calibrated parameters
    - image file comparing flow and pressure of each measured element after calibration
    - image file comparing Pressure_Head_Demand for each element
    - image file comparing Flow_Head_Pressure for each element
    - image file comparing before and after calibration percetange change of the parameters
    - xlsx file comparing before and after calibration values of the parameters
    - xlsx file saving signals of flow and pressure
"""
# Load required packages
import wntr
import neat
from spotpy.objectivefunctions import nashsutcliffe as nse
from statistics import mean
import pandas as pnd
import shutil
import os
import matplotlib.pyplot as plt
from text_neat import copy_text_neat
from dijkstra import dijkstra

class calibration:
    def __init__(self, network, net_name, Real_data, Elements,
                 calib_roughness, calib_minorloss, calib_demand, calib_leak, calib_valves_coeff, calib_valve_status, calib_demand_pattern, 
                 generations, population, loops= 1, hidden_layers = 0, disc_mass = 10, disc_energy = 4,
                 percent_closure_valves = 20, rang_valve_coeff = [0.2,10.0], rang_minorloss = [0.1,2.0], rang_leak = [0,1], rang_demand = [0.5,1.5]):
        # Create inp network to store results
        resul_ex = net_name+"_Calibrated.inp"
        shutil.copy(network, resul_ex)
        
        # Create results folder
        u = network.split("/")
        u[-1] = "Results"
        results = "/".join(u)
        if not os.path.exists(results):
            os.makedirs(results)
            
        #Load the network in a wntr variable
        wn_ex = wntr.network.WaterNetworkModel(resul_ex)
            
        #Units & discretization on the mass and energy parameters
        equivalence = {'LPS':1000,'LPM':60000,'CMH':3600,'CMD':86400}
        unit_factor = equivalence[wn_ex.options.hydraulic.inpfile_units]
        discretization_mass = [1.0]*disc_mass
        discretization_energy = [1.0]*disc_energy
        
        # Define vectors for each calibrated parameter to save the final results
        if calib_roughness:
            final_r = []
        if calib_minorloss:
            final_ml = []
        if calib_demand:
            final_dem = []
        if calib_leak:
            final_leak = []
        if calib_valves_coeff:
            final_valcoeff = []
        if calib_valve_status:
            final_valstat = []
            
        # Define roughness ranges
        if calib_roughness or calib_minorloss:
            # Roughness coefficients
            coeff_roug = {"DCI":{"C1":10.533,"C2":1.523,"C3":1.084},
                        "GCI":{"C1":22.760,"C2":-4.180,"C3":1.2},
                        "PEH":{"C1":0.302,"C2":0.274,"C3":0.993},
                        "PVC":{"C1":0.302,"C2":0.274,"C3":0.993},
                        "STE":{"C1":74.681,"C2":-5.436,"C3":1.311},
                        "FCE":{"C1":7.627,"C2":-0.171,"C3":1.211},
                        "CON":{"C1":7.627,"C2":-0.171,"C3":1.211}}
            # Initial roughness for each material
            ini_roug = {"DCI":0.26,
                       "GCI":0.8,
                       "PEH":0.007,
                       "PVC":0.0015,
                       "STE":0.045,
                       "FCE":0.3,
                       "CON":0.6}
            
            #Round approximation for each material
            round_roug = {"DCI":2,
                       "GCI":1,
                       "PEH":3,
                       "PVC":4,
                       "STE":3,
                       "FCE":1,
                       "CON":1}
            range_roug = {}
            materials = coeff_roug.keys()
            
            #Calculate the roughness of each material based on the equation described in the manual.
            for material in materials:
                range_roug[material] = {}
                for i in range(0,50,1):
                    if i == 0:
                        range_roug[material][str(i)] = round(ini_roug[material],round_roug[material])
                    else:
                        range_roug[material][str(i)] = round(ini_roug[material]*((coeff_roug[material]["C1"]*((i/50)**2))+(coeff_roug[material]["C2"]*(i/50))+(coeff_roug[material]["C3"])),round_roug[material])
            # Group the roughness by material
            material_pipes = []
            for pipe_name, pipe in wn_ex.pipes():
                material_pipes.append(pipe.tag[0:6])
            # Define which unique materials are contained in the network
            real_materials = []
            [real_materials.append(matter) for matter in material_pipes if matter not in real_materials]
        #Save initial demands for each node
        if calib_demand:
            init_demands = pnd.DataFrame(index = wn_ex.junction_name_list, columns = ["Demand"])
            for node_name, node in wn_ex.junctions():
                init_demands["Demand"][node_name] = node.demand_timeseries_list[0].base_value
        # Define which valves do not cause an error due to their independent closure
        if calib_valve_status or calib_valves_coeff:
            # Load valve test network in a wntr parameter
            wn_valve = wntr.network.WaterNetworkModel(resul_ex)
            sim_valve = wntr.sim.EpanetSimulator(wn_valve)
            # Define time duration as 0
            wn_valve.options.time.duration = 0
            # Open all the valves
            for valve_name, valve in wn_valve.valves():
                valve.initial_status = wntr.network.base.LinkStatus["Active"]
            # Create vectors to save which valves are suitable to manipulate in the calibration process
            suitable_valves = []
            no_suitable_valves = []
            # Test each valve independently 
            for valve_name, valve in wn_valve.valves():
                # Select TCV valves
                if valve.valve_type == "TCV":
                    # Change status to closed
                    valve.initial_status = wntr.network.base.LinkStatus["Closed"]
                    # Verify the convergence criteria for the simulated network
                    try:
                        results = sim_valve.run_sim(convergence_error=True)
                    # Save valve as no suitable in case of error
                    except:
                        print("Handled error")
                        no_suitable_valves.append(valve_name)
                    # Save valve as suitable in case of no error
                    else:
                        local_dir = os.path.dirname(__file__)
                        report = os.path.join(local_dir,"temp.rpt")
                        file1 = open(report, "r")
                        readfile = file1.read()
                        if 'disconnected' in readfile: 
                            print("Disconnected")
                            no_suitable_valves.append(valve_name)
                        else:
                            suitable_valves.append(valve_name)
                    # Open the tested valve
                    finally:
                        valve.initial_status = wntr.network.base.LinkStatus["Active"]
            # Load suitable open and closed valves for plot purposes
            # Define all elevation parameters as none for plot purposes
            for z in wn_valve.node_name_list:
                wn_valve.get_node(z).elevation = None
            # Load valve information
            for valve_name, valve in wn_valve.valves():
                if valve_name in suitable_valves:
                    wn_valve.get_node(valve.start_node).elevation = 0
                elif valve_name in no_suitable_valves:
                    wn_valve.get_node(valve.start_node).elevation = 100
            # Plot the valve information
            wntr.graphics.plot_network(wn_valve,title="Suitable open valves", node_attribute='elevation',node_colorbar_label='No Suitable(Red)\nSuitable(Blue)',node_cmap='jet', node_size=10)
            plt.savefig("Suitable valves",dpi=1000,orientation="landscape")  
        # Comparison before calibration
        # Simulate the network before calibration
        sim = wntr.sim.EpanetSimulator(wn_ex)
        results = sim.run_sim()
        # Compare each element with the input data and sum the fitness function
        d = 0
        div = len(Real_data)/len(Elements["Type"])
        fitness = []
        for i, j in Elements.iterrows():
            #Compare links "T" (tubes)
            if j["Type"] == "T":
                caudal1 = results.link['flowrate'].loc[:,j["ID"]].tolist()
                caudal1 = [caudal1 * unit_factor for caudal1 in caudal1]
                fitness.append(nse(Real_data[int(d*div):int(((d+1)*div))],caudal1))
                d += 1
            #Compare junctions "N" (nodes)
            if j["Type"] == "N":
                presion = results.node['pressure'].loc[:,j["ID"]].tolist()
                fitness.append(nse(Real_data[int(d*div):int(((d+1)*div))],presion))
                d += 1
        #Print initial fitness before calibration
        print(sum(fitness))
        # Implement the pattern calibration methodology based on the nearest measured element
        if calib_demand_pattern:
            # Initialize the measured elements and element types
            element = []
            elemtype = []
            for n, ty in Elements.iterrows():
                element.append(ty["ID"])
                elemtype.append(ty["Type"])
            # Find the nearest measured element to each node
            nearest_elem =dijkstra(network, element, elemtype).run()
            # Calculate the normalized pattern to each measured element and save it in norm_pattern
            norm_pattern = []
            # Normalize each pattern of each measured element
            d = 0
            div = len(Real_data)/len(Elements["Type"])
            for i, j in Elements.iterrows():
                # Normalize each flow pattern of each measured element
                if j["Type"] == "T":
                    norm_pattern.append([x/(mean(Real_data[int(d*div):int(((d+1)*div))])) for x in Real_data[int(d*div):int(((d+1)*div))]])
                    d += 1
                # Normalize and inverse each pressure pattern of each measured element
                if j["Type"] == "N":
                    norm_pattern.append([1+(1-(x/(mean(Real_data[int(d*div):int(((d+1)*div))])))) for x in Real_data[int(d*div):int(((d+1)*div))]])
                    d += 1
            # Delete previous patterns
            # Delete the patterns in each node
            for node_name, node in wn_ex.junctions():
                del node.demand_timeseries_list[:]
            # Save new network without patterns
            wn_ex.write_inpfile(wn_ex.name)
            # Reload network without previous patterns
            wn_ex = wntr.network.WaterNetworkModel(resul_ex)
            # Exclude reservoir patterns
            exclude_patterns = []
            for name, x in wn_ex.reservoirs():
                exclude_patterns.append(wn_ex.get_node(name).head_pattern_name)
            # Delete patterns from the pattern registry
            for k in wn_ex.pattern_name_list:
                if k not in exclude_patterns:
                    wn_ex.remove_pattern(str(k))
            # Create new patterns based on the measured elements
            for i in range(len(norm_pattern)):
                wn_ex.add_pattern('pat'+str(i), norm_pattern[i])
            # Assign each pattern to each node based on the nearest measured element
            for t in range(len(wn_ex.junction_name_list)):
                wn_ex.get_node(wn_ex.junction_name_list[t]).demand_timeseries_list[0].pattern_name = 'pat'+str(nearest_elem[t])
            # Rewrite the network
            wn_ex.write_inpfile(wn_ex.name)  
        # Initialize the number of loops for the mass-energy calibration process
        loop = 1
        while loop <= loops:
            ### Start NEAT Mass analysis: Calibrate leaks, demand and valve status
            if calib_demand or calib_leak or calib_valve_status:
                # Load wntr network to perform mass calibration
                wn_ex = wntr.network.WaterNetworkModel(resul_ex)
                # Define the number of steps taken for each calibration range in the mass calibration process
                number_steps_mass = len(discretization_mass)
                # Define the number of outputs based on the calibrated parameters
                num_outputs_m = 0
                if calib_demand:
                    num_outputs_m += wn_ex.num_junctions
                if calib_leak:
                    num_outputs_m += wn_ex.num_junctions
                if calib_valve_status:
                    num_outputs_m += len(suitable_valves)
                # Create a copy of the parameter file for the neat methodology in the mass calibration process
                copy_text_neat("m",population,hidden_layers,number_steps_mass,num_outputs_m)
                # Define objective function for the mass calibration process
                def eval_genomes_m(genomes, config):
                    #Test the performance of each individual
                    for genome_id, genome in genomes:
                        # Create ANN for each individual
                        net = neat.nn.FeedForwardNetwork.create(genome, config)
                        # Calculate the outputs
                        output = net.activate(discretization_mass)
                        # Create count variable m to define the number of element
                        m = 0
                        # Load the network file
                        wn_ex = wntr.network.WaterNetworkModel(resul_ex)
                        # Save the demand outputs
                        if calib_demand:
                            for node_name, node in wn_ex.junctions():
                                node.demand_timeseries_list[0].base_value = (((output[m]+1)*(init_demands["Demand"][node_name]*(rang_demand[1]-rang_demand[0])))/2)+init_demands["Demand"][node_name]*rang_demand[0]
                                m += 1
                        # Save the leaks outputs
                        if calib_leak:
                            for node_name, node in wn_ex.junctions():
                                if output[m]>0.9:
                                    node.emitter_coefficient = round(((((output[m]+1)*((rang_leak[1]-rang_leak[0])/unit_factor))/2)+(rang_leak[0]/unit_factor)),4)
                                m += 1
                        # Save the valve status outputs
                        if calib_valve_status:
                            for valve_name, valve in wn_ex.valves():
                                if (valve_name in suitable_valves) and (valve.valve_type == "TCV"):                                  
                                    if output[m] >= int((100-percent_closure_valves)/100):
                                        valve.initial_status = wntr.network.base.LinkStatus["Closed"]
                                    else:
                                        valve.initial_status = wntr.network.base.LinkStatus["Active"]
                                    m += 1  
                        #Simulate the network with the new outputs
                        sim = wntr.sim.EpanetSimulator(wn_ex)
                        # Verify the convergence criteria for the simulated network
                        try:
                            results = sim.run_sim(convergence_error=True)
                        #If a error is raised the fitness is increases to discard the individual in the next generation
                        except:
                            print("Handled error")
                            genome.fitness = -99999
                        # Else calculate the fitness based on the procedure 
                        else:
                            # Verify if any element is disconnected
                            local_dir = os.path.dirname(__file__)
                            report = os.path.join(local_dir,"temp.rpt")
                            # opening a text file
                            file1 = open(report, "r")
                            # read file content
                            readfile = file1.read()
                            if 'disconnected' in readfile: 
                                print("Disconnected")
                                genome.fitness = -99999
                            # Compare each element with the input data and sum the fitness function
                            else: 
                                d = 0
                                div = len(Real_data)/len(Elements["Type"])
                                fitness = []
                                for i, j in Elements.iterrows():
                                    #Compare links "T" (tubes)
                                    if j["Type"] == "T":
                                        caudal1 = results.link['flowrate'].loc[:,j["ID"]].tolist()
                                        caudal1 = [caudal1 * unit_factor for caudal1 in caudal1]
                                        fitness.append(nse(Real_data[int(d*div):int(((d+1)*div))],caudal1))
                                        d +=1
                                    #Compare junctions "N" (nodes)
                                    if j["Type"] == "N":
                                        presion = results.node['pressure'].loc[:,j["ID"]].tolist()
                                        fitness.append(nse(Real_data[int(d*div):int(((d+1)*div))],presion))
                                        d +=1
                                # Sum fitness function
                                genome.fitness = sum(fitness)
                            file1.close()
                # Create a function to execute the mass calibration
                def run_m(config_file):
                    # Define the NEAT configuration
                    config = neat.Config(neat.DefaultGenome, neat.DefaultReproduction,
                                         neat.DefaultSpeciesSet, neat.DefaultStagnation,
                                         config_file)
                    # Load population configuration
                    p_m = neat.Population(config)
                    # Create temporal checkpoint path to reload best fitness individual to the next process
                    t_m = network.split("/")
                    # Temporal checkpoint based on loop and number of generations
                    t_m[-1] = "checkpoint-m-"+str(int(int(generations)*(loop-1))-(1*(loop-1)))
                    temp_check_m = "/".join(t_m)
                    # Load previous checkpoint if exist
                    if os.path.exists(temp_check_m):
                        p_m = neat.Checkpointer.restore_checkpoint("checkpoint-m-"+str(int(int(generations)*(loop-1))-(1*(loop-1))))
                        print("Restored checkpoint "+t_m[-1])
                    # Add reporter
                    p_m.add_reporter(neat.StdOutReporter(True))
                    stats = neat.StatisticsReporter()
                    p_m.add_reporter(stats)
                    # Create checkpoint based on loop and number of generations
                    p_m.add_reporter(neat.Checkpointer(generation_interval=int(int(generations)*loop-(1*(loop-1))),time_interval_seconds=None,filename_prefix="checkpoint-m-"))
                    # Execute fitness function of each individual and save best fitness individual 
                    winner = p_m.run(eval_genomes_m, int(generations))
                    # Save best fitness ANN
                    winner_net = neat.nn.FeedForwardNetwork.create(winner, config)
                    # Save outputs
                    output = winner_net.activate(discretization_mass)
                    # Load network to save best fitness individual parameters
                    wn_ex = wntr.network.WaterNetworkModel(resul_ex)
                    # Create counter based on the number of parameter
                    n = 0
                    # Save best demand information
                    if calib_demand:
                        for node_name, node in wn_ex.junctions():
                            node.demand_timeseries_list[0].base_value = (((output[n]+1)*(init_demands["Demand"][node_name]*(rang_demand[1]-rang_demand[0])))/2)+init_demands["Demand"][node_name]*rang_demand[0]
                            n += 1
                            final_dem.append(node.demand_timeseries_list[0].base_value)
                    # Save best leak information
                    if calib_leak:
                        for node_name, node in wn_ex.junctions():
                            if output[n]>0.9:
                                node.emitter_coefficient = round(((((output[n]+1)*((rang_leak[1]-rang_leak[0])/unit_factor))/2)+(rang_leak[0]/unit_factor)),4)
                            n += 1
                            final_leak.append(node.emitter_coefficient)
                    # Save best valve status information
                    if calib_valve_status:
                        for valve_name, valve in wn_ex.valves():
                            if (valve_name in suitable_valves) and (valve.valve_type == "TCV"):                                  
                                if output[n] >= int((100-percent_closure_valves)/100):
                                    valve.initial_status = wntr.network.base.LinkStatus["Closed"]
                                else:
                                    valve.initial_status = wntr.network.base.LinkStatus["Active"]
                                n += 1
                                final_valstat.append(valve.initial_status)
                    # rewrite best fitness information in the network file
                    wn_ex.write_inpfile(wn_ex.name)       
                # Execute run mass function
                local_dir = os.path.dirname(__file__)
                config_path_m = os.path.join(local_dir, 'temp_neat_m.txt')
                run_m(config_path_m)
            ### End NEAT Mass analysis
            ### Start NEAT Energy analysis: Calibrate roughness, minor losses and valve coefficients
            if calib_roughness or calib_minorloss or calib_valves_coeff:
                # Load wntr network to perform energy calibration
                wn_ex = wntr.network.WaterNetworkModel(resul_ex)
                # Define the number of steps taken for each calibration range in the energy calibration process
                number_steps_energy = len(discretization_energy)
                # Define the number of outputs based on the calibrated parameters
                num_outputs_e = 0
                if calib_roughness:
                    num_outputs_e += len(real_materials)
                if calib_minorloss:
                    num_outputs_e += len(real_materials)
                if calib_valves_coeff:
                    num_outputs_e += len(suitable_valves)
                # Create a copy of the parameter file for the neat methodology in the energy calibration process
                copy_text_neat("e",population,hidden_layers,number_steps_energy,num_outputs_e)
                # Define objective function for the energy calibration process
                def eval_genomes_e(genomes, config):
                    #Test the performance of each individual
                    for genome_id, genome in genomes:
                        # Create ANN for each individual
                        net = neat.nn.FeedForwardNetwork.create(genome, config)
                        # Calculate the outputs
                        output = net.activate(discretization_energy)
                        # Create count variable m to define the number of output
                        m = 0
                        # Load the network file
                        wn_ex = wntr.network.WaterNetworkModel(resul_ex)
                        # Save the roughness outputs
                        if calib_roughness:
                            for pipe_name, pipe in wn_ex.pipes():
                                if int(pipe.tag[4:6]) < 10:
                                    pipe.roughness = round((((output[m+int(real_materials.index(pipe.tag[0:6]))]+1)*(range_roug[pipe.tag[0:3]][str(20)]-range_roug[pipe.tag[0:3]][str(0)]))/2)+range_roug[pipe.tag[0:3]][str(0)],round_roug[pipe.tag[0:3]])
                                else:
                                    pipe.roughness = round((((output[m+int(real_materials.index(pipe.tag[0:6]))]+1)*(range_roug[pipe.tag[0:3]][str(int(pipe.tag[4:6])+10)]-range_roug[pipe.tag[0:3]][str(int(pipe.tag[4:6])-10)]))/2)+range_roug[pipe.tag[0:3]][str(int(pipe.tag[4:6])-10)],round_roug[pipe.tag[0:3]])
                            m += len(real_materials)
                        # Save the minor losses outputs
                        if calib_minorloss:
                            for pipe_name, pipe in wn_ex.pipes():
                                pipe.minor_loss = round((pipe.length*((((output[m+int(real_materials.index(pipe.tag[0:6]))]+1)*(rang_minorloss[1]-rang_minorloss[0]))/2)+rang_minorloss[0])),2)
                            m += len(real_materials)
                        # Save the valve coefficient outputs
                        if calib_valves_coeff:
                            for valve_name, valve in wn_ex.valves():
                                if (valve_name in suitable_valves) and (valve.valve_type == "TCV"):
                                    valve.initial_setting = round((((output[m]+1)*(rang_valve_coeff[1]-rang_valve_coeff[0]))/2)+rang_valve_coeff[0],0)
                                    m += 1
                        # Simulate the network with the new outputs    
                        sim = wntr.sim.EpanetSimulator(wn_ex)
                        # Verify the convergence criteria for the simulated network
                        try:
                            results = sim.run_sim(convergence_error=True)
                        #If a error is raised the fitness is increases to discard the individual in the next generation
                        except:
                            print("Handled error")
                            genome.fitness = -99999
                        # Else calculate the fitness based on the procedure 
                        else:
                            # Verify if any element is disconnected
                            local_dir = os.path.dirname(__file__)
                            report = os.path.join(local_dir,"temp.rpt")
                            # opening a text file
                            file1 = open(report, "r")
                            # read file content
                            readfile = file1.read()
                            # checking condition for string found or not
                            if 'disconnected' in readfile: 
                                print("Disconnected")
                                genome.fitness = -99999
                            # Compare each element with the input data and sum the fitness function
                            else: 
                                d = 0
                                div = len(Real_data)/len(Elements["Type"])
                                fitness = []
                                for i, j in Elements.iterrows():
                                    #Compare links "T" (tubes)
                                    if j["Type"] == "T":
                                        caudal1 = results.link['flowrate'].loc[:,j["ID"]].tolist()
                                        caudal1 = [caudal1 * unit_factor for caudal1 in caudal1]
                                        fitness.append(nse(Real_data[int(d*div):int(((d+1)*div))],caudal1))
                                        d +=1
                                    #Compare junctions "N" (nodes)
                                    if j["Type"] == "N":
                                        presion = results.node['pressure'].loc[:,j["ID"]].tolist()
                                        fitness.append(nse(Real_data[int(d*div):int(((d+1)*div))],presion))
                                        d +=1
                                # Sum fitness function
                                genome.fitness = sum(fitness)
                            file1.close()
                # Create a function to execute the energy calibration
                def run_e(config_file):
                    # Define the NEAT configuration
                    config = neat.Config(neat.DefaultGenome, neat.DefaultReproduction,
                                         neat.DefaultSpeciesSet, neat.DefaultStagnation,
                                         config_file)
                    # Load population configuration
                    p_e = neat.Population(config)
                    # Create temporal checkpoint path to reload best fitness individual to the next process
                    t_e = network.split("/")
                    # Temporal checkpoint based on loop and number of generations
                    t_e[-1] = "checkpoint-e-"+str(int(int(generations)*(loop-1))-1)
                    temp_check_e = "/".join(t_e)
                    # Load previous checkpoint if exists
                    if os.path.exists(temp_check_e):
                        p_e = neat.Checkpointer.restore_checkpoint("checkpoint-e-"+str(int(int(generations)*(loop-1))-1))
                        print("Restored checkpoint e "+t_e[-1])
                    # Add reporter
                    p_e.add_reporter(neat.StdOutReporter(True))
                    stats = neat.StatisticsReporter()
                    p_e.add_reporter(stats)
                    # Create checkpoint based on loop and number of generations
                    p_e.add_reporter(neat.Checkpointer(generation_interval=int(int(generations)*loop-(loop-1)),time_interval_seconds=None,filename_prefix="checkpoint-e-"))
                    # Execute fitness function of each individual and save best fitness individual 
                    winner = p_e.run(eval_genomes_e, int(generations))
                    # Save best fitness ANN
                    winner_net = neat.nn.FeedForwardNetwork.create(winner, config)
                    # Save outputs
                    output = winner_net.activate(discretization_energy)
                    # Load network to save best fitness individual parameters
                    wn_ex = wntr.network.WaterNetworkModel(resul_ex)
                    # Create counter based on the number of parameter
                    n = 0
                    # Save best roughness information
                    if calib_roughness:
                        for pipe_name, pipe in wn_ex.pipes():
                            if int(pipe.tag[4:6]) < 10:
                                pipe.roughness = round((((output[n+int(real_materials.index(pipe.tag[0:6]))]+1)*(range_roug[pipe.tag[0:3]][str(20)]-range_roug[pipe.tag[0:3]][str(0)]))/2)+range_roug[pipe.tag[0:3]][str(0)],round_roug[pipe.tag[0:3]])
                            else:
                                pipe.roughness = round((((output[n+int(real_materials.index(pipe.tag[0:6]))]+1)*(range_roug[pipe.tag[0:3]][str(int(pipe.tag[4:6])+10)]-range_roug[pipe.tag[0:3]][str(int(pipe.tag[4:6])-10)]))/2)+range_roug[pipe.tag[0:3]][str(int(pipe.tag[4:6])-10)],round_roug[pipe.tag[0:3]])
                            final_r.append(pipe.roughness)
                        n += len(real_materials)
                    # Save best minor losses information
                    if calib_minorloss:
                        for pipe_name, pipe in wn_ex.pipes():
                            pipe.minor_loss = round((pipe.length*((((output[n+int(real_materials.index(pipe.tag[0:6]))]+1)*(rang_minorloss[1]-rang_minorloss[0]))/2)+rang_minorloss[0])),2)
                            final_ml.append(pipe.minor_loss)
                        n += len(real_materials)
                    # Save best valve coefficient information        
                    if calib_valves_coeff:
                        for valve_name, valve in wn_ex.valves():
                            if (valve_name in suitable_valves) and (valve.valve_type == "TCV"):
                                valve.initial_setting = round((((output[n]+1)*(rang_valve_coeff[1]-rang_valve_coeff[0]))/2)+rang_valve_coeff[0],0)
                                n += 1
                                final_valcoeff.append(str(valve.initial_setting))
                    # rewrite best fitness information in the network file
                    wn_ex.write_inpfile(wn_ex.name)
                # Execute run energy function
                local_dir = os.path.dirname(__file__)
                config_path_e = os.path.join(local_dir, 'temp_neat_e.txt')
                run_e(config_path_e)
            # Print "Process #n finished!"
            print("Process #"+str(loop)+" finished!")
            # Count +1 loop
            loop += 1
            #End NEAT Energy analysis
        print("Saving results...")
        #Execute simulation final results
        wnia = wntr.network.WaterNetworkModel(resul_ex)
        sim2 = wntr.sim.EpanetSimulator(wnia)
        resex = sim2.run_sim()
        # Exexute simulation initial results
        wn_init = wntr.network.WaterNetworkModel(network)
        #Plot flow and pressure Comparison after calibration
        # Counter o for each element calibrated
        o = 0
        div = len(Real_data)/len(Elements["Type"])
        for i, j in Elements.iterrows():
            # Plot pressure Comparison for each node
            if j["Type"] == "N":     
                pinr = pnd.DataFrame(columns = ["Real", "Calibrated"])
                pinr["Calibrated"]= resex.node['pressure'].loc[:,j["ID"]]
                pinr["Real"]= Real_data[int(o*div):int(((o+1)*div))]
                pinr.index /= 3600 # convert time to hours
                lab = pinr.plot(kind="line",figsize=(8,4),grid = True)
                lab.set_xlabel("Time [hr]",fontsize=10)
                lab.set_ylabel("Pressure [m]",fontsize=10)
                lab.set_title("Pressure_node_"+j["ID"])
                plt.savefig('Results/'+net_name+"_PressureComparison_Element_"+str(o)+".png",dpi=800)
                plt.close()
                o += 1
            # Plot flow Comparison for each link
            if j["Type"] == "T":
                #Pipes
                qinr = pnd.DataFrame(columns = ["Real", "Calibrated"])
                qinr["Calibrated"]= resex.link['flowrate'].loc[:,j["ID"]]
                qinr["Calibrated"]= qinr["Calibrated"]*unit_factor
                qinr["Real"] = Real_data[int(o*div):int(((o+1)*div))]
                qinr.index /= 3600 # convert time to hours
                lab = qinr.plot(kind="line",figsize=(8,4),grid = True)
                lab.set_xlabel("Time [hr]",fontsize=10)
                lab.set_ylabel("Flow ["+str(wn_ex.options.hydraulic.inpfile_units)+"]",fontsize=10)
                lab.set_title("Flow_pipe_"+j["ID"])
                plt.savefig('Results/'+net_name+"_FlowComparison_Element_"+str(o)+".png",dpi=800)
                plt.close()
                o +=1
        # Plot Comparison figures
        # Calculate nearest path to each reservoir
        element = []
        elemtype = []
        for name, x in wn_ex.reservoirs():
            element.append(name)
            elemtype.append('N')
        nearest_elem = dijkstra(resul_ex, element, elemtype).run()
        # Counter o for each element calibrated
        o = 0
        for i, j in Elements.iterrows():
            # Plot Comparison_Pressure_Head_Demand for each element
            if j["Type"] == "N":     
                pinr = pnd.DataFrame(columns = ["Real", "Calibrated"])
                pinr["Calibrated"]= resex.node['pressure'].loc[:,j["ID"]]
                pinr["Real"]= Real_data[int(o*div):int(((o+1)*div))]
                pinr.index /= 3600 # convert time to hours
                # Create subplot figure
                fig, axes = plt.subplots(nrows=3, ncols=1, constrained_layout=True)
                # Plot Pressure
                lab = pinr.plot(kind="line",figsize=(8,4),grid = True, ax=axes[0])
                lab.set_xlabel("Time [hr]",fontsize=10)
                lab.set_ylabel("Pressure [m]",fontsize=10)
                lab.set_title("Pressure_node_"+j["ID"])
                # Plot Head in the nearest reservoir
                near_head = pnd.DataFrame(columns = ["Head"])
                near_head["Head"]= resex.node['head'].loc[:,str(wn_ex.reservoir_name_list[nearest_elem[wn_ex.node_name_list.index(j["ID"])]])]
                near_head.index /= 3600 # convert time to hours
                lab = near_head.plot(kind="line",figsize=(8,4),grid = True, ax=axes[1])
                lab.set_xlabel("Time [hr]",fontsize=10)
                lab.set_ylabel("Head [m]",fontsize=10)
                lab.set_title("Head_reservoir_"+j["ID"])
                # Plot Demand in the node
                near_demand = pnd.DataFrame(columns = ["Demand"])
                near_demand["Demand"]= resex.node['demand'].loc[:,j["ID"]]
                near_demand.index /= 3600 # convert time to hours
                lab = near_demand.plot(kind="line",figsize=(8,4),grid = True, ax=axes[2])
                lab.set_xlabel("Time [hr]",fontsize=10)
                lab.set_ylabel("Demand ["+str(wn_ex.options.hydraulic.inpfile_units)+"]",fontsize=10)
                lab.set_title("Demand_"+j["ID"])
                # Save figure
                fig.savefig('Results/'+net_name+"_Comparison_Pressure_Head_Demand_Element_"+str(o)+".png",dpi=1000)
                plt.close()
                o += 1
            # Plot Comparison_Flow_Head_Pressure for each element
            if j["Type"] == "T":
                qinr = pnd.DataFrame(columns = ["Real", "Calibrated"])
                qinr["Calibrated"]= resex.link['flowrate'].loc[:,j["ID"]]
                qinr["Calibrated"]= qinr["Calibrated"]*unit_factor
                qinr["Real"] = Real_data[int(o*div):int(((o+1)*div))]
                qinr.index /= 3600 # convert time to hours
                # Create subplot figure
                fig, axes = plt.subplots(nrows=3, ncols=1, constrained_layout=True)
                # Plot Flow
                lab = qinr.plot(kind="line",figsize=(8,4),grid = True, ax=axes[0])
                lab.set_xlabel("Time [hr]",fontsize=10)
                lab.set_ylabel("Flow [CMH]",fontsize=10)
                lab.set_title("Flow_pipe_"+j["ID"])
                # Plot Head in the nearest reservoir
                near_head = pnd.DataFrame(columns = ["Head"])
                near_head["Head"]= resex.node['head'].loc[:,str(wn_ex.reservoir_name_list[nearest_elem[wn_ex.node_name_list.index(str(wn_ex.get_link(j["ID"]).start_node))]])]
                near_head.index /= 3600 # convert time to hours
                lab = near_head.plot(kind="line",figsize=(8,4),grid = True, ax=axes[1])
                lab.set_xlabel("Time [hr]",fontsize=10)
                lab.set_ylabel("Head [m]",fontsize=10)
                lab.set_title("Head_reservoir_"+j["ID"])
                # Plot Pressure in the nearest node
                near_pressure = pnd.DataFrame(columns = ["Pressure"])
                near_pressure["Pressure"]= resex.node['pressure'].loc[:,str(wn_ex.get_link(j["ID"]).start_node)]
                near_pressure.index /= 3600 # convert time to hours
                lab = near_pressure.plot(kind="line",figsize=(8,4),grid = True, ax=axes[2])
                lab.set_xlabel("Time [hr]",fontsize=10)
                lab.set_ylabel("Pressure [CMH]",fontsize=10)
                lab.set_title("Pressure_"+str(wn_ex.get_link(j["ID"]).start_node))
                # Save figure
                fig.savefig('Results/'+net_name+"_Comparison_Flow_Head_Pressure_Element_"+str(o)+".png",dpi=1000)
                plt.close()
                o +=1
        # Save Comparison before and after calibration of the parameters
        # Save Comparison before and after calibration of the demands
        if calib_demand:
            #Load before and after calibration networks
            wn_init = wntr.network.WaterNetworkModel(network)
            wn_calib = wntr.network.WaterNetworkModel(resul_ex)
            c = 0
            # Load before and after calibration information for plot purposes
            for node_name, node in wn_init.junctions():
                if node.demand_timeseries_list[0].base_value == 0:
                    node.elevation = None
                else:
                    node.elevation = ((wn_calib.get_node(node_name).demand_timeseries_list[0].base_value-node.demand_timeseries_list[0].base_value)/node.demand_timeseries_list[0].base_value)
                c += 1
            # Verify possible raised errors plotting the figure
            try:
                wntr.graphics.plot_network(wn_init,title= "Comparison before and after calibration of the demands", node_attribute='elevation',node_size=15, node_cmap="Accent_r")
            except:
                print("Base demand Image Error")
                pass
            else:
                plt.savefig('Results/'+net_name+"Comparison_Demand.png",dpi=800,orientation="landscape")
        # Save Comparison before and after calibration of the leaks
        if calib_leak:
            #Load before and after calibration networks
            wn_init = wntr.network.WaterNetworkModel(network)
            wn_calib = wntr.network.WaterNetworkModel(resul_ex)
            # Load before and after calibration information for plot purposes
            c = 0
            for node_name, node in wn_calib.junctions():
                if wn_calib.get_node(node_name).emitter_coefficient == None:
                    node.elevation = None
                else:
                    node.elevation = wn_calib.get_node(node_name).emitter_coefficient*3600
                    c += 1
            # Verify possible raised errors plotting the figure
            try:
                wntr.graphics.plot_network(wn_calib,title= "Comparison before and after calibration of the leaks", node_attribute='elevation',node_size=50, node_cmap="rainbow")
            except:
                print("Emitter Coeff Image Error")
                pass
            else:
                plt.savefig('Results/'+net_name+"Comparison_Leak.png",dpi=800,orientation="landscape") 
        # Save Comparison before and after calibration of the roughness
        if calib_roughness:  
            #Load before and after calibration networks
            wn_init = wntr.network.WaterNetworkModel(network)
            wn_calib = wntr.network.WaterNetworkModel(resul_ex)
            # Load before and after calibration information for plot purposes
            c = 0
            for pipes_name, pipe in wn_init.pipes():
                if pipe.roughness == wn_calib.get_link(pipes_name).roughness: 
                    pipe.lenght = None
                else:
                    pipe.length = ((wn_calib.get_link(pipes_name).roughness-pipe.roughness)/pipe.roughness)*100
                c += 1
            # Verify possible raised errors plotting the figure
            try:
                wntr.graphics.plot_network(wn_init,title="Comparison before and after calibration of the roughness", link_attribute='length', node_size=0.1, link_width=2)
            except:
                print("Roughness Image Error")
                pass
            else:
                plt.savefig('Results/'+net_name+"Comparison_Roguhness.png",dpi=800,orientation="landscape")
        # Save Comparison before and after calibration of the minor loss
        if calib_minorloss: 
            #Load before and after calibration networks
            wn_init = wntr.network.WaterNetworkModel(network)
            wn_calib = wntr.network.WaterNetworkModel(resul_ex)
            # Load before and after calibration information for plot purposes
            c = 0
            for pipes_name, pipe in wn_init.pipes():
                if wn_calib.get_link(pipes_name).minor_loss == pipe.minor_loss or wn_calib.get_link(pipes_name).minor_loss == 0:
                    pipe.length = None
                else:
                    pipe.length = ((wn_calib.get_link(pipes_name).minor_loss))
                c += 1
            # Verify possible raised errors plotting the figure
            try:
                wntr.graphics.plot_network(wn_init,title="Comparison before and after calibration of the minor loss", link_attribute='length', node_size=0.1, link_width=2)
            except:
                print("Minor Loss Image Error")
                pass
            else:
                plt.savefig('Results/'+net_name+"Comparison_Minorloss.png",dpi=800,orientation="landscape")
        # Save calibrated information in a xlsx file
        # Save roughness information in DataFrame
        if calib_roughness:
            # Save initial roughness
            init_r = []
            for pipe_name, pipe in wn_init.pipes():
                init_r.append(pipe.roughness)
            rough_name = pnd.DataFrame(wn_init.pipe_name_list, columns=["Pipe Name"])
            roug_init = pnd.DataFrame(init_r, columns=["R_init [mm]"])
            roug_sim = pnd.DataFrame(final_r, columns=["R_calib [mm]"])
            
        # Save minor loss information in DataFrame
        if calib_minorloss:
            # Save initial minor loss
            init_ml = []
            for pipe_name, pipe in wn_init.pipes():
                init_ml.append(pipe.minor_loss)
            km_name = pnd.DataFrame(wn_init.pipe_name_list, columns=["Pipe Name"])
            km_init = pnd.DataFrame(init_ml, columns=["Km_init [-]"])
            km_sim = pnd.DataFrame(final_ml, columns=["Km_calib [-]"])
            
        # Save demand information in DataFrame
        if calib_demand:
            init_dem = []
            # Save initial demand information
            for node_name, node in wn_init.junctions():
                init_dem.append(node.demand_timeseries_list[0].base_value)
            dem_name = pnd.DataFrame(wn_init.node_name_list, columns=["Node Name"])
            dem_init = pnd.DataFrame(init_dem, columns=["Dem_init ["+str(wn_ex.options.hydraulic.inpfile_units)+"]"])
            dem_sim = pnd.DataFrame(final_dem, columns=["Dem_calib ["+str(wn_ex.options.hydraulic.inpfile_units)+"]"])
            
        # Save leak information in DataFrame
        if calib_leak:
            init_leak = []
            # Save initial leak information
            for node_name, node in wn_init.junctions():
                final_leak.append(node.emitter_coefficient)
            leak_name = pnd.DataFrame(wn_init.node_name_list, columns=["Node Name"])
            leak_init = pnd.DataFrame(init_leak, columns=["Leak_init [-]"])
            leak_sim = pnd.DataFrame(final_leak, columns=["Leak_calib [-]"])
            
        # Save valve coefficient information in DataFrame
        if calib_valves_coeff:
            init_valcoeff = []
            # Save initial valve coeff information
            for valve_name, valve in wn_init.valves():
                if (valve_name in suitable_valves) and (valve.valve_type == "TCV"):
                    init_valcoeff.append(valve.initial_setting)
            valcoeff_name = pnd.DataFrame(suitable_valves, columns=["Valve Name"])
            valcoeff_init = pnd.DataFrame(init_valcoeff, columns=["Val_coeff_init [-]"])
            valcoeff_sim = pnd.DataFrame(final_valcoeff, columns=["Val_coeff_calib [-]"])
            
        # Save valve status information in DataFrame
        if calib_valve_status:
            init_valstat = []
            # Save initial valve status information
            for valve_name, valve in wn_init.valves():
                if (valve_name in suitable_valves) and (valve.valve_type == "TCV"):                                  
                    init_valstat.append(str(valve.initial_status))
            valstat_name = pnd.DataFrame(suitable_valves, columns=["Valve Name"])
            valstat_init = pnd.DataFrame(init_valstat, columns=["Val_Stat_init [-]"])
            valstat_sim = pnd.DataFrame(final_valstat, columns=["Val_Stat_calib [-]"])
        
        # Save information in xlsx using i as column counter
        i = 0
        with pnd.ExcelWriter('Results/'+net_name+"_Comparison_Element_Results.xlsx") as writer:  
            # Save roughness
            if calib_roughness:
                rough_name.to_excel(writer, startcol = i, index = False)
                i += 1
                roug_init.to_excel(writer, startcol = i, index = False)
                i += 1
                roug_sim
                roug_sim.to_excel(writer, startcol = i, index = False)
                i += 2
            # Save minor loss
            if calib_minorloss:
                km_name.to_excel(writer, startcol = i, index = False)
                i += 1
                km_init.to_excel(writer, startcol = i, index = False)
                i += 1
                km_sim.to_excel(writer, startcol = i, index = False)
                i += 2
            # Save demand
            if calib_demand:
                dem_name.to_excel(writer, startcol = i, index = False)
                i += 1
                dem_init.to_excel(writer, startcol = i, index = False)
                i += 1
                dem_sim.to_excel(writer, startcol = i, index = False)
                i += 2
            # Save leaks
            if calib_leak:
                leak_name.to_excel(writer, startcol = i, index = False)
                i += 1
                leak_init.to_excel(writer, startcol = i, index = False)
                i += 1
                leak_sim.to_excel(writer, startcol = i, index = False)
                i += 2
            # Save valves coefficient
            if calib_valves_coeff:
                valcoeff_name.to_excel(writer, startcol = i, index = False)
                i += 1
                valcoeff_init.to_excel(writer, startcol = i, index = False)
                i += 1
                valcoeff_sim.to_excel(writer, startcol = i, index = False)
                i += 2
            # Save valve status
            if calib_valve_status:
                valstat_name.to_excel(writer, startcol = i, index = False)
                i += 1
                valstat_init.to_excel(writer, startcol = i, index = False)
                i += 1
                valstat_sim.to_excel(writer, startcol = i, index = False)
                i += 2
        # Save signals of flow and pressure in a xlsx file
        # presure xlsx vector
        p_xls_c = []
        # flow xlsx vector
        q_xls_c = []
        # Organize information to save
        div = len(Real_data)/len(Elements["Type"])
        for i, j in Elements.iterrows():
        #Create xlsx file
            if j["Type"] == "N":     
                p_xls_c.append("N_"+j["ID"]) 
            if j["Type"] == "T":
                q_xls_c.append("T_"+j["ID"]) 
        # Save data in Dataframe
        p_xls = pnd.DataFrame(columns=p_xls_c)
        q_xls = pnd.DataFrame(columns=q_xls_c)
        # Recalculate de flow units using the unit factor
        for i, j in Elements.iterrows():
            if j["Type"] == "T":
                caudal1 = resex.link['flowrate'].loc[:,j["ID"]].tolist()
                caudal1 = [caudal1 * unit_factor for caudal1 in caudal1]
                q_xls["T_"+j["ID"]] = caudal1
            if j["Type"] == "N":     
                p_xls["N_"+j["ID"]] = resex.node['pressure'].loc[:,j["ID"]].tolist()
        # Save signals in xlsx files
        with pnd.ExcelWriter('Results/'+net_name+"_Flow_Pressure_Signals.xlsx") as writer:
            p_xls.to_excel(writer, sheet_name="Pressure")
            q_xls.to_excel(writer, sheet_name="Flow")
        
        print("Calibration process successfully")