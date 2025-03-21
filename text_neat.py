# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 10:44:52 2021

@author: Cristian Camilo Gómez Cortés

### Description:
Create a copy of the parameter neat file 'config-feedforward.txt' and restructured with the inputs provided.

### The input parameters are the following:
    - prefix: (str) - prefix name provided to identify if the file corresponce to a mass or energy calibration
    - pop_size: (int) - population size of each generation. Recomended: 100
    - num_hidden: (int) - number of nodes in the hidden layers of initial artificial neural network. Recomended: 0
    - num_inputs:(int) - Number of input parameters
    - num_outputs: (int) - Number of output parameters, dependin on the number of variables calibrated.
### The output parameters are the following:
    - New .txt file with the modified parameters
"""
# Load required packages
import shutil
import os

class copy_text_neat:
    def __init__(self ,prefix ,pop_size, num_hidden, num_inputs, num_outputs):
        search1 = "reset_on_extinction   = False"
        search3 = "reset_on_extinction   = True"
        search2 = "# network parameters"
        local_dir = os.path.dirname(__file__)
        config_path = os.path.join(local_dir, 'config-feedforward_'+str(prefix)+'.txt')
        temp = os.path.join(local_dir, 'temp_neat_'+prefix+'.txt')
        shutil.copy(config_path, temp)
        with open(temp, "r") as read_file:
            data = read_file.readlines()
        num = 0
        for line in data:
            line = line.rstrip()
            if search1 in line or search3 in line:
                data[num-1] = "pop_size              = "+str(pop_size)+"\n"
            if search2 in line:
                data[num+1] = "num_hidden              = "+str(num_hidden)+"\n"
                data[num+2] = "num_inputs              = "+str(num_inputs)+"\n"
                data[num+3] = "num_outputs             = "+str(num_outputs)+"\n"
            num += 1
        with open(temp, 'w') as file:
            file.writelines( data )