# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 15:21:42 2022

@author: Cristian Camilo Gómez Cortés

### Description:
Dijkstra’s algorithm is an algorithm for finding the shortest paths between nodes in a graph from a given list of elements.

### The input parameters are the following:
    - network: (str) - File path to the .inp file of the current network.
    - element: (list[str]) - list containing the unic id (nmae) of each element.
    - elemtype: (list[str]) - list containing the types of each element, either N for nodes or T for links(tubes). 
### The output parameters are the following:
    - Image file with the nearest elements
    - (nearest_elem) List containing the nearest element of each node in the network
"""
# Load required packages
import heapq
import numpy as np
from numpy import Inf
import wntr
import math
import matplotlib.pyplot as plt

class dijkstra:
    #Initialize the elements of the function
    def __init__(self, network, element, elemtype):
        self.network = network
        self.element = element
        self.elemtype = elemtype
        #self.flow = flow
    #Execute the Dijkstra's algoritm
    def dijkstras(self, graph, root):
        n = len(graph)
        # set up "inf" distances
        dist = [Inf for _ in range(n)]
        # set up root distance
        dist[root] = 0
        # set up visited node list
        visited = [False for _ in range(n)]
        # set up priority queue
        pq = [(0, root)]
        # while there are nodes to process
        while len(pq) > 0:
            # get the root, discard current distance
            _, u = heapq.heappop(pq)
            # if the node is visited, skip
            if visited[u]:
                continue
            # set the node to visited
            visited[u] = True
            # check the distance and node and distance
            for v, l in graph[u]:
                # if the current node's distance + distance to the node we're visiting
                # is less than the distance of the node we're visiting on file
                # replace that distance and push the node we're visiting into the priority queue
                if dist[u] + l < dist[v]:
                    dist[v] = dist[u] + l
                    heapq.heappush(pq, (dist[v], v))
        return dist
    # Organize the results based on the topological information of the network and plot it.
    def run(self):
        self.wn = wntr.network.WaterNetworkModel(self.network)

        self.edges = self.wn.node_name_list
        self.graph = {}
        #Create the graph for the network information
        for i in self.edges:
            self.graph[self.edges.index(i)] = []
        
        for pipe_name, pipe in self.wn.links():
            if pipe.link_type == 'Valve':
                self.graph[self.edges.index(pipe.start_node_name)].append((self.edges.index(pipe.end_node_name),0))
                self.graph[self.edges.index(pipe.end_node_name)].append((self.edges.index(pipe.start_node_name),0))
            if pipe.link_type == 'Pipe':
                self.graph[self.edges.index(pipe.start_node_name)].append((self.edges.index(pipe.end_node_name),math.ceil(pipe.length/100)))
                self.graph[self.edges.index(pipe.end_node_name)].append((self.edges.index(pipe.start_node_name),math.ceil(pipe.length/100)))
                
        seed = []
        for k in range(len(self.element)):
            if self.elemtype[k] == "T":
                seed.append(self.wn.get_link(self.element[k]).start_node_name)
            if self.elemtype[k] == "N":
                seed.append(self.element[k])
 
        dijk = []
        # Run the Dijkstra's algoritm
        for r in seed:
            dijk.append(self.dijkstras(self.graph, self.edges.index(r)))
        dijk = np.transpose(dijk)
        nearest_elem = []
        # Organize the outputs
        for l in range(len(dijk)):
            nearest_elem.append(dijk[l].argmin())
        # Save the information to plot it
        for z in range(len(self.edges)):
            self.wn.get_node(self.wn.node_name_list[z]).elevation = nearest_elem[z]
        # Plot the results
        wntr.graphics.plot_network(self.wn,title="Nearest nodes", node_attribute='elevation',node_colorbar_label='Zone', node_size=10)
        plt.savefig("Nearest nodes.png",dpi=300,orientation="landscape")
        # Return 
        return nearest_elem


    