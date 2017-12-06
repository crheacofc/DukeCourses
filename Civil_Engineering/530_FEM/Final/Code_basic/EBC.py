#Application of EBC to F and K matrices
import numpy as np
from read_mesh import read_in


def Apply_EBC(K,F,mesh,BC_nodes,BC_vals):
    #first we need to see which nodes are not in the boundary!
    nodes = []
    print(BC_vals)
    for i in range(len(mesh)):
        nodes.append(i)
    Nodes_wout_IC = []
    for i in range(len(nodes)):
        if (nodes[i] in BC_nodes):
            pass
        else:
            Nodes_wout_IC.append(i)
    #print(Nodes_wout_IC)
    #print(len(K.T))
    #now to fix F and K matrices
    K_corrected = np.zeros((len(Nodes_wout_IC), len(Nodes_wout_IC)))
    F_corrected = np.zeros((len(Nodes_wout_IC),1))
    print(F_corrected)
    for i in range(len(Nodes_wout_IC)):
        current_node_without = Nodes_wout_IC[i]
        for j in range(len(BC_nodes)):
            #print(len(BC_nodes))
            current_node_with = BC_nodes[j]
            #print(current_node_with,current_node_without)
            #now take away BC_val at node with bc from affected nodes
            if j==1: #if the first one affected
                F_corrected[i] = F[current_node_without] - BC_vals[j]*K[current_node_without,current_node_with]
            else: #all others based on what we already have for F_corrected just calculated
                F_corrected[i] -= BC_vals[j]*K[current_node_without,current_node_with]
        #end j
    #end i





    for k in range(len(Nodes_wout_IC)):
        row = Nodes_wout_IC[k]
        for c in range(len(Nodes_wout_IC)):
            col = Nodes_wout_IC[c]
            K_corrected[k,c] = K[row,col]

    return K_corrected, F_corrected, Nodes_wout_IC
