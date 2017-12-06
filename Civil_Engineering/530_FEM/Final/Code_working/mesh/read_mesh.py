#This file will handle the reading of mesh!
import numpy as np
def read_in(mesh,connect,sides):
    # This opens a handle to the file, in 'r' read mode
    file_handle = open(mesh, 'r')
    # Read in all the lines of your file into a list of lines
    lines_list = file_handle.readlines()
    # Extract dimensions from first line. Cast values to integers from strings.
    x, y = (float(val)-1 for val in lines_list[0].split())
    # Do a double-nested list comprehension to get the rest of the data into your matrix
    coordinates = [[float(val) for val in line.split()] for line in lines_list[:]]
    #in the form of [x,y]
    #Lets grab the connectivity matrix! form of [n1,n2,n3,n4] for each element! So row - element, col - node
    # This opens a handle to the file, in 'r' read mode
    file_handle = open(connect, 'r')
    # Read in all the lines of your file into a list of lines
    lines_list = file_handle.readlines()
    # Extract dimensions from first line. Cast values to integers from strings.
   # N1,N2,N3,N4 = (int(val)-1 for val in lines_list[0].split())
    # Do a double-nested list comprehension to get the rest of the data into your matrix
    connectivity_matrix = [[int(val)-1 for val in line.split()] for line in lines_list[:]]
    #to access row, con_mat[i] and element con_mat[i][j]

    #read sides
    # Lets grab the connectivity matrix! form of [n1,n2,n3,n4] for each element! So row - element, col - node
    # This opens a handle to the file, in 'r' read mode
    file_handle = open(sides, 'r')
    # Read in all the lines of your file into a list of lines
    lines_list = file_handle.readlines()
    # Extract dimensions from first line. Cast values to integers from strings.
    #S1,S2,S3,S4 = (int(val)-1 for val in lines_list[0].split())
    # Do a double-nested list comprehension to get the rest of the data into your matrix
    sides_mat = [[int(val) - 1 for val in line.split()] for line in lines_list[:]]
    top = []; left = []; bottom = []; right = []
    for i in range(len(sides_mat)):
        top.append(sides_mat[i][0])
        left.append(sides_mat[i][1])
        bottom.append(sides_mat[i][2])
        right.append(sides_mat[i][3])


    return coordinates,connectivity_matrix,top,left,bottom,right

'''
mesh = '../mesh/mesh_nodes_1x1.txt'
connect = '../mesh/mesh_connectivity_1x1.txt'
sides = '../mesh/sides_1x1.txt'
read_in(mesh,connect,sides)
'''
