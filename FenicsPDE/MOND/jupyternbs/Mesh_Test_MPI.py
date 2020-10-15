#!/usr/bin/env python
# coding: utf-8

# In[1]:


from dolfin import *

#Importing MPI for parallel computing
# from mpi4py import MPI

# #Importing the PETSc module for parallel use
# from petsc4py import PETSc

#importing mshr for all mesh functions
import mshr as mshr

# Use SymPy to compute f from the manufactured solution u
import sympy as sym

#Option to avoid printing redundant information from each core when running the code in parallel from
#a python (.py) script obtained from the jupyter notebook.
# parallel_run = True

#MPI communicator
comm = MPI.comm_world

#Rank of each process (its ID essentially)
rank = MPI.rank(comm)

#Total number of processes
number_processes = MPI.size(comm)

print(f'This is process {rank} out of {number_processes-1}')

if number_processes <2:

    #Increasing the width of the notebook (visual difference only)
    from IPython.core.display import display, HTML
    display(HTML("<style>.container { width:100% !important; }</style>"))
    
    #have to define where to put plots BEFORE importing matplotlib
    get_ipython().run_line_magic('matplotlib', 'notebook')

#Importing matplotlib to plot the results
from matplotlib import pyplot as plt

#Importing numpy to work with arrays
import numpy as np

#Importing time to compute how long each segment takes
import time

#importing regex to change every instance of radius_tot so we change the ones in the C++ code
#at the same time too
import re

#varname gives the name of the variable as a string
from varname import varname

#Needed to use the 3D scatter
from mpl_toolkits.mplot3d import Axes3D

#Importing the decimal package to be able to specify arbitrary accuracy, needed e.g. when
#calculating the jacobian for the lensing
from decimal import *

#Importing all quantities, constants etc used in the calculations
from MONDquantities import *

#Importing all classes I created
from MONDclasses import *

#Importing the functions I made from the MONDfunctions file
from MONDfunctions import *

#Importing all expressions for weak forms, initial guesses/BCs and sources
from MONDexpressions import *

#Needed if want to use the adapt function for mesh refinement, see:
#https://fenicsproject.org/qa/6719/using-adapt-on-a-meshfunction-looking-for-a-working-example/
#If using 'plaza' instead of 'plaza_with_parent_facets', it's faster by about 30%!
parameters["refinement_algorithm"] = "plaza"

#ParMETIS is optimised to partition the mesh in parallel. SCOTCH is the default I used for serial
#ParMETIS: http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview. Unfortunately Dolfin wasnt
#built with ParMETIS so i cant use it.
parameters['mesh_partitioner'] = 'SCOTCH'

#Setting compiler parameters.
#Optimisation
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True

#Ghost mode for when using MPI. Each process gets ghost vertices for the part of the domain it does not
#own. Have to set to 'none' instead or I get Error 'Unable to create BoundaryMesh with ghost cells.'
parameters['ghost_mode'] = 'none'


# In[2]:


# info(parameters,True)


# In[3]:


## starting time of whole PDE solver
starting_time = time.time()

#starting an empty list to contain all of the run_time objects to plot later
section_times = []

print('Starting mesh generation...\n')
mesh_generation_start = time.time()

#Making mesh from function defined above
mesh = make_spherical_mesh(domain_size, mesh_resolution)

mesh_generation_end = time.time()
mesh_generation_time = run_time(mesh_generation_end - mesh_generation_start, 'Mesh Generation')
section_times.append(mesh_generation_time)
print('Mesh generated in {} s \n'.format(mesh_generation_time.time))

#Setting the MPI communicator for the mesh (doesnt seem to do anything right now)
# mesh.mpi_comm = comm

# print(f'Process {rank} says marmamamam')


# In[4]:


#For all the points to be within a given radius, each coordinate must be smaller than
#radius_population/sqrt(3)
random_max_distance = radius_population/sqrt(3)

#Setting a given seed so we can always have the same random numbers for now
np.random.seed(1)

#We want a mean of 0 so center of mass is in center, and the same standard deviation as the gaussian
#pulse. This means we sample from the same distribution as the smooth one, and have the same mean.
#This is exactly what we want to compare coarse and smooth distributions
mu, sigma = 0, stand_dev
random_coordinates_x = np.random.normal(mu, sigma, source_number)
random_coordinates_y = np.random.normal(mu, sigma, source_number)

#If we want all source to be in the same plane, we set the z axis to be 0 for all of them. Otherwise,
#random as above
if coplanar_sources == True:
    
    random_coordinates_z = np.zeros((source_number, 1)).ravel()

else:
    
    random_coordinates_z = np.random.normal(mu, sigma, source_number)

#If we dont need Gaussian, defining a source_number*3 array of random numbers between 0 and 1 and
#multiplying by the radius just defined so all points are inside a sphere of radius_tot.
#Subtracting 0.5 so #we're sampling equally from the positive and negative instead of from 0 to 1
# random_coordinates = random_max_distance * (np.random.rand(source_number, 3)-0.5)

# Uncomment for test case with two equal masses on the xy plane at a given distance
# their_distance = 3

# # random_coordinates[0][0] = -domain_size/their_distance
if central_mass:

    random_coordinates_x[0] = 0
    random_coordinates_y[0] = 0
    random_coordinates_z[0] = 0

# random_coordinates[1][0] = domain_size/their_distance
# random_coordinates[1][1] = 0
# random_coordinates[1][2] = 0

#Overall array containing all coordinates. If over-writing the random position, this has to go after it,
#otherwise the c++ array for the source sets the wrong position!
random_coordinates = np.array((random_coordinates_x, random_coordinates_y, random_coordinates_z))
random_coordinates = np.transpose(random_coordinates)

#Obtaining the center of each source as a list of points
source_centers = [Point(random_coordinates_x[i], random_coordinates_y[i], random_coordinates_z[i]) for i in range(source_number)]

len(source_centers)


# In[5]:


print('Starting mesh refinement...\n')
mesh_refine_start = time.time()


#Finding the #cells in each mesh, so if the ID of the collision is larger than the # cells, we are sure
#that submesh doesnt contain that point.
number_cells = mesh.num_cells()

print(f'The submesh of process {rank} has {number_cells} cells')

cell_containing = intersect(mesh, source_centers[0]).intersected_cells()

print(f'Process {rank}\'s cell_containing is: {cell_containing}')

#We only want to refine the submesh if it contains the point. If it doesnt contain the point, the list
#cell_containing will be empty! We can check this by the length of the list .If the list has length
#0 that means that submesh doesn't include the point! 
#The reason the function was giving problems to begin with is that we were trying to index that list
#that had zero elements in it, so the 0th element was already too large an index and was out of range!

how_many = refine_times

#Starting # cells before we refine, to compute growth factor
starting_cells = mesh.num_cells()

for i in range(how_many):

    #Declaring Boolean Mesh Function to individuate cell containing point
    contain_function = MeshFunction("bool", mesh, 3)

    #Setting function to False everywhere
    contain_function.set_all(False)

    #Initial number of cells before refinement
    initial_cells = mesh.num_cells()

    #List comprehension containing the cell IDs for the cells containing a source
    #The if statement make sure that the intersect_list is only populated for points
    #that are inside the mesh. For parallel, where the mesh is split, this is necessary!
    intersect_list = [intersect(mesh, source).intersected_cells() for source in source_centers]

    #Setting the cell function contain_function to true for each cell containing a source
    for cell_index in intersect_list:

        #For MPI, the intersect_list might be empty in case the point is not inside the 
        #submesh! So we first need to check that there's an element present. If we don't
        #we'll get an error about index out of range cause index 0 is out of range for an
        #empty list!
        if not len(cell_index) == 0: 

            contain_function[cell_index[0]] = True

    #Refining the mesh only for cells that contain a source
    mesh = refine(mesh, contain_function, redistribute=True)    

    #Final # cells after refinement
    final_cells = mesh.num_cells()

    partial_growth_factor = final_cells/initial_cells

    print(('Iteration {} of {}: The Cell number went from {} to {}, up by a factor {}\n'
          .format(i+1, how_many, initial_cells, final_cells, partial_growth_factor)))

#ratio between # cells at beginning and end of refinement
total_growth_factor = final_cells/starting_cells

print('Cell number went up by a factor {}\n'.format(total_growth_factor))

# cell_containing

mesh_refine_end = time.time()
mesh_refine_time = run_time(mesh_refine_end - mesh_refine_start, 'Mesh Refinement')
print('Mesh refined in {} s \n'.format(mesh_refine_time.time))


# In[6]:


#Checking if computing the collision with the bounding box tree is faster than checking if point is
#inside a cell

mesh


# In[7]:


#See this link: https://fenicsproject.org/qa/12431/reading-hdf5-file-with-fenics-2016-2/

# meshname = 'parallel_mesh'

# f = HDF5File(mesh.mpi_comm(), meshname+".hdf5", 'w')
# f.write(mesh, meshname)


# In[8]:


# mesh = Mesh()
# f = HDF5File(MPI.comm_world, meshname+"hdf5", 'r')
# f.read(mesh, meshname, False)

