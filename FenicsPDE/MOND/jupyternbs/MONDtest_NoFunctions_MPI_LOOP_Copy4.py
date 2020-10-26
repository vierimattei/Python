#!/usr/bin/env python
# coding: utf-8

# # Defining all the constants required (copying them from Matlab and adjusting syntax)

# In[ ]:


from dolfin import *

#pandas is needed to import the cluster database
import pandas as pd

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

#Importing tempfile to save numpy arrays from the main script so we can get them back and plot them
#interatively rather than saving a pdf or png!
from tempfile import TemporaryFile

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
#If using 'plaza' instead of 'plaza_with_parent_facets', it's faster by about 30%! Also, I get the
#'*** Warning: Cannot calculate parent facets if redistributing cells'. So for MPI no need to use
#with parent facets!
parameters["refinement_algorithm"] = "plaza"

#Setting compiler parameters.
#Optimisation
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True

#Nonzero initial guess for the Krylov solver. Doesnt seem to make a difference for nonlinear problems
# parameters['krylov_solver']['nonzero_initial_guess'] = True

#Ghost mode for when using MPI. Each process gets ghost vertices for the part of the domain it does not
#own. Have to set to 'none' instead or I get Error 'Unable to create BoundaryMesh with ghost cells.'
parameters['ghost_mode'] = 'none'

#Start of overall for loop over the parameters of the database

#I deleted the original table by mistake, but luckily I had the dataframe open in adifferent notebook
#so I saved it as a pickle. Now all I need to do to import it is read the pickle!
df = pd.read_pickle('cluster_data_pickle.pkl')

#Putting each column in its own list to loop over. Scaling rho_0 by 10^22, and rc by kp
cluster_name = df.loc[:, 'name']
cluster_rc = df.loc[:, 'r_c_frame']*kp
cluster_rho0 = df.loc[:, 'rho_0_frame']*10**(-22)
cluster_beta = df.loc[:, 'beta_frame']

for i in np.arange(66,88,1):
    
    #Setting the parameters for the beta distribution that we want to compute with
    #The domain size is the Abell radius for all the clusters 
    domain_size = radius_abell
    
    #All the other quantities come from the cluster database for the respective index
    r_c = cluster_rc[i]
    rho_0 = cluster_rho0[i]
    beta = cluster_beta[i]
    
    #Making a dummy variable going from 0 to the domain_size to get the total mass by integrating the beta
    #gas density over the whole domain
    r_integration = np.linspace(0,domain_size,10000)

    mgb = (np.trapz(rho_0/(1+(r_integration/r_c)**2)**(3*beta/2)*
                          4*pi*r_integration**2, x = r_integration))
    
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

    print(f'The mesh of process {rank} has {mesh.num_cells()} cells')

    ## Defining coordinates for some test mass distributions

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

    # # random_coordinates[0][0] = -domain_size/their_distance
    if central_mass:

        random_coordinates_x[0] = 0
        random_coordinates_y[0] = 0
        random_coordinates_z[0] = 0

    #Overall array containing all coordinates. If over-writing the random position, this has to go after it,
    #otherwise the c++ array for the source sets the wrong position!
    random_coordinates = np.array((random_coordinates_x, random_coordinates_y, random_coordinates_z))
    random_coordinates = np.transpose(random_coordinates)

    #Obtaining the center of each source as a list of points
    source_centers = [Point(random_coordinates_x[i], random_coordinates_y[i], random_coordinates_z[i]) for i in range(source_number)]

    #For the current case in which all sources have the same mass, we simply divide by #sources
    center_of_mass_x = random_coordinates[:,0].sum()/source_number
    center_of_mass_y = random_coordinates[:,1].sum()/source_number
    center_of_mass_z = random_coordinates[:,2].sum()/source_number

    #Overall center of mass
    center_of_mass = [center_of_mass_x, center_of_mass_y, center_of_mass_z]

    # center_of_mass
    print(f'Process {rank} about to refine')

    print('Starting mesh refinement...\n')
    mesh_refine_start = time.time()
    new_mesh = more_modified_refinement(mesh, source_centers, refine_times)
    # new_mesh = local_refinement(mesh, source_centers, radius_refine, refine_times, technique = 'ring')
    mesh_refine_end = time.time()
    mesh_refine_time = run_time(mesh_refine_end - mesh_refine_start, 'Mesh Refinement')
    section_times.append(mesh_refine_time)
    print('Mesh refined in {} s \n'.format(mesh_refine_time.time))

    mesh = new_mesh

    # Gathering all the data from the mesh AFTER having done the mesh refinement and defined the mesh for plotting

    print('Rearranging mesh data\n')
    rearrange_start = time.time()

    V, vertex_number, x_coords, y_coords, z_coords, r_coords, sorting_index, x_sorted, y_sorted, z_sorted, r_sorted = rearrange_mesh_data(mesh, center_of_mass, degree_PDE)

    #To be able to gather the coordinate arrays with MPI, the coordinates need to be C_contiguous
    x_coords, y_coords, z_coords, r_coords = [np.ascontiguousarray(coord_array) for coord_array in [x_coords, y_coords, z_coords, r_coords]]

    rearrange_end = time.time()
    rearrange_time = run_time(rearrange_end - rearrange_start, 'Mesh data rearrange')
    section_times.append(rearrange_time)
    print('Mesh data rearranged in {} s \n'.format(rearrange_time.time))

    # Defining a few BVP from combinations we use often. Naming scheme: 'weak form_source'

    #BVPs for a discrete dirac mass distribution, for Newton and MOND with/out interpolations
    newton_dirac = BVP(F_Newton, u_Newton, f_multiple_dirac, 'Newton, discrete dirac')
    mond_deep_dirac = BVP(F_MOND_deep, u_displaced_cpp, f_multiple_dirac, 'Deep MOND, discrete dirac')
    mond_simple_dirac = BVP(F_MOND_simple, u_displaced_cpp, f_multiple_dirac, 'Simple MOND, discrete dirac')
    mond_standard_dirac = BVP(F_MOND_standard, u_displaced_cpp, f_multiple_dirac, 'Standard MOND, discrete dirac')

    #BVPs for a discrete gauss mass distribution.
    newton_gauss = BVP(F_Newton, u_Newton, f_multiple_gauss, 'Newton, discrete gauss')
    mond_deep_gauss = BVP(F_MOND_deep, u_displaced_cpp, f_multiple_gauss, 'Deep MOND, discrete gauss')
    mond_simple_gauss = BVP(F_MOND_simple, u_displaced_cpp, f_multiple_gauss, 'Simple MOND, discrete gauss')
    mond_standard_gauss = BVP(F_MOND_standard, u_displaced_cpp, f_multiple_gauss, 'Standard MOND, discrete gauss')

    #BVPs for a continuous distribution, for Newton and MOND with/out interpolations
    newton_continuous = BVP(F_Newton, u_Newton, f_exponent_test, 'Newton, continuous gauss')
    mond_deep_continuous = BVP(F_MOND_deep, u_displaced_cpp, f_exponent_test, 'Deep MOND, continuous gauss')
    mond_simple_continuous = BVP(F_MOND_simple, u_displaced_cpp, f_exponent_test, 'Simple MOND, continuous gauss')
    mond_standard_continuous = BVP(F_MOND_standard, u_displaced_cpp, f_exponent_test, 'Standard MOND, continuous gauss')

    #BVPs for a three parameter beta distribution
    newton_beta = BVP(F_Newton, u_Newton, f_gas_three_beta, 'Newton, beta')
    mond_deep_beta = BVP(F_MOND_deep, u_sphere_cpp, f_gas_three_beta, 'Deep MOND, beta')
    mond_simple_beta = BVP(F_MOND_simple, u_sphere_cpp, f_gas_three_beta, 'Simple MOND, beta')
    mond_standard_beta = BVP(F_MOND_standard, u_sphere_cpp, f_gas_three_beta, 'Standard MOND, beta')
    
    # Trying an alternative method for assigning values inside the c++ expressions by using exec, to avoid the limit on eval!

    #Defining a function for the boundary. Since we only have one BC for the whole boundary, we
    #can make a simple function that returns true for each value on the boundary
    #the on_boundary built-in function takes each point in domain and returns true if on boundary
    def boundary(x, on_boundary):
        return on_boundary

    radius_tot=r_c

    def solve_PDE(the_BVP):
        '''Function takes in a BVP object, which defines the weak form, initial guess/BC and source for
        a PDE, and computes its solution'''

        ## starting time of PDE solver
        solver_start = time.time()
        print('Starting PDE Solver...\n')

        #Defining the source term here, cause the make_source_string function creates a string that 
        #evaluate the expression for a variable called 'source'
        source = the_BVP.source

        #Evaluating the source term obtained from the make_source_string function
        f = eval(make_source_string(source_number))

        #Declaring the expression for the initial guess
        u = (Expression(the_BVP.initial_guess,
        degree = degree_PDE, a0 = a0, ms = ms,mgb = mgb, G = G,  ly = ly, kp = kp, radius_tot = radius_tot,
        volume_out = volume_out, center_of_mass_x = center_of_mass_x,
        center_of_mass_y = center_of_mass_y, center_of_mass_z = center_of_mass_z,
        source_number = source_number, source_mass = source_mass))

        #Declaring the expression for the boundary condition with displaced CM (center of mass)
        boundary_CM = u

        #Declaring the boundary condition. It takes three arguments: function space, value of BC, 
        #section of the boundary (in our case the whole boundary).
        bc = DirichletBC(V, boundary_CM, boundary)

        #Defining the variational problem
        #u is the solution. for linear problems, we'd have to define it as TrialFunction, but for 
        #non-linear we define it as Function directly
        u = interpolate(u, V)

        #defining the test function
        v = TestFunction(V)

        #defining the weak form to be solved
        F = eval(the_BVP.weak_form)

        #Computing the solution for normal deep MOND
        (solve(F == 0, u, bc, solver_parameters={"newton_solver":{"relative_tolerance":1e-6},
                                                 "newton_solver":{"maximum_iterations":200}}))

        solver_end = time.time()
        solver_time = run_time(solver_end - solver_start, 'PDE Solver')
        section_times.append(solver_time)

        print('PDE solved in {}\n'.format(solver_time.time))

        return u, f

    #Waiting for each process to have completed before moving on to solve the PDE
    # MPI.barrier(comm)

    #Defined the quantity BVP_to_solve in the MONDquantities file as a string, so to use it we need to 
    #evaluate it with eval.
    u, f = solve_PDE(eval(BVP_to_solve))

    data_collection_start = time.time()
    print('Collecting data from PDE...\n')

    if plotting_option == True:
        mesh = mesh_for_plots
        V_plot, vertex_number, x_coords, y_coords, z_coords, r_coords, sorting_index, x_sorted, y_sorted, z_sorted, r_sorted = rearrange_mesh_data(mesh, center_of_mass)
        u_plot = interpolate(u, V_plot)
        f_plot = interpolate(f, V_plot)
        u = u_plot
        #Calling the rearrange_mesh_data function to get coordinates and order them based on the 
        #distance from the center of mass

    #The value of the function at each vertex of the mesh is stored in a np array. Its order
    #corresponds to the otder of the mesh.coordinates() values
    potential = u.compute_vertex_values()

    #The value of the source at each vertex of the mesh
    source = 1/(4*pi*G)*f.compute_vertex_values(mesh)

    #Getting the degree from the scalar function space V from the PDE
    degree = V.ufl_element().degree()

    #Laplacian of the solution to get back the scaled mass distribution
    lap = div(grad(u))

    apparent_mass_project = project(lap, V)

    #Projecting the acceleration onto a vector space is expensive, so don't do it unless needed
    #If not needed, set acceleration to 0 everywhere
    if acceleration_needed:

        #To obtain the values for the acceleration, we need to define a new function space, since the 
        #gradient is a vector function is the function space for the PDE is a scalar function space
        W = VectorFunctionSpace(mesh, 'P', degree)

        #Projecting (similar to interpolating) the grad(u) field onto W, gives a function
        acceleration_project = project(grad(u), W)

        #The result of project is n*3,1 np.array, with 3 (x,y,z) values for each of the n vertices
        acceleration = acceleration_project.compute_vertex_values()

        #reshaping the array to split the x,y,z components into their own column each
        acceleration = np.reshape(acceleration, (3, int(acceleration.shape[0]/3)))

    else:

        acceleration = np.zeros((3, len(potential)))


    acceleration_x = acceleration[0]
    acceleration_y = acceleration[1]
    acceleration_z = acceleration[2]

    #Finding the magnitude of the acceleration
    acceleration_magnitude = np.linalg.norm(acceleration, axis=0)

    #Sorting the potential, acceleration and source according to thr r of the vertex they pertain to
    potential_sorted = potential[sorting_index]
    acceleration_magnitude_sorted = acceleration_magnitude[sorting_index]
    source_sorted = source[sorting_index]

    data_collection_end = time.time()
    data_collection_time = run_time(data_collection_end - data_collection_start, 'Data Collection')
    section_times.append(data_collection_time)
    print('Data collected in {} s\n'.format(data_collection_time.time))

    # Calculating the Laplacian of the potential to obtain the apparent dark matter distribution.

    # #Projecting the divergence above onto the same scalar function space as the potential
    # apparent_mass_project = project(apparent_mass_divergence, V)

    integral = assemble(lap*dx)

    #Gathering the values of the mass distribution 
    apparent_mass_distribution = 1/(4*pi*G)*apparent_mass_project.compute_vertex_values()

    #Sorting the mass distribution values
    apparent_mass_distribution_sorted = apparent_mass_distribution[sorting_index]

    (integral/(4*pi*G))/mgb

    # Gathering the potential and coordinate numpy array onto process 0 to have the full solution.

    #First, we need to know how many vertices we have in total in the full mesh to preallocate the array
    #for both the potential and the coordinates. We do this with the MPI reduce operation MPI_SUM
    print(f'Process {rank}: potential has {len(potential)} elements.')

    #We need the total #vertices as an int to define an array. Calling MPI.sum with communicator and
    #value to be summed from each process
    total_mesh_vertices = int(MPI.sum(comm, len(potential)))

    if rank == 0:

        print(f'Process {rank}: the overall potential has {total_mesh_vertices} elements.')

    #Now we can gather all values of the potential and coordinates. First, we define arrays to hold the
    #result, the size of the total potential on process 0

    #Have to initialise the receving buffer for the potential to None on all processes or we get an error
    potential_total = None
    x_coords_total = None
    y_coords_total = None
    z_coords_total = None
    r_coords_total = None
    source_total = None
    apparent_mass_total = None

    if rank == 0:

        #There is a problem with the receive buffer not being big enough. A simple fix for now is to 
        #multiply its size by 1.5, then we can remove all the trailing zeros
        receiver_size = int(1.5*total_mesh_vertices)

        potential_total = np.empty(receiver_size, dtype = type(potential[0]))
        x_coords_total = np.empty(receiver_size, dtype = type(x_coords[0]))
        y_coords_total = np.empty(receiver_size, dtype = type(x_coords[0]))
        z_coords_total = np.empty(receiver_size, dtype = type(x_coords[0]))
        r_coords_total = np.empty(receiver_size, dtype = type(x_coords[0]))
        source_total = np.empty(receiver_size, dtype = type(potential[0]))
        apparent_mass_total = np.empty(receiver_size, dtype = type(potential[0]))

    #IMPORTANT: Have to use Gatherv, not Gather, or it won't work!
    comm.Gatherv(potential, potential_total, root = 0)
    comm.Gatherv(x_coords, x_coords_total, root = 0)
    comm.Gatherv(y_coords, y_coords_total, root = 0)
    comm.Gatherv(z_coords, z_coords_total, root = 0)
    comm.Gatherv(r_coords, r_coords_total, root = 0)
    comm.Gatherv(source, source_total, root = 0)
    comm.Gatherv(apparent_mass_distribution, apparent_mass_total, root = 0)

    #Now we want to sort as usual, now for the total potential and based on the overall r coordinates

    if rank == 0:

        #Storing the index to sort according to the total r
        sorting_index_total = r_coords_total.argsort()

        #Sorting all total quantities
        r_total_sorted = r_coords_total[sorting_index_total]

        x_total_sorted = x_coords_total[sorting_index_total]

        y_total_sorted = y_coords_total[sorting_index_total]

        z_total_sorted = z_coords_total[sorting_index_total]

        potential_total_sorted = potential_total[sorting_index_total]

        source_total_sorted = source_total[sorting_index_total]

        apparent_mass_total_sorted = apparent_mass_total[sorting_index_total]

        #Finding the zero elements in the sorted r array, and removing them. We do this by only keeping the
        #Finding the indices for which r is larger than the smallest r on process 0, divided by 10**5 just to
        #make sure. There should definitely not be any mesh points with distances smaller than that!
        total_nonzero_indices = (r_total_sorted > r_sorted[0]/(10**5))

        #Taking the non-padding components of radius, potential, source and mass distribution
        x_total_sorted = x_total_sorted[total_nonzero_indices]
        y_total_sorted = y_total_sorted[total_nonzero_indices]
        z_total_sorted = z_total_sorted[total_nonzero_indices]
        r_total_sorted = r_total_sorted[total_nonzero_indices]
        potential_total_sorted = potential_total_sorted[total_nonzero_indices]
        source_total_sorted = source_total_sorted[total_nonzero_indices]
        apparent_mass_total_sorted = apparent_mass_total_sorted[total_nonzero_indices]

        #The dark matter is the difference between the apparent and source masses
        dark_mass_total_sorted = apparent_mass_total_sorted - source_total_sorted

        #Saving all these numpy arrays so we can plot them again in Python, instead of just having a saved figure
        #that is not interactive!

        #Saving all the quantities to the respective files
        np.save(f'database_results/{cluster_name[i]}/potential_{cluster_name[i]}.npy', potential_total_sorted)
        np.save(f'database_results/{cluster_name[i]}/source_{cluster_name[i]}.npy', source_total_sorted)
        np.save(f'database_results/{cluster_name[i]}/apparent_{cluster_name[i]}.npy', apparent_mass_total_sorted)
        np.save(f'database_results/{cluster_name[i]}/dark_mass_{cluster_name[i]}.npy', dark_mass_total_sorted)
        np.save(f'database_results/{cluster_name[i]}/x_sorted_{cluster_name[i]}.npy', x_total_sorted)
        np.save(f'database_results/{cluster_name[i]}/y_sorted_{cluster_name[i]}.npy', y_total_sorted)
        np.save(f'database_results/{cluster_name[i]}/z_sorted_{cluster_name[i]}.npy', z_total_sorted)
        np.save(f'database_results/{cluster_name[i]}/r_sorted_{cluster_name[i]}.npy', r_total_sorted)

    if rank == 0:

        print(f'potential_total has length {len(potential_total)}')

        potential_total_no_zeros = potential_total[np.nonzero(potential_total)]

        print(f'potential_total_no_zeros has length {len(potential_total_no_zeros)}')


        print(f'x_coords has length: {len(x_coords)}')
        print(f'x_coords_total has length: {len(x_coords_total)}')
        # x_coords_total
        x_total_no_zeros = x_coords_total[np.nonzero(x_coords_total)]
        print(f'x_total_no_zeros has length: {len(x_total_no_zeros)}')

    # Other instance of main solver to either compare solutions or explore parameter space etc.

    ## First, we compare the three interpolation functions (deep, simple, standard) for some different mass distributions.

    #Lists of same source, different weak form.
    BVP_dirac_list = [newton_dirac, mond_deep_dirac, mond_simple_dirac, mond_standard_dirac]
    BVP_gauss_list = [newton_gauss, mond_deep_gauss, mond_simple_gauss, mond_standard_gauss]
    BVP_continuous_list = [newton_continuous, mond_deep_continuous, mond_simple_continuous, mond_standard_continuous]

    #list of same weak form, different source.
    BVP_deep_list = [mond_deep_dirac, mond_deep_gauss, mond_deep_continuous]
    BVP_simple_list = [mond_simple_dirac, mond_simple_gauss, mond_simple_continuous]
    BVP_standard_list = [mond_standard_dirac, mond_standard_gauss, mond_standard_continuous]
    BVP_newton_list = [newton_dirac, newton_gauss, newton_continuous]

    if make_comparison:

        #Running the compare function
        discrete_list = compare_solutions(BVP_dirac_list, stand_dev, 'stand_dev', 1, '\sigma = ', 'Mpc')


# In[ ]:


# potential_total_sorted = np.load(f'database_results/{cluster_name[0]}/potential_{cluster_name[0]}.npy')
# source_total_sorted = np.load(f'database_results/{cluster_name[0]}/source_{cluster_name[0]}.npy')
# apparent_mass_total_sorted = np.load(f'database_results/{cluster_name[0]}/apparent_{cluster_name[0]}.npy')
# dark_mass_total_sorted = np.load(f'database_results/{cluster_name[0]}/dark_mass_{cluster_name[0]}.npy')
# x_total_sorted = np.load(f'database_results/{cluster_name[0]}/x_sorted_{cluster_name[0]}.npy')
# y_total_sorted = np.load(f'database_results/{cluster_name[0]}/y_sorted_{cluster_name[0]}.npy')
# z_total_sorted = np.load(f'database_results/{cluster_name[0]}/z_sorted_{cluster_name[0]}.npy')
# r_total_sorted = np.load(f'database_results/{cluster_name[0]}/r_sorted_{cluster_name[0]}.npy')


# In[ ]:


# #Analytic potential for a MOND homogeneous sphere. Adding MOND potential at boundary for the offset
# potential_sphere_MOND = (np.heaviside(r_total_sorted - radius_tot, 0.5)*sqrt(G*mgb*a0)*np.log(r_total_sorted) +
# (np.heaviside(radius_tot - r_total_sorted, 0.5))*(4/3*sqrt(pi/3*a0*G*mgb/volume_out)*np.power(r_total_sorted,3/2)+
# sqrt(G*mgb*a0)*ln(radius_tot) - 4/3*sqrt(pi/3*a0*G*mgb/volume_out)*radius_tot**(3/2)))

# #Potential for a homogeneous sphere in Newton
# potential_sphere_Newton = (np.heaviside(r_total_sorted - radius_tot, 0.5)*(-G*mgb/r_total_sorted) +
# (np.heaviside(radius_tot - r_total_sorted, 0.5))*G*mgb/(2*radius_tot**3)*(r_total_sorted**2-
# 3*radius_tot**2)+sqrt(G*mgb*a0)*ln(domain_size))


# plt.figure()

# plt.plot(r_total_sorted, potential_total_sorted)
# plt.plot(r_total_sorted, potential_sphere_Newton, linestyle = '--')


# In[ ]:


# plt.figure()

# plt.plot(r_total_sorted, source_total_sorted)


# In[1]:


#array to check all quantities are finite
# is_finite = np.zeros((8,))

# #Checking if all arrays have finite values only
# for i, element in enumerate([potential_total_sorted, source_total_sorted, apparent_mass_total_sorted,
#                 dark_mass_total_sorted, x_total_sorted, y_total_sorted, z_total_sorted,r_total_sorted]):
    
#     is_finite[i] = np.isfinite(element).any()
    
# print(f'The finite arrays: {is_finite}')


# In[2]:




