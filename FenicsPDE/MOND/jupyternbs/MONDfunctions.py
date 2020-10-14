from dolfin import *

#importing mshr for all mesh functions
import mshr as mshr

# Use SymPy to compute f from the manufactured solution u
import sympy as sym

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

#Needed if want to use the adapt function for mesh refinement, see:
#https://fenicsproject.org/qa/6719/using-adapt-on-a-meshfunction-looking-for-a-working-example/
parameters["refinement_algorithm"] = "plaza_with_parent_facets"
#Importing the quantities files since there are a bunch of default values in these functions 

from MONDquantities import *

from MONDclasses import *

# # Function to define the initial mesh

# In[3]:

def make_spherical_mesh(size, resolution):

    #Defining the domain for the mesh using the Sphere function from mshr
    domain = mshr.Sphere(origin, size)

    #Meshing the sphere generated with resolution 10 (cells per dimension)
    initial_mesh = mshr.generate_mesh(domain, resolution)
    
    return initial_mesh

def inside_circle (mesh, radius):
    '''Function to mark all elements that have their center farther than a certain
    radius.'''
    
    #Defining the origin based on the geometric dimension of the mesh (plane 2 volume 3)
    origin = Point(np.zeros((mesh.geometric_dimension(),1)))
        
    #List comprehension with a single True corresponding to the cell containing the point
    is_in_list = [cell.distance(origin) < radius for cell in cells(mesh)]
    
    inside_cell_index = np.nonzero(is_in_list)
    
    #Converting list to a np array so it's faster (and has 0 and 1 so can be sorted)
    is_in_numpy = np.fromiter(is_in_list, float, mesh.num_cells())
    
    #Indices at which we have non-zero in the numpy array are the cell indices we need
    is_outside_index = np.nonzero(is_in_numpy)

    return inside_cell_index


# ## Trying to create mesh from scratch using the MeshEditor

# ## Next, using the inside_circle function to exclude cells outside a given radius

# In[6]:


def circle_grid_mesh(mesh, radius, show_mesh = True):
    
    #Getting the MeshFunction and numpy array containing whether each cell satisfies the function
    inside_cell_index = inside_circle(mesh, radius)

    #Getting the vertex indices only for the cells that are inside the required radius
    inside_cells_vertices_old = mesh.cells()[inside_cell_index]

    #Number of remaining cells, given by the number of rows in the array (shape[0])
    inside_cell_number = inside_cells_vertices_old.shape[0]

    #Finding the indices of all the vertices included in the remaining cells (only need each index
    #once, hence we use unique)
    inside_vertex_old = np.unique(inside_cells_vertices_old)

    #Getting the coordinates of each vertex of cells contained in the radius
    inside_vertex_coordinates = mesh.coordinates()[inside_vertex_old]

    #Number of remaining vertices (same way as for cells)
    inside_vertex_number = inside_vertex_old.shape[0]

    #New indices for the vertices, fom 0 to #vertices
    inside_vertex_new = np.arange(inside_vertex_number)

    #Replacing the old vector indices by the new vector indices in the cell array
    #Approach from StackExchange: https://stackoverflow.com/questions/55949809/
    #efficiently-replace-elements-in-array-based-on-dictionary-numpy-python
    #Making np. zero array with as many elements as the initial # vertices
    mapping_old_new = np.zeros(inside_vertex_old.max()+1,dtype=inside_vertex_new.dtype)

    #Replacing the zero at the position of each old vertex, with the corresponding new vertex!
    mapping_old_new[inside_vertex_old] = inside_vertex_new

    #Now indexing each element of the mapping by an element of the cells_vertices array. So each
    #old vertex is replaced by a new vertex! Super efficient and fast
    inside_cells_vertices_new = mapping_old_new[inside_cells_vertices_old]
    
    #Now we have a list of cells and vertices that we want to put into the mesh
    #Declaring an empty object from the MeshEditor class
    editor = MeshEditor()

    #Declaring an empty mesh to use in the mesh editor
    test_mesh = Mesh()

    #Opening the line mesh in the mesh editor, with given topological and geometric dimensions
    editor.open(test_mesh, mesh.cell_name(), mesh.topology().dim(), mesh.geometric_dimension())

    #Initialising the mesh with the # vertices corresponding to the cells inside the radius
    editor.init_vertices(inside_vertex_number)

    #Initialising an amount of cells equal to those inside the given radius
    editor.init_cells(inside_cell_number)

    #Adding all the vertices from cells inside radius
    for vertex_index, vertex in enumerate(inside_vertex_coordinates):
        editor.add_vertex(vertex_index, vertex)

    #Giving list of cells with the new indexing
    for cell_index, cell_vertices in enumerate(inside_cells_vertices_new):
        editor.add_cell(cell_index, cell_vertices)

    #Closing the mesh editor
    editor.close()
    
    #Plotting mesh if necessary
    if show_mesh:

        plt.figure()
        plot(test_mesh, color = 'w')
        
    #Returning the mesh object
    return test_mesh


# In[7]:


def make_cube_mesh(side_length, resolution):
    
    #Defining extremes of main diagonal to obtain side_length*3 cube
    box_edge_low = Point(-side_length*np.array([1,1,1]))
    box_edge_high = Point(side_length*np.array([1,1,1]))
    
    #Using dolfin builtin mesh, defining edges and resolution for all 3 axes
    cubic_mesh = BoxMesh(box_edge_low, box_edge_high, resolution, resolution, resolution)
    
    #Object returned by BoxMesh is not a dolfin mesh, so we define it as such to be able to use it
    cubic_mesh = Mesh(cubic_mesh)
    
    return cubic_mesh

def more_modified_refinement (mesh, location, how_many):
    '''Function to refine mesh locally, only refining cells containing one of the points described by the
    list of Dolfin Points location,a Point for each source. 
    '''
    
    #Starting # cells before we refine, to compute growth factor
    starting_cells = mesh.num_cells()
    
    #How_many gives the amount of time the mesh should be refined
    if how_many > 0:
    
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
            intersect_list = [intersect(mesh, source).intersected_cells() for source in location if mesh.bounding_box_tree().compute_first_entity_collision(source) < mesh.num_cells()]
            
            #Setting the cell function contain_function to true for each cell containing a source
            for cell_index in intersect_list:
                
                contain_function[cell_index[0]] = True
            
            #Refining the mesh only for cells that contain a source
            mesh = refine(mesh, contain_function)    
            
            #Final # cells after refinement
            final_cells = mesh.num_cells()
            
            partial_growth_factor = final_cells/initial_cells
            
            print(('Iteration {} of {}: The Cell number went from {} to {}, up by a factor {}\n'
                  .format(i+1, how_many, initial_cells, final_cells, partial_growth_factor)))
    
        #ratio between # cells at beginning and end of refinement
        total_growth_factor = final_cells/starting_cells

        print('Cell number went up by a factor {}\n'.format(total_growth_factor))
    
    #returning the refined mesh
    return mesh


# # Using local mesh refinement after implementing a mesh function to mark individual cells to be refined based on some condition (page 187 of book)

# In[14]:


def local_refinement (mesh, location, radius, how_many, technique = 'inside'):
    '''Function to refine mesh locally, based on the distance from a given point
    '''
    
    #Starting # cells before we refine, to compute growth factor
    starting_cells = mesh.num_cells()
    
    if how_many > 0:
    
        for i in range(how_many):
        
            for j, source in enumerate(location): 
                
                print(f'Source {j+1} of {source_number}', end="\r", flush=True)
                
                #IMPORTANT: Need to delcare (and set to False) the cell_to_refine function
                #for each source, as the mesh is different after each iteration!
                #Initial mesh cell count
                initial_cells = mesh.num_cells()

                #We want to refine based on cell location, so we want to mark cells within a given radius
                #We create a MeshFunction containing a bool for each cell (3 stands for the topological
                #dimension, which is vertex->1, face->2, cell->3), on the given mesh
                cell_to_refine = MeshFunction("bool", mesh, 3)

                #Initialising all marker to false
                cell_to_refine.set_all(False)

                #cells(mesh_trial) is an iterator, so we can loop over it directly, without needing to know
                #anything about the cells
                for k, cell in enumerate(cells(mesh)):
                    
                    if (k+1)%1000 == 0:
                        print(f'Cells done: {int(k+1)} out of {mesh.num_cells()}', end="\r", flush=True)
                    
                    #With the close technique, we refine all cells with centers within a given radius
                    if str(technique) == 'close':

                        #looking at the center of the cell
                        p = cell.midpoint()

                        #if the cell is within the required radius, we set the corresponding marker to True
                        if p.distance(source) < radius:

                            cell_to_refine[cell] = True
                        
                    elif str(technique) == 'ring':

                        #looking at the center of the cell
                        p = cell.midpoint()

                        #if the cell is within the required radius, we set the corresponding marker to True
                        if (radius*0.9 < p.distance(source) and  (p.distance(source) < radius*1.1)):

                            cell_to_refine[cell] = True    
                    
                    elif technique == 'inside': 

                        if cell.contains(source):

                            cell_to_refine[cell] = True

                #Refining the mesh only where the markers are True, so inside the desired radius 
                mesh = refine(mesh, cell_to_refine)

                final_cells = mesh.num_cells()

                partial_growth_factor = final_cells/initial_cells
            
            print(('Iteration {} of {}: The Cell number went from {} to {}, up by a factor {}\n'
                  .format(i+1, how_many, initial_cells, final_cells, partial_growth_factor)))
    
        #ratio between # cells at beginning and end of refinement
        total_growth_factor = final_cells/starting_cells

        print('Cell number went up by a factor {}\n'.format(total_growth_factor))
    
    #returning the refined mesh
    return mesh


def rearrange_mesh_data(mesh, C_O_M = origin, degree = 1):
    V = FunctionSpace(mesh, 'CG', degree)

    #Finding the total number of vertices in the mesh
    vertex_number = mesh.num_vertices()

    #Storing x, y, z coordinates of each point of the mesh in a separate numpy array
    x_coords = mesh.coordinates()[:,0]
    y_coords = mesh.coordinates()[:,1]
    z_coords = mesh.coordinates()[:,2]

    #using the numpy linalg.norm function to get the radius(norm) of each vertex
    r_coords = np.linalg.norm(mesh.coordinates(), axis=1)
    
    r_coords_CM = (np.sqrt(np.square(x_coords-C_O_M [0]) +
    np.square(y_coords-C_O_M [1]) + np.square(z_coords-C_O_M [2])))
    
    #Storing the index to sort according to r
    sorting_index = r_coords_CM.argsort()

    #Sorted radial distances
    r_sorted = r_coords_CM[sorting_index]

    #Sorted coordinates
    x_sorted = x_coords[sorting_index]
    y_sorted = y_coords[sorting_index]
    z_sorted = z_coords[sorting_index]

    return V, vertex_number, x_coords, y_coords, z_coords, r_coords, sorting_index, x_sorted, y_sorted, z_sorted, r_sorted

# # Function to make string for discrete delta peak distribution of source_number #masses.

def make_discrete_dirac(how_many_sources):
    '''Function to generate a c++ string giving a discrete distribution of delta peaks (modelled as
    spheres of uniform mass with a very small radius) at arbitrary points. The number of sources is
    given by how_many_sources. The positions have to be specified explicitly, one by one, which is
    achieved by the make_source_string function'''

    #Starting the string with an open bracket so we can have the whole conditional statement between
    #brackets
    discrete_string = '('

    #Looping over the nummber of sources to generate a conditional statement for each
    for i in range(how_many_sources):

        #For all loop iterations apart from the last one, we add the condition for a sphere
        #at the respective location, and add an OR(||) at the end to chain it to the next
        if i != how_many_sources - 1:

            discrete_string = discrete_string + ('((pow((x[0]-position_x_{}),2)+pow((x[1]-position_y_{}),2)'
            '+pow((x[2]-position_z_{}),2)) <=pow(radius_tot, 2)) ||'.format(i,i,i))

        #If it's the last iteration, we don't add an OR but a conditional delimiter (?), and add the
        #body of the loop containing if/not the condition is satisfied
        else:

            discrete_string = discrete_string + ('((pow((x[0]-position_x_{}),2)+pow((x[1]-position_y_{}),2)'
            '+pow((x[2]-position_z_{}),2)) <=pow(radius_tot, 2)))? '.format(i,i,i))

            discrete_string = discrete_string + '4*pi*G*source_mass/volume_out : 0'
        
    return discrete_string

# # Function to make a discrete distribution, but using Gaussians instead of delta peaks.

# In[20]:


def make_discrete_gauss(how_many_sources):
    '''Function to generate a c++ string giving a discrete distribution of delta peaks (modelled as
    spheres of uniform mass with a very small radius) at arbitrary points. The number of sources is
    given by how_many_sources. The positions have to be specified explicitly, one by one, which is
    achieved by the make_source_string function'''

    #Starting the string with an open bracket so we can have the whole conditional statement between
    #brackets
    discrete_string = ('(4*pi*G*pow(2*pi,-1.5)*source_mass/pow(stand_dev_peak,3)*(')

    #Looping over the nummber of sources to generate a conditional statement for each
    for i in range(how_many_sources):

        #For all loop iterations apart from the last one, we add the condition for a sphere
        #at the respective location, and add an OR(||) at the end to chain it to the next
        if i != how_many_sources - 1:

            discrete_string = discrete_string + ('exp(-(pow((x[0]-position_x_{}),2) +'
            ' pow((x[1]-position_y_{}),2) +'
            ' pow((x[2]-position_z_{}),2))/(2*(pow(stand_dev_peak,2))))+'.format(i,i,i))

        #If it's the last iteration, we don't add an OR but a conditional delimiter (?), and add the
        #body of the loop containing if/not the condition is satisfied
        else:

            discrete_string = discrete_string + ('exp(-(pow((x[0]-position_x_{}),2) +'
            ' pow((x[1]-position_y_{}),2) +'
            ' pow((x[2]-position_z_{}),2))/(2*(pow(stand_dev_peak,2))))))'.format(i,i,i))
        
    return discrete_string

# # Function to take a c++ string defining multiple discrete masses and defining all quantities inside the Expression, and all source positions one by one.

# In[27]:


def make_source_string(how_many_sources):
    '''Function to take an initial c++ expression, add all the required constants, and then adding
    all the coordinates for each source we want to add from the random_coordinates arrays'''
    
    #Defining the initial string that we want to add the variable source location definitions to
    executable_string = ('(Expression(source,'
    'degree = degree_PDE, a0 = a0, ms = ms,mgb = mgb, G = G,  ly = ly, kp = kp, radius_tot = radius_tot,'
    'volume_out = volume_out, stand_dev = stand_dev, stand_dev_peak=stand_dev_peak, p=p,'
    'source_number = source_number, source_mass = source_mass')

    #Looping over each source and generating a string assigning the value of the source position
    #in x,y,z, directions to the corresponding variable
    for i in range(how_many_sources):

        executable_string = (executable_string + ',position_x_{} = random_coordinates_x[{}],'
        'position_y_{} = random_coordinates_y[{}], position_z_{} = random_coordinates_z[{}]'
        .format(i,i,i,i,i,i))

    #Adding two closing brackets to close the expression
    executable_string = executable_string+'))'
    
#     executable_string = eval(executable_string)
    
    return executable_string


def solve_PDE(the_BVP):
    '''Function takes in a BVP object, which defines the weak form, initial guess/BC and source for
    a PDE, and computes its solution'''
    
    ## starting time of PDE solver
    solver_start = time.time()
    print('Starting PDE Solver...\n')

    #defining the x,y,z coordinates from the coordinate array in sympy
#     x, y, z = sym.symbols('x[0], x[1], x[2]')

    #VERY IMPORTANT: If using sympy, use sym. in front of every mathematical operator, or the sym. and UFL (the
    #mathematical syntax used in fenics) collide and an error about UFL conversion appears
    
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


def plot_format(plot_to_format, dimension, radial, scale = 6):
    '''Format the 1D, 2D and 3D plots. Takes as input the dimension of the graph (1D,2D,3D),
    whether the plot is radial in 1D to set the x axis limits, and the scale of the system,
    indicating the power of parsec: 6 is Mpc, 3 is Kpc and so on.'''
    
    #Setting a different unit depending on the scale, to be added to the corresponding axis
    if scale == 6:
        dim = ' (Mpc)'
    else:
        dim = ' (kpc)'
    
    plt.legend()
    plt.tight_layout()
    #The grid can be given for a subplot, all the other things must be given as the plot itself
    #So the input to the function can actually be a subplot now
    plot_to_format.grid()
    #*args gives an arbitrary # of inputs. we take the 1st one to define if the plot is radial
    #if it is (=1), we only need to plot from 0, not -domainsize
    if radial == 0:
        plt.xlabel('x'+ dim)
        (plt.xticks(np.linspace(-domain_size,domain_size,11), 1/(10**(scale-3)*kp)
                    *np.linspace(-domain_size,domain_size,11)))
    else:
        plt.xlabel('r'+ dim)
        plt.xticks(np.linspace(0,domain_size,11), 1/(10**(scale-3)*kp)*np.linspace(0,domain_size,11))
    if dimension == 2:
        plt.ylabel('y'+ dim)
        plt.yticks(np.linspace(-domain_size,domain_size,11),1/(10**(scale-3)*kp)
                    *np.linspace(-domain_size,domain_size,11))
    elif dimension == 3:
#         plt.zlabel('z'+ dim)
        plt.zticks(np.linspace(0,domain_size,11),1/(10**(scale-3)*kp)
                    *np.linspace(-domain_size,domain_size,11))
    else:
        None


# ### Defining a function to add useful annotations to radial graphs

# In[34]:


def plot_annotations(plot_to_annotate):
    #Adding a vertical line to denote the end of the mass distribution
    plt.axvline(x = radius_tot, color = 'k', linestyle = '-.', linewidth = 1, label = 'Mass Density')

    #and another to signal where the transition between Newton and MOND occurs
    plt.axvline(x = sqrt(G*mgb/a0), color ='k',linestyle ='--',linewidth = 1,label = 'Transition')
    
  
def line_integral(integrand, A, B, n, show = False): 
    '''Integrate u over segment [A, B] partitioned into n elements'''
    
    #This code comes from the Fenics Forum. Changed it to suit what I need to do. Link:
    #https://fenicsproject.org/qa/13863/integrate-along-axes-lines/
    #Assert statements are used to check expressions are of the correct form. If the assert 
    #statement would return False, then the program throws an error
    #value_rank has to be 0 because we can only integrate over scalar Function Spaces!
    assert u.value_rank() == 0
    assert len(A) == len(B) > 1 and np.linalg.norm(A-B) > 0
    assert n > 0
    
    #Defining an np array starting at A, getting to B in n steps
    mesh_points = [A + t*(B-A) for t in np.linspace(0, 1, n+1)]
    
    #Uncomment to check the mesh points are the correct ones
#     print(f'The mesh points are {mesh_points}\n')
    
    #the topological and geometric dimensions of the points(topological is always one for a
    #point), geometric depends on the #coordinates, here will be 3 as we're in 3D
    tdim, gdim = 1, len(A)
    
    #Declaring an empty mesh
    line_mesh = Mesh()
    
    #Declaring an empty object from the MeshEditor class
    editor = MeshEditor()
    
    #Opening the line mesh in the mesh editor, with given topological and geometric dimensions
    #Had too add 'interval' for the degree, as it is required in the newest version
    editor.open(line_mesh, 'interval', tdim, gdim)
    
    #Initialising the line mesh with the #vertices declared in the function definition
    editor.init_vertices(n+1)
    
    #Initialising n cells (1 less than the vertices, one cell for every two adjacent vertices)
    editor.init_cells(n)
    
    #Adding the declared vertices to the MeshEditor
    for vi, v in enumerate(mesh_points): editor.add_vertex(vi, v)
    
    #Declaring each cell to be made of the two current adjacent vertices
    for ci in range(n): editor.add_cell(ci, np.array([ci, ci+1], dtype='uintp'))
    
    #Closing the mesh editor
    editor.close()

    #Setting up the function space from the function we want to plot
    elm = integrand.function_space().ufl_element()
    
    #Finding family and degree properties of the potential ('P' and 1 usually)
    family = elm.family()
    degree = elm.degree()
    
    #Defining a function space for the function on the 1-D line mesh defined above
    V = FunctionSpace(line_mesh, family, degree)
    
    #Interpolating the potential over the function space for the 1-D line mesh
    v = interpolate(integrand, V)
    
    #Performing the integration on the defined line mesh
    integral = assemble(v*dx)
    
    #Set show = True to check the line is being plotted correctly
    
    if show:
        
        #Adding figure 
        figure = plt.figure()
        
        #projection='3d' needed to specify a 3D scatter plot
        integral_scatter = figure.add_subplot(111, projection='3d')
        
        #Plotting the mesh
        plot(mesh, alpha = 0.5, color = 'w')
        
        #Scattering the line_mesh points as blue crosses
        integral_scatter.scatter(mesh_points, 0, 0, marker = '+', color = 'b', s = 0.5)
        
        #Scattering the sources
        integral_scatter.scatter(random_coordinates[:,0],random_coordinates[:,1], 
        random_coordinates[:,2], marker = 'o', color = 'r', s=1)
    
    #Returning the integral of the interpolated potential over the line mesh declared
    return integral


def sum_individual_contributions(mesh, origin, coordinates):
    '''Summing individual contributions of each source linearly'''
    
    #Initialising the sum for the potentials to be 0, then add each contribution to it from each source
    potential_off_center = 0
    
    #Setting the distances to be from the origin so we have the same coordinates for all individual potential
    #and we can coherently add them together
    V, vertex_number, x_coords, y_coords, z_coords, r_coords, sorting_index, x_sorted, y_sorted, z_sorted, r_sorted = rearrange_mesh_data(mesh, origin)
    
    #Assigning each coordinate from the array given as input
    coordinates_x = coordinates[:,0]
    coordinates_y = coordinates[:,1]
    coordinates_z = coordinates[:,2]
    
    #Summing over each source
    for i in range(len(coordinates_x)):
    
        #Defining the potential with respect to a displaced source, but with the usual coordinates
        potential_off_center += (sqrt(G*source_mass*a0)*1/2*np.log((x_coords-coordinates_x[i])**2+
                    (y_coords-coordinates_y[i])**2 + (z_coords-coordinates_z[i])**2))
    
    return potential_off_center

def radial_dist_hist(r_coords, mesh, density, bins):
    
    #Having independent y axes so we don't need to adjust the function to compare the scaling 
    fig, histo = plt.subplots()
    color = 'tab:blue'
    histo.set_ylabel('radial distance density', color=color)
    histo.tick_params(axis='y', labelcolor=color)
    
    #Looking at how the points are distributed radially. If they are uniform, their density
    #should increase with r^2 due to sphere surface
    #plotting histogram of point density radially
    to_hist = histo.hist(r_coords, density=density, bins=bins, label = 'distribution', color = color)
    
    #Adding the value of each bin on top of the bin
    for i in range(bins):
        plt.text(to_hist[1][i],to_hist[0][i],str(int(to_hist[0][i])))
    
    #Putting the grid with the histogram values as those are the interesting quantities
    plt.grid()
    
    #Second y axis
    quadratic = histo.twinx()
    color = 'tab:orange'
    quadratic.set_ylabel('Quadratic', color=color)
    quadratic.tick_params(axis='density', labelcolor=color)
    
    #plotting a cubic relation, scaled by the max element^3 to be of order unity
    quadratic.plot(r_coords, np.power(r_coords,2), label = 'quadratic', color = color)
    
    #Forcing the lower limit in both plots to 0 so there's no offset
    plt.ylim(bottom=0)
    
    plt.tight_layout()
    plt.legend()
    
   
def plot_mesh(mesh, center_of_mass, degree, figure, vertex_fraction, color, show_mesh = False, alpha = 0.3):
    '''
    Plotting the first vertex_number/vertex_fraction points of the mesh, based on radial
    distance. Need to input a figure name of format figure = plt.figure() for the subplot
    to work. Through this, we can embed this plot into other plots. Adding optional mesh and
    alpha inputs.
    '''
    
    V, vertex_number, x_coords, y_coords, z_coords, r_coords, sorting_index, x_sorted, y_sorted, z_sorted, r_sorted = rearrange_mesh_data(mesh, center_of_mass, degree)

    
    vertex_number = mesh.num_vertices()
    
    #Getting the number of vertices required as an intege of the fraction given
    how_many = int(vertex_number/vertex_fraction)
    
    #projection='3d' needed to specify a 3D scatter plot
    mesh_scatter = figure.add_subplot(111, projection='3d')
    
    #plotting the total/vertex_fraction closest vertices to the origin
    #s gives the size of the dots, multiplying it by vertex_fraction so when we have less dots
    #we can make them more visible
    (mesh_scatter.scatter(x_sorted[0:how_many], y_sorted[0:how_many], z_sorted[0:how_many],
                         marker = '.', s=1*vertex_fraction, c = color))
    
    #Adding optional to show the mesh at the same time
    if show_mesh:
        plot(mesh, alpha = alpha, color = 'w')
    else:
        pass
    
    
def slice_mesh(amount, mesh, center_of_mass, degree, height = 0, portion = True, values = True, *args):
    '''Selecting only points of the mesh that are close to the xy plane. First, sorting points
    according to their z coordinate, then selecting either amount # of them or a portion of the
    total, dpeending on 'portion'. Optionally output a function args[0] at those points
    '''
    
    
    V, vertex_number, x_coords, y_coords, z_coords, r_coords, sorting_index, x_sorted, y_sorted, z_sorted, r_sorted = rearrange_mesh_data(mesh, center_of_mass, degree)
    
    #Getting abs(z), distance from xy plane, and getting the index that sorts this new array
    xy_plane_distance = np.abs(z_coords - height)
    
    #index of points sorted by distance from xy plane (absolute value of z coordinate)
    z_sorting_index = xy_plane_distance.argsort()
    
    #x,y,z coordinates of point sorted by distance from xy plane
    z_xy_plane = z_coords[z_sorting_index]
    x_xy_plane = x_coords[z_sorting_index]
    y_xy_plane = y_coords[z_sorting_index]
    
    #If using the portion option, we take all points below a cutoff given by domain_size/amount
    if portion == True:
        
        #counting number of points with xy_plane distance > domain_size/amount
        amount = np.count_nonzero(xy_plane_distance < int(domain_size/amount))
    
    #If not using the portion option, the slice has num_vertices/amount total # points
    else:
        
        amount = int(mesh.num_vertices()/amount)
                
    #The slice has amount # points
    x_xy_plane = x_xy_plane[0:amount]
    y_xy_plane = y_xy_plane[0:amount]
    z_xy_plane = z_xy_plane[0:amount]

    if values == True:
        function_xy_plane = args[0][z_sorting_index]
        function_xy_plane = function_xy_plane[0:amount]
        
        return x_xy_plane, y_xy_plane, z_xy_plane, function_xy_plane, amount
    
    else:
        return x_xy_plane, y_xy_plane, z_xy_plane, amount

    
def plot_mesh_slice(amount, figure, mesh, center_of_mass, degree, source_coordinates, height = 0, portion = True, show_mesh = True, cross_section = False, alpha = 0.3, show_source = True):
    '''
    Plotting a slice of points on the mesh close to the xy axis, so for small z.
    Need to create a figure before calling the function. Optional showing mesh on top
    of scatter, True by default, and transparency of mesh set by default to 0.3
    '''
    
    #Calling the slice_mesh function to slice the mesh before plotting
    x_xy_plane, y_xy_plane, z_xy_plane, amount = (slice_mesh(amount, mesh, center_of_mass, degree, height = height, portion = portion, values = False))
    
    #projection='3d' needed to specify a 3D scatter plot
    mesh_xy_axis = figure.add_subplot(111, projection='3d')
    
    vertex_number = mesh.num_vertices()
    
    #plotting the total/vertex_fraction closest vertices to the origin
    #s gives the size of the dots, multiplying it by vertex_fraction so when we have less dots
    #we can make them more visible. Having vertex_number/(amount+1) cause on MPI it could be that
    #a given mesh slice has amount=0!
    (mesh_xy_axis.scatter(x_xy_plane, y_xy_plane, z_xy_plane,
                         marker = '.', s=vertex_number/(amount+1), c = z_xy_plane))
    
    #Plotting the source(s) on top of the mesh points, to check if we refine mesh correctly
    if show_source == True:
        
        #Again adding 1 to amount when we divide by it or might divide by 0 when using MPI if 
        #a submesh has amount = 0
        (mesh_xy_axis.scatter(source_coordinates[:,0],source_coordinates[:,1],
        source_coordinates[:,2], marker = 'o', s=vertex_number/((amount+1)*5), c = 'r', label='Source'))
            
    
    #Projecting each point on the xy plane, at z=0 (zs=0)
    mesh_xy_axis.plot(x_xy_plane, y_xy_plane, 'b+', zdir='z', zs=height, label = 'xy-plane projection')
    
    mesh_xy_axis.legend()
    
    #Adding optional to show the mesh at the same time
    if show_mesh:
        plot(mesh, alpha = alpha, color = 'w')
    else:
        pass
    
    #Option to plot the cross section, color coded by distance from xy plane
    if cross_section:
        plt.figure()
        plt.scatter(x_xy_plane, y_xy_plane, c = z_xy_plane, marker = '+')
        
        if show_source:
            (plt.scatter(source_coordinates[:,0],source_coordinates[:,1], marker = 'o',
            color = 'r', s=vertex_number/((amount+1)*5)))
        
    else:
        pass

    
# ## Defining a function to plot a trisurf graph of a slice of the domain, here the xy axis

# In[57]:


def trisurf_function_slice(figure, function, amount, height, mesh, center_of_mass, degree, slices = 50, high_low = 'low', project = True):
    '''Plot a trisurf along a plane (currently only the xy-plane) of a given function. The slice
    of points defining the xy plane is obtained from the slice_mesh function. Then, project
    the triangulated surface onto each axis if optional arg project = True. high_low determines
    whether to project the xy plane contour above or below surface
    '''
    
    #Obtaining the x,y coordinates and functions to plot from slice_mesh
    x_xy_plane, y_xy_plane, z_xy_plane, function_xy_plane, amount = (slice_mesh(amount,
                              mesh, center_of_mass, degree,  height, True, True, function))
    
    #Adding a subplot to the figure input
    plane_trisurf = figure.add_subplot(111, projection='3d')
    
    #Plotting the trisurf (triangulated surface) on the plane
    #Only do this for x,y,z >3 or throws error. This is necessary as when using MPI certain submeshes
    #might not have the needed # elements. They're all the same length, so checking one is sufficient
    if (len(x_xy_plane) > 3):
        plane_trisurf.plot_trisurf(x_xy_plane, y_xy_plane, function_xy_plane, cmap = 'jet')

        #Setting the x and y limits to leave some space for the contours to be clearer
        x_bottom, x_top = plt.xlim(-domain_size*1.2, domain_size*1.2) 
        y_bottom, y_top = plt.ylim(-domain_size*1.2, domain_size*1.2) 

        #Setting the z limit conditionally, based on the function being increasing or decreasing
        if high_low == 'high':
            z_limit = function.max()

        else:
            z_limit = function.min()

        #If project = True, project contours on each axis
        if project:

            #Projecting contours of the surface in each cartesian direction using zdir and offset
            #for direction and offset from the axes origins respectively
            plane_trisurf.tricontourf(x_xy_plane, y_xy_plane, function_xy_plane, slices, zdir='z', offset=z_limit, cmap = 'jet')
            plane_trisurf.tricontourf(x_xy_plane, y_xy_plane, function_xy_plane, slices, zdir='x', offset=y_bottom, cmap = 'jet')
            plane_trisurf.tricontourf(x_xy_plane, y_xy_plane, function_xy_plane, slices, zdir='y', offset=x_top, cmap = 'jet')

        
# ## Plotting contour lines of the potential, so we can do that for different values of z and see the whole domain.

# In[61]:


def tricontour_function_slice(figure, function, center_of_mass, amount, levels, height, mesh, degree):
    '''Plot a trisurf along a plane (currently only the xy-plane) of a given function. The slice
    of points defining the xy plane is obtained from the slice_mesh function
    '''
    
    #Obtaining the x,y coordinates and functions to plot from slice_mesh
    x_xy_plane, y_xy_plane, z_xy_plane, function_xy_plane, amount = (slice_mesh(amount,
                              mesh, center_of_mass, degree, height, True, True, function))
    
    #Adding a subplot to the figure input. Not adding the projection = 3d option so we have it 
    #all in one plane and can plot multiple for e.g. different values of z
    plane_tricontour = figure.add_subplot(111)
    
    #Defining the levels of the contour to be between the max and min value of the function
    #in 'levels' homogeneous steps
    contour_levels = np.linspace(function.min(), function.max(), levels)
    
    #Plotting the contour
    #Only do this for x,y,z >3 or throws error. This is necessary as when using MPI certain submeshes
    #might not have the needed # elements. They're all the same length, so checking one is sufficient
    if (len(x_xy_plane) > 3):
        plane_tricontour.tricontour(x_xy_plane, y_xy_plane, function_xy_plane, contour_levels, cmap = 'jet')
        plane_tricontour.set_aspect('equal', 'box')
    
    
def contour_3D_slices(figure, function, amount, height, slices = 50, high_low = 'low'):
    
    #Adding a subplot to the figure input
    plane_trisurf = figure.add_subplot(111, projection='3d')

    #Setting the x and y limits to leave some space for the contours to be clearer
    x_bottom, x_top = plt.xlim(-domain_size*1.2, domain_size*1.2) 
    y_bottom, y_top = plt.ylim(-domain_size*1.2, domain_size*1.2) 
#     z_bottom, z_top = Axes3D.zlim(-domain_size*1.2, domain_size*1.2) 
    
    #Setting the z limit conditionally, based on the function being increasing or decreasing
    if high_low == 'high':
        z_limit = function.max()
    
    else:
        z_limit = function.min()
    
    #Obtaining a range of n values for both the function and the domain (from min to max)
    function_to_loop = np.linspace(function.min(), function.max(), 11)
    height_to_loop = np.linspace(-domain_size*0.9, domain_size*0.9, 11)
    
    #Stacking the two ranges together so we can loop over them, and transpose so that each
    #element of the array is a pair of corresponding function value-domain location
    total_to_loop = np.vstack((function_to_loop, height_to_loop))
    total_to_loop = np.transpose(total_to_loop)
    
    #Looping over the function-domain slice pairs
    for elevations in total_to_loop:
        
        #IMPORTANT! For this loop to work, we need to have the same #function points for each
        #slice. So we need to set the portion option from slice_mesh to False, so it takes
        #the same #points every time, and not whatever number is within the given range
        #Obtaining the x,y coordinates and functions to plot from slice_mesh
        x_xy_plane, y_xy_plane, z_xy_plane, function_xy_plane, amount = (slice_mesh(amount,
                                                mesh, elevations[1], False, True, function))

        #Projecting contours of the surface in each cartesian direction using zdir and offset
        #for direction and offset from the axes origins respectively
        plane_trisurf.tricontourf(x_xy_plane, y_xy_plane, function_xy_plane, slices, zdir='z', offset=elevations[0],cmap = 'jet')
    #     plane_trisurf.tricontourf(x_xy_plane, y_xy_plane, function_xy_plane, slices, zdir='x', offset=y_bottom, cmap = 'jet')
    #     plane_trisurf.tricontourf(x_xy_plane, y_xy_plane, function_xy_plane, slices, zdir='y', offset=x_top, cmap = 'jet')