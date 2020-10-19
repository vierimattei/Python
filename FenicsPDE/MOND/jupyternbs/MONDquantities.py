from dolfin import *

#importing mshr for all mesh functions
import mshr as mshr

c = 2.998*10**8
G = 6.674*10**(-11) 
a0 = 1.2*10**(-10) 
#Mass of the sun
ms = 1.989*10**30 
mgd = 10**12*ms 
ly = 9.461*10**15 
kp = 3261.56*ly 
kpm = 20*kp 
au = 0.0000158*ly 
lst = 200*au 
lsu = 50*au 
lg = 52850*ly 
lgu = 40*10**3*ly 
lgb = 3*kp 
rs = 0.00465*au 
#Mass of the milky way
mgb = 1.1*6.6*10**10*ms
#Mass of the virgo cluster (inlcudes dark matter)
mass_virgo = 10**15*ms
#Mass of the coma cluster (without dark matter, from the Brownstein 2006 paper). This is not the same 
#as the gas mass calculated Newtonianly! For that, we need to take the 3-beta model and use its 
#distribution directly!
mass_coma_MOND = 3.81*10**14*ms
#Gas mass in the coma cluster
mass_coma_gas = 1.13*(10**14)*ms
#Radius of the coma cluster
radius_coma = 1954*kp

#Here replacing the mass of the galaxy with that of virgo, since we're working on the cluster scale
mgb = mass_coma_gas

#Setting domain_size to the radius of the Virgo Cluster, 2.3 Mpc
domain_size = radius_coma
#Smallest size of a galaxy in Virgo has a radius = r_Virgo/250, so that's the resolution we need
radius_tot = domain_size/250
#Origin at (0,0,0) for the mesh
origin = Point(0,0,0)
radius_refine = radius_tot
volume_out = 4/3*pi*(radius_tot**3)
#Standard deviation for coarse mass distribution and location of peaks. domain_size/3 corresponds to have
#99.7% of the mass distribution inside the domain. Dividing that by some other factor makes it less 
#likely that the masses will be outside the domain and throw an error
stand_dev = domain_size/3/2
#Standard deviation for the gaussian peaks themselves, Radius tot/3 so 99.7% of the mass is inside the 
#equivalent dirac delta made with a uniform sphere
stand_dev_peak = radius_tot/3*100
mesh_resolution = 21
#Coefficient for GEA changing the potential based on how spherically symmetric the mass distribution is
c_2 = -1.8
#Coefficient for GEA giving the magnitude of the K^3/2 term in the Lagrangian, which determines the
#interpolation function
beta = 6/sqrt(2+c_2)
#Resolution of the uniform grid onto which we interpolate our results to have nicer plots
plot_resolution = mesh_resolution
#Size of the mesh for plotting. Should be bigger than the normal one or some points might be outside its domain
mesh_plot_size = domain_size*0.8
refine_times = 10
p = 1*kp
source_number = 1
source_mass = mgb/source_number
radius_population = domain_size/2
#Degree of the functionspace we want to solve the PDE on
#IMPORTANT: Optimal degree = 3, increasing it to 4 does not make the computation more accurate!
degree_PDE = 2
#Values for the three parameter beta model
beta = 0.654
#Characteristic radius
r_c = 242.3*kp
#Density. The power 10^-31 is the power needed for the integrated density to have the correct value!
#The units given in the paper dont add up!
rho_0 = 0.06*(10**(-31))
#BVP To be solved: we use a string so we dont need to define the object in here and instead we
#evaluate it in the main code
BVP_to_solve = 'mond_standard_beta'
#IMPORTANT!!! It might be that interpolating on a linear space makes all derivatives disappear:
#https://fenicsproject.org/qa/9893/simple-question-about-derivation/
#However, I should be fine cause I do the derivative, then interpolate it on a linear space and do it again
#So I never take a second derivative on the same space technically
#Interpolating all the lensing parameters onto new function spaces is expensive. Only do if needed
lensing_interpolations = False
#giving the option to use a finer, regular mesh and all functions interpolated on it for all plotting
plotting_option = False
#Deciding if we want the sources to all be in the same plane (can make it easier to get good contours, or
#for thin lens approximation for lensing)
coplanar_sources = True
#Option to make the comparison between two solutions after the whole code has run
make_comparison = False
#Option to have the first mass of the distribution placed in the origin
central_mass = True
#Option to calculate acceleration by interpolating the gradient of the solution on a vector space. 
#The operation is very expensive, so if it's not strictly necessary it's better to not do it at all
acceleration_needed = False
#Option to plot 3D graphs
plot_3D_graphs = False
#Option to avoid printing redundant information from each core when running the code in parallel from
#a python (.py) script obtained from the jupyter notebook.
# parallel_run = True