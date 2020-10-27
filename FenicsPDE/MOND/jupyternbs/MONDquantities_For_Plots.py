from dolfin import *

#importing mshr for all mesh functions
import mshr as mshr

c = 2.998*10**8
G = 6.674*10**(-11) 

#Value of a0. Currently changed it so the value of the Hubble constant when we take it as H0=2pi*a0/c
#Is close, but below the most recent estimate of the Hubble constant. We do this cause for the cluster
#database we need to make a pick for the numerical value of H0, and setting it w.r.t. a0 gives a 
#coherent picutre for MOND on cluster scales.
a0 = 1.14*10**(-10) 
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
radius_tot = domain_size/200
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
stand_dev_peak = radius_tot/3*10
mesh_resolution = 21
#Coefficient for GEA changing the potential based on how spherically symmetric the mass distribution is
c_2 = -1.8
#Coefficient for GEA giving the magnitude of the K^3/2 term in the Lagrangian, which determines the
#interpolation function
beta_GEA = 6/sqrt(2+c_2)
#Resolution of the uniform grid onto which we interpolate our results to have nicer plots
plot_resolution = mesh_resolution
#Size of the mesh for plotting. Should be bigger than the normal one or some points might be outside its domain
mesh_plot_size = domain_size*0.8
refine_times = 0
p = 1*kp
source_number = 1
source_mass = mgb/source_number
radius_population = domain_size/2
#Degree of the functionspace we want to solve the PDE on
#IMPORTANT: Optimal degree = 3, increasing it to 4 does not make the computation more accurate!
degree_PDE = 1
#Values for the three parameter beta model
beta = 0.654
#Characteristic radius
r_c = 242.3*kp
#Density for the beta distribution. In paper it's given as gram/cm^3, but I'm using everything in terms
#of meters and kilograms, and 1g/cm^3 = 1000kg/m^3, so need to multiply by 1000. from 10^-25 to 10^-22
#so this should be in terms of kg/m^3
rho_0 = 0.06*(10**(-22))
#All the measurements use a reduced Hubble constant. As the Hubble constant can be expressed in terms
#of Miglrom's constant and speed of light, we express it as that. I modified a0 = 1.14 * 10**-10 so 
#that H0 is now very close to the latest (upper) estimate. Now a0 and H0 are consistent.
H0 = 2*pi*(a0)/c
#The Hubble constant used is h50, given as H0/50 km/s/Mpc. We need to multiply most of the quantities in
#the table by this to obtain 
h50 = H0/(50*10**3/(1000*kp))
#The size of the domain is the size over which the measurements were carried out. Technically, we should
#only consider the mass integrated up to r200, where the gas is supposed to have 200 times the critical
#mass density of the universe. However, the tail should not make a big difference (can check and change
#the beta density expression if the assumption does not hold). 
#On the other hand, the total gravitational mass of all systems is calculated within the same radius,
#the so called Abell radius. We use this as the domain size for all of the clusters.
radius_abell = 3/h50*1000*kp
#BVP To be solved: we use a string so we dont need to define the object in here and instead we
#evaluate it in the main code
BVP_to_solve = 'mond_deep_dirac'
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
#Plots evaluating performance for various parameters in the simulation
plots_for_thesis = True