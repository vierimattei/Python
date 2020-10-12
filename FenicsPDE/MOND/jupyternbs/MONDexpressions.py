#Need to import the functions as we use one to create the source string for multiple sources
from MONDfunctions import * 


# # C++ Expressions for the sources

# In[21]:


#VERY IMPORTANT: Each constant used in declaring an expression in c++ code needs to be
#declared inside Expression, after the degree! If not, the c++ compiler will throw an error!
#IMPORTANT!: Cannot use the log function cause that name is apparently already taken by some
#other function, not the logarithm! So using log10 and rescaling by 0.43429.
#IMPORTANT!! When changing the radius_tot in the above script is doesn't change in here cause
#it has to be defined independently for all of these C++ expressions! If you don't change it 
#inside here for both u and f, the solution will be wrong!

#Sphere of constant density
#Using separate '...' on each line so I can write the C++ code over multiple lines, less messy
f_sphere_cpp = ('(pow(x[0],2) + pow(x[1],2) + pow(x[2],2))<=pow(radius_tot, 2) ?'
'4*pi*G*mgb/volume_out : 0')

# Dirac Delta function using conditional for each coordinate (doesnt work)
f_dirac_coord_cond_cpp = ('(abs(x[0]-x_close) <= x_close && '
                                   'abs(x[1]-y_close) <= y_close && '
                                   'abs(x[2]-z_close) <= z_close)  ? '
                                   '4*pi*G*mgb: 0')

#Dirac delta using very small radius (doesnt work)
f_dirac_radius_cond_cpp = ('pow(x[0],2) + pow(x[1],2) + pow(x[2],2) <= 2*pow(radius_close,2)'
                           ' ? mgb : 0')

#Dirac Delta using definition (as suggested here https://fenicsproject.org/qa/7941/applying-
#dirac-function-using-fenics-pointsource-function/)
f_dirac_analytic1_cpp = ('4*pi*G*mgb*eps/pi*1/(pow(x[0],2) + pow(x[1],2) + pow(x[2],2) + pow(eps,2))')

#Defining a Gaussian pulse with a standard deviation radiustot/3, so 99.7% of mass is inside
#radius tot and we can compare with the analytic solution for a sphere
f_gauss_cpp = ('4*pi*G*mgb*(1/pow((stand_dev*sqrt(2*pi)),3))*exp(-(1/2)*((pow(x[0]-x_close ,2)'
               '+ pow(x[1] - y_close,2) + pow(x[2] - z_close,2))/(pow(stand_dev,2))))')

#Defining the source and initial guess for an isothermal distribution
f_isothermal_cpp = ('4*pi*G*3*mgb/(4*pi)*pow((pow(p, 1/2)/(pow(p, 3/2) +'
                    'pow((pow(x[0],2) + pow(x[1],2) + pow(x[2],2)), 3/4))), 3)')

#Testing the exponential function in the source cause it doesnt seem to work, THIS WORKS!
#APART FROM HAVING SLIGHTLY MORE MASS THAN IT SHOULD!
f_exponent_test = ('4*pi*G*pow(2*pi,-1.5)*source_mass/pow(stand_dev,3)*exp(-(pow(x[0],2) + pow(x[1],2) + pow(x[2],2))/(2*(pow(stand_dev,2))))')

#Trying to implement a for loop inside the Expression statement
f_displaced_cpp = ('(pow((x[0]-position_x),2)+pow((x[1]-position_y),2)+pow((x[2]-position_z),2))'
                '<=pow(radius_tot, 2) ?'
                '4*pi*G*source_mass/volume_out : 0')

#Calling the function to make the cpp code for the discrete distribution and assigning the output
#to the string_cpp variable
f_multiple_dirac = make_discrete_dirac(source_number)

#Source expression for multiple sources using a for loop over the source locations defining
#the string itself, as in the cell above
f_multiple_gauss = make_discrete_gauss(source_number)


# # C++ Expressions for the Initial Guesses/BCs

# In[22]:


#Solution for a sphere of uniform density, including the solution inside and outside the sphere
u_sphere_cpp = ('(pow(x[0],2) + pow(x[1],2) + pow(x[2],2))>=pow(radius_tot, 2) ?'
'sqrt(G*mgb*a0)*1/2*2.3026*log10(pow(x[0],2)+pow(x[1],2)+pow(x[2],2)) :'
'(4/3*sqrt(pi/3*a0*G*mgb/volume_out)*3/2*pow((pow(x[0],2)+pow(x[1],2)+pow(x[2],2)),3/4)+'
'sqrt(G*mgb*a0)*2.3026*log10(radius_tot)-4/3*sqrt(pi/3*a0*G*mgb/volume_out)*'
'3/2*pow(radius_tot,3/2))')

#Dirac delta in the origin
u_dirac_cpp = ('sqrt(G*mgb*a0)*1/2*2.3026*log10(pow(x[0],2)+pow(x[1],2)+pow(x[2],2))')

#Solution for a mass distribution follwoing an isothermal profile
u_isothermal_cpp = ('2/3*sqrt(G*mgb*a0/6)*2.3206*log10(1 + pow((pow(x[0],2) + '
                    'pow(x[1],2) + pow(x[2],2)), 3/4) / pow(p, 3/2))')

#Boundary condition for displaced source/center of mass
u_displaced_cpp = ('sqrt(G*source_mass*a0)*1/2*2.3026*log10(pow((x[0] - center_of_mass_x),2)'
                          '+pow((x[1] - center_of_mass_y),2)+pow((x[2] - center_of_mass_z),2))')

#Boundary condition for Newtonian Gravity
u_Newton = ('-G*mgb/sqrt(pow((x[0] - center_of_mass_x),2)'
                          '+pow((x[1] - center_of_mass_y),2)+pow((x[2] - center_of_mass_z),2))')


# # C++ Expressions for the weak form of a specific PDE, e.g. different interpolation of GEA modification. Defined as strings to be evaluated with eval cause Python doesnt know what u is yet

# In[23]:


#Weak for for Newton
F_Newton = 'inner(grad(u), grad(v))*dx + f*v*dx'

#Weak form for deep MOND
F_MOND_deep = 'inner(sqrt(inner(grad(u), grad(u))) * grad(u)/a0, grad(v))*dx + f*v*dx'

#Weak form for the simple interpolation function
F_MOND_simple = ('inner(sqrt(inner(grad(u), grad(u)))/(sqrt(inner(grad(u), grad(u)))+a0)'
              '* grad(u), grad(v))*dx + f*v*dx')
    
#Weak form for the standard interpolation function
F_MOND_standard = ('inner(sqrt(inner(grad(u), grad(u)))/sqrt((inner(grad(u), grad(u))+a0**2))'
              '* grad(u), grad(v))*dx + f*v*dx')

#Weak form for the exponential interpolation function. Doesnt seem to work but it behaves almost exactly
#the same as the standard inteprolation function, so it's not so bad.
F_MOND_exponential  = ('inner(1-exp(-sqrt((inner(grad(u), grad(u))))/a0)'
              '* grad(u), grad(v))*dx + f*v*dx')

# #Defining the r unit vector in terms of cartesian coordinates for use in weak form of GEA
r_unit_1 = Expression('x[0]/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])', degree = 1)
r_unit_2 = Expression('x[1]/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])', degree = 1)
r_unit_3 = Expression('x[2]/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])', degree = 1)

#Weak form for GEA with spherical symmetry component. Multiplying the distribution by 6/beta as for GEA
#we use the acceleration scale M=6a_0, so a_0-> 6a_0 and we need to divide by beta
F_MOND_GEA = ('(inner(sqrt(2*inner(grad(u_GEA), grad(u_GEA)) - c_2*((u_GEA.dx(0)*r_unit_1)+'
        '(u_GEA.dx(1)*r_unit_2)+(u_GEA.dx(2)*r_unit_3))**2)* grad(u_GEA), grad(v))*dx + 6/beta*f*v*dx')

