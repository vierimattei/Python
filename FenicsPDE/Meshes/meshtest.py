# -*- coding: utf-8 -*-
#importing the meshr module and the fenics module
from fenics import *
from mshr import *

#importing matplotlib to plot results 
from matplotlib import pyplot as plt

#importing timeit to time length of mesh creation
import time

#importing numpy to work with arrays
import numpy as np

#importing random to distribute masses on a mesh
import random

#need math to use sine and cosine
from math import sin,cos

# Built-in meshes
#defining the radius of a circle
big_radius = 10
#defining center of the circle
big_center = Point(0,0)
#Defining a Circle (case sensitive!), taking input the center and radius
big_circle = Circle(big_center, big_radius)

#declaring the domain using the built in command 
domain = big_circle

#radius of subdomains
sub_radius = 0.5

#number of masses
sub_num = 10

#angles equally spaced apart
angles_N = [i*2*pi/sub_num for i in range(sub_num)]

#centers are now around
sub_center_x = [cos(angle) for angle in angles_N]
sub_center_y = [sin(angle) for angle in angles_N]

#grouping x and y coordinates to feed to Circle
sub_coords = np.vstack((sub_center_x, sub_center_y))
sub_coords = np.transpose(sub_coords)

#Points of origin of each mass. each element of sub_coords is pair of coordinates. We make a point from each
#then make a circle around that origin with the chosen radius
sub_coords_points = [Circle(Point(coords),sub_radius) for coords in sub_coords]

#for each coordinate and circle defined above we make a subdomain. enumerate return the index i of the coords,
#so we use it to give the subdomain index. need to add 1 cause cannot use subdomain 0
for i, mass in enumerate(sub_coords_points):
    domain.set_subdomain(i+1, mass)

mesh = generate_mesh(domain, 10)
#mesh
