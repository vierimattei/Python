# -*- coding: utf-8 -*-

#Importing pandas for data acquisition/manipulation
import pandas as pd

#Importing re to use regular expressions
import re

#Importing astropy to use constants and coordinates, and units to do conversions
#and the WMAP cause we need the Hubble constant
from astropy import constants as const
from astropy import units as u
from astropy.cosmology import WMAP9 as cosmo

def make_quantity(value, unit):
    """Taking as input a numerical value and a unit and making a quantity value
    through the astropy module
    """
    
    #Checking if the type of the input is a frame series
    if isinstance(value, pd.Series):
        
        #cannot build quantity from frame series, only from numpy array.
        #series.values returns numpy array with the same elements as series
        value = value.values
    else:
        #otherwise use value directly
        value = value
    
    #build quantity from value*unit
    quantity = value*unit
    
    return quantity


def convert_unit(value, unit1, unit2):
    """Takes in a scalar, sequence or numpy array, and two astropy units. Converts
    the input from unit1 to unit2 and outputs its magnitude"""
    
    quantity = make_quantity(value, unit1)
        
    #building an astropy quantity from the input value, then converting it to a 
    #different unit and taking its magnitude
    
    quantity = quantity.to(unit2).value
    
    return quantity, str(unit2)


def update_series_name(frame, series, unit):
    """Takes as input a series from dataframe, renames it with unit to which it
    was converted in the convert_unit function"""
    
#    #Checking if the type of the input is a frame series
    if isinstance(frame[series], pd.Series):

        pattern = re.compile(r'\(\D+\)')
        new_name = re.sub(pattern, "("+unit+")", frame[series].name)
        
        #to rename the series it is best to rename the column in the dataframe
        #directly, or the series property changes but not the label in the frame
        #inplace=True means it changes it in the current frame, not a new one
        frame.rename(columns={series:new_name}, inplace=True)
    
    else:
        print("Input is not a series")
        

def convert_vel_dist(velocity, distance_unit, time_unit):
    """Takes as input a recession velocity in a given unit and converts it to
    a distance using Hubble's law with H(0) from astropy.cosmo
    """
    
    velocity_unit = distance_unit/time_unit
    
    #making the velocity magnitude into a quantity
    velocity = make_quantity(velocity, velocity_unit)
    #obtaining distance from Hubble's lsw with H(0) the present day value of H
    distance = velocity / cosmo.H(0)
    
    #splitting the quantity into its value and distance for output
    distance_magnitude = distance.value
    distance_unit = str(distance.unit)
    

    return distance_magnitude, distance_unit