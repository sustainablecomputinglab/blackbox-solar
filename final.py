from __future__ import division
import os
import datetime as dt

try:
    from importlib import reload
except ImportError:
    try:
        from imp import reload
    except ImportError:
        pass

import csv
import sys
import math
import pytz
import datetime
import matplotlib
import numpy as np
import pandas as pd
import sys , getopt
from math import sqrt
from tzwhere import tzwhere
from sklearn import metrics
from pvlib import atmosphere
import matplotlib.pyplot as plt
from pvlib.location import Location
from sklearn.metrics import mean_squared_error
from pvlib import clearsky, atmosphere, solarposition
from pvlib.tools import datetime_to_djd, djd_to_datetime

def get_solarposition(time, latitude, longitude, altitude=None, pressure=None, temperature=10, **kwargs):

    #calculate the altitude and pressure based on default vales
    if altitude is None and pressure is None:
        altitude = 0
        pressure = 101325
    elif altitude is None:
        altitude = atmosphere.pres2alt(pressure)
    elif pressure is None:
        pressure = atmosphere.alt2pres(altitude)
    
    #converting the time into DatetimeIndex 
    time = pd.DatetimeIndex([time,])

    #Normalizing the time
    t=time.tz_localize(tz = 'UTC')
    spa_python(t, latitude, longitude, altitude, pressure, temperature, **kwargs)
        
def _spa_python_import():

    from pvlib import spa
    os.environ['PVLIB_USE_NUMBA'] = '0'
    spa = reload(spa)
    del os.environ['PVLIB_USE_NUMBA']

    return spa

def spa_python(time, latitude, longitude, altitude=0, pressure=101325, temperature=10, delta_t=67.0, atmos_refract=None, numthreads=4, **kwargs):

    global result
    lati = latitude
    longi = longitude
    elevi = altitude
    pressure = pressure / 100 
    atmos_refract = atmos_refract or 0.5667
    
    #calculate the timezone of the given latitude and longitude
    tz = tzwhere.tzwhere()
    timezone_str = tz.tzNameAt(lati, longi)
    timezone_str

    global timezone
    timezone = pytz.timezone(timezone_str)

    #conerting the time into local time based on timezone
    ltime = time.tz_convert(timezone)

    #converting the time into unixtime
    unixtime = np.array(ltime.astype(np.int64) / 10**9)

    spa = _spa_python_import()
    delta_t = delta_t or spa.calculate_deltat(ltime.year, ltime.month)

    #calculate azimuth and zenith
    app_zenith, zenith, app_elevation, elevation, azimuth, eot = \
        spa.solar_position(unixtime, lati, longi, elevi, pressure, temperature,
                           delta_t, atmos_refract, numthreads)

    result = pd.DataFrame({'apparent_zenith': app_zenith, 
                           'zenith': zenith,
                           'apparent_elevation': app_elevation,
                           'elevation': elevation, 
                           'azimuth': azimuth,
                           'equation_of_time': eot}
                          )
    return result

#converting the arguments to float
latitude = float(sys.argv[1])
longitude = float(sys.argv[2])

get_solarposition('12/30/2015  23:00:00', latitude, longitude)


#reading the data.csv file into a pandas dataframe
df = pd.read_csv(sys.argv[3], delimiter=',')

#obtaining the date and time in a pandas dataframe
time_df = df['Date & Time']

#making the Date & Time column as index column 
df.set_index(['Date & Time'], inplace=True)

#copying the Solar Generation values into a pandas dataframe
af = df['Generation [kW]']

#obatin the start and end time 
length = len(df) - 1 
end_time = time_df[0]
start_time = time_df[length]

#obtaining the clear sky solar generation for the given latitude and longitude for the given start and end time
place = Location(latitude, longitude, timezone , 0 , 'easthampton')
times = pd.DatetimeIndex(start=start_time, end=end_time,  freq='1H', tz=place.tz)
cs = place.get_clearsky(times)

#copying the ghi values into pandas dataframe
bf = cs['ghi']

#inverting the pandas dataframe
bf = bf.iloc[::-1]

#obtaing the ambient temperatue into pandas dataframe
tf = df['Ambient Temp [°c]']


#conerting the azimuth and zenith angles from degrees to radians
azimuth = result['azimuth'] * 0.0175
zenith = result['zenith'] * 0.0175


def parameter_K(clear_sky, gen_data, z, surface_tilt, surface_azimuth):

    k = z

    solar_data = clear_sky * k * ( math.cos(1.5708-zenith) * math.sin(surface_tilt) * math.cos(surface_azimuth - azimuth) + math.sin(1.5708-zenith) * math.cos(surface_tilt))

    #varying the values of k to obtain the tight upper bound of the solar data using binary search  
    if (solar_data > clear_sky):
        k = k/2
        solar_data = clear_sky * k * ( math.cos(1.5708-zenith) * math.sin(surface_tilt) * math.cos(surface_azimuth - azimuth) + math.sin(1.5708-zenith) * math.cos(surface_tilt))
        
        #untill calculated data is equal to solar data divide the k by half 
        while (solar_data > gen_data and gen_data > 0):
            k = k/2
            solar_data = clear_sky * k * ( math.cos(1.5708-zenith) * math.sin(surface_tilt) * math.cos(surface_azimuth - azimuth) + math.sin(1.5708-zenith) * math.cos(surface_tilt))      

    
    if (solar_data < clear_sky):
        k = k/2
        solar_data = clear_sky * k * ( math.cos(1.5708-zenith) * math.sin(surface_tilt) * math.cos(surface_azimuth - azimuth) + math.sin(1.5708-zenith) * math.cos(surface_tilt)) 
        
        while (solar_data > gen_data):
            #untill calculated data is equal to solar data divide the k by half
            k = k/2
            solar_data = clear_sky * k * ( math.cos(1.5708-zenith) * math.sin(surface_tilt) * math.cos(surface_azimuth - azimuth) + math.sin(1.5708-zenith) * math.cos(surface_tilt))       
    
        while (solar_data < gen_data and solar_data > 0 ):
            #untill calculated data is lesser than solar data multiply the k by 2
            k = k*2
            solar_data = clear_sky * k * ( math.cos(1.5708-zenith) * math.sin(surface_tilt) * math.cos(surface_azimuth - azimuth) + math.sin(1.5708-zenith) * math.cos(surface_tilt)) 

    return k

def parameter_azimuth(clear_sky, gen_data,new_k,surface_tilt, surface_azimuth):

    solar_data = 0
    perfect_surface_azimuth = surface_azimuth
    new_surface_azimuth = surface_azimuth
    new_solar_data = 0

    #vary the value of azimuth from 0 degrees to 360 degrees to obtain the one which gives tightest bound on data
    for i in np.arange(0, 6.28319, 0.08727):
        
        #copy of older azimuth and its data values
        surface_azimuth = i
        old_solar_data = new_solar_data
        
        old_surface_azimuth = new_surface_azimuth
        
        solar_data = clear_sky * new_k * ( math.cos(1.5708-zenith) * math.sin(surface_tilt) * math.cos(surface_azimuth - azimuth) + math.sin(1.5708-zenith) * math.cos(surface_tilt))
        
        #copy of new azimuth and its data values
        new_solar_data = solar_data
        new_surface_azimuth = surface_azimuth

        #discard the values if it is outside the bound i.e above clear sky and solar generation
        if (solar_data < gen_data):
            continue
        elif (solar_data > clear_sky):
            continue
        else:
            #inside the bound obtain the one which gives tightest bound on data
            if(new_solar_data < old_solar_data):
                perfect_surface_azimuth = old_surface_azimuth
            else:
                perfect_surface_azimuth = new_surface_azimuth
            
    return perfect_surface_azimuth

def parameter_tilt(clear_sky, gen_data,new_k,surface_tilt, surface_azimuth):

    solar_data = 0
    perfect_surface_tilt = surface_tilt
    new_surface_tilt = surface_tilt
    new_solar_data = 0
    
    #vary the value of tilt from 0 degrees to 90 degrees to obtain the one which gives tightest bound on data
    for i in np.arange(0, 1.5708, 0.08727):
        
        surface_azimuth = i
        
        #copy of older tilt and its data values
        old_solar_data = new_solar_data
        
        old_surface_tilt = new_surface_tilt
        
        solar_data = clear_sky * new_k * ( math.cos(1.5708-zenith) * math.sin(surface_tilt) * math.cos(surface_azimuth - azimuth) + math.sin(1.5708-zenith) * math.cos(surface_tilt))
        
        #copy of new tilt and its data values
        new_solar_data = solar_data
        new_surface_tilt = surface_tilt

        #discard the values if it is outside the bound i.e above clear sky and solar generation
        if (solar_data < gen_data):
            continue
        elif (solar_data > clear_sky):
            continue
        else:
            #inside the bound obtain the one which gives tightest bound on data
            if(new_solar_data < old_solar_data):
                perfect_surface_tilt = new_surface_tilt
            else:
                perfect_surface_tilt = old_surface_tilt
            
    return perfect_surface_tilt

def parameters (z, surface_azimuth, surface_tilt):
    gen_data = []
    data = []
    rmse = []
    k_list = []
    sa_list = []
    st_list = []
    clean_data = []
    for x,y in zip(bf,af):

        #obtain all the solar generation values in a list
        gen_data.append(y)

        #call the parameter_k function to get the optimised value of k
        optimised_k = parameter_K(x,y,z,surface_tilt, surface_azimuth)
        k_list.append(optimised_k)

        new_k = optimised_k
        old_surface_azimuth = surface_azimuth
       
        #call the parameter_azimuth function to get the optimised value of surface azimuth(orientaion)
        optimised_surface_azimuth = parameter_azimuth(x,y,new_k,surface_tilt, old_surface_azimuth)
        sa_list.append(optimised_surface_azimuth)

        new_surface_azimuth = optimised_surface_azimuth
        old_surface_tilt = surface_tilt

        #call the parameter_tilt function to get the optimised value of surface tilt
        optimised_surface_tilt = parameter_tilt(x,y,new_k,old_surface_tilt, new_surface_azimuth)
        st_list.append(optimised_surface_tilt)

        

        #calculate the solar data for the optimised values of k, orientaion and tilt
        solar_data = x * optimised_k * ( math.cos(1.5708-zenith) * math.sin(optimised_surface_tilt) * math.cos(optimised_surface_azimuth - azimuth) + math.sin(1.5708-zenith) * math.cos(optimised_surface_tilt))
        data.append(solar_data)

        #discard the values if it is 0 as they will affect the rmse value
        if (solar_data != 0 and solar_data != -0.0):
            clean_data.append(solar_data)
        
        #calculate the rmse for non zero values    
        rmse_data = np.sqrt(metrics.mean_squared_error(gen_data, data))
        rmse.append(rmse_data)
    
    clean_data_length = len(clean_data) - 1

    data_clean_first_index = data.index(clean_data[0])
    data_clean_last_index = data.index(clean_data[clean_data_length])
    
    minimum_rmse = min(rmse[data_clean_first_index : data_clean_last_index])
    index_rmse = rmse.index(minimum_rmse)

    #obtain the values of k, orientaion and tilt values for which we get minimum rmse
    next_k = k_list[index_rmse]
    next_surface_azimuth = sa_list[index_rmse]
    next_surface_tilt = st_list[index_rmse]
    
    return (next_k, next_surface_azimuth, next_surface_tilt, data)

#calling the parameter function to get the optimised values of k, orientaion and tilt values
next1 = [100, 3.1415, 0.737409, None]
final =  parameters(next1[0], next1[1], next1[2])

#obtain the index value of calculated data which has the minimum differnce with the solar generation data 
s_data = final[3]
s_clean_data = []
for x in s_data:
    if ( x != 0.0 and x != -0.0 ):
        s_clean_data.append(x)
        
s_clean_data_length = len(s_clean_data) - 1

s_data_clean_first_index = s_data.index(s_clean_data[0]) 
s_data_clean_last_index = s_data.index(s_clean_data[s_clean_data_length])
s_af = af[s_data_clean_first_index:s_data_clean_last_index]
    
diff_list = []
for d,g in zip(s_clean_data , s_af):
    if(d > g):
        diff = d - g
        diff_list.append(diff)
    else:
        diff = g - d
        diff_list.append(diff)
   

minimum_diff = min(diff_list)
minimum_diff_index = diff_list.index(minimum_diff) + s_data_clean_first_index
#print(diff_list)

def parameter_k_adjust(clear_sky, gen_data,optimised_k, optimised_surface_azimuth, optimised_surface_tilt, Tbase, Tair):
    
    global c
    c_list = []
    azimuth = 2.05819
    zenith = 1.639965
    solar_data = 0
    perfect_k_adjust = 0
    new_k_adjust = 0
    new_solar_data = 0
    
    #calculating the value of c for wich it eives the tightest upper bound
    for i in np.arange(0, 0.5, 0.001):
        
        c = i
        old_solar_data = new_solar_data
        
        old_k_adjust = new_k_adjust
        
        #calculate the adjusted k for the varying c
        k_adjust = optimised_k * (1 + c * ( Tbase - Tair))

        solar_data = clear_sky * k_adjust * ( math.cos(1.5708-zenith) * math.sin(optimised_surface_tilt) * math.cos(optimised_surface_azimuth - azimuth) + math.sin(1.5708-zenith) * math.cos(optimised_surface_tilt))
        
        new_solar_data = solar_data
        new_k_adjust = k_adjust

        if (solar_data < gen_data):
            continue
        elif (solar_data > clear_sky):
            continue
        else:
            if(new_solar_data < old_solar_data):
                perfect_k_adjust = new_k_adjust
            else:
                perfect_k_adjust = old_k_adjust
    
    return (perfect_k_adjust)

def parameters_k_adjusted (k, surface_azimuth, surface_tilt):
    gen_data = []
    data = []
    rmse = []
    k_list = []
    sa_list = []
    st_list = []
    clean_data = []
    Tair = []
    Tair = tf

    Tbase = Tair[minimum_diff_index]
    for x,y,t in zip(bf,af,Tair):
        
        gen_data.append(y)
        
        Tair = t

        #obtain the temperature adjusted value of k
        adjusted_k = parameter_k_adjust(x,y,k, surface_azimuth,surface_tilt, Tbase, Tair)
        k_list.append(adjusted_k)
        
        #calculate the solar data based on the adjusted k 
        solar_data = x * adjusted_k * ( math.cos(1.5708-zenith) * math.sin(surface_tilt) * math.cos(surface_azimuth - azimuth) + math.sin(1.5708-zenith) * math.cos(surface_tilt))
        data.append(solar_data)

        if (solar_data != 0 and solar_data != -0.0):
            clean_data.append(solar_data)
        
        #calculate the rmse for non zero values    
        rmse_data = np.sqrt(metrics.mean_squared_error(gen_data, data))
        rmse.append(rmse_data)
    
    clean_data_length = len(clean_data) - 1

    data_clean_first_index = data.index(clean_data[0])
    data_clean_last_index = data.index(clean_data[clean_data_length])
    
    minimum_rmse = min(rmse)
    index_rmse = rmse.index(minimum_rmse) + data_clean_first_index

    #obtain the values of k, orientaion and tilt values for which we get minimum rmse
    next_k = k_list[index_rmse]
   
    return (next_k, surface_azimuth, surface_tilt,data,c)

para = parameters_k_adjusted (final[0],final[1],final[2])

#the value of k
k = para[0]
print("The value of K is :" , k)

#the value of orientation in degress by multiplying it with 57.296
orientation = para[1] * 57.296
print("The value of orientation is :" , orientation)

#the value of tilt in degress by multiplying it with 57.296
surface_tilt = para[2] * 57.296
print("The value of surface_tilt is :" , surface_tilt)

#the value of c
c = para[4]
print("The value of c is :" , c)
