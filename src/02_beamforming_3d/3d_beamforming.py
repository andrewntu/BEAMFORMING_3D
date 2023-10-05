#!/usr/bin/env python
import cmath
import math
from math import radians, cos, sin, asin, sqrt, pi, acos, tan, atan, exp
import obspy
from obspy.core import read
import obspy.signal
from obspy.signal.filter import envelope
import glob
import sys
import numpy as np
import os
from numba import njit
from numba import jit
import time
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy
import pymap3d as pm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib

#==== Needed functions ====#
# 1. create the grids 
# 2. single beamforming
# 3. calculate distance and degree transformation

@jit
def get_dist(rel_x, rel_y, rel_z, x_grid, y_grid, z_grid):
    #=== Combine grid location with station location ===#
    all_dist = np.zeros((len(x_grid),len(y_grid),len(z_grid),len(rel_x)))
    all_time_shift = np.zeros((len(x_grid),len(y_grid),len(z_grid),len(rel_x)))
    for ss in range(0,len(rel_x)):
        for xx in range(0,len(x_grid)):
            for yy in range(0,len(y_grid)):
                for zz in range(0,len(z_grid)):
                    loc_x = x_grid[xx] - rel_x[ss];
                    loc_y = y_grid[yy] - rel_y[ss];
                    loc_z = z_grid[zz] + rel_z[ss];
                    all_dist[xx,yy,zz,ss] = np.sqrt((loc_x)**2 + (loc_y)**2 + (loc_z)**2)/1000
                    all_time_shift[xx,yy,zz,ss] = all_dist[xx,yy,zz,ss]/vel_mod # slowness * distance = shift_time
    return all_dist, all_time_shift                    

def BF_shift_and_stacking(data, waveform, omega, time_shift):
    save_amp = np.zeros((len(x_grid), len(y_grid), len(z_grid), len(waveform), len(omega)))
    save_env = np.zeros((len(x_grid), len(y_grid), len(z_grid), len(omega)))
    result_env = np.zeros((len(x_grid), len(y_grid), len(z_grid)))
    for xx in range(0,len(x_grid)):
        for yy in range(0,len(y_grid)):
            for zz in range(0,len(z_grid)):
                amp = np.zeros((len(omega)))
                amp = np.array(amp,dtype=complex)
                for i in range(0,len(waveform)):
                    #print(time_shift[xx, yy, zz, i])
                    dtcs_total = complex(time_shift[xx, yy, zz, i])
                    test_omega = omega*dtcs_total
                    shift_waveform = waveform[i,:]*np.exp(-1j*test_omega)
                    tdata_shift_temp = np.real(np.fft.irfft(shift_waveform[0:int(len(omega)/2)+1]))
#                    tdata_shift_temp = np.real(np.fft.irfft(shift_waveform))
#                    print(shift_waveform.shape)
                    '''
                    plt.figure()
                    plt.plot(tdata_shift_temp, 'k-')
                    plt.plot(data[i,:],'r-')
#                    plt.plot(waveform[i,:],'g-')
                    plt.show()
#                    '''
                    amp = amp + waveform[i,:]*np.exp(-1j*test_omega)
                    save_amp[xx,yy,zz,i,:] = tdata_shift_temp
                tdata_stack_temp = np.fft.irfft(amp[0:int(len(omega)/2)+1])
                tdata_stack = np.real(tdata_stack_temp)/len(waveform)
                envelope_amp = scipy.signal.hilbert(tdata_stack)
                result_env[xx,yy,zz] = round(np.max(np.abs(envelope_amp)),6)
#                save_amp[xx,yy,zz,:] = tdata_stack
                save_env[xx,yy,zz,:] = np.abs(envelope_amp)
                '''
                plt.figure()
                plt.plot(tdata_stack, 'k-')
                plt.plot(np.abs(envelope_amp),'r-')
                plt.show()
#                '''
#                sys.exit()
#                print("=========================")
    return result_env, save_amp, save_env

@jit 
def cal_rel_time_shift(rel_x, rel_y, rel_z, x_grid, y_grid, z_grid, shift_time, idx):
    rel_shift_time = np.zeros(shift_time.shape)
    for ss in range(0,len(rel_x)):
        rel_shift_time[:,:,:,ss] = shift_time[:,:,:,idx] - shift_time[:,:,:,ss]
    return rel_shift_time
        
def hanning(data):
    new_data = np.zeros((len(data)))
    for i in range(0,len(data)):
            new_data[i] = data[i]*(0.5-0.5*cos(2*pi*i/len(data)-1))
    return new_data

def fft_parameter(dd):
    n = next_power_of_2(len(dd[0,:]))
    n = int(n)
    freq = np.ones(n)
    for i in range(0,len(freq)):
        freq[i] = sample*freq[i]*i/n
    omega = 2*pi*freq
    omega = np.array(omega,dtype=complex)
    temp2 = np.zeros((len(dd),n))
    fdata = np.array(temp2,dtype=complex)
    for i in range(0,len(dd)):
        fdata[i,:] = np.fft.fft(dd[i,:],n)
    return fdata,omega

def next_power_of_2(x):  
    n = 2**(x - 1).bit_length()
    return n

def cal_snr(data,dist,cut_left,cut_right):
        snr_follow = []
        center = int(len(data[0,:])/2)
        for i in range(0,len(data)):
                left = int(sample*dist[i]/max_v)
                right = int(sample*dist[i]/min_v)
                snr_fol = np.max(np.abs(data[i,center+left:]))/np.sqrt(np.mean(data[i,noi_b:]**2))
                snr_follow.append(snr_fol)
        return snr_follow

if __name__ == "__main__":
    global lon_r, lat_r, lon_s, lat_s, event_lat, event_lon, sample, size, envelope, min_v, max_v, freq, vel_mod, center
    global noi_b, sig_b, sig_e, center_pts, x_size, y_size, z_size
    if len(sys.argv[:]) != 6:
            print('proper usage: python ./src/02_beamforming_3d/3d_beamforming.py stalst infold csta freq stage')
            sys.exit()
    tStart = time.time()
    n = 0

#=== Parameter ===#
    SaveDirectory = os.getcwd()
    in_stalst = str(sys.argv[1])
    in_fold = str(sys.argv[2])
    csta = str(sys.argv[3])
    freq = str(sys.argv[4])
    stage = str(sys.argv[5])
    #=== value from Cros et al. (2011) & Vandemeulebrouck et al. (2013). ===#
    vel_mod = 0.5
    #=== Specify the cutting window ===#
    l_win = -2
    r_win = 2


#=== Read in/Set up the directory ===#
    if not os.path.isdir('Figure'):
        os.mkdir('Figure')
    if not os.path.isdir('Output'):
        os.mkdir('Output')
    dir0 = 'Output/freq_'+freq
    if not os.path.isdir(dir0):
        os.mkdir(dir0)
    f_out = open(dir0+'_3d_beam.txt','a+')

    sta = []
    data = []
    stlon = []
    stlat = []
    stel = []
    evlon = []
    evlat = []
    dist = []
    trace = []
    staname = []

#=== Load Data ===#
    #f = open('Output/freq_'+str(freq)+'/'+csta+'_'+str(freq)+'_all.txt','a+')
    f = open('Output/freq_'+str(freq)+'/'+csta+'_'+str(freq)+'.txt','a+')
    csta_idx = 0
    stalst = pd.read_table(in_stalst, header = None, sep = ' ')
    for r in range(0,len(stalst)):
        rec = str(stalst.iloc[r,0])
        pair = csta+'-'+rec
        ccor = glob.glob(SaveDirectory+'/'+in_fold+'/'+csta+'/'+freq+'hz/'+pair+'.ZZ.stage.'+stage+'.sac.norm')
        if len(ccor) == 1:
            st = read(ccor[0],debug_header=True)
            sta.append(rec)
            tr = st[0]
            sample = tr.stats.sampling_rate
            trace.append(tr)
            #f.writelines(pair+' '+str(tr.stats.sac.dist)+'\n')
            center_pts = int(len(tr.data)/2)
            normal = np.max(np.abs(tr.data))
            tr.data[0:int(l_win * sample) + center_pts] = 0
            tr.data[int(r_win * sample) + center_pts:] = 0
            #data.append(tr.data[int(l_win * sample) + center_pts : int(r_win*sample) + center_pts]/normal)
            data.append(tr.data/normal)
            stlon.append(tr.stats.sac.stlo)
            stlat.append(tr.stats.sac.stla)
            stel.append(tr.stats.sac.stel)
            staname.append(rec)
            evlon.append(tr.stats.sac.evlo)
            evlat.append(tr.stats.sac.evla)
            event_lon = tr.stats.sac.evlo
            event_lat = tr.stats.sac.evla
            baz = tr.stats.sac.baz
            '''
            plt.figure()
            #plt.plot(tr.data[int(l_win * sample) + center_pts : int(r_win*sample) + center_pts]/normal)
            plt.plot(tr.data)
            plt.show()
#            '''
            if str(stalst.iloc[r,0]) == rec:
                lon_r = stalst.iloc[r,1]
                lat_r = stalst.iloc[r,2]
                lon_s = tr.stats.sac.evlo
                lat_s = tr.stats.sac.evla
            if str(stalst.iloc[r,0]) == csta:   
                csta_idx = r
#    f.close()
#    sys.exit() 
    print("event loc:",csta,lon_s,lat_s, csta_idx)
    #sys.exit()
    '''
    sig_b = int((dist_stack/vel_mod)*sample + center_pts)
    sig_e = int((dist_stack/vel_mod)*sample + center_pts)
    noi_b = sig_e
    '''
    data = np.array(data)
    evlon = np.array(evlon)
    evlat = np.array(evlat)
    stlon = np.array(stlon)
    stlat = np.array(stlat)
    stel = np.array(stel)
    staname = np.array(staname)

#=== Create Grid ===# !!! The Old Faithful cone is the referenced center.
    center = [-110.828211897200987, 44.460437153479248, 2240]
    #center = [-110.829643, 44.464347, 2240]
    # !!! The grid size in meters is 70*70*20
    '''
    x_range = 2000
    y_range = 2000
    z_range = 2000
    x_size = 50
    y_size = 50
    z_size = 50
    '''
    x_range = 200
    y_range = 200
    z_range = 100
    x_size = 5
    y_size = 5
    z_size = 2
    x_grid = np.arange(-1*(x_range/2), (x_range/2) + x_size, x_size)
    y_grid = np.arange(-1*(y_range/2), (y_range/2) + y_size, y_size)
    z_grid = np.arange( 0, z_range + z_size, z_size)

    x_grid_1 = np.arange(-1*(x_range/2) - x_size, (x_range/2) + x_size, x_size)
    y_grid_1 = np.arange(-1*(y_range/2) - y_size, (y_range/2) + y_size, y_size)
    z_grid_1 = np.arange( 0 - z_size, z_range + z_size, z_size)
#    print(x_grid_1)

#=== Tranfer Time Series to Frequency domain ===#
    fdata,omega = fft_parameter(data)

#=== Calculating the distance array ===#
    rel_x, rel_y, rel_z = pm.geodetic2enu(stlat, stlon, stel, center[1], center[0], center[2], ell = pm.Ellipsoid.from_name('grs80'))
    dist, time_shift = get_dist(rel_x, rel_y, rel_z, x_grid, y_grid, z_grid)
    rel_time_shift = cal_rel_time_shift(rel_x, rel_y, rel_z, x_grid, y_grid, z_grid, time_shift, csta_idx)
    #print(rel_x[0], rel_y[0], rel_z[0], x_grid[0], y_grid[0], z_grid[0])
    #print(dist, time_shift)

    '''
    #=== check distance in each station ===#
    for xx in range(0,len(x_grid)):
        for yy in range(0,len(y_grid)):
            for zz in range(0,len(z_grid)):
                for i in range(0,len(stlon)):
                    print(rel_x[i], x_grid[xx], rel_y[i], y_grid[yy], rel_z[i], z_grid[zz], dist[xx,yy,zz,i], time_shift[xx,yy,zz,i],staname[i])

    sys.exit()
    '''

#=== Shift-and-stack Beamforming ===#
    print('Beamforming Start !')
    output_env, save_amp, save_env= BF_shift_and_stacking(data, fdata, omega, rel_time_shift)
    best_grid_loc = np.unravel_index(output_env.argmax(), output_env.shape) 
    print(best_grid_loc)
    #print(output_env[best_grid_loc[0], best_grid_loc[1], best_grid_loc[2]] == output_env.max())
    print(output_env[best_grid_loc[0], best_grid_loc[1], best_grid_loc[2]])
    output_env_norm = output_env/(output_env[best_grid_loc[0], best_grid_loc[1], best_grid_loc[2]])
    print(rel_time_shift[best_grid_loc[0], best_grid_loc[1], best_grid_loc[2]])
    print(time_shift[best_grid_loc[0], best_grid_loc[1], best_grid_loc[2]])
    print(dist[best_grid_loc[0], best_grid_loc[1], best_grid_loc[2]])
#=== Part 1: check all the shifted waveform ===#
    '''
    plt.figure()
    for i in range(0,len(stlon)):
        plt.plot(save_amp[best_grid_loc[0], best_grid_loc[1], best_grid_loc[2], i, :], label = stalst.iloc[i,0])
    plt.legend()
    plt.show()
    sys.exit()
#    '''

#=== Part 2: check the raw/shifted waveform ===#
    '''
    plt.figure()
    for i in range(0,len(stlon)):
        plt.plot(save_amp[best_grid_loc[0], best_grid_loc[1], best_grid_loc[2], i, :], 'r-', label = str(stalst.iloc[i,0])+' shifted')
        plt.plot(data[i,:], 'k-', label = str(stalst.iloc[i,0])+' raw')
        plt.legend()
        plt.show()
    sys.exit()
#    '''

    '''
    for xx in range(0,len(x_grid)):
        for yy in range(0,len(y_grid)):
            print("X = ", x_grid[xx], " Y= " , y_grid[yy], output_env[xx,yy,:])
    '''
    f.writelines(stage+' '+str(x_grid[best_grid_loc[0]])+' '+str(y_grid[best_grid_loc[1]])+' '+str(z_grid[best_grid_loc[2]])+'\n')
    f.close()

#=== Plot 2d surface hit count map ===#
#    '''
    X, Y = np.meshgrid(x_grid, y_grid)
    plt.subplot(121)
    plt.imshow(output_env_norm[:,:,0], interpolation = 'nearest', cmap=cm.rainbow, origin = 'lower', vmin = 0.0, vmax = 1.0)  
    plt.subplot(122)
    plt.imshow(output_env_norm[:,:,10], interpolation = 'nearest', cmap=cm.rainbow, origin = 'lower', vmin = 0.0, vmax = 1.0)  
#    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()
#    sys.exit()
#    '''


#=== Plot 3d hit count map ===#
    X, Y, Z = np.meshgrid(x_grid_1, y_grid_1, z_grid_1)
    cmap = plt.get_cmap("rainbow")
    norm = plt.Normalize(0,1)
    norm_c = matplotlib.colors.Normalize(vmin=0, vmax=1)
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.voxels(X, Y, Z, output_env_norm, facecolors=cmap(norm(output_env_norm)), alpha = 0.5)
    ax.set_xlabel('X-direction')
    ax.set_ylabel('Y-direction')
    plt.plot(x_grid[best_grid_loc[0]], y_grid[best_grid_loc[1]], z_grid[best_grid_loc[2]], 'r+')
#    plt.scatter(rel_x,rel_y, s=80, marker = '^', color = 'red')
    '''
    cube = ax.scatter(X, Y, Z, zdir='z', c=output_env_norm, cmap=plt.cm.rainbow)
    fig.colorbar(cube, shrink=0.6, aspect=5)
    '''
    '''
    surf = ax.plot_surface(X,Y,output_env_norm[:,:,0],cmap=cm.rainbow)  
    fig.colorbar(surf, shrink=0.5, aspect=5)
    '''
    m = cm.ScalarMappable(cmap=cm.rainbow, norm = norm_c)
    m.set_array([])
    plt.colorbar(m)
    ax.invert_zaxis()
    plt.show()
    sys.exit()

    fmt='%1.4f', '%1.4f', '%1.5f'
    np.savetxt(str(SaveDirectory)+'/'+dir0+'/1st_'+csta+'_'+rec+'_coarse', amp_c, delimiter=' ', fmt=fmt , newline='\n', header='', encoding=None)

    tEnd = time.time()
    print("It cost %f sec" % (tEnd - tStart))
