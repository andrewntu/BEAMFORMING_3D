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
import scipy


#==== Needed functions ====#
# 1. create the grids 
# 2. single beamforming
# 3. calculate distance and degree transformation

@jit
def geodistance(lon1,lat1,lon2,lat2):
        lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
        dlon=lon2-lon1
        dlat=lat2-lat1
        a=sin(dlat/2)**2+cos(lat1)*cos(lat2)*sin(dlon/2)**2 
        dis=2*asin(sqrt(a))*6371
        return dis
@jit
def theta(lon_r,lat_r,sta_lon,sta_lat,eq_lon,eq_lat,degree_tran):
#        print(lon_r,lat_r,sta_lon,sta_lat,eq_lon,eq_lat,degree_tran)
        x2 = (eq_lon - lon_r)*101.751561277
        y2 = (eq_lat - lat_r)*110.760
        #=== Degree rotation ===#
        new_c1 = x2*math.cos(degree_tran*math.pi/180)-y2*math.sin(degree_tran*math.pi/180)
        new_l1 = x2*math.sin(degree_tran*math.pi/180)+y2*math.cos(degree_tran*math.pi/180)
        #=== Find new location ===#
        c1 = new_c1/101.751561277 + lon_r
        l1 = new_l1/110.760 + lat_r
        #=== Calculate distance and azimuth ===#
        dist_eq_sta  = geodistance(lon_r,lat_r,c1,l1)
        x1 = (c1 - lon_r)*101.751561277
        y1 = (l1 - lat_r)*110.760
        x2 = (sta_lon - lon_r)*101.751561277
        y2 = (sta_lat - lat_r)*110.760        
        dist = (x1*x2+y1*y2)/dist_eq_sta
        return dist

def matrix_beamforming(waveform,omega,coarse_slow_map):
        result_sum = np.zeros((int(len(coarse_slow_map)/len(waveform)),3))
        result_env = np.zeros((int(len(coarse_slow_map)/len(waveform)),len(omega)))
        index = 0
        for i in range(0,len(coarse_slow_map),len(waveform)):
                amp = np.zeros((len(omega)))
                amp = np.array(amp,dtype=complex)
                rescale_slow_map = coarse_slow_map[i:i+len(waveform),2]
#                print(len(rescale_slow_map))
                for x in range(0,len(waveform)):
                        test_omega = omega*rescale_slow_map[x]
                        amp = amp + waveform[x,:]*np.exp(-1j*test_omega)
                tdata_stack_temp = np.fft.irfft(amp[0:int(len(omega)/2)+1])
                tdata_stack = np.real(tdata_stack_temp)#/len(waveform)
                envelope_amp = scipy.signal.hilbert(tdata_stack)
                envelope_amp = np.abs(envelope_amp)
                result_sum[index,0] = coarse_slow_map[i,0]
                result_sum[index,1] = coarse_slow_map[i,1]
                result_sum[index,2] = round(np.max(np.abs(envelope_amp)),8)
                result_env[index,:] = envelope_amp[:]
                '''
                plt.plot(tdata_stack)
                plt.plot(envelope_amp,'r-')
                plt.show()
                '''
                index += 1
        return result_sum, result_env


def stacking(data,omega,evlon,evlat,stlon,stlat,slowness_s,slowness_r,cut_window,baz):
        shift_waveform = []
        waveform,omega = fft_parameter(data)
        amp = np.zeros((len(omega)))
        amp = np.array(amp,dtype=complex)
        for x in range(0,len(waveform)):
                new_waveform = np.zeros((len(omega)))
                new_waveform = np.array(new_waveform,dtype=complex)
                distance_s = theta(lon_s,lat_s,evlon[x],evlat[x],lon_r,lat_r,0)
                distance_r = theta(lon_r,lat_r,stlon[x],stlat[x],lon_s,lat_s,0)
                dtcs_s = round(distance_s,5)*round(slowness_s,5)
                dtcs_r = round(distance_r,5)*round(slowness_r,5)
                dtcs_total = dtcs_s + dtcs_r
                dtcs_total = complex(dtcs_total)
                for i in range(0,len(waveform[x,:])):
                        amp[i] = amp[i]+waveform[x,i]*cmath.exp(-1j*(omega[i]*dtcs_total))
                        new_waveform[i] = waveform[x,i]*cmath.exp(-1j*(omega[i]*dtcs_total))
                tdata_stack_temp = np.fft.irfft(amp[0:int(len(omega)/2)+1])
                tdata_stack = np.real(tdata_stack_temp)/len(waveform)
                tdata_shift_temp = np.fft.irfft(new_waveform[0:int(len(omega)/2)+1])
                tdata_shift = np.real(tdata_shift_temp)/len(waveform)
                shift_waveform.append(tdata_shift)
        return tdata_stack,shift_waveform


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

@jit
def construct_slow_matrix(waveform, matrix, x_grid, y_grid, z_grid):
        inx = 0
        for i in range(0,len(x_grid)):
            for x in range(0,len(y_grid)):
                for y in range(0,len(z_grid)):
                    for z in range(0,len(waveform)):
                                    #distance_r = theta(lon_r,lat_r,stlon[z],stlat[z],lon_s,lat_s,0)
                                    distance_r = theta(lon_r,lat_r,stlon[z],stlat[z],center[0] + x_grid[i]/1000, center[0] + y_grid[x]/1000,0)
                                    dtcs_r = round(distance_r,12)*round(vel_mod,3); 
                                    matrix[i,0] = x_grid[i]
                                    matrix[x,1] = x_grid[x]
                                    matrix[y,2] = z_grid[y]
                                    matrix[z,3] = dtcs_r
                                    inx += 1
 #                                   print(lon_r,lat_r,stlon[z],stlat[z], distance_r)
#                    sys.exit()                                    
        return matrix

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
    global lon_r, lat_r, lon_s, lat_s, event_lat, event_lon, deg, coarse_slowness, fine_slowness, sample, size, envelope, min_v, max_v, freq, vel_mod, center
    global noi_b, sig_b, sig_e, center_pts, x_size, y_size, z_size
    if len(sys.argv[:]) != 5:
            print('proper usage: python ./src/02_beamforming_3d/3d_beamforming.py csta rec freq stage')
            sys.exit()
    tStart = time.time()
    n = 0

#=== Parameter ===#
    SaveDirectory = os.getcwd()
    stalst = str(sys.argv[1])
    csta = str(sys.argv[2])
    freq = str(sys.argv[3])
    stage = str(sys.argv[4])
    # value from Cros et al. (2011) & Vandemeulebrouck et al. (2013).
    vel_mod = 0.13
    degree = np.arange(0, 360.5, 0.5)

#=== Read in/Set up the directory ===#
    print('1-5hz_30m/'+csta+'/'+rsta)
    if os.path.isfile('1-5hz_30m/'+csta+'/'+rsta):
        Receiver = pd.read_table('1-5hz_30m/'+csta+'/'+rsta,header = None,sep = ' ')
    else:
        print("No selected Receiver!!")
        sys.exit()
    if len(Receiver) < 4:
        print("Not enough receiver station!!")
        sys.exit()
    if csta == rsta:
        sys.exit()
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

#=== Load Data ===#
    f = open('Output/freq_'+str(freq)+'/'+csta+'_'+rsta+'_'+str(freq)+'.txt','w+')
    for r in range(0,len(Receiver)):
        rec = str(Receiver.iloc[r,0])
        pair = csta+'-'+rec
        ccor = glob.glob(SaveDirectory+'/STACK_ERUPTION/'+csta+'/'+freq+'hz/'+pair+'.ZZ.stage.'+stage+'.sac.norm')
        if len(ccor) == 1:
                st = read(ccor[0],debug_header=True)
                sta.append(rec)
                tr = st[0]
                sample = tr.stats.sampling_rate
                center_pts = int(len(tr.data)/2)
                trace.append(tr)
                f.writelines(pair+' '+str(tr.stats.sac.dist)+'\n')
                normal = np.max(np.abs(tr.data[center_pts:]))
                tr.data[:center_pts] = 0.0
                data.append(tr.data/normal)
                stlon.append(tr.stats.sac.stlo)
                stlat.append(tr.stats.sac.stla)
                stel.append(tr.stats.sac.stel)
                dist.append(tr.stats.sac.dist)
                evlon.append(tr.stats.sac.evlo)
                evlat.append(tr.stats.sac.evla)
                event_lon = tr.stats.sac.evlo
                event_lat = tr.stats.sac.evla
                baz = tr.stats.sac.baz
                if str(Receiver.iloc[r,0]) == rec:
                    lon_r = Receiver.iloc[r,1]
                    lat_r = Receiver.iloc[r,2]
                    lon_s = tr.stats.sac.evlo
                    lat_s = tr.stats.sac.evla
    f.close()
    
    print("center station loc:",rsta,lon_r,lat_r)
    print("event loc:",csta,lon_s,lat_s)
    dist_stack = geodistance(lon_s, lat_s, lon_r, lat_s)
    sig_b = int((dist_stack/vel_mod)*sample + center_pts)
    sig_e = int((dist_stack/vel_mod)*sample + center_pts)
    noi_b = sig_e
    data = np.array(data)
    evlon = np.array(evlon)
    evlat = np.array(evlat)
    stlon = np.array(stlon)
    stlat = np.array(stlat)
    dist = np.array(dist)

#=== Create Grid ===# !!! The Old Faithful cone is the referenced center.
    #old_faithful = [-110.828211897200987, 44.460437153479248, 2240]
    center = [-110.828211897200987, 44.460437153479248, 2240]
    # !!! The grid size in meters is 70*70*20
    x_range = 70
    y_range = 70
    z_range = 20
    x_size = 5
    y_size = 5
    z_size = 2
    x_grid = np.arange(-1*(x_range/2), (x_range/2) + x_size, x_size)
    y_grid = np.arange(-1*(y_range/2), (y_range/2) + y_size, y_size)
    z_grid = np.arange( 0, z_range + z_size, z_size)

#=== Slowness range and searching increment ===#
    coarse_slow_map = np.zeros((len(x_grid)*len(y_grid)*len(z_grid)*len(data),4))

#=== Tranfer Time Series to Frequency domain ===#
    fdata,omega = fft_parameter(data)

#=== First Coarse Grid Search ===#
    print('Beamforming Start !')
    coarse_slow_map = construct_slow_matrix(data,coarse_slow_map,x_grid, y_grid, z_grid)
    print(coarse_slow_map)
    sys.exit()
    cut_win = int(1000)
    amp_c, env_c  = matrix_beamforming(fdata,omega,coarse_slow_map)

#=== First Fine Grid Search ===#

    max_r_index = int(np.argmax(amp_c[:,2])*len(data))
    fine_low_s =  coarse_slow_map[max_r_index,0] - 0.05
    fine_high_s = coarse_slow_map[max_r_index,0] + 0.05
    fine_low_r =  coarse_slow_map[max_r_index,1] - 0.05
    fine_high_r = coarse_slow_map[max_r_index,1] + 0.05
    if fine_low_s < 0.0:
            fine_low_s = 0.0
    if fine_low_r < 0.0:
            fine_low_r = 0.0
    fine_slowness_s_f = np.arange(fine_low_s,fine_high_s,float(fine_slow_step))
    fine_slowness_r_f = np.arange(fine_low_r,fine_high_r,float(fine_slow_step))
    fine_slow_map = np.zeros((len(fine_slowness_s_f)*len(fine_slowness_r_f)*len(data),3))
    fine_slow_map = construct_slow_matrix(data,fine_slow_map,fine_slowness_s_f,fine_slowness_r_f)
    amp, env = matrix_beamforming(fdata,omega,fine_slow_map)
    best_slowness_s = amp[np.argmax(amp[:,2]),0]
    best_slowness_r = amp[np.argmax(amp[:,2]),1]
    env[np.argmax(np.abs(amp[:,2])),:] = env[np.argmax(np.abs(amp[:,2])),:]/np.max(env[np.argmax(np.abs(amp[:,2])),:])

    re_data_2 = np.array(re_data_2)
    amp_data_2 = np.array(amp_data_2)
    fmt='%1.4f', '%1.4f', '%1.5f'
    fund_amp,fund_stack = stacking(data,omega,evlon,evlat,stlon,stlat,best_slowness_s,best_slowness_r,cut_win,0)
    fund_stack = np.array(fund_stack)
    np.savetxt(str(SaveDirectory)+'/'+dir0+'/1st_'+csta+'_'+rec+'_coarse', amp_c, delimiter=' ', fmt=fmt , newline='\n', header='', encoding=None)
    idx_1 = 0
    idx_2 = 0
    f_fund_env = open('Output/freq_'+freq+'/1st_'+csta+'_'+rec+'_env','w+')
    f_fund_slow = open('Output/freq_'+freq+'/1st_'+csta+'_'+rec+'_slow','w+')
    f_fund_shift = open('Output/freq_'+freq+'/1st_'+csta+'_'+rec+'_shift','w+')
    f_fund = open('Output/freq_'+freq+'/1st_'+csta+'_'+rec+'_wave','w+')
    for i in range(0,cut_win):
            sec = round(i/sample,3)
            f_fund_env.writelines(str(sec)+' '+str(env[np.argmax(np.abs(amp[:,2])),center_pts+i]/np.max(env[np.argmax(np.abs(amp[:,2])),:]))+'\n')
            f_fund.writelines(str(sec)+' '+str(fund_amp[center_pts+i]/np.max(np.abs(fund_amp)))+'\n')
    for i in range(0,len(trace_fund)):
            trace_fund[i].data = fund_stack[i,:]
            trace_fund[i].write('bp_'+str(freq)+'/'+csta+'_'+rec+'_'+sta_src[i]+'_'+sta[i]+'_shift_fund.SAC',format='SAC')
            f_fund_shift.writelines('bp_'+str(freq)+'/'+csta+'_'+rec+'_'+sta_src[i]+'_'+sta[i]+'_shift_fund.SAC'+' '+str(dist[i])+'\n')
    for i in range(0,len(coarse_slowness_r)):
            for x in range(0,1000):
                    best_time = round(x/sample,3)
                    value = env_c[i,center_pts+x]
                    f_fund_slow.writelines(str(best_time)+' '+str(coarse_slowness_r[i])+' '+str(value)+'\n')
    max_signal = np.max(np.abs(fund_amp[sig_b:sig_e]))
    #noise = np.sqrt(np.mean(fund_amp[noi_b:]**2)+np.mean(fund_amp[center_pts:sig_b]**2))
    noise = np.sqrt(np.mean(fund_amp[noi_b:]**2))
    snr = round(max_signal/noise,3)
    if float(freq)>=4:
        f_out.writelines(csta+' '+str(lon_s)+' '+str(lat_s)+' '+rec+' '+str(lon_r)+' '+str(lat_r)+' '+str(best_slowness_s)+' '+str(best_slowness_r)+' '+str(snr)+'\n')
    elif float(freq)<4:
        f_sus_fund.writelines(csta+' '+str(lon_s)+' '+str(lat_s)+' '+rec+' '+str(lon_r)+' '+str(lat_r)+' '+str(best_slowness_s)+' '+str(best_slowness_r)+' '+str(snr)+'\n')

    f_out.close()
    f_sus_fund.close()
    f_fund_env.close()
    f_fund_slow.close()        
    f_fund_shift.close()
    f_fund.close()
    print("Only have 1 mode !! Abort!!")
    sys.exit()

# Tranfer Time Series to Frequency domain====

    fdata_2,omega = fft_parameter(data_2)
    coarse_slow_map = np.zeros((len(coarse_slowness_s)*len(coarse_slowness_r)*len(data_2),3))
    coarse_slow_map = construct_slow_matrix(data_2,coarse_slow_map,coarse_slowness_s,coarse_slowness_r)

# Second Coarse Grid Search==================

    print('Second Beamforming Start !')
    amp_2_c, env_2_c = matrix_beamforming(fdata_2,omega,coarse_slow_map)

# Second Fine Grid Search====================

    max_r_index_2 = int(np.argmax(amp_2_c[:,2])*len(data_2)) #int(np.argmax(amp_c[:,2])*len(data))
    fine_low_s =  coarse_slow_map[max_r_index_2,0] - 0.05
    fine_high_s = coarse_slow_map[max_r_index_2,0] + 0.05
    fine_low_r =  coarse_slow_map[max_r_index_2,1] - 0.05
    fine_high_r = coarse_slow_map[max_r_index_2,1] + 0.05
    if fine_low_s < 0.0:
            fine_low_s = 0.0
    if fine_low_r < 0.0:
            fine_low_r = 0.0
    fine_slowness_s_s = np.arange(fine_low_s,fine_high_s,float(fine_slow_step))
    fine_slowness_r_s = np.arange(fine_low_r,fine_high_r,float(fine_slow_step))
    fine_slow_map_2 = np.zeros((len(fine_slowness_s_s)*len(fine_slowness_r_s)*len(data),3))
    fine_slow_map_2 = construct_slow_matrix(data_2,fine_slow_map_2,fine_slowness_s_s,fine_slowness_r_s)
    amp_2, env_2 = matrix_beamforming(fdata_2, omega, fine_slow_map_2)
    best_slowness_s_2 = amp_2[np.argmax(amp_2[:,2]),0]
    best_slowness_r_2 = amp_2[np.argmax(amp_2[:,2]),1]
    env_2[np.argmax(np.abs(amp_2[:,2])),:] = env_2[np.argmax(np.abs(amp_2[:,2])),:]/np.max(env_2[np.argmax(np.abs(amp_2[:,2])),:])
#        print(best_slowness_s,best_slowness_r,best_slowness_s_2,best_slowness_r_2)
#        sys.exit()

# Output with normalization=========
    f_fund_env = open('Output/freq_'+freq+'/1st_'+csta+'_'+rec+'_env','w+')
    f_high_env = open('Output/freq_'+freq+'/2nd_'+csta+'_'+rec+'_env','w+')
    f_fund_slow = open('Output/freq_'+freq+'/1st_'+csta+'_'+rec+'_slow','w+')
    f_high_slow = open('Output/freq_'+freq+'/2nd_'+csta+'_'+rec+'_slow','w+')
    f_fund_shift = open('Output/freq_'+freq+'/1st_'+csta+'_'+rec+'_shift','w+')
    f_high_shift = open('Output/freq_'+freq+'/2nd_'+csta+'_'+rec+'_shift','w+')
    idx_1 = 0
    idx_2 = 0
    idx_3 = 0
    idx_4 = 0
#        print(np.argmax(np.abs(env[np.argmax(amp[:,2]),:])), np.argmax(np.abs(env_2[np.argmax(amp_2[:,2]),:])), best_slowness_r, best_slowness_r_2)
    if np.argmax(np.abs(env[np.argmax(amp[:,2]),:])) < np.argmax(np.abs(env_2[np.argmax(amp_2[:,2]),:])) and best_slowness_r < best_slowness_r_2 and abs(best_slowness_r - best_slowness_r_2) >= 0.1: # 1 is higher moode
            env_f = env_2[np.argmax(amp_2[:,2]),:]
            env_h = env[np.argmax(amp[:,2]),:]
            slow_f_s = best_slowness_s_2
            slow_f_r = best_slowness_r_2
            slow_h_s = best_slowness_s
            slow_h_r = best_slowness_r
            data_f = data_2
            data_h = data
            np.savetxt(str(SaveDirectory)+'/'+dir0+'/1st_'+csta+'_'+rec+'_coarse', amp_2_c, delimiter=' ', fmt=fmt , newline='\n', header='', encoding=None)
            np.savetxt(str(SaveDirectory)+'/'+dir0+'/2nd_'+csta+'_'+rec+'_coarse', amp_c, delimiter=' ', fmt=fmt , newline='\n', header='', encoding=None)
            for i in range(0,len(coarse_slowness_r)):
                    for x in range(0,1000):
                            best_time = round(x/sample,3)
                            value_2 = env_2_c[i,center_pts+x]
                            f_fund_slow.writelines(str(best_time)+' '+str(coarse_slowness_r[i])+' '+str(value_2)+'\n')
                            best_time = round(x/sample,3)
                            value = env_c[i,center_pts+x]
                            f_high_slow.writelines(str(best_time)+' '+str(coarse_slowness_r[i])+' '+str(value)+'\n')

            max_signal_fund = np.max(env_2[np.argmax(np.abs(amp_2[:,2])),center_pts:])
            #noise_fund = np.sqrt(np.mean(env_2[np.argmax(np.abs(amp_2[:,2])),center_pts:sig_b]**2)+np.mean(env_2[np.argmax(np.abs(amp[:,2])),noi_b:]**2))
            noise_fund = np.sqrt(np.mean(env_2[np.argmax(np.abs(amp_2[:,2])),noi_b:]**2))
            snr_fund = round(max_signal_fund/noise_fund,3)

            max_signal_high = np.max(env[np.argmax(np.abs(amp[:,2])),center_pts:])
            #noise_high = np.sqrt(np.mean(env[np.argmax(np.abs(amp[:,2])),center_pts:sig_b]**2)+np.mean(env[np.argmax(np.abs(amp[:,2])),noi_b:]**2))
            noise_high = np.sqrt(np.mean(env[np.argmax(np.abs(amp[:,2])),noi_b:]**2))
            snr_high = round(max_signal_high/noise_high,3)

            f_out.writelines(csta+' '+str(lon_s)+' '+str(lat_s)+' '+rec+' '+str(lon_r)+' '+str(lat_r)+' '+str(best_slowness_s_2)+' '+str(best_slowness_r_2)+' '+str(snr_fund)+'\n')
            f_all_high.writelines(csta+' '+str(lon_s)+' '+str(lat_s)+' '+rec+' '+str(lon_r)+' '+str(lat_r)+' '+str(best_slowness_s)+' '+str(best_slowness_r)+' '+str(snr_high)+'\n')
            f_out.close()
            f_all_high.close()
            f_fund_slow.close()
            f_high_slow.close()
    elif np.argmax(np.abs(env_2[np.argmax(amp_2[:,2]),:])) < np.argmax(np.abs(env[np.argmax(amp[:,2]),:])) and best_slowness_r_2 < best_slowness_r and abs(best_slowness_r - best_slowness_r_2) >= 0.1: # 2 is higher moode
            env_f = env[np.argmax(amp[:,2]),:]
            env_h = env_2[np.argmax(amp_2[:,2]),:]
            slow_f_s = best_slowness_s
            slow_f_r = best_slowness_r
            slow_h_s = best_slowness_s_2
            slow_h_r = best_slowness_r_2
            data_f = data
            data_h = data_2
            np.savetxt(str(SaveDirectory)+'/'+dir0+'/1st_'+csta+'_'+rec+'_coarse', amp_c, delimiter=' ', fmt=fmt , newline='\n', header='', encoding=None)
            np.savetxt(str(SaveDirectory)+'/'+dir0+'/2nd_'+csta+'_'+rec+'_coarse', amp_2_c, delimiter=' ', fmt=fmt , newline='\n', header='', encoding=None)
            for i in range(0,len(coarse_slowness_r)):
                    for x in range(0,1000):
                            best_time = round(x/sample,3)
                            value_2 = env_2_c[i,center_pts+x]
                            f_high_slow.writelines(str(best_time)+' '+str(coarse_slowness_r[i])+' '+str(value_2)+'\n')
                            best_time = round(x/sample,3)
                            value = env_c[i,center_pts+x]
                            f_fund_slow.writelines(str(best_time)+' '+str(coarse_slowness_r[i])+' '+str(value)+'\n')

            max_signal_fund = np.max(env[np.argmax(np.abs(amp[:,2])),center_pts:])
            #noise_fund = np.sqrt(np.mean(env[np.argmax(np.abs(amp[:,2])),center_pts:sig_b]**2)+np.mean(env[np.argmax(np.abs(amp[:,2])),noi_b:]**2))
            noise_fund = np.sqrt(np.mean(env[np.argmax(np.abs(amp[:,2])),noi_b:]**2))
            snr_fund = round(max_signal_fund/noise_fund,3)

            max_signal_high = np.max(env_2[np.argmax(np.abs(amp_2[:,2])),center_pts:])
            #noise_high = np.sqrt(np.mean(env_2[np.argmax(np.abs(amp_2[:,2])),center_pts:sig_b]**2)+np.mean(env_2[np.argmax(np.abs(amp[:,2])),noi_b:]**2))
            noise_high = np.sqrt(np.mean(env_2[np.argmax(np.abs(amp_2[:,2])),noi_b:]**2))
            snr_high = round(max_signal_high/noise_high,3)

            f_out.writelines(csta+' '+str(lon_s)+' '+str(lat_s)+' '+rec+' '+str(lon_r)+' '+str(lat_r)+' '+str(best_slowness_s)+' '+str(best_slowness_r)+' '+str(snr_fund)+'\n')
            f_all_high.writelines(csta+' '+str(lon_s)+' '+str(lat_s)+' '+rec+' '+str(lon_r)+' '+str(lat_r)+' '+str(best_slowness_s_2)+' '+str(best_slowness_r_2)+' '+str(snr_high)+'\n')
    
            f_out.close()
            f_all_high.close()
            f_fund_slow.close()
            f_high_slow.close()
    else:
            print("Still only have good 1 mode!!")
            fund_amp,fund_stack = stacking(data,omega,evlon,evlat,stlon,stlat,best_slowness_s,best_slowness_r,cut_win,0)
            fund_stack = np.array(fund_stack)
            np.savetxt(str(SaveDirectory)+'/'+dir0+'/1st_'+csta+'_'+rec+'_coarse', amp_c, delimiter=' ', fmt=fmt , newline='\n', header='', encoding=None)
            idx_1 = 0
            idx_2 = 0
            f_fund_env = open('Output/freq_'+freq+'/1st_'+csta+'_'+rec+'_env','w+')
            f_fund_slow = open('Output/freq_'+freq+'/1st_'+csta+'_'+rec+'_slow','w+')
            f_fund_shift = open('Output/freq_'+freq+'/1st_'+csta+'_'+rec+'_shift','w+')
            f_fund = open('Output/freq_'+freq+'/1st_'+csta+'_'+rec+'_wave','w+')
            for i in range(0,cut_win):
                    sec = round(i/sample,3)
                    f_fund_env.writelines(str(sec)+' '+str(env[np.argmax(np.abs(amp[:,2])),center_pts+i]/np.max(env[np.argmax(np.abs(amp[:,2])),:]))+'\n')
                    f_fund.writelines(str(sec)+' '+str(fund_amp[center_pts+i]/np.max(np.abs(fund_amp)))+'\n')
            for i in range(0,len(trace_fund)):
                    trace_fund[i].data = fund_stack[i,:]
                    trace_fund[i].write('bp_'+str(freq)+'/'+csta+'_'+rec+'_'+sta_src[i]+'_'+sta[i]+'_shift_fund.SAC',format='SAC')
                    f_fund_shift.writelines('bp_'+str(freq)+'/'+csta+'_'+rec+'_'+sta_src[i]+'_'+sta[i]+'_shift_fund.SAC'+' '+str(dist[i])+'\n')
            for i in range(0,len(coarse_slowness_r)):
                    for x in range(0,1000):
                            best_time = round(x/sample,3)
                            value = env_c[i,center_pts+x]
                            f_fund_slow.writelines(str(best_time)+' '+str(coarse_slowness_r[i])+' '+str(value)+'\n')
            max_signal = np.max(env[np.argmax(np.abs(amp[:,2])),center_pts:])
            #noise = np.sqrt(np.mean(env[np.argmax(np.abs(amp[:,2])),noi_b:]**2)+np.mean(env[np.argmax(np.abs(amp[:,2])),center_pts:sig_b]**2))
            noise = np.sqrt(np.mean(env[np.argmax(np.abs(amp[:,2])),noi_b:]**2))
            snr = round(max_signal/noise,3)
            f_out.writelines(csta+' '+str(lon_s)+' '+str(lat_s)+' '+rec+' '+str(lon_r)+' '+str(lat_r)+' '+str(best_slowness_s)+' '+str(best_slowness_r)+' '+str(snr)+'\n')

            f_out.close()
            sys.exit()

    for i in range(0,cut_win):
            sec = round(i/sample,3)
            f_high_env.writelines(str(sec)+' '+str(env_h[center_pts+i]/np.max(env_h))+'\n')
            f_fund_env.writelines(str(sec)+' '+str(env_f[center_pts+i]/np.max(env_f))+'\n')
    f_high_env.close()
    f_fund_env.close()        

    high_amp,high_stack = stacking(data_h,omega,evlon,evlat,stlon,stlat,slow_h_s,slow_h_r,cut_win,0)
    fund_amp,fund_stack = stacking(data_f,omega,evlon,evlat,stlon,stlat,slow_f_s,slow_f_r,cut_win,0)
    high_stack = np.array(high_stack)
    fund_stack = np.array(fund_stack)

    for i in range(0,len(trace_high)):
            trace_high[i].data = high_stack[i,:]
            trace_fund[i].data = fund_stack[i,:]
            trace_high[i].write('bp_'+str(freq)+'/'+csta+'_'+rec+'_'+sta_src[i]+'_'+sta[i]+'_shift_high.SAC',format='SAC')
            trace_fund[i].write('bp_'+str(freq)+'/'+csta+'_'+rec+'_'+sta_src[i]+'_'+sta[i]+'_shift_fund.SAC',format='SAC')
            f_high_shift.writelines('bp_'+str(freq)+'/'+csta+'_'+rec+'_'+sta_src[i]+'_'+sta[i]+'_shift_high.SAC'+' '+str(dist[i])+'\n')        
            f_fund_shift.writelines('bp_'+str(freq)+'/'+csta+'_'+rec+'_'+sta_src[i]+'_'+sta[i]+'_shift_fund.SAC'+' '+str(dist[i])+'\n')                
    f_high_shift.close()                
    f_fund_shift.close()                
    
    f_high = open('Output/freq_'+freq+'/2nd_'+csta+'_'+rec+'_wave','w+')
    f_fund = open('Output/freq_'+freq+'/1st_'+csta+'_'+rec+'_wave','w+')

    for i in range(0,cut_win):
            sec = round(i/sample,3)
            f_high.writelines(str(sec)+' '+str(high_amp[center_pts+i]/np.max(np.abs(high_amp)))+'\n')
            f_fund.writelines(str(sec)+' '+str(fund_amp[center_pts+i]/np.max(np.abs(fund_amp)))+'\n')
    f_high.close()
    f_fund.close()


    tEnd = time.time()
    print("It cost %f sec" % (tEnd - tStart))
