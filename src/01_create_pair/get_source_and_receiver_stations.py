import math
from math import radians, cos, sin, asin, sqrt
import numpy as np
import pandas as pd
import os
import sys
def geodistance(lng1,lat1,lng2,lat2):
        lng1, lat1, lng2, lat2 = map(radians, [lng1, lat1, lng2, lat2])
        dlon=lng2-lng1
        dlat=lat2-lat1
        a=sin(dlat/2)**2+cos(lat1)*cos(lat2)*sin(dlon/2)**2
        dis=2*asin(sqrt(a))*6371
        return dis
'''
if len(sys.argv[:]) != 2:
    print('proper usage: python src/01_create_pair/get_source_and_receiver_stations.py stalst')
    sys.exit()
#'''

#=== Parameters Setting ===#
period = []
freq = ["1-5","10-20"]
wavel_ref = 0.3 
wavelength = 0.01


#=== Read-in station list ===#
#in_stalst = sys.argv[1]
in_stalst = 'UGB21_FALL-Locs_204.txt'
sta_list = pd.read_table(in_stalst,sep = ' ', header = None)
#print(sta_list)
stn = sta_list.iloc[:,0]
stlon = sta_list.iloc[:,3]
stlat = sta_list.iloc[:,2]
stlon = np.array(stlon)
stlat = np.array(stlat)
bc = pd.read_table(in_stalst,sep = ' ', header = None)
bcst = bc.iloc[:,0]
bclo = bc.iloc[:,3]
bcla = bc.iloc[:,2]
bclo = np.array(bclo)
bcla = np.array(bcla)
#sys.exit()

#=== Main Code ===#

refer_vel  = 3
#wavel = per*refer_vel
#dir0 = 'db_period_'+str(per)
dir0 = '1-5hz_30m'
if not os.path.isdir(dir0):
        os.mkdir(dir0)
for i in range(0,len(bclo)):
        loc_s_stn = str(bcst[i])
        loc_s_lon = bclo[i]
        loc_s_lat = bcla[i]        
        print(loc_s_stn, loc_s_lon, loc_s_lat)
        #sys.exit()
        dir1 = str(dir0)+'/'+loc_s_stn
        if not os.path.isdir(dir1):
                os.mkdir(dir1)
        for j in range(0,len(bclo)):
                loc_r_stn = str(bcst[j])
                loc_r_lon = bclo[j]
                loc_r_lat = bcla[j]
                #if round(geodistance(loc_s_lon,loc_s_lat,loc_r_lon,loc_r_lat),3) < 3.0*wavel:
                #print(round(geodistance(loc_s_lon,loc_s_lat,loc_r_lon,loc_r_lat),5))
                if round(geodistance(loc_s_lon,loc_s_lat,loc_r_lon,loc_r_lat),8) < 1.5*wavelength:
                        continue
                f2 = open(dir1+'/'+loc_r_stn,'w+')
                for k in range(len(stn)):
                        distance = round(geodistance(stlon[k],stlat[k],loc_r_lon,loc_r_lat),8)
                        if distance <= wavelength:
                                f2.writelines(str(stn[k])+' '+str(stlon[k])+' '+str(stlat[k])+' '+str(distance)+'\n')
                f2.close()
