import numpy as np
import obspy
from obspy.core import read
import sys
import os
import glob
import pandas as pd


#'''
if len(sys.argv[:]) != 6:
    print('proper usage: python src/extract_max_value.py csta stage freq')
    sys.exit()
#'''


csta = str(sys.argv[1])
stage = str(sys.argv[2])
freq = str(sys.argv[3])
infold = str(sys.argv[4])
stalst = str(sys.argv[5])

station_list = pd.read_table(stalst, header = None, sep = ' ')

f = open(csta+'_time_amp_stage.'+stage+'.'+freq+'hz.txt','w+')
#for wave in glob.glob(infold+'/'+csta+'/'+freq+'hz/'+csta+'-*.ZZ.stage.'+stage+'.sac.norm'):
for i in range(0, len(station_list)):
    rec = station_list.iloc[i,0]
    pair = csta+'-'+str(rec)
    wave = glob.glob(infold+'/'+csta+'/'+freq+'hz/'+pair+'.ZZ.stage.'+stage+'.sac.norm')
    if len(wave) < 1:
        continue
    st = read(wave[0])
    tr = st[0]
    sample = round(1/tr.stats.sampling_rate,2)
    t = np.arange(-1*sample*int(len(tr.data)/2), sample*int(len(tr.data)/2) + sample, sample)
    min_arg = np.argmin(tr.data)
    max_arg = np.argmax(tr.data)
    min_val = np.min(tr.data)
    max_val = np.max(tr.data)
    if abs(min_val) > max_val:
        final_arg = min_arg
        final_amp = min_val
    else:
        final_arg = max_arg
        final_amp = max_val

    f.writelines(pair+' '+str(tr.stats.sac.stlo)+' '+str(tr.stats.sac.stla)+' '+str(round(t[final_arg],2))+' '+str(final_amp)+'\n')
        
f.close()
