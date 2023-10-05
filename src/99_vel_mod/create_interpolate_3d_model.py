import numpy as np
import os
import sys

def interploate_model(in_vel, wanted_depth, inc_depth, out_vel):
    for i in range(0, len(wanted_depth)):
        for x in range(0, len(in_vel)-1):
#            print(wanted_depth[i], in_vel[x,0], in_vel[x,1], in_vel[x+1,0], in_vel[x+1,1])
            if wanted_depth[i] >= in_vel[x,0] and wanted_depth[i] <= in_vel[x+1, 0]:
                #print(wanted_depth[i], in_vel[x,0], in_vel[x,1], in_vel[x+1,0], in_vel[x+1,1])
                out_vel[i] = in_vel[x, 1] + (wanted_depth[i]-in_vel[x,0])*(in_vel[x+1,1] - in_vel[x,1])/(in_vel[x+1,0] - in_vel[x,0])
    return out_vel

in_mod = './Vel_mod/1D_Cros'
in_vel = np.loadtxt(in_mod, delimiter = ' ')
inc = 0.01
wanted_depth = np.arange(0, 0.2 + inc, 0.01)
out_vel = np.zeros((len(wanted_depth)))
out_vel = interploate_model(in_vel, wanted_depth, inc, out_vel)

print(out_vel)


