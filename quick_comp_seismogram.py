import matplotlib.pyplot as plt
import numpy as np

f1='/home/leon/code/axisem/SOLVER/ak135_iso_10s_mtr_vsv/Data_Postprocessing/SEISMOGRAMS/JJJJ_II_disp_post_mij_conv0010_E.dat'
f2='/home/leon/code/axisem/SOLVER/ak135_iso_10s_mtr/Data_Postprocessing/SEISMOGRAMS/JJJJ_II_disp_post_mij_conv0010_E.dat'

inArr1 = np.loadtxt(f1)
inArr2 = np.loadtxt(f2)
t1= inArr1[:,0]
t2=inArr2[:,0]
s1 = inArr1[:,1]
s2 = inArr2[:,1]

plt.plot(t1, s1, 'k-', lw=3, label='dvsv=-0.1')
plt.plot(t2,s2, 'r--', lw=3, label='ref')
plt.legend()
plt.show()