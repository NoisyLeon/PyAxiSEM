import symdata, waveform_func

dset1=symdata.axisemASDF('ak135_iso_10s_mrr.asdf')
dset1.readtxt(workingdir='/home/leon/code/axisem/SOLVER/ak135_iso_10s_mrr')
dset1.read_event('/home/leon/code/axisem/SOLVER/ak135_iso_10s_mrr/inparam_source')

dset2=symdata.axisemASDF('ak135_radial_vsv_10s_40_exp.asdf')
# dset2.readtxt(workingdir='/home/leon/code/axisem/SOLVER/ak135_radial_vsv_10s_40_exp')
# dset2.read_event('/home/leon/code/axisem/SOLVER/ak135_radial_vsv_10s_40_exp/inparam_source')
# 
st1=dset1.get_stream(staid='II.AAAA')
st2=dset2.get_stream(staid='II.IIII')

wpair = waveform_func.STpair(st1=st1, st2=st2)
wpair.compare(fmin=0.01, fmax=0.1, vmin=2.5, vmax=5.0)
