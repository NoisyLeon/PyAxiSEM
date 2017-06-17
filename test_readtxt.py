import symdata

dset=symdata.axisemASDF('ak135_iso_10s_mpr.asdf')
dset.readtxt(workingdir='/home/leon/code/axisem/SOLVER/ak135_iso_10s_mpr')
dset.read_event('/home/leon/code/axisem/SOLVER/ak135_iso_10s_mpr/inparam_source')


dset.write2sac(staid='II.AAAA')
tr=dset.get_trace(staid='II.AAAA')
st=dset.get_stream(staid='II.AAAA')