import stations

SLst = stations.StaLst()
# SLst.homo_stalst(minlat=0., Nlat=90, minlon=0., Nlon=1, dlat=1., dlon=1.)
SLst.read('/scratch/summit/life9360/axisem/SOLVER/STATIONS')
SLst.write('./')