import vmodel

# hetLst = vmodel.heteroLst(anisotype='radial')
# hetLst.add_radial(dr=1., dthe=0.1, rmin=6351., dvsv=0.1)
# hetLst.write('/scratch/summit/life9360/axisem/SOLVER/discr_radial_test.sph')

# hetLst = vmodel.heteroLst(anisotype='iso')
# hetLst.add_iso(dr=1., dthe=0.1, rmin=6351., dvs=0.1)
# hetLst.write('/scratch/summit/life9360/axisem/SOLVER/discr_iso.sph')

hetLst = vmodel.heteroLst(anisotype='iso', rel=True)
hetLst.add_iso(dr=1., dthe=0.1, rmin=6351., dvs=-10, dvp=-10, drho=0.)
hetLst.write('/scratch/summit/life9360/axisem/SOLVER/discr_iso_rel.sph')