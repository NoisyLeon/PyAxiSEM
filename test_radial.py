import vmodel

# hetLst = vmodel.heteroLst(anisotype='radial')
# hetLst.add_radial(dr=1., dthe=0.1, rmin=6351., dvsv=0.1)
# hetLst.write('/scratch/summit/life9360/axisem/SOLVER/discr_radial_test.sph')

# hetLst = vmodel.heteroLst(anisotype='iso')
# hetLst.add_iso(dr=1., dthe=0.1, rmin=6351., dvs=0.1)
# hetLst.write('/scratch/summit/life9360/axisem/SOLVER/discr_iso.sph')

# hetLst = vmodel.heteroLst(anisotype='radial', rel=True)
# # hetLst.add_iso(dr=10., dthe=1., rmin=6351., dvs=-10, dvp=-10, drho=0.)
# hetLst.add_radial(dr=1., dthe=1., rmin=6201., dvsv=-30.)
# hetLst.write('/home/leon/code/axisem/SOLVER/discr_radial_vsv_rel.sph')

hetLst = vmodel.heteroLst(anisotype='radial', rel=True)
hetLst.add_radial(dr=1., dthe=1., rmin=6201., dvsv=-40.)
hetLst.write('/home/leon/code/axisem/SOLVER/discr_radial_vsv_rel.sph')

# 
hetLst = vmodel.heteroLst(anisotype='hex', rel=True)
hetLst.add_hex(dr=1., dthe=1., rmin=6201., dvsv=-40., ftheta=0.)

# hetLst.add_hex_ned(dr=1., dthe=1., rmin=6201., dvsv=-30., ftheta=180.)
# hetLst.write('/home/leon/code/axisem/SOLVER/discr_hex_ned_vsv_rel.sph')
hetLst.write('/home/leon/code/axisem/SOLVER/discr_hex_vsv_rel.sph')