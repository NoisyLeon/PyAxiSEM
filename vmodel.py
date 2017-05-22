import numpy as np
# import pysat

class hetero(object):
    def __init__(self, radius, theta, dvpv, dvsv, drho, dvph=None, dvsh=None, deta=None, ftheta=None, fphi=None):
        self.radius = radius
        self.theta  = theta
        self.dvpv    = dvpv
        self.dvph    = dvph
        self.dvsv    = dvsv
        self.dvsh    = dvsh
        self.drho   = drho
        self.deta   = deta
        self.ftheta = ftheta
        self.fphi   = fphi
        

class heteroLst(object):
    """
    An object to handle a list of ahetero object for SW4
    """
    def __init__(self, heteros=None, rel=True, anisotype='hex'):
        self.heteros = []
        if isinstance(heteros, hetero):
            heteros  = [heteros]
        if heteros:
            self.heteros.extend(heteros)
        self.rel        = rel
        self.anisotype  = anisotype
        self.Nhetero    = len(self.heteros)

    def __add__(self, other):
        """
        Add two heteroLst with self += other.
        """
        if isinstance(other, hetero):
            other = heteroLst([other])
        if not isinstance(other, heteroLst):
            raise TypeError
        heteros = self.heteros + other.heteros
        return self.__class__(heteros=heteros)

    def __len__(self):
        """
        Return the number of Traces in the heteroLst object.
        """
        return len(self.heteros)

    def __getitem__(self, index):
        """
        __getitem__ method of heteroLst objects.
        :return: hetero objects
        """
        if isinstance(index, slice):
            return self.__class__(heteros=self.heteros.__getitem__(index))
        else:
            return self.heteros.__getitem__(index)

    def append(self, inhetero):
        """
        Append a single hetero object to the current heteroLst object.
        """
        if isinstance(inhetero, hetero):
            self.heteros.append(inhetero)
        else:
            msg = 'Append only supports a single Blockmodel object as an argument.'
            raise TypeError(msg)
        return self
    
    def add(self, radius, theta, dvpv, dvsv, drho, dvph=None, dvsh=None, deta=None, ftheta=None, fphi=None):
        
        if (dvph==None or dvsh==None or deta==None) and ( self.anisotype != 'iso'):
            raise ValueError('For anisotropy, dvph, dvsh, deta needs to be specified!')
        if ftheta == None: ftheta = theta
        if (ftheta==None or fphi==None ) and self.anisotype == 'hex':
            raise ValueError('For anisotropy, fast axis orientation needs to be specified!')
        self.heteros.append(hetero(radius=radius, theta=theta, dvpv=dvpv, dvsv=dvsv, drho=drho, dvph=dvph, dvsh=dvsh,
                                deta=deta, ftheta=ftheta, fphi=fphi))
        self.Nhetero += 1
        return
    
    def add_radial(self, dr, dthe, themin=0., themax=180., rmin=0., rmax=6371., dvpv=0., dvsv=0., drho=0., dvph=0., dvsh=0., deta=0.):
        if self.anisotype !='radial':
            raise ValueError('Anisotropic type is NOT radial !')
        if dvpv==0. and dvsv==0. and dvph==0. and dvsh==0.:
            raise ValueError('At least one of the velocity parameter need to be non-zero !')
        if not self.rel and (dvpv==0. or dvsv==0. or dvph==0. or dvsh==0. or drho==0. or deta==0.):
            raise ValueError('Absolute value option! All velocity parameters need to be specified!')
        rArr    = np.arange(rmin, rmax+dr, dr)
        thetaArr= np.arange(themin, themax+dthe, dthe)
        for r in rArr:
            for theta in thetaArr:
                if r <0. or r > 6371.:
                    print 'Warning: r =', r
                    continue
                if theta <0. or theta > 180.:
                    print 'Warning: theta =', theta
                    continue
                self.add(radius=r, theta=theta, dvpv=dvpv, dvsv=dvsv, drho=drho, dvph=dvph, dvsh=dvsh, deta=deta)
        return
    
    def add_iso(self, dr, dthe, themin=0., themax=180., rmin=0., rmax=6371., dvs=0., dvp=0., drho=0.):
        if self.anisotype !='iso':
            raise ValueError('Anisotropic type is NOT iso !')
        if dvp==0. and dvs==0. and drho==0.:
            raise ValueError('At least one of the velocity parameter need to be non-zero !')
        if not self.rel and (dvp==0. or dvs==0. or drho==0.):
            raise ValueError('Absolute value option! All velocity parameters need to be specified!')
        rArr    = np.arange(rmin, rmax+dr, dr)
        thetaArr= np.arange(themin, themax+dthe, dthe)
        for r in rArr:
            for theta in thetaArr:
                if r <0. or r > 6371.:
                    print 'Warning: r =', r
                    continue
                if theta <0. or theta > 180.:
                    print 'Warning: theta =', theta
                    continue
                self.add(radius=r, theta=theta, dvpv=dvp, dvsv=dvs, drho=drho)
        return
    
    # def add_aniso(self, rmin, rmax, dr, themin, themax, dthe, dip, strike):
        
    def write(self, outfname):
        with open(outfname, 'wb') as f:
            f.writelines('%d \n' %self.Nhetero)
            for het in self.heteros:
                if self.anisotype=='radial': tempstr = '%E %E %E %E %E %E %E %E \n' %(het.radius, het.theta,
                              het.dvpv, het.dvsv, het.dvph, het.dvsh, het.drho, het.deta)
                elif self.anisotype=='iso': tempstr = '%E %E %E %E %E\n' %(het.radius, het.theta, het.dvpv, het.dvsv, het.drho)
                f.writelines(tempstr)
        return
    
    
    