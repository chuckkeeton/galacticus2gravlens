############################################################
"""
Code to take semi-analytic models from Galacticus and produce input
files for gravlens.

Developed by Sean Brennan and Chuck Keeton, drawing on earlier code
by Andrew Benson and Anthony Pullen.

This version 2019/3/27.
"""
############################################################

import numpy as np
from astropy.constants import G, c
from astropy.cosmology import FlatLambdaCDM

############################################################
"""
Cosmological quantities.  Assumes flat LCDM cosmology.

Required inputs:
  zlens = lens redshift
  zsrc = source redshift
Optional inputs:
  H0 = Hubble constant in units of km/s/Mpc [default 70]
  Om0 = density parameter [default 0.3]

Class parameters:
  zlens
  zsrc
  D_l = angular diameter distance to lens
  D_s = angular diameter distance to source
  D_ls = angular diameter distance from lens to source
  Sigma_crit = critical surface mass density for lensing, in Msun/Mpc^2
  as_kpc = arcseconds per kpc
  as_Mpc = arcseconds per Mpc
"""

class cosmology:
    def __init__(self, zlens, zsrc, H0=70, Om0=0.3):
        # assume a flat LCDM cosmology
        cosmo = FlatLambdaCDM(H0=H0,Om0=Om0)
        # compute cosmology quantities
        self.zlens = zlens
        self.zsrc = zsrc
        self.D_l = cosmo.angular_diameter_distance(zlens) 
        self.D_s = cosmo.angular_diameter_distance(zsrc) 
        self.D_ls = cosmo.angular_diameter_distance_z1z2(zlens,zsrc)
        # critical density for lensing, in units of Msun/Mpc^2
        E_cr = ((c**2)/(4.0*np.pi*G))*(self.D_s/(self.D_ls*self.D_l))
        self.Sigma_crit = E_cr.to('Msun/Mpc^2')
        # arcsec per kpc and Mpc
        self.as_kpc = cosmo.arcsec_per_kpc_proper(zlens)
        self.as_Mpc = self.as_kpc.to('arcsec/Mpc')

############################################################
"""
Take an array of 3d positions and compute a projection.

Required input:
  xyzarr[n,3] = array of 3d positions
Optional input:
  angles[2] = polar angles [theta,phi]; if not specified, will be set randomly

Output:
  xyarr[n,2] = array of projected 2d positions
"""

def project_positions(xyzarr, angles=[]):

    if len(angles)<2:
        # pick random angles on unit sphere
        u = np.random.uniform(0,1)
        v = np.random.uniform(0,1)    
        phi = 2*np.pi*u
        theta = np.arccos(2*v - 1)
    else:
        theta = angles[0]
        phi = angles[1]

    # construct vectors to handle the projection
    sinTheta = np.sin(theta)
    cosTheta = np.cos(theta)
    sinPhi = np.sin(phi)
    cosPhi = np.cos(phi)
    nhat = np.array([sinTheta*cosPhi,sinTheta*sinPhi,cosTheta])
    that = np.array([cosTheta*cosPhi,cosTheta*sinPhi,-sinTheta])
    phat = np.array([-sinPhi,cosPhi,0.0])

    # do the projection
    rp = xyzarr - np.outer(np.dot(xyzarr,nhat),nhat)
    rpx = np.dot(rp,phat)
    rpy = np.dot(rp,that)

    return np.transpose(np.vstack((rpx,rpy)))

############################################################
"""
Process Galacticus data into a set of lensing parameters.

Required inputs:
  galacticus_data = structure read from Galacticus file
  cosmo = instance of cosmology class
Optional inputs:
  npop = desired number of realizations of the clump population;
         npop<=0 means use all [default]
  angles[2] = polar angles [theta,phi]; if not specified, will be set randomly
  mcut = clumps above mcut (in Msun) will not be included [default 1.0e15]

Returns:
  alllens[n,5] = array of lens parameters [kappa_s, x, y, r_s (arcsec), tau]
                 for all clumps in all realizations
"""

def galacticus2gravlens(galacticus_data, cosmo, npop=0, angles=[], mcut=1.0e15):

    # extract key quantities from galacticus data
    dat1 = galacticus_data[u'Output1']
    ndata = dat1[u'nodeData']
    mass1 = ndata[u'satelliteBoundMass']
    niso_all = ndata[u'nodeIsIsolated']
    rs1 = ndata[u'darkMatterProfileScale']
    rv1 = ndata[u'nodeVirialRadius']
    massbasic1 = ndata[u'basicMass']
    mtcount = dat1[u'mergerTreeCount']
    mtweight = dat1[u'mergerTreeWeight']
    x1 = ndata[u'satellitePositionX']
    y1 = ndata[u'satellitePositionY']
    z1 = ndata[u'satellitePositionZ']
    r1 = np.transpose(np.vstack((x1,y1,z1)))

    # truncation scaling with mass (where does this come from?)
    xm = np.loadtxt('root_trunc.dat')
    a = xm[:,0]
    xm = xm[:,1]

    # initialize array to hold results
    alllens = []

    # loop over merger trees
    if (npop<=0)or(npop>len(mtcount)): npop = len(mtcount)
    for i in range(npop):

        # find starting index for this population
        if i==0:
            ic = 0
        else:
            ic = np.sum(mtcount[:i])

        # number of halos in this population
        nhalo = mtcount[i]

        # positions
        x = x1[ic:ic+nhalo]
        y = y1[ic:ic+nhalo]
        z = z1[ic:ic+nhalo]
        r = np.transpose(np.vstack((x,y,z)))
        # projection
        rproj = project_positions(r,angles=angles)

        # physical quantities: scale radii, mass, etc.
        rs = rs1[ic:ic+nhalo]
        rv = rv1[ic:ic+nhalo]
        cm = rv/rs
        massbasic = massbasic1[ic:ic+nhalo]
        mass = mass1[ic:ic+nhalo]
        rhos = massbasic/(4.0*np.pi*rs**3)/(np.log(1.0+cm)-cm/(1.0+cm))
        M0 = 4.0*np.pi*rs**3*rhos

        # lensing strength
        kappas = rhos*rs/cosmo.Sigma_crit.value  

        # exclude the main halo by using only subhalos that are not isolated
        niso = niso_all[ic:ic+nhalo]
        niso_i = np.where(niso==0)
        r = r[niso_i]
        rproj = rproj[niso_i]
        rs = rs[niso_i]
        mass = mass[niso_i]
        M0 = M0[niso_i]
        kappas = kappas[niso_i]

        # handle truncation
        ai = mass/M0
        tau = np.interp(ai,a,xm,left=-1,right=-1)

        # limit to cases where tau is positive
        tau_pos_i = np.where(tau>0.0)
        r = r[tau_pos_i]
        rproj = rproj[tau_pos_i]
        rs = rs[tau_pos_i]
        mass = mass[tau_pos_i]
        M0 = M0[tau_pos_i]
        kappas = kappas[tau_pos_i]
        tau = tau[tau_pos_i]

        # impose mass cut
        mcut_i = np.where(mass<mcut)
        r = r[mcut_i]
        rproj = rproj[mcut_i]
        rs = rs[mcut_i]
        mass = mass[mcut_i]
        M0 = M0[mcut_i]
        kappas = kappas[mcut_i]
        tau = tau[mcut_i]

        # convert from Mpc to arcsec
        rproj_as = rproj*cosmo.as_Mpc.value
        rs_as = rs*cosmo.as_Mpc.value

        # assemble lensing parameters for this population
        tmplens = np.transpose(np.vstack((kappas,rproj_as[:,0],rproj_as[:,1],rs_as,tau)))
        print('Population',i,'Nclump',len(tmplens))
        alllens.append(tmplens)

    return alllens

############################################################
"""
Write lens parameters formatted as a gravlens startup file.

Inputs:
  outbase = base name for output files
  lensparams[n,5] = array of lens parameters from galacticus2gravlens

Results:
  A set of files to be used with the gravlens "setlens" command.
  The structure of the file names is: outbaseN.start
"""

def write_gravlens(outbase, lensparams):
    for ipop in range(len(lensparams)):
        f = open(outbase+str(ipop)+'.start','w')
        f.write(str(len(lensparams[ipop]))+' 1\n')
        for p in lensparams[ipop]:
            f.write('tnfw3 {} {} {} 0 0 0 0 {} {} 1\n'.format(p[0],p[1],p[2],p[3],p[4]))
        for p in lensparams[ipop]:
            f.write('0 0 0 0 0 0 0 0 0 0\n')
        f.close()

############################################################

