"""
In this example, we process just the first two realizations
but examine 10 random projections for each.
"""

import numpy as np
import h5py
import galacticus2gravlens as g2g

# set cosmology
cosmo = g2g.cosmology(0.5,1.0)
print('cosmology:')
print('Dl,Ds,Dls=',cosmo.D_l,cosmo.D_s,cosmo.D_ls)
print('Sigma_crit= {:e}'.format(cosmo.Sigma_crit))

# read galacticus file
galfile = '../pspec/galacticus_mw19200_fsh.hdf5'
galtmp = h5py.File(galfile,'r')
galdat = galtmp[u'Outputs']

# compute lensing parameters for a single population with
# multiple projections
for i in range(10):
    lens = g2g.galacticus2gravlens(galdat,cosmo,npop=2)
    # write parameters for gravlens
    g2g.write_gravlens('example2/proj'+str(i)+'pop',lens)

