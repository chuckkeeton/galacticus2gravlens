This repository contains code to take semi-analytic models of galaxies from the [Galacticus code](https://bitbucket.org/galacticusdev/galacticus/wiki/Home) and produce input files for gravitational lensing calculations with the [gravlens package](http://www.physics.rutgers.edu/~keeton/gravlens/2012WS).

Developed by Sean Brennan and Chuck Keeton, drawing on earlier code by Andrew Benson and Anthony Pullen.

Last updated 2019/3/27.

### Code

The main python code is [galacticus2gravlens.py](galacticus2gravlens.py). It uses data in [root_trunc.dat](root_trunc.dat).

#### Example 1

Here is an example that uses the python code to process the first 100 realizations using the default projection ([example1.py](example1.py)).

```
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

# compute lensing parameters for default projection
lens0 = g2g.galacticus2gravlens(galdat,cosmo,npop=100,angles=[0,0])
# write parameters for gravlens
g2g.write_gravlens('example1/pop',lens0)
```

Here is a bare bones gravlens input file that uses a few of the realizations to generate convergence maps ([example1.in](example1.in)).

```
set gridflag = 0

setlens example1/pop0.start
plotkappa example1/pop0.fits 3 -60.0 60.0 200 -60.0 60.0 200

setlens example1/pop1.start
plotkappa example1/pop1.fits 3 -60.0 60.0 200 -60.0 60.0 200

setlens example1/pop2.start
plotkappa example1/pop2.fits 3 -60.0 60.0 200 -60.0 60.0 200

quit
```

#### Example 2

Another example uses the python code to process just the first two realizations but with 10 random projections for each: [example2.py](example2.py) and [example2.in](example2.in)
