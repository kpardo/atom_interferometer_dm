import numpy as np
import astropy.units as u
import astropy.constants as const

u.set_enabled_equivalencies(u.mass_energy())

## conversion factors
lp = (1.*u.GeV*0.197e-15*u.m)**(-1) ## convert m to GeV
tp = (const.c*lp).to((u.GeV*u.s)**(-1)) ## convert s to GeV

## DM constants
Omx = 0.05 ## fraction of dark matter in X.
rhox= Omx*(0.4*u.GeV/((0.01*u.m)**3*lp**3)).to(u.GeV**4) ## local DM density
vesc = (600*u.km/u.s*1./const.c).to('') ## esc vel of MW
vdm = (230*u.km/u.s*1./const.c).to('') ## local DM vel.
