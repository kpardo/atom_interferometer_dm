import numpy as np
import astropy.units as u
import astropy.constants as const

u.set_enabled_equivalencies(u.mass_energy())

lp = (1.*u.GeV*0.197e-15*u.m)**(-1)
tp = (const.c*lp).to((u.GeV*u.s)**(-1))
rhox=(0.4*u.GeV/((0.01*u.m)**3*lp**3)).to(u.GeV**4)
vesc = (600*u.km/u.s*1./const.c).to('')
vdm = (230*u.km/u.s*1./const.c).to('')
