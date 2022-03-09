import numpy as np
import astropy.units as u
import astropy.constants as const
from dataclasses import dataclass
from scipy.integrate import quad_vec
from scipy.special import spherical_jn

from atom_interferometer_dm.const import tp, lp, rhox, vesc, vdm
from atom_interferometer_dm.experiments import GDM
from atom_interferometer_dm.rate import *

u.set_enabled_equivalencies(u.mass_energy())
gdm = GDM() ## to set default experiment


def noise(gammavis = 0.5, ex=gdm):
    return (4*(gammavis**(-1)-1)/ex.Nmeas)**(1./2.)

def bkg(gammavis=0.5):
    return np.log(gammavis**(-1))

def cs_limit(mx, ex=gdm, medtype='light', mphi=None):
    return ((noise(ex=ex)+bkg())/rate(mx, ex=ex, medtype=medtype, mphi=mphi))*u.cm**2

def cs_to_axion(cs, deltamn=1.*u.GeV):
    cs = (cs*lp**2)
    return ((deltamn**2/(256*np.pi*cs))**(1./4.)).to(u.GeV)

def cs_to_gb(cs, ex=gdm):
    mN=ex.mT/ex.N
    cs = (cs*lp**2)
    return ((16*np.pi*cs/3*mN**2)**(1./4.)*3.).to('')
