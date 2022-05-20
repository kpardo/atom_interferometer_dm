import numpy as np
import astropy.units as u
import astropy.constants as const
from dataclasses import dataclass
from scipy.integrate import quad_vec
from scipy.special import spherical_jn

from aidm.const import tp, lp, rhox, vesc, vdm
from aidm.experiments import GDM
from aidm.rate import *

u.set_enabled_equivalencies(u.mass_energy())
gdm = GDM() ## to set default experiment


def cs_limit(mx, ex=gdm, medtype='light',mphi=None, phase=False):
    ave={'v':0.5, 'phi':0.}
    std={'v':0.1*0.5, 'phi':1.*ex.phimin}
    if phase:
        return (std['phi'])*1./rate(mx,ex=ex, medtype=medtype, mphi=mphi, phase=phase)*u.cm**2
    else:
        max_V = ave['v'] - std['v'] ## e.g., 1 sigma lower than average.
        stds = std['v']/ave['v']
        # smin = -1.*np.log(max_V)
        return ((stds+(-1.*np.log(ave['v'])))/ex.Nmeas**(0.5))*1./rate(mx,ex=ex, medtype=medtype, mphi=mphi, phase=phase)*u.cm**2

def noise_mod(gammavis = 0.5, ex=gdm):
    return (4*(gammavis**(-1)-1)/ex.Nmeas)**(1./2.)

def bkg_mod(gammavis=0.5, ex=gdm):
    return ex.etabkg*np.log(gammavis**(-1))

def cs_limit_mod(mx, ex=gdm, medtype='light', mphi=None, phase=False):
    return ((noise_mod(ex=ex)+bkg_mod(ex=ex))/(ex.etadm*rate(mx, ex=ex, medtype=medtype, mphi=mphi, phase=phase)))*u.cm**2

def cs_to_axion(cs, deltamn=1.*u.GeV):
    cs = (cs*lp**2)
    return ((deltamn**2/(256*np.pi*cs))**(1./4.)).to(u.GeV)

def cs_to_gb(cs, ex=gdm):
    mN=ex.mT/ex.N
    cs = (cs*lp**2)
    return ((16*np.pi*cs/3*mN**2)**(1./4.)*3.).to('')
