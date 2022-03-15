import numpy as np
import astropy.units as u
import astropy.constants as const
from dataclasses import dataclass
from scipy.integrate import quad_vec
from scipy.special import spherical_jn

from aidm.const import tp, lp, rhox, vesc, vdm
from aidm.experiments import GDM
from aidm.dm_velocity import *

u.set_enabled_equivalencies(u.mass_energy())
gdm = GDM() ## to set default experiment


def rate_prefac(mx, ex=gdm, cs=1.*u.cm**2,vdm=vdm, N0=N0(vdm=vdm, vesc=vesc)):
    cs = (cs*lp**2).to(u.GeV**(-2))
    T = (ex.Texp*tp).to(u.GeV**(-1))
    return (T*ex.N*np.pi*cs*rhox*vdm**2/(mx**3*N0)).to(u.MeV**(-2))

def light_rate_prefac(mx, ex=gdm, cs=1.*u.cm**2,vdm=vdm, N0=2.4582e-9):
    cs = (cs*lp**2).to(u.GeV**(-2))
    T = (ex.Texp*tp).to(u.GeV**(-1))
    return (T*ex.N*np.pi*cs*rhox*vdm**6*mx/N0).to(u.MeV**(2))

def formfac(x):
    return 3.*spherical_jn(1,x)/x

def helmformfac(q, ex=gdm, s=(1.e-15*u.m*lp).to(u.MeV**(-1))):
    expfactor = np.exp(-q**2/s.value**2)
    ra = ex.A**(1./3.)*(1.2e-15*u.m*lp).to(u.MeV**(-1))
    return 3.*spherical_jn(1,q*ra.value)/(q*ra.value)*expfactor


def rate_integral(mx, ex=gdm):
    formfacexp = lambda q: (1+ex.A*helmformfac(q, ex=ex)**2+ex.N*(formfac(q*ex.r.value))**2)
    integrand = lambda q: q*(1-np.sin(q*ex.deltax.value)/(q*ex.deltax.value))*formfacexp(q)*expfac(q, mx)
    return quad_vec(integrand, ex.qmin.value, np.inf)[0]*u.MeV**2

def light_rate_integral(mx, ex=gdm, mphi=None):
    if mphi == None:
        mphi = 1.e-5*mx
    formfacexp = lambda q: (1+ex.A*helmformfac(q, ex=ex)**2+ex.N*(formfac(q*ex.r.value))**2)
    integrand = lambda q: q/(q**2 + mphi.value**2)**2*(1-np.sin(q*ex.deltax.value)/(q*ex.deltax.value))*formfacexp(q)*expfac(q, mx)
    return quad_vec(integrand, ex.qmin.value, np.inf)[0]*u.MeV**(-2)

def rate(mx, ex=gdm, medtype = 'light', mphi=None):
    if medtype == 'heavy':
        rp = rate_prefac(mx, ex=ex)
        ri = rate_integral(mx, ex=ex)

    elif medtype == 'light':
        rp = light_rate_prefac(mx, ex=ex)
        ri = light_rate_integral(mx, ex=ex, mphi=mphi)
    return rp*ri
