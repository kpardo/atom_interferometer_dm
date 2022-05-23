import numpy as np
import astropy.units as u
import astropy.constants as const
from dataclasses import dataclass
from scipy.integrate import quad_vec, simpson
from scipy.special import spherical_jn, erf, erfc

from aidm.const import tp, lp, rhox, vesc, vdm, ve
from aidm.experiments import GDM
from aidm.dm_velocity import *

import warnings
# warnings.filterwarnings("ignore", category=RuntimeWarning)

u.set_enabled_equivalencies(u.mass_energy())
gdm = GDM() ## to set default experiment


def rate_prefac(mx, ex=gdm, cs=1.*u.cm**2,vdm=vdm, N0=N0(vdm=vdm, vesc=vesc)):
    cs = (cs*lp**2).to(u.GeV**(-2))
    T = (ex.Texp*tp).to(u.GeV**(-1))
    return (T*ex.N*np.pi*cs*rhox*vdm**2/(mx**3*N0)).to(u.MeV**(-2))

def light_rate_prefac(mx, ex=gdm, mphi=None, cs=1.*u.cm**2,vdm=vdm, N0=N0(vdm=vdm, vesc=vesc)):
    if mphi == None:
        mphi = 1.e-5*mx
    cs = (cs*lp**2).to(u.GeV**(-2))
    T = (ex.Texp*tp).to(u.GeV**(-1))
    return (T*ex.N*np.pi*cs*rhox*vdm**2*((mx*vdm)**2+mphi**2)**2 / (mx**3*N0)).to(u.MeV**(2))

def formfac(x):
    return 3.*spherical_jn(1,x)/x

def helmformfac(q, ex=gdm, s=(0.9e-15*u.m*lp).to(u.MeV**(-1))):
    expfactor = np.exp(-0.5*q**2*s.value**2)
    ra = ex.A**(1./3.)*(1.2e-15*u.m*lp).to(u.MeV**(-1))
    return 3.*spherical_jn(1,q*ra.value)/(q*ra.value)*expfactor

def phase_integral(q, mx, ex=gdm, vdm=vdm, ve=ve, vesc=vesc, real=True):
    ## first check vmin:
    # vmin = vesc - 1/q*np.abs(q*ve + q**2/(2*mx))
    # try:
    #     if vmin < 0:
    #         return 0.
    # except ValueError:
    #     if all(vmin) < 0:
    #         return np.zeros((len(vmin)))
    ve = ve.value
    vdm = vdm.value
    vesc = vesc.value
    prefac = -0.25j*np.sqrt(np.pi)*vdm/ve*np.exp(-q**2*ex.deltax.value*(2j*ve+mx*vdm**2*ex.deltax.value)/(4*mx*ve**2))
    term1 = -1 + erf((q/mx-2*ve-1.j*q*vdm**2*ex.deltax.value/ve)/(2*vdm))
    term2mult = np.exp(1j*q**2*ex.deltax.value/(mx*ve))
    term2a = -1.*erf((q/mx-2*ve+1.j*q*vdm**2*ex.deltax.value/ve)/(2*vdm))
    term2b = 1.*erf((q/mx+2*ve+1.j*q*vdm**2*ex.deltax.value/ve)/(2*vdm))
    term2 = term2mult*(term2a + term2b)
    term3 = erfc((q/mx+2*ve-1.j*q*vdm**2*ex.deltax.value/ve)/(2*vdm))
    if real:
        res = np.real(prefac*(term1 + term2 + term3))
    if real == False:
        res = np.imag(prefac*(term1 + term2 + term3))
    # res[vmin < 0] = 0.
    try:
        # print(res)
        checknan = np.isnan(res)
        if checknan:
            return 0.
    except:
        res[np.isnan(res)] = 0.
    return res

def rate_integral(mx, ex=gdm, phase=False, exactphase=False):
    formfacexp = lambda q: (1+ex.A*helmformfac(q, ex=ex)**2+ex.N*(formfac(q*ex.r.value))**2)
    if phase:
        integrand = lambda q: q*phase_integral(q, mx.value, ex=ex)*formfacexp(q)
        if exactphase:
            return quad_vec(integrand, 0., 1./ex.r.value+1.e10)[0]*u.MeV**(-2)
        else:
            ### this is good to ~1% for GDM, but doesn't fit all other exp.
            # qs = np.logspace(np.log10(1./ex.deltax.value)-4.5, np.log10(1./ex.deltax.value)+2, 30000)
            ### this is good to <1% for GDM and fits others. But takes a bit longer.
            qs = np.logspace(np.log10(1./ex.deltax.value)-10, np.log10(1./ex.deltax.value)+10, 100000)
            return simpson(integrand(qs), x=qs)*u.MeV**(-2)
    else:
        integrand = lambda q: q*(1.-np.sin(q*ex.deltax.value)/(q*ex.deltax.value))*formfacexp(q)*expfac(q, mx)
    return quad_vec(integrand, 0., np.inf)[0]*u.MeV**2

def light_rate_integral(mx, ex=gdm, mphi=None, phase=False, exactphase=False):
    if mphi == None:
        mphi = 1.e-5*mx
    formfacexp = lambda q: (1+ex.A*helmformfac(q, ex=ex)**2+ex.N*(formfac(q*ex.r.value))**2)
    if phase:
        integrand = lambda q: q/(q**2 + mphi.value**2)**2*phase_integral(q, mx.value, ex=ex)*formfacexp(q)
        if exactphase:
            return quad_vec(integrand, 0., 1./ex.r.value+1.e10)[0]*u.MeV**(-2)
        else:
            ### this is good to ~1% for GDM, but doesn't fit all other exp.
            # qs = np.logspace(np.log10(1./ex.deltax.value)-4.5, np.log10(1./ex.deltax.value)+2, 30000)
            ### this is good to <1% for GDM and fits others. But takes a bit longer.
            qs = np.logspace(np.log10(1./ex.deltax.value)-10, np.log10(1./ex.deltax.value)+10, 100000)
            return simpson(integrand(qs), x=qs)*u.MeV**(-2)
    else:
        integrand = lambda q: q/(q**2 + mphi.value**2)**2*(1.-np.sin(q*ex.deltax.value)/(q*ex.deltax.value))*formfacexp(q)*expfac(q, mx)
        return quad_vec(integrand, 1.e-40, np.inf)[0]*u.MeV**(-2)


def rate(mx, ex=gdm, medtype = 'light', mphi=None, phase=False):
    if medtype == 'heavy':
        rp = rate_prefac(mx, ex=ex)
        ri = rate_integral(mx, ex=ex, phase=phase)

    elif medtype == 'light':
        rp = light_rate_prefac(mx, ex=ex, mphi=mphi)
        ri = light_rate_integral(mx, ex=ex, mphi=mphi, phase=phase)
    return rp*ri
