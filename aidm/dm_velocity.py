import numpy as np
import astropy.units as u
import astropy.constants as const
from dataclasses import dataclass
from scipy.special import erf

from aidm.const import tp, lp, rhox, vesc, vdm

u.set_enabled_equivalencies(u.mass_energy())

def vmin(q, mx, vesc=vesc, theta=0.):
    vq = q/(2.*mx*np.cos(theta))
    try:
        vq[vq>vesc] = vesc
    except TypeError:
        if vq > vesc:
            vq = vesc
    return vq

def expfac(q, mx, vdm=vdm, vesc=vesc):
    return np.exp(-vmin(q, mx.value)**2/vdm**2)-np.exp(-vesc**2/vdm**2)

def N0(vdm=vdm, vesc=vesc):
    prefac = np.pi**(1.5)*vdm**3
    term1 = erf(vesc/vdm)
    term2 = -2./np.sqrt(np.pi)*vesc/vdm*np.exp(-vesc**2/vdm**2)
    return prefac*(term1+term2)
