import numpy as np
import astropy.units as u
import astropy.constants as const
from dataclasses import dataclass

from atom_interferometer_dm.const import tp, lp, rhox, vesc, vdm

u.set_enabled_equivalencies(u.mass_energy())

def vmin(q, mx, vesc=vesc):
    vq = q/(2.*mx)
    vq[vq>vesc] = vesc
    return vq

def expfac(q, mx, vdm=vdm, vesc=vesc):
    return np.exp(-vmin(q, mx.value)**2/vdm**2)-np.exp(-vesc**2/vdm**2)
