import numpy as np
import astropy.units as u
import astropy.constants as const
from dataclasses import dataclass
from aidm.const import tp, lp

u.set_enabled_equivalencies(u.mass_energy())


@dataclass
class Background:
    name: str
    energies: np.ndarray
    flux: np.ndarray


@dataclass
class Photons:
    name: 'photons'


@dataclass
class CosmicRays:
    name: 'crs'


@dataclass
class Neutrinos:
    name: 'nus'


@dataclass
class Dust:
    name: 'dust'
