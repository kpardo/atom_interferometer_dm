import numpy as np
import astropy.units as u
import astropy.constants as const
from dataclasses import dataclass
from aidm.const import tp, lp

u.set_enabled_equivalencies(u.mass_energy())


@dataclass
class Experiment:
    name: str
    N: float
    A: int
    amin: int
    r: float
    deltax: float
    Nmeas: float
    Texp: float
    phimin: float
    qmin: float = 1.e-20*u.MeV

    def __post_init__(self):
        self.mT = (const.m_p*self.N).to(u.MeV)

    #def get_qmin(self):
    #    return (0.5*self.mT*self.amin*self.Texp).to(u.MeV)

@dataclass
class GDM(Experiment):
    name: str = "GDM"
    N: float = 8.7e9
    A: float = 87
    amin: float = (2.2e-15*u.m/u.s**2*1./const.c).to(u.Hz)
    r: float = (1.e-3*u.m*lp).to(u.MeV**(-1))
    deltax: float = (25*u.m*lp).to(u.MeV**(-1))
    Texp: float = 20.*u.s
    Nmeas: int = int((1.*u.yr/Texp).to(''))
    etadm: float = 1.
    etabkg: float = 1.
    phimin: float = 1.4e-10

@dataclass
class BECCAL(Experiment):
    name: str="BECCAL"
    N: float = 8.7e7
    A: float = 87
    amin: float = (3.e-12*u.m/u.s**2*1./const.c).to(u.Hz)
    r: float = (1.5e-4*u.m*lp).to(u.MeV**(-1))
    deltax: float = (3.e-3*u.m*lp).to(u.MeV**(-1))
    Texp: float = 2.6*u.s
    Nmeas: int = int((1.*u.yr/Texp).to(''))
    etadm: float = 0.5
    etabkg: float = 0.001
    phimin: float = 8.6e-6

@dataclass
class MAQRO(Experiment):
    name: str='MAQRO'
    N: float = 1.e10
    A: float = 60
    amin: float = (1.e-13*u.m/u.s**2*1./const.c).to(u.Hz)
    r: float = (120.e-9*u.m*lp).to(u.MeV**(-1))
    deltax: float = (100.e-9*u.m*lp).to(u.MeV**(-1))
    Texp: float = 100.*u.s
    Nmeas: int = int((1.*u.yr/Texp).to(''))
    etadm: float = 1.0
    etabkg: float = 1.0
    phimin: float = 2.5e-10

@dataclass
class Pino(Experiment):
    name: str = "Pino"
    N: float = 2.2e13
    A: float = 93
    amin: float = (2.5e-12*u.m/u.s**2*1./const.c).to(u.Hz)
    r: float = (1.e-6*u.m*lp).to(u.MeV**(-1))
    deltax: float = (290.e-9*u.m*lp).to(u.MeV**(-1))
    Texp: float = 0.483*u.s
    Nmeas: int = int((1.*u.yr/Texp).to(''))
    etadm: float = 0.5
    etabkg: float = 0.001
    phimin: float = 4.0e-4
