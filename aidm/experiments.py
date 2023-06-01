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
    qmin: float = 1.e-30*u.MeV

    def __post_init__(self):
        self.mT = (const.m_p*self.N).to(u.MeV)
        self.phimin = self.get_phimin()

    def get_phimin(self):
        return (0.25*const.m_p.to(u.MeV)*self.A*self.deltax*self.amin*self.Texp).to('')


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


@dataclass
class BECCAL(Experiment):
    name: str = "BECCAL"
    N: float = 8.7e7
    A: float = 87
    amin: float = (3.e-12*u.m/u.s**2*1./const.c).to(u.Hz)
    r: float = (1.5e-4*u.m*lp).to(u.MeV**(-1))
    deltax: float = (3.e-3*u.m*lp).to(u.MeV**(-1))
    Texp: float = 2.6*u.s
    Nmeas: int = int((1.*u.yr/Texp).to(''))
    etadm: float = 0.5
    etabkg: float = 0.001


@dataclass
class MAQRO(Experiment):
    name: str = 'MAQRO'
    N: float = 1.e10
    A: float = 60
    amin: float = (1.e-13*u.m/u.s**2*1./const.c).to(u.Hz)
    r: float = (120.e-9*u.m*lp).to(u.MeV**(-1))
    deltax: float = (100.e-9*u.m*lp).to(u.MeV**(-1))
    Texp: float = 100.*u.s
    Nmeas: int = int((1.*u.yr/Texp).to(''))
    etadm: float = 1.0
    etabkg: float = 1.0


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


@dataclass
class Stanford(Experiment):
    # drop tower described in Asenbaum + (2020) PRL
    name: str = 'Stanford'
    N: float = 87*4.e6
    A: float = 87
    amin: float = (1.372e-10*u.m/u.s**2*1./const.c).to(u.Hz)
    r: float = (200.e-6*u.m*lp).to(u.MeV**(-1))
    deltax: float = (0.067*u.m*lp).to(u.MeV**(-1))
    Texp: float = 1.910*u.s
    Nmeas: int = int((1.*u.yr/Texp).to(''))


@dataclass
class MAGIS(Experiment):
    # using http://arxiv.org/abs/1812.00482
    # and https://arxiv.org/abs/2104.02835
    # FIX -- this doesn't have enough info!
    pass


@dataclass
class AION10(Experiment):
    # using goal #s from https://arxiv.org/abs/1911.11755 Tab 1
    name: str = 'AION-10'
    N: float = 88*1.e5  # FIX
    A: float = 88.  # use Sr-88 or Sr-87
    amin: float = (1.*u.m/u.s**2*1./const.c).to(u.Hz)  # FIX
    r: float = (1*u.m*lp).to(u.MeV**(-1))  # FIX
    deltax: float = (10.*u.m*lp).to(u.MeV**(-1))  # FIX
    Texp: float = 1.4*u.s
    Nmeas: int = int((1.*u.yr/Texp).to(''))


@dataclass
class AION100(Experiment):
    # using goal #s from https://arxiv.org/abs/1911.11755 Tab 1
    name: str = 'AION-100'
    N: float = 88*1.e7  # FIX - assume 100x more than 10.
    A: float = 88.  # use Sr-88 or Sr-87
    amin: float = (1.*u.m/u.s**2*1./const.c).to(u.Hz)  # FIX
    r: float = (1*u.m*lp).to(u.MeV**(-1))  # FIX
    deltax: float = (10.*u.m*lp).to(u.MeV**(-1))  # FIX
    Texp: float = 1.4*u.s
    Nmeas: int = int((1.*u.yr/Texp).to(''))


@dataclass
class AIONkm(Experiment):
    # from https://arxiv.org/abs/1911.11755 Tab 1
    name: str = 'AION-km'
    N: float = 88*1.e8  # FIX -- assume 10x more than 100.
    A: float = 88.  # use Sr-88 or Sr-87
    amin: float = (1.*u.m/u.s**2*1./const.c).to(u.Hz)  # FIX
    r: float = (1*u.m*lp).to(u.MeV**(-1))  # FIX
    deltax: float = (10.*u.m*lp).to(u.MeV**(-1))  # FIX
    Texp: float = 5.*u.s
    Nmeas: int = int((1.*u.yr/Texp).to(''))


@dataclass
class AEDGE(Experiment):
    # from https://arxiv.org/abs/1908.00802 Tab 1
    name: str = 'AEDGE'
    N: float = 87*1.e11  # FIX -- assume 10x more than 1 km.
    A: float = 87.  # use Sr-87 (696 nm clock transition)
    amin: float = (1.*u.m/u.s**2*1./const.c).to(u.Hz)  # FIX
    r: float = (1*u.m*lp).to(u.MeV**(-1))  # FIX
    deltax: float = (0.9*u.m*lp).to(u.MeV**(-1))  # FIX
    Texp: float = 600.*u.s
    Nmeas: int = int((1.*u.yr/Texp).to(''))
