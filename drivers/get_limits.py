import numpy as np
import astropy.units as u
import astropy.constants as const
from dataclasses import dataclass
from scipy.integrate import quad_vec
from scipy.special import spherical_jn
from astropy.table import Table
import datetime

from atom_interferometer_dm.const import tp, lp, rhox, vesc, vdm
import atom_interferometer_dm.experiments as x
from atom_interferometer_dm.cross_sections import cs_limit

exps = ['GDM', 'MAQRO', 'Pino', 'BECCAL']
mxs = np.logspace(-10, 4, 100)*u.MeV
mphiratios = 1.e-5


def get_lim(ex, qmin=True, mphi_ratio = mphiratios, medtype='light'):
    print(datetime.datetime.now())
    print(f'Getting limits for {ex} with qmin={qmin}')
    ## Define exp. w/ qmin
    class_ = getattr(x, ex)
    exp = class_()
    if qmin == False:
        exp.qmin = 1.e-20*u.MeV

    ## Get limits
    try:
        lims = [cs_limit(mxs, medtype=medtype, mphi=mpr*mxs) for mpr in mphi_ratio]
    except TypeError:
        lims = cs_limit(mxs, medtype=medtype, mphi=mphi_ratio*mxs)

    ## Save to file.
    if medtype=='light':
        fn = f'../results/limits/{ex}_{qmin}_{medtype[0]}_{int(-1.*np.log10(mphi_ratio))}.dat'
    else:
        fn = f'../results/limits/{ex}_{qmin}_{medtype[0]}.dat'
    t = Table([mxs, lims],
        names=('mx', 'sigma'))
    t.write(fn, format='ascii.csv')
    print(f'Saved to file {fn}')
    return 0

## Run
for ex in exps:
    get_lim(ex)
