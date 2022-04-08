import numpy as np
import astropy.units as u
import astropy.constants as const
from dataclasses import dataclass
from scipy.integrate import quad_vec
from scipy.special import spherical_jn
from astropy.table import Table, hstack
import datetime

from aidm.const import tp, lp, rhox, vesc, vdm
import aidm.experiments as x
from aidm.cross_sections import cs_limit

exps = ['MAQRO', 'GDM', 'Pino', 'BECCAL']
mxs = np.logspace(-3.5, 3.5, 1000)*u.MeV
mphiratios = [1.e-10, 1.e-7, 1.e-5, 1.e-4, 1.e-3, 1.e-2]
qmin = False
med = 'heavy'
phase = False

def get_lim(ex, qmin=True, mphi_ratio = mphiratios, medtype='light',
    phase = phase):
    if phase:
        print(f'{datetime.datetime.now()}: Getting phase limits for {ex} with qmin={qmin}')
    else:
        print(f'{datetime.datetime.now()}: Getting decoherence limits for {ex} with qmin={qmin}')
    ## Define exp. w/ qmin
    class_ = getattr(x, ex)
    exp = class_()
    if qmin == False:
        exp.qmin = 1.e-30*u.MeV

    ## Get limits
    try:
        lims = [cs_limit(mxs, ex=exp, medtype=medtype, mphi=mpr*mxs, phase=phase) for mpr in mphi_ratio]
    except TypeError:
        lims = cs_limit(mxs, ex=exp, medtype=medtype, mphi=mphi_ratio*mxs, phase=phase)

    ## Save to file.
    if medtype == 'heavy':
        fn = f'../results/limits/{ex}_{qmin}_{medtype[0]}'
    else:
        try:
            len(mphi_ratio)
            fn = f'../results/limits/{ex}_{qmin}_{medtype[0]}'
        except TypeError:
            fn = f'../results/limits/{ex}_{qmin}_{medtype[0]}_{int(-1.*np.log10(mphi_ratio))}'

    if phase == True:
        fn += "_phase.dat"
    else:
        fn += ".dat"

    tm = Table()
    tm['mx'] = mxs
    try:
        limT = Table([l for l in lims])
    except TypeError:
        limT = Table([lims])
    t = hstack([tm, limT])
    if medtype == 'heavy':
        names = [f'sigma']
        old = [f'col0']
    else:
        try:
            names = [f'sigma_{int(-1.*np.log10(mpr))}' for mpr in mphi_ratio]
            old = [f'col{i}' for i in range(len(mphi_ratio))]
        except TypeError:
            names = [f'sigma_{int(-1.*np.log10(mphi_ratio))}']
            old = ['col0']
    t.rename_columns(names=old, new_names=names)
    t.write(fn, format='ascii.ecsv', overwrite=True)
    print(f'Best cross section is {np.min(lims)} at mx = {mxs[np.argmin(lims)]}.')
    print(f'{datetime.datetime.now()}: Saved to file {fn}')
    return 0

## Run
for ex in exps:
    get_lim(ex, medtype=med, qmin=qmin)
    # if get_noqmin:
    #     get_lim(ex, qmin=False, medtype=med)
