'''
this script gets the cross section limits for given inputs.
inputs: mediator type, phase

example:
python get_limits.py 'light' True
'''

import numpy as np
import astropy.units as u
import astropy.constants as const
from dataclasses import dataclass
from scipy.integrate import quad_vec
from scipy.special import spherical_jn
from astropy.table import Table, hstack
import datetime
import sys

from aidm.const import tp, lp, rhox, vesc, vdm
import aidm.experiments as x
from aidm.cross_sections import cs_limit, cs_limit_mod

exps = ['GDM', 'MAQRO', 'Pino', 'BECCAL']
mxs = np.logspace(-6.5, 3.5, 1000)*u.MeV

try:
    med = sys.argv[1]
    phase = bool(sys.argv[2])
except:
    print("Oops you didn't give an input! Running default.")
    med = 'light'
    phase = False

if med == 'light':
    mphiratios = [1.e-10, 1.e-9, 1.e-7, 1.e-6, 1.e-5, 1.e-4, 1.e-3, 1.e-2]
elif med == 'fixed_light':
    mphiratios = (1.*u.eV).to(u.MeV).value
else:
    mphiratios = 1.

def get_lim(ex, mphi_ratio = mphiratios, medtype='light',
    phase = phase):

    if phase:
        print(f'{datetime.datetime.now()}: Getting {medtype} mediator phase limits for {ex}')
    else:
        print(f'{datetime.datetime.now()}: Getting {medtype} mediator decoherence limits for {ex}')
    ## Define exp. w/ qmin
    class_ = getattr(x, ex)
    exp = class_()

    ## Get limits
    if (phase) and (med == 'light'):
        lims = [[cs_limit(mx, ex=exp, medtype=medtype, mphi=mpr*mx, phase=phase) for mx in mxs] for mpr in mphi_ratio]
    elif (phase) and (med == 'heavy'):
        lims = [cs_limit(mx, ex=exp, medtype=medtype, phase=phase) for mx in mxs]
    elif (phase) and (med == 'fixed_light'):
        lims = [cs_limit(mx, ex=exp, medtype='light', mphi=mphiratios*u.MeV, phase=phase) for mx in mxs]
    elif (~phase) and (med == 'fixed_light'):
        lims = cs_limit(mxs, ex=exp, medtype='light', mphi=mphiratios*u.MeV, phase=phase)
    else:
        try:
            lims = [cs_limit(mxs, ex=exp, medtype=medtype, mphi=mpr*mxs, phase=phase) for mpr in mphi_ratio]
        except TypeError:
            lims = cs_limit(mxs, ex=exp, medtype=medtype, mphi=mphi_ratio*mxs, phase=phase)
    ## Save to file.
    if medtype == 'heavy':
        fn = f'../results/limits/{ex}_{medtype[0]}'
    else:
        try:
            len(mphi_ratio)
            fn = f'../results/limits/{ex}_{medtype[0]}'
        except TypeError:
            fn = f'../results/limits/{ex}_{medtype[0]}_{int(-1.*np.log10(mphi_ratio))}'

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
    print(f'{datetime.datetime.now()}: Saved to file {fn}')
    try:
        print(f'Best cross section is {np.min(lims)} at mx = {mxs[np.argmin(lims)]}.')
    except:
        return 0
    return 0

## Run
for ex in exps:
    get_lim(ex, medtype=med, phase=phase)
