'''
this script gets the cross section limits for given inputs.
inputs: mediator type, phase, N2_only

example:
python get_limits.py 'light' False False
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

# exps = ['GDM', 'MAQRO', 'Pino', 'BECCAL', 'Stanford']
exps = ['MAQRO', 'Pino', 'BECCAL', 'Stanford', 'GDM']
# exps = ['GDM', 'BECCAL', 'Stanford']
# mxs = np.logspace(-6.5, 3.5, 1000)*u.MeV
mxs = np.logspace(-5, 3, 1000)*u.MeV

try:
    med = sys.argv[1]
    phase = eval(sys.argv[2])
    N2_only = eval(sys.argv[3])
    assert isinstance(phase, bool)
    assert isinstance(N2_only, bool)
except:
    print("Oops you didn't give a valid input! Running default.")
    med = 'light'
    phase = False
    N2_only = False

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
    # if (phase) and (med == 'light'):
    if (med == 'light'):
        lims = [[cs_limit(mx, ex=exp, medtype=medtype, mphi=mpr*mx, phase=phase, N2_only=N2_only) for mx in mxs] for mpr in mphi_ratio]
    elif (phase) and (med == 'heavy'):
        if ex in ['BECCAL', 'GDM', 'Stanford']:
            ## for numerical issues
            mxs_ = np.logspace(-5, 3, 5000)*u.MeV
        else:
            mxs_ = np.logspace(-6.5, 3.5, 1000)*u.MeV
        lims = [cs_limit(mx, ex=exp, medtype=medtype, phase=phase, N2_only=N2_only) for mx in mxs_]
    elif (~phase) and (med == 'heavy'):
        if ex in ['BECCAL', 'GDM', 'Stanford']:
            ## for numerical issues
            mxs_ = np.logspace(-5, 3, 5000)*u.MeV
        else:
            mxs_ = np.logspace(-6.5, 3.5, 1000)*u.MeV
        lims = cs_limit(mxs_, ex=exp, medtype=medtype, phase=phase)
    elif (phase) and (med == 'fixed_light'):
        lims = [cs_limit(mx, ex=exp, medtype='light', mphi=mphiratios*u.MeV, phase=phase, N2_only=N2_only) for mx in mxs]
    elif (~phase) and (med == 'fixed_light'):
        lims = cs_limit(mxs, ex=exp, medtype='light', mphi=mphiratios*u.MeV, phase=phase, N2_only=N2_only)
    else:
        try:
            lims = [cs_limit(mxs, ex=exp, medtype=medtype, mphi=mpr*mxs, phase=phase, N2_only=N2_only) for mpr in mphi_ratio]
        except TypeError:
            lims = cs_limit(mxs, ex=exp, medtype=medtype, mphi=mphi_ratio*mxs, phase=phase, N2_only=N2_only)
    ## Save to file.
    if medtype == 'heavy':
        fn = f'../results/limits/{ex}_{medtype[0]}'
    else:
        try:
            len(mphi_ratio)
            fn = f'../results/limits/{ex}_{medtype[0]}'
        except TypeError:
            fn = f'../results/limits/{ex}_{medtype[0]}_{int(-1.*np.log10(mphi_ratio))}'
    if N2_only:
        fn +="_N2only"

    if phase == True:
        fn += "_phase.dat"
    else:
        fn += ".dat"

    tm = Table()
    tm['mx'] = mxs
    if (med == 'heavy') and (ex in ['BECCAL', 'GDM', 'Stanford']):
        tm['mx'] = np.logspace(-5, 3, 5000)*u.MeV
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
