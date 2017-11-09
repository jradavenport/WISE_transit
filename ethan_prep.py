from GetData import WISE_LC
import pandas as pd
import numpy as np
import os
from glob import glob

# colnames = ['designation', 'ra', 'dec', 'sigra', 'sigdec', 'sigradec', 'glon', 'glat', 'elon', 'elat', 'wx', 'wy',
#             'w1mpro', 'w1sigmpro', 'w1snr', 'w1rchi2', 'w2mpro', 'w2sigmpro', 'w2snr', 'w2rchi2',
#             'w3mpro', 'w3sigmpro', 'w3snr', 'w3rchi2', 'w4mpro', 'w4sigmpro', 'w4snr', 'w4rchi2',
#             'rchi2', 'nb', 'na', 'w1sat', 'w2sat', 'w3sat', 'w4sat', 'satnum', 'ra_pm', 'dec_pm',
#             'sigra_pm', 'sigdec_pm', 'sigradec_pm', 'pmra', 'sigpmra', 'pmdec', 'sigpmdec', 'w1rchi2_pm',
#             'w2rchi2_pm', 'w3rchi2_pm', 'w4rchi2_pm', 'rchi2_pm', 'pmcode', 'cc_flags', 'ext_flg', 'var_flg',
#             'ph_qual', 'det_bit', 'moon_lev', 'w1nm', 'w1m', 'w2nm', 'w2m', 'w3nm', 'w3m', 'w4nm', 'w4m',
#             'w1cov', 'w2cov', 'w3cov', 'w4cov', 'tmass_key', 'r_2mass', 'pa_2mass', 'n_2mass', 'j_m_2mass',
#             'j_msig_2mass', 'h_m_2mass', 'h_msig_2mass', 'k_m_2mass', 'k_msig_2mass', 'dist', 'angle']


if not os.path.exists('lc'):
    os.makedirs('lc')

data = glob('data/*neowiser_p1bs_psd*')

Nlimit = 9000

for k in range(len(data)):
    df4 = pd.read_csv(data[k])

    ok4 = (df4['ph_qual'].str[0] == 'A') & (df4['nb'] == 1) & (df4['cc_flags'].str[0:2] == '00') & (df4['w1rchi2'] < 5)

    if (sum(ok4) >= Nlimit):
        med_mag = np.nanmedian(df4['w1mpro'][ok4])
        w1flux = 10**((df4['w1mpro'][ok4] - med_mag) / (-2.5))
        w1fluxerr = np.abs(df4['w1sigmpro'][ok4] * np.log(10) / (-2.5) * w1flux)

        # write simple 3-col file for Ethan
        df_tmp = pd.DataFrame(data={'time':df4['mjd'][ok4], 'flux':w1flux, 'err':w1fluxerr})

        df_tmp.to_csv('lc/' + data[k][5:-4] + '.dat', columns=('time','flux','err'), index=False, sep=' ')


print()
print('> package and send to Ethan')