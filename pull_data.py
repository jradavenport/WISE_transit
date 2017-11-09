from GetData import WISE_LC
import pandas as pd
import numpy as np
import time

colnames = ['designation', 'ra', 'dec', 'sigra', 'sigdec', 'sigradec', 'glon', 'glat', 'elon', 'elat', 'wx', 'wy',
            'w1mpro', 'w1sigmpro', 'w1snr', 'w1rchi2', 'w2mpro', 'w2sigmpro', 'w2snr', 'w2rchi2',
            'w3mpro', 'w3sigmpro', 'w3snr', 'w3rchi2', 'w4mpro', 'w4sigmpro', 'w4snr', 'w4rchi2',
            'rchi2', 'nb', 'na', 'w1sat', 'w2sat', 'w3sat', 'w4sat', 'satnum', 'ra_pm', 'dec_pm',
            'sigra_pm', 'sigdec_pm', 'sigradec_pm', 'pmra', 'sigpmra', 'pmdec', 'sigpmdec', 'w1rchi2_pm',
            'w2rchi2_pm', 'w3rchi2_pm', 'w4rchi2_pm', 'rchi2_pm', 'pmcode', 'cc_flags', 'ext_flg', 'var_flg',
            'ph_qual', 'det_bit', 'moon_lev', 'w1nm', 'w1m', 'w2nm', 'w2m', 'w3nm', 'w3m', 'w4nm', 'w4m',
            'w1cov', 'w2cov', 'w3cov', 'w4cov', 'tmass_key', 'r_2mass', 'pa_2mass', 'n_2mass', 'j_m_2mass',
            'j_msig_2mass', 'h_m_2mass', 'h_msig_2mass', 'k_m_2mass', 'k_msig_2mass', 'dist', 'angle']

NEP = pd.read_table('NEP_sources.txt.gz', delim_whitespace=True, skiprows=180, names=colnames, na_values='null')
SEP = pd.read_table('SEP_sources.txt.gz', delim_whitespace=True, skiprows=180, names=colnames, na_values='null')

#bright, likely main sequnence objects, with continuous coverage
okN = ((NEP['w1mpro'] < 14.5) & (NEP['w1mpro'] > 7.8) & (NEP['dist'] <= 23.5*60) &
       (NEP['j_m_2mass'] - NEP['w1mpro'] < 1.2) & (NEP['j_m_2mass'] - NEP['w1mpro'] > 0.75) &
       (NEP['w1mpro'] - NEP['w2mpro'] > -0.03))

okS = ((SEP['w1mpro'] < 14.5) & (SEP['w1mpro'] > 7.8) & (SEP['dist'] <= 23.5*60) &
       (SEP['j_m_2mass'] - SEP['w1mpro'] < 1.2) & (SEP['j_m_2mass'] - SEP['w1mpro'] > 0.75) &
       (SEP['w1mpro'] - SEP['w2mpro'] > -0.03))

# print(sum(okN), sum(okS))

# PULL NORTH DATA
print('> getting data for ' + str(sum(okN)) + ' NEP sources')
for k in range(sum(okN)):
    print(k, 'WISE ' + NEP['designation'][okN].values[k])
    WISE_LC('WISE ' + NEP['designation'][okN].values[k])
    time.sleep(2) # put a 2second sleep in, to try and not anger the sys-admins

print()

# PULL SOUTH DATA
print('> getting data for ' + str(sum(okS)) + ' SEP sources')
for k in range(sum(okS)):
    print(k,' WISE ' + SEP['designation'][okS].values[k])
    WISE_LC('WISE ' + SEP['designation'][okS].values[k])
    time.sleep(2) # put a 2second sleep in, to try and not anger the sys-admins

print()
print('> happy hunting')
