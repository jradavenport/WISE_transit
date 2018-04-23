import pandas as pd
import numpy as np
import os
from glob import glob


def ethan_prep(Nlimit = 9000, tfix=45):
    '''
    Nlimit=9000, need at least this many visits to run
    tfix=45min, the minimum time separation to have between epochs, a hack fix
    '''
    if not os.path.exists('lc'):
        os.makedirs('lc')

    # Only grab the NeoWISE data for now, since we want to do the search on this
    # best-sampled data, and then verify with the original WISE data.
    data = glob('data/*neowiser_p1bs_psd*')

    for k in range(len(data)):
        print(data[k])
        # only running the NeoWISE stuff Will use original WISE as verification for any candidate
        df4 = pd.read_csv(data[k])
        ok4 = np.where((df4['ph_qual'].str[0] == 'A') & (df4['nb'] == 1) &
                       (df4['cc_flags'].astype('str').str[0] == '0') & (df4['w1rchi2'] < 5)  & (df4['qual_frame'] > 8))[0]

        if (sum(ok4) >= Nlimit):
            med_mag = np.nanmedian(df4['w1mpro'][ok4])
            w1flux = 10**((df4['w1mpro'].values[ok4] - med_mag) / (-2.5))
            w1fluxerr = np.abs(df4['w1sigmpro'].values[ok4] * np.log(10) / (-2.5) * w1flux)

            # fix the multiple exposures within a few seconds problem...
            # this is *probably* due to a cosmic ray being cleaned out, and a second epoch being generated
            # it SHOULD be correctable via flag cuts, BUT astroquery won't allow you to use all the flags...
            # SO, executive decision for now: just pick the *latter epoch* when this occurs
            dtime = df4['mjd'].values[ok4][:-1] - df4['mjd'].values[ok4][1:]
            tok = np.where(( dtime*24.*60. < (-1*tfix)  ))

            # write simple 3-col file for Ethan
            df_tmp = pd.DataFrame(data={'time':df4['mjd'].values[ok4][tok], 'flux':w1flux[tok], 'err':w1fluxerr[tok]})

            df_tmp.to_csv('lc/' + data[k][5:-4] + '.dat', columns=('time','flux','err'), index=False, sep=' ')


if __name__ == "__main__":
    # import sys
    ethan_prep()
