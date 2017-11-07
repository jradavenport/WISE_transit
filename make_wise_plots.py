import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astroquery.irsa import Irsa # had to install astroquery w/ pip
# from astroquery.simbad import Simbad

matplotlib.rcParams.update({'font.size':18})
matplotlib.rcParams.update({'font.family':'serif'})

Irsa.ROW_LIMIT = 0

# a list of targets that astroquery can resolve
# need list of WISE targets within 47 arcmin of the Ecliptic Poles
# Simbad defines the NEP = 18 00 00.000 +66 33 38.55
#   SEP = 06 00 00.000 -66 33 38.55
targets = ['WISE J060224.34-661926.1', 'WISE J180214.15+661150.5', 'WISE J180022.32+663315.1',
           'WISE J055625.16-662924.6', 'WISE J175946.15+663746.7', 'WISE J060051.07-664439.4',
           'WISE J060113.59-662905.1', 'WISE J055552.81-662920.9', 'WISE J060105.40-663425.1',
           'WISE J060132.35-660815.0', 'WISE J180201.72+663739.0']

# the WISE tables to search
cats = ['neowiser_p1bs_psd', 'allsky_4band_p1bs_psd', 'allsky_3band_p1bs_psd', 'allsky_2band_p1bs_psd']

for obj in targets:
    print('Running '+ str(obj))

    table1 = Irsa.query_region(obj, catalog=cats[0], spatial='Cone', width='2 arcsec')
    table2 = Irsa.query_region(obj, catalog=cats[1], spatial='Cone', width='2 arcsec')
    table3 = Irsa.query_region(obj, catalog=cats[2], spatial='Cone', width='2 arcsec')
    table4 = Irsa.query_region(obj, catalog=cats[3], spatial='Cone', width='2 arcsec')

    table1.sort('mjd')
    table2.sort('mjd')
    table3.sort('mjd')
    table4.sort('mjd')

    df1 = table1.to_pandas()
    df2 = table2.to_pandas()
    df3 = table3.to_pandas()
    df4 = table4.to_pandas()

    print(' found '+str(len(df1))+' + '+ str(len(df2))+' + '+ str(len(df3))+' + '+ str(len(df4))+' visits')

    # dump to CSV files for possible later analysis
    df1.to_csv('data/' + obj + cats[0] + '.csv')
    df2.to_csv('data/' + obj + cats[1] + '.csv')
    df3.to_csv('data/' + obj + cats[2] + '.csv')
    df4.to_csv('data/' + obj + cats[3] + '.csv')


    #### NEED TO ADD QUALITY CUTS
    # can't add this to the latter 3 surveys...  (df1['qual_frame'] > 8)
    ok1 = (df1['ph_qual'].str[0] == 'A') & (df1['nb'] == 1) & (df1['cc_flags'].str[0:2] == '00') & (df1['qual_frame'] > 8) & (df1['w1rchi2'] < 5)
    ok2 = (df2['ph_qual'].str[0] == 'A') & (df2['nb'] == 1) & (df2['cc_flags'].str[0:2] == '00') & (df2['w1rchi2'] < 5)
    ok3 = (df3['ph_qual'].str[0] == 'A') & (df3['nb'] == 1) & (df3['cc_flags'].str[0:2] == '00') & (df3['w1rchi2'] < 5)
    ok4 = (df4['ph_qual'].str[0] == 'A') & (df4['nb'] == 1) & (df4['cc_flags'].str[0:2] == '00') & (df4['w1rchi2'] < 5)


    # colors
    colors = ['#1f77b4', '#ff7f0e', '#c5b0d5', '#d62728']

    ## make 3 basic figures:
    # 1) W1 light curve
    plt.figure(figsize=(13,8))
    plt.errorbar(df1['mjd'][ok1], df1['w1mpro'][ok1], yerr=df1['w1sigmpro'][ok1],
                 marker='o', linestyle='none', alpha=0.25, color=colors[0])
    plt.errorbar(df2['mjd'][ok2], df2['w1mpro'][ok2], yerr=df2['w1sigmpro'][ok2],
                 marker='o', linestyle='none', alpha=0.25, color=colors[1])
    plt.errorbar(df3['mjd'][ok3], df3['w1mpro'][ok3], yerr=df3['w1sigmpro'][ok3],
                 marker='o', linestyle='none', alpha=0.25, color=colors[2])
    plt.errorbar(df4['mjd'][ok4], df4['w1mpro'][ok4], yerr=df4['w1sigmpro'][ok4],
                 marker='o', linestyle='none', alpha=0.25, color=colors[3])
    plt.ylabel('W1 (mag)')
    plt.gca().invert_yaxis()
    plt.title(obj)
    plt.savefig('img/'+obj + '_W1.png', dpi=150, bbox_inches='tight', pad_inches=0.25)
    # plt.show()
    plt.close()



    # 2) W1-W2 color light curve
    plt.figure(figsize=(13,8))
    plt.errorbar(df1['mjd'][ok1], df1['w1mpro'][ok1] - df1['w2mpro'][ok1],
                 yerr=np.sqrt(df1['w1sigmpro'][ok1]**2 + df1['w2sigmpro'][ok1]**2),
                 marker='o', linestyle='none', alpha=0.25, color=colors[0])
    plt.errorbar(df2['mjd'][ok2], df2['w1mpro'][ok2] - df2['w2mpro'][ok2],
                 yerr=np.sqrt(df2['w1sigmpro'][ok2]**2 + df2['w2sigmpro'][ok2]**2),
                 marker='o', linestyle='none', alpha=0.25, color=colors[1])
    plt.errorbar(df3['mjd'][ok3], df3['w1mpro'][ok3] - df3['w2mpro'][ok3],
                 yerr=np.sqrt(df3['w1sigmpro'][ok3]**2 + df3['w2sigmpro'][ok3]**2),
                 marker='o', linestyle='none', alpha=0.25, color=colors[2])
    plt.errorbar(df4['mjd'][ok4], df4['w1mpro'][ok4] - df4['w2mpro'][ok4],
                 yerr=np.sqrt(df4['w1sigmpro'][ok4]**2 + df4['w2sigmpro'][ok4]**2),
                 marker='o', linestyle='none', alpha=0.25, color=colors[3])
    plt.xlabel('MJD (days)')
    plt.ylabel('W1-W2 (mag)')
    plt.title(obj)
    plt.savefig('img/'+obj + '_W1W2.png', dpi=150, bbox_inches='tight', pad_inches=0.25)
    # plt.show()
    plt.close()



    # 3) CMD
    # plt.figure(figsize=(8,8))
    # plt.errorbar(df1['w1mpro'] - df1['w2mpro'], df1['w1mpro'],
    #              xerr=np.sqrt(df1['w1sigmpro']**2 + df1['w2sigmpro']**2), yerr=df1['w1sigmpro'],
    #              marker='o', linestyle='none', alpha=0.25, color='#1f77b4')
    # plt.errorbar(df2['w1mpro_ep'] - df2['w2mpro_ep'], df2['w1mpro_ep'],
    #              xerr=np.sqrt(df2['w1sigmpro_ep']**2 + df2['w2sigmpro_ep']**2 ), yerr=df2['w1sigmpro_ep'],
    #              marker='o', linestyle='none', alpha=0.25, color='#ff7f0e')
    #
    # plt.ylabel('W1 (mag)')
    # plt.xlabel('W1-W2 (mag)')
    # plt.gca().invert_yaxis()
    # plt.savefig('img/'+obj + '_cmd.png', dpi=150, bbox_inches='tight', pad_inches=0.25)
    # plt.close()



    # bonus: RA,Dec to make sure not a blend, etc
    # plt.figure(figsize=(8, 8))
    # plt.scatter(df1['ra'], df1['dec'],
    #              marker='o', alpha=0.25, color='#1f77b4')
    # plt.scatter(df2['ra'], df2['dec'],
    #              marker='o', alpha=0.25, color='#ff7f0e')
    # plt.xlabel('RA (deg)')
    # plt.ylabel('Dec (deg)')
    # plt.savefig('img/' + obj + '_radec.png', dpi=150, bbox_inches='tight', pad_inches=0.25)
    # plt.close()


    '''

    ## now make zoom-in light curves for every independent visit
    # find breaks in data (should be ~6mo unless near the poles)

    dlim = 20 # days, limit of time for a window to be split at

    dt1 = df1['mjd'].values[1:] - df1['mjd'].values[0:-1]
    #bk1 = np.append(np.append([-1], np.where((dt1 > dlim))[0]), len(df1))
    bk1 = np.where((dt1 > dlim))[0]
    R1 = np.append(bk1+1, [len(df1)+1])
    L1 = np.append([0], bk1+1)

    dt2 = df2['mjd'].values[1:] - df2['mjd'].values[0:-1]
    #bk2 = np.append(np.append([-1], np.where((dt2 > dlim))[0]), len(df2))
    bk2 = np.where((dt2 > dlim))[0]
    R2 = np.append(bk2+1, [len(df2)+1])
    L2 = np.append([0], bk2+1)


    # now step through the chunks of continuous data, make a light curve for each
    if len(df2)>0:
        for k in range(len(R2)):
            tmin = np.nanmin(df2['mjd'].values[L2[k]:R2[k]])
            tmin_s = str(tmin)
            plt.figure(figsize=(9,5))
            plt.errorbar(df2['mjd'].values[L2[k]:R2[k]] - tmin, df2['w1mpro_ep'].values[L2[k]:R2[k]],
                         yerr=df2['w1sigmpro_ep'].values[L2[k]:R2[k]],
                         marker='o', linestyle='none', alpha=0.25, color='#ff7f0e')
            plt.xlabel('MJD - '+tmin_s+' (days)')
            plt.ylabel('W1 (mag)')
            plt.gca().invert_yaxis()
            plt.savefig('img/'+obj + '_visit'+str(k+1)+'.png', dpi=150, bbox_inches='tight', pad_inches=0.25)
            plt.close()

    k0 = len(R2) # counter to keep using for visit number

    if len(df1)>0:
        for k in range(len(R1)):
            tmin = np.nanmin(df1['mjd'].values[L1[k]:R1[k]])
            tmin_s = str(tmin)
            plt.figure(figsize=(9,5))
            plt.errorbar(df1['mjd'].values[L1[k]:R1[k]] - tmin, df1['w1mpro'].values[L1[k]:R1[k]],
                         yerr=df1['w1sigmpro'].values[L1[k]:R1[k]],
                         marker='o', linestyle='none', alpha=0.25, color='#1f77b4')
            plt.xlabel('MJD - '+tmin_s+' (days)')
            plt.ylabel('W1 (mag)')
            plt.gca().invert_yaxis()
            plt.savefig('img/'+obj + '_visit'+str(k+k0+1)+'.png', dpi=150, bbox_inches='tight', pad_inches=0.25)
            plt.close()

    '''
