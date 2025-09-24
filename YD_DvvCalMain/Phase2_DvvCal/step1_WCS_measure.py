from xwt import xwt
from os.path import join, basename
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from copy import deepcopy
from scipy.optimize import curve_fit
import os
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import time
from obspy.signal.regression import linear_regression

def define_qc(sub_Weight):
    non_zero_ratio = np.count_nonzero(sub_Weight) / sub_Weight.size * 100
    return non_zero_ratio > 25

def define_weight(WXamp, WCoh, WXdt):
    Weight = []
    log_WXamp = np.log(WXamp+np.finfo(float).eps)
    max_val = np.max(np.abs(log_WXamp))
    Weight = log_WXamp / max_val
    # print(Weight)
    Weight[WCoh<0.75] = 0
    Weight[np.abs(WXdt)>0.3] = 0
    return Weight

def xwt_measure(trace_ref, trace_current):
    fs, ns, nt, vpo, freqmin, freqmax, nptsfreq = 50, 3, 0.25, 10, 2, 8, 200
    WXamp, WXspec, WXangle, WCoh, WXdt, freqs, coi = xwt(trace_ref, trace_current, fs, ns, nt, vpo, freqmin, freqmax, nptsfreq)
    return WXamp, WCoh, WXdt, freqs

def dvov_cal_main(trace_ref, trace_current):
    '''
    Parameter:
    lag_time
    t1, t2, t3, t4
    '''
    WXamp, WCoh, WXdt, freqs = xwt_measure(trace_ref, trace_current)
    lag_time = np.linspace(-30, 30, 3001)
    Weight = define_weight(WXamp, WCoh, WXdt)
    # Weighted_WXdt = WXdt * Weight
    pos_begin_coda, neg_begin_coda = 2, -2
    length_coda = []

    for i in range(len(freqs)):
        end_time = (1/freqs[i]) * 20
        if end_time > 60:
            length_coda.append(60)
        else:
            length_coda.append(end_time)
    
    # Coda wave Window selection
    length_coda = np.array(length_coda)
    
    t1, t2, t3, t4 = neg_begin_coda - length_coda[i], neg_begin_coda, pos_begin_coda, pos_begin_coda + length_coda[i]

    mask = ((lag_time >= t1) & (lag_time <= t2)) | ((lag_time >= t3) & (lag_time <= t4))


    dvovs = []
    for i in range(len(Weight)):
        sub_WXdt = WXdt[i][mask]
        sub_Weight = Weight[i][mask]
        none_zero_mask = sub_Weight > 0
        Used_WXdt = sub_WXdt[none_zero_mask]
        Used_Weight = sub_Weight[none_zero_mask]
        Used_lag = lag_time[mask][none_zero_mask]

        qc = define_qc(sub_Weight)

        if qc == True:
            ## Apply the linear regression: -dt/t = dv/v
            # popt, pcov = curve_fit(linear_dvov, Used_lag, Used_WXdt)
            
            p0, p1 = linear_regression(Used_lag, Used_WXdt, weights=Used_Weight)
            dvovs.append(p0 * -1)
            # plt.scatter(Used_lag, Used_WXdt, marker="o", color="red")
            # plt.plot(Used_lag, Used_lag*popt[0]+popt[1], lw=2, ls="--", color="black")
            # plt.show()
        else:
            dvovs.append(np.nan)
    return np.array(dvovs), freqs

def worker(trace_ref, trace_cur):
    return dvov_cal_main(trace_ref, trace_cur)


project_dir = ".."
scf_dir = join(project_dir, "data", "scorr_npz")
dvov_dir = join(project_dir, "output", "scorr_dvov")
scs = ['Wiener_en', 'Wiener_ez', 'Wiener_nz']
used_cpu = 3

if __name__ == '__main__':
    start = time.time()

    scf_files = glob(join(scf_dir, "*"))
    for scf_file in scf_files[:]:
        prefix = basename(scf_file)
        ntnm, stnm = prefix.split('.')[0], prefix.split('.')[1]
        index = os.path.exists(join(dvov_dir, ntnm, stnm))
        if index == False:
            scf_data = np.load(scf_file)
            for sc in scs:
                scf = scf_data[sc]
                trace_ref = scf.mean(axis=0)
                print("-"*60)
                print(f"{ntnm} {stnm} {sc} Measure Start, Total {len(scf)} traces")
                with ProcessPoolExecutor(max_workers=used_cpu) as executor:
                    dvov_func = partial(worker, trace_ref)
                    results = list(executor.map(dvov_func, scf))

                Dvovs = np.array([r[0] for r in results])
                freqs = results[0][1]

                os.makedirs(join(dvov_dir, ntnm, stnm), exist_ok=True)
                dvov_file = f"{sc}.dvov.npz"
                np.savez(join(dvov_dir, ntnm, stnm, dvov_file), Dvovs=Dvovs, freqs=freqs)
                print(f"{ntnm} {stnm} {sc} Measure End, Total {len(scf)} traces")
                print("-"*60)

        else:
            print(f"{ntnm} {stnm} Have Been Calculated!")

    end = time.time()

    print(end - start)