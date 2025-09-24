import os
from glob import glob
from os.path import join, basename
import numpy as np
import datetime

def get_dates():
    start_time = datetime.datetime(2022, 7, 31, 6, 0, 0)
    end_time = datetime.datetime(2022, 8, 19, 18, 0, 0)

    delta = datetime.timedelta(hours=1)

    dates = []
    current = start_time
    while current < end_time:
        dates.append(current)
        current += delta
    return dates

project_dir = ".."
# scorr_dvov = join(project_dir, "scorr_dvov_qc")
scorr_dvov = join(project_dir, "output", "scorr_dvov")

nets = glob(join(scorr_dvov, "*"))
for net in nets[:]:
    ntnm = basename(net)
    stas = glob(join(net, "*"))
    for sta in stas[:]:
        stnm = basename(sta)
        dvov_files = glob(join(sta, "*"))
        Dvovs = []
        for dvov_file in dvov_files:
            dvov_data = np.load(dvov_file)
            # print(dvov_file)
            sub_dvov, freqs = dvov_data['Dvovs'], dvov_data['freqs']
            Dvovs.append(sub_dvov)
        
        Dvovs = np.array(Dvovs)
        mean_Dvovs = np.nanmean(Dvovs, axis=0)
        mean_dvov = np.nanmean(mean_Dvovs, axis=1)
        dates = get_dates()
        np.savez(join(sta, f"{ntnm}.{stnm}.dvov.npz"), Dvovs=mean_Dvovs, mean_dvov=mean_dvov, dates=dates, freqs=freqs)
        print(f"{ntnm} {stnm}  Over!")
