import os
from glob import glob
from os.path import join, basename
import numpy as np
from obspy import UTCDateTime
import h5py
from scipy.linalg import svd
from scipy.signal import wiener

def wiener_filter(scf, SVD_N):
    U, S, V = svd(scf)
    C = np.zeros((len(U), len(V)))

    for i in range(SVD_N):
        si = S[i]
        ui = U[:, i].reshape(-1, 1)
        vi = V[i, :].reshape(1, -1)

        sub_C = si * np.dot(ui, vi)
        sub_C_filtered = wiener(sub_C, (3, 3))

        C += sub_C_filtered
    C = wiener(C, (3, 3))
    return C

def get_SCF(sta, Begin_Date, End_Date):
    Begin_Date = UTCDateTime("2022-07-31T06:00:00")
    End_Date = UTCDateTime("2022-08-19T18:00:00")

    Total_Hour = int((End_Date.timestamp - Begin_Date.timestamp)/3600)  
    en, ez, nz = [], [], []
    Dates = []
    for i in range(Total_Hour):
        ind = str(UTCDateTime(Begin_Date) + i * 3600)[0:19]
        scffile = glob(join(sta, f"*{ind}*"))[0]
        scf = h5py.File(scffile, "r")
        sub_en = scf['Csub_en'][:].mean(axis=0)
        sub_ez = scf['Csub_ez'][:].mean(axis=0)
        sub_nz = scf['Csub_nz'][:].mean(axis=0)
        if np.isnan(sub_en[0]) == True:
            print(f"Error Value {ind}")
        en.append(sub_en)
        ez.append(sub_ez)
        nz.append(sub_nz)
        Dates.append(ind)
    return en, ez, nz, Dates

wave_dir = "../data/scorr"
out_dir = "../data/scorr_npz"
os.makedirs(out_dir, exist_ok=True)
Begin_Date = UTCDateTime("2022-07-31T06:00:00")
End_Date = UTCDateTime("2022-08-19T18:00:00")
Total_Hour = int((End_Date.timestamp - Begin_Date.timestamp)/3600)

nets = glob(join(wave_dir, "*"))
for net in nets[:]:
    ntnm = basename(net)
    stas = glob(join(net, "*"))
    for sta in stas[:]:
        stnm = basename(sta)

        npz_file = f"{ntnm}.{stnm}.npz"
        if os.path.exists(join(out_dir, npz_file)) == False:
            print("-"*10)
            print(f"{ntnm} {stnm} Begin Processing!")

            scf_sta = join(wave_dir, ntnm, stnm)

            # Get the SC from the h5 files
            en, ez, nz, Dates = get_SCF(scf_sta, Begin_Date, End_Date)

            # Apply the Wiener Filter
            Wiener_en = wiener_filter(en, SVD_N=10)
            Wiener_ez = wiener_filter(ez, SVD_N=10)
            Wiener_nz = wiener_filter(nz, SVD_N=10)

            npz_file = f"{ntnm}.{stnm}.npz"

            np.savez(join(out_dir, npz_file), en=en, ez=ez, nz=nz, Dates=Dates, Wiener_ez=Wiener_ez, Wiener_en=Wiener_en, Wiener_nz=Wiener_nz)
            
            print("-"*10)
            print(f"{ntnm} {stnm} End Processing!")
        else:
            print(f"{ntnm} {stnm} Exist!")





