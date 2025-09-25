using Distributed
# Choose the nums of multi.... Consider the RAM Pleas
# addprocs(10)

# @everywhere begin

using SeisNoise, SeisIO, Glob
using HDF5
function pre_process_fft(s, fs, cc_len, cc_step, freqmin, freqmax)
    process_raw!(s, fs)
    R = RawData.([s], cc_len, cc_step)
    detrend!.(R)
    taper!.(R)
    bandpass!.(R, freqmin, freqmax)
    whiten!.(R, freqmin, freqmax)
    # whiten!.(R, freqmin, freqmax)
    onebit!.(R)
    F = rfft.(R)
    # whiten!.(F, freqmin, freqmax)
    return F
end

function scorr_cal(F_e, F_n, F_z, maxlag)
    Csub_en = correlate(F_e[1], F_n[1], maxlag, corr_type="CC")
    Csub_ez = correlate(F_e[1], F_z[1], maxlag, corr_type="CC")
    Csub_nz = correlate(F_n[1], F_z[1], maxlag, corr_type="CC")
    return Csub_en, Csub_ez, Csub_nz
end

function get_scorr_cal_list(data_dir, out_dir, fs, cc_len, cc_step, freqmin, freqmax, maxlag)
    scorr_para_list_list = []
    nets = glob("*", data_dir)
    for net in nets
        stas = glob("*", net)
        ntnm = basename(net)
        for sta in stas
            stnm = basename(sta)
            if !isdir(joinpath(out_dir, ntnm, stnm))
                mkpath(joinpath(out_dir, ntnm, stnm))
            end

            efiles = glob("*BHE*", sta)

            for efile in efiles
                file = basename(efile)
                file_info = split(file, ".")
                # print(efile)
                nfile = join([file_info[1], file_info[2], file_info[3], file_info[4], "BHN", "sac"], ".")
                zfile = join([file_info[1], file_info[2], file_info[3], file_info[4], "BHZ", "sac"], ".")
                nfile = joinpath(data_dir, ntnm, stnm, nfile)
                zfile = joinpath(data_dir, ntnm, stnm, zfile)
                h5_file = join([ntnm, stnm, file_info[3], "h5"], ".")
                o_h5_file = joinpath(out_dir, ntnm, stnm, h5_file)
                scorr_para_list = [efile, nfile, zfile, fs, cc_len, cc_step, freqmin, freqmax, maxlag, o_h5_file]
                push!(scorr_para_list_list, scorr_para_list)
            end
        end
    end
    return scorr_para_list_list
end

function scorr_cal_main(param)
    efile, nfile, zfile = param[1], param[2], param[3]
    fs, cc_len, cc_step = param[4], param[5], param[6]
    freqmin, freqmax, maxlag = param[7], param[8], param[9]
    o_h5_file = param[10]
    s_e = read_data("sac", efile)
    s_n = read_data("sac", nfile)
    s_z = read_data("sac", zfile)
    F_e = pre_process_fft(s_e, fs, cc_len, cc_step, freqmin, freqmax)
    F_n = pre_process_fft(s_n, fs, cc_len, cc_step, freqmin, freqmax)
    F_z = pre_process_fft(s_z, fs, cc_len, cc_step, freqmin, freqmax)
    Csub_en, Csub_ez, Csub_nz = scorr_cal(F_e, F_n, F_z, maxlag)
    h5open(o_h5_file, "w") do file
        write(file, "Csub_en", Csub_en.corr)
        write(file, "Csub_ez", Csub_ez.corr)
        write(file, "Csub_nz", Csub_nz.corr)
    end
    print(o_h5_file, " Over !", "\n")
end

fs = 50.0
cc_len, cc_step = 1800, 900
freqmin, freqmax = 2, 8
maxlag = 60

data_dir = "../data/hourlydaya"
out_dir = "../data/scorr"
scorr_para_list_list = get_scorr_cal_list(data_dir, out_dir, fs, cc_len, cc_step, freqmin, freqmax, maxlag)

start_time = time()

for scorr_para_list in scorr_para_list_list
    scorr_cal_main(scorr_para_list)
end

end_time = time()

print("Total Minutes Use For Scorr Cal: ", (end_time - start_time) / 60, "Min")
