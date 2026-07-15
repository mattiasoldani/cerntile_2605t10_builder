# analyse the FCC TileCal beamtest QC data and build common lists of synced events
# by M. Soldani, 2026

''' compilation:
python3 sync_lists.py
'''

# import imports
import numpy as np
import time

# set settings
pathdatamaster = "../data/"
pathqc = pathdatamaster + "QC/"
pathin = pathdatamaster + "input/"
pathout = pathdatamaster + "output/sync_lists/"
listname = "FERS_runList_beams" # FERS_runList_total / FERS_runList_beams
fers_bd = [0 , 1]
digi_ch = {0:6 , 1:7}
fers_ch = 0

# open run list
runs = np.genfromtxt(pathin+"%s.txt"%listname, dtype=int, comments="#")

####################
# loop over all runs, search files existing in QC output data, build common list of events

out_dict = {}
for run in runs:
    fers_list = {}
    
    try:
        for bd in fers_bd:
            fers_listname = "merged_run%d_digi%d_fers%d_hg%d_events.txt" % (run, digi_ch[bd], bd, fers_ch)
            fers_list[bd] = np.genfromtxt(pathqc+"run%d/"%run+fers_listname, usecols=(2, 19), unpack=True)

    except:
        print("Run%d not (fully) found, skipping..." % run)
        print("---")
        continue

    print("Run%d found --> building common list..." % run)
    t0 = time.time()

    ls_digi_id = list(set(np.concatenate([fers_list[i][0] for i in fers_list.keys()]).astype(int)))
    
    out_dict[run] = {}
    
    for idigi in ls_digi_id:
        out_dict[run][idigi] = [np.nan for i in fers_bd]
        for bd in fers_list:
            if (idigi in fers_list[bd][0]): 
                out_dict[run][idigi][bd] = int(fers_list[bd][1][fers_list[bd][0]==idigi][0])

    with open(pathout+"%d.txt"%run, "w") as f:
        for i in out_dict[run].keys():
            out_str = "%d" % i
            for ibd, bd in enumerate(fers_bd):
                out_str += ",-1" if bool(np.isnan(out_dict[run][i][ibd])) else ",%d" % out_dict[run][i][ibd]
            f.write(out_str+"\n")

    t1 = time.time()
    print("Done! In %d seconds." % (t1-t0))
    print("---")

####################

print("All done!")
