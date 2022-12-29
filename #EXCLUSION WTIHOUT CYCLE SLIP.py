#EXCLUSION WTIHOUT CYCLE SLIP

import numpy as np
import math
import sys
import os
import datetime
import glob
import time

path_tec="C:\\"
START=time.time()
threshold=1.000
eps=0.01

year4 = 2017
year2 = year4 % 100
#
#
#
days = [240,257]
start_time=[20.9500,21.9500]

n_file=[0 for i in range(len(start_time))]

start=time.time()

for h in range(len(days)):
    day=days[h]
    s_time=start_time[h]

    first_day_path=path_tec+"data\\"+str(year4)+"\\"+str(day)
    sat_folders=glob.glob(first_day_path+"/*")

    for sat_folder in sat_folders:
        print(sat_folder)
        sat=sat_folder[-2:]
        print(sat)
        files=glob.glob(sat_folder+"/*")
        for file in files:
            rec=file[-8:-4]
            print("sat {s:02} , rec {r:04} finished : {t:.3f}".format(s=sat,r=rec,t=time.time()-start))
            with open(file,"r") as f:
                fline=[s.strip() for s in f.readlines()]
            
            before_cycle_slip=18
            has_valid_time_data=False
            
            for i in range(18,len(fline)):
                if "# cycle slip" in fline[i]:
                    before_cycle_slip=i
                else:
                    tmp_time=float(fline[i].split()[0])
                    if abs(tmp_time-s_time)<eps:
                        has_valid_time_data=True
                        break
            
            if not has_valid_time_data:
                continue
            
            n_file[h]+=1
            
            os.makedirs(path_tec+"selected_data\\"+str(year4)+"\\"+str(day)+"\\"+"G{s:02}".format(s=sat),exist_ok=True)
            output_file=path_tec+"selected_data\\"+str(year4)+"\\"+str(day)+"\\"+"G{s:02}".format(s=sat)+"\\"+"{r:04}".format(r=rec)+".txt"

            with open(output_file,"w") as o_f:
                for i in range(17):
                    print(fline[i],file=o_f)
                o_f.write("# SELECTED TIME : {t:2.4f}\n".format(t=s_time))

            idx=before_cycle_slip+1
            
            while True:
                if "# cycle slip" in fline[idx]:
                    break
                if idx==len(fline)-1:
                    with open(output_file,"a") as o_f:
                        print(fline[idx],file=o_f)
                        break
                else:
                    with open(output_file,"a") as o_f:
                        print(fline[idx],file=o_f)
                        idx+=1

print("number of datas : ",n_file)
