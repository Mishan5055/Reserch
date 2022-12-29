# TOMO continuous days connect
import numpy as np
import math
import sys
import os
import datetime
import glob
import time

path_tec = "D:\\"
START=time.time()
threshold=1.000

year4 = 2017
year2 = year4 % 100
#
#
#
days = [257,258]
# number of record

path_tec="C:\\"

first_day=days[0]
second_day=days[1]

first_day_path=path_tec+"data\\"+str(year4)+"\\"+str(first_day)
sat_folders=glob.glob(first_day_path+"/*")

for sat_folder in sat_folders:
    first_files=glob.glob(sat_folder+"/*")
    for first_file in first_files:
        print(first_file,"{:.3f}".format(time.time()-START))
        first_file_splited=first_file.split(str(first_day),1)
        second_file=first_file_splited[0]+str(days[1])+first_file_splited[1]
        
        with open(first_file,"r") as f_f:
            f_lines=[s.strip() for s in f_f.readlines()]
        
        with open(second_file,"r") as s_f:
            s_lines=[s.strip() for s in s_f.readlines()]
        # sslines=[line.split() for line in s_lines]
        for k in range(1,33):
            os.makedirs("C:\\data\\{y:04}\\{s1:03}-{s2:03}\\G{sat:02}\\".format(y=year4,s1=first_day,s2=second_day,sat=k),exist_ok=True)
        output_file=path_tec+"data\\"+str(year4)+"\\"+str(first_day)+"-"+str(second_day)+first_file_splited[1]

        flag=True
        last_obs_bool=True
        with open(output_file,"w") as o_f:
            for i in range(len(f_lines)):
                print(f_lines[i],file=o_f)
            if len(f_lines)>0:
                last_obs=float(f_lines[len(f_lines)-1].split()[1])
            else:
                last_obs_bool=False
            for i in range(18,len(s_lines)):
                if "#" in s_lines[i]:
                    print(s_lines[i],file=o_f)
                else:
                    obs_time=float(s_lines[i][:6])+24.0
                    obs_data=s_lines[i][7:]
                    if flag:
                        flag=False
                        
                        first_obs=float(s_lines[i].split()[1])
                        if last_obs_bool and abs(first_obs-last_obs)>threshold:
                            print("# cycle slip detected",file=o_f)
                    print("{:.4f}".format(obs_time),obs_data,file=o_f)

