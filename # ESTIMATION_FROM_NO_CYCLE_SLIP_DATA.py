# ESTIMATION_FROM_NO_CYCLE_SLIP_DATA

import numpy as np
import math
import sys
import os
import datetime
import glob
import time


def AIC(est_data,real_data,deg):
    L=len(est_data)
    first=L*(math.log(2*math.pi)+1)
    sigma=0.0
    for i in range(L):
        sigma+=(est_data[i]-real_data[i])**2
    sigma/=L
    if sigma==0.0:
        second=sys.float_info.min
    else:
        second=L*math.log(sigma**2)
    third=2*(deg+2)
    #print(first,second,third)
    return first+second+third

def polyfit(t_data,x_data,mu_data,deg):
    N=len(t_data)
    A=np.full((deg+1,deg+1),0.0,dtype="float64")
    # x=np.full((deg+1,1),0.0,dtype="float128")
    b=np.full((deg+1,1),0.0,dtype="float64")
    
    for i in range(deg+1):
        r_tmp=0.0
        for k in range(N):
            r_tmp+=mu_data[k]*x_data[k]*np.power(t_data[k],i)
        b[i]=r_tmp
        for j in range(deg+1):
            tmp=0.0
            for k in range(N):
                tmp+=mu_data[k]*np.power(t_data[k],i+j)
            A[i,j]=tmp
    
    x=np.linalg.solve(A,b)
    return x

def poly2data(t_data,a) -> np.ndarray:
    L=len(t_data)
    M=len(a)
    ans=np.full((L),0.0,dtype="float64")
    for k in range(L):
        tmp=0.0
        for l in range(M):
            tmp+=a[l]*np.power(t_data[k],l)
        ans[k]=tmp
    return ans

def best_deg_selector(t_data,x_data,mu_data,m_deg,M_deg):
    L=len(t_data)
    best_deg=-1
    best_AIC=10000000.0
    best_est_data=[]
    est_data_list=[]
    best_cor=np.ndarray
    for deg in range(m_deg,M_deg):
        cor=polyfit(t_data,x_data,mu_data,deg)
        est_data=poly2data(t_data,cor)
        deg_AIC=AIC(est_data,x_data,deg)
        if deg_AIC<best_AIC:
            best_cor=cor
            best_AIC=deg_AIC
            best_deg=deg
            best_est_data=est_data
    return best_deg,best_AIC,best_est_data,best_cor

path_tec="C:\\"
START=time.time()
threshold=1.000
eps=0.01
m_deg=1
M_deg=12

year4 = 2017
year2 = year4 % 100
#
#
#
days = [240,257]

exclusive_window_start=[20.967,21.950]
exclusive_window_end=[21.200,22.267]
n_file=[0 for i in range(len(exclusive_window_start))]
margin=0.4

for h in range(len(days)):
    day=days[h]
    w_start=exclusive_window_start[h]
    w_end=exclusive_window_end[h]
    
    path=path_tec+"selected_data\\"+str(year4)+"\\"+str(day)
    sat_folders=glob.glob(path+"/*")
    
    for sat_folder in sat_folders:
        sat=sat_folder[-2:]
        obs_files=glob.glob(sat_folder+"/*")
        for file in obs_files:
            rec=file[-8:-4]
            print("sat {s:02} , rec {r:04} finished : {t:.3f}".format(s=sat,r=rec,t=time.time()-START))
            with open(file,"r") as f:
                lines=[s.strip() for s in f.readlines()]
                
            slines=[line.split() for line in lines]
            
            t_data=[]
            stec_data=[]
            sat_x=[]
            sat_y=[]
            sat_z=[]
            rec_x=[]
            rec_y=[]
            rec_z=[]
            
            
            start_idx=0
            
            for i in range(len(lines)):
                if not "#" in lines[i]:
                    start_idx=i
                    break
                    
            # print(start_idx,len(lines))
            
            for i in range(start_idx,len(lines),1):
                t_data.append(float(slines[i][0]))
                stec_data.append(float(slines[i][1]))
                sat_x.append(float(slines[i][2]))
                sat_y.append(float(slines[i][3]))
                sat_z.append(float(slines[i][4]))
                rec_x.append(float(slines[i][5]))
                rec_y.append(float(slines[i][6]))
                rec_z.append(float(slines[i][7]))
            
            # print(len(t_data),t_data[0],t_data[len(t_data)-1])
            
            if t_data[0]>w_start-margin or t_data[len(t_data)-1]<w_end+margin:
                continue
            
            n_file[h]+=1
            L=len(t_data)
            mu=[1.0 for i in range(L)]
            for i in range(L):
                if w_start<t_data[i]<w_end:
                    mu[i]=0.0
                else:
                    mu[i]=1.0
            
            start_time=t_data[0]
            
            for i in range(L):
                t_data[i]-=start_time
            
            best_deg,best_AIC,best_est_data,best_cor=best_deg_selector(t_data,stec_data,mu,m_deg,M_deg)
            
            for i in range(L):
                t_data[i]+=start_time
            
            os.makedirs(path_tec+"estimated_data\\"+str(year4)+"\\"+str(day)+"\\"+"G{s:02}".format(s=sat),exist_ok=True)
            output_file=path_tec+"estimated_data\\"+str(year4)+"\\"+str(day)+"\\"+"G{s:02}".format(s=sat)+"\\"+"{r:04}".format(r=rec)+".txt"
            
            with open(output_file,"w") as o_f:
                for i in range(18):
                    print(lines[i],file=o_f)
                for i in range(L):
                    o_f.write("{t:.4f} {s:.5f} {o_x:.4f} {o_y:.4f} {o_z:.4f} {r_x:.4f} {r_y:.4f} {r_z:.4f}\n".format(
                        t=t_data[i],s=best_est_data[i]-stec_data[i],o_x=sat_x[i],o_y=sat_y[i],o_z=sat_z[i],r_x=rec_x[i],r_y=rec_y[i],r_z=rec_z[i]
                    ))
            

print("number of datas : ",n_file)