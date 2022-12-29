import numpy as np
import math
import sys
import os
import datetime
import glob
import time

maxrec=1500
maxsat = 40
threshold=1.0

vel_light = 2.9979*pow(10, 8)
f_l1 = 1575420000
f_l2 = 1227600000

year4 = 2016
year2 = year4 % 100
days = [105,107]

def xyztoll(xyz):
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    a = 6378137.0
    f = 1.0/298.257223563
    b = a*(1.0-f)
    p = math.sqrt(x*x+y*y)
    theta = np.arctan2(z*a, p*b)
    e0 = (a*a-b*b)/(a*a)
    e1 = (a*a-b*b)/(b*b)
    phi = np.arctan2((z+e1*b*(np.sin(theta)**3)), (p-e0*a*(np.cos(theta)**3)))
    lam = np.arctan2(y, x)
    lat = np.rad2deg(phi)
    lon = np.rad2deg(lam)
    return [lat, lon]


def kepler(dmk, e):
    thres = pow(10, -14)
    niteration = 0
    ek = dmk
    diff = ek+e*math.sin(ek)-dmk
    while abs(diff) > thres:
        diff = ek+e*math.sin(ek)-dmk
        partial = 1-e*math.cos(ek)
        ek = ek-diff/partial
        niteration += 1
        if niteration > 100:
            print("The calculation was terminated because the iteration of Newton's method in the Kep1er function exceeded 100.")
            break
    return ek


def gtxyz(time1, time0, ele):

    # data format of navigation file
    #       0:IODE  1:Crs  2:delta-n  3:m0
    #       4:Cuc   5:e    6:Cus      7:root-a
    #       8:Toe   9:Cic 10:Omega   11:Cis
    #      12:i0   13:Crc 14:omega   15:OmegaDot
    #      16:iDot 17-27: not used

    GM = 3.986005*pow(10, 14)
    omega_dot_e = 7.292115*pow(10, -5)
    # (1)
    a = ele[7]**2
    # (2)
    dnzero = math.sqrt(GM*pow(1/a, 3))
    # (3)
    tk = (time1-time0)*60.0*60.0
    # (4)
    dn = dnzero+ele[2]
    # (5)
    dmk = ele[3]+dn*tk
    # (6)
    ek = kepler(dmk, ele[5])
    # (7)
    cosvk = (math.cos(ek)-ele[5])/(1.0-ele[5]*math.cos(ek))
    sinvk = math.sqrt(1.0-ele[5]*ele[5])*math.sin(ek)/(1.0-ele[5]*math.cos(ek))
    vk = math.atan2(sinvk, cosvk)
    # (8)
    phik = vk+ele[14]
    # (9)
    delta_uk = ele[6]*math.sin(2.0*phik)+ele[4]*math.cos(2.0*phik)
    uk = phik+delta_uk
    # (10)
    delta_rk = ele[1]*math.sin(2.0*phik)+ele[13]*math.cos(2.0*phik)
    rk = a*(1.0-ele[5]*math.cos(ek))+delta_rk
    # (11)
    delta_dik = ele[11]*math.sin(2.0*phik)+ele[9]*math.cos(2.0*phik)
    dik = ele[12]+delta_dik+ele[16]*tk
    # (12)
    xdashk = rk*math.cos(uk)
    ydashk = rk*math.sin(uk)
    # (13)
    omegak = ele[10]+(ele[15]-omega_dot_e)*tk-ele[8]*omega_dot_e
    # (14)
    dx = xdashk*math.cos(omegak)-ydashk*math.cos(dik)*math.sin(omegak)
    dy = xdashk*math.sin(omegak)+ydashk*math.cos(dik)*math.cos(omegak)
    dz = ydashk*math.sin(dik)
    return [dx, dy, dz]

def load_navi_file_v3_02(file): # for GPS
    idsat = np.zeros(maxrec,dtype=int)
    tiempo = np.zeros(maxrec)
    oel = np.zeros((maxrec,28))
    i1stsat = [[] for i in range(maxsat)]
    with open(file,"r") as f_n:
        # Header
        while True:
            s_line=f_n.readline()
            if "END OF HEADER" in s_line:
                break
            if not s_line:
                break

        # Navigation
        # krec is record index
        krec=0
        s_line=f_n.readline()
        while True:
            idsat[krec] = int(s_line[1:3])
            ih = int(s_line[15:17])
            imin = int(s_line[18:20])
            sec = float(s_line[21:23])
            # UTC[h]
            tiempo[krec] = ih+imin/60.0+sec/3600.0
            s_line = f_n.readline()
            if not s_line:
                break
            nav_flag = 0
            for irec in range(7):
                ii = irec*4
                for j in range(4):
                    m_str = s_line[4+19*j:4+15+19*j].strip()
                    e_str = s_line[20+19*j:20+3+19*j].strip()
                    if not m_str == "" and not e_str == "":
                        oel[krec][ii+j] = float(m_str)*pow(10,int(e_str))
                        nav_flag += 1
                    else:
                        oel[krec][ii+j] = 0.0
                s_line = f_n.readline()
            if nav_flag >=26:
                krec += 1
            if not s_line:
                break
        mrec = krec

    #
    # finding 1st appearance of satellites
    # i1stsat[isat] ... index which No.[isat] satellite's record start
    #

    for krec in range(mrec):
        i1stsat[idsat[krec]].append(krec)
    return idsat,tiempo,oel,i1stsat

def load_navi_file_v2_11(file): # for GPS
    oel = np.zeros((maxrec,28)) # storing navigation data
    idsat = np.zeros(maxrec,dtype=int)
    tiempo = np.zeros(maxrec)
    i1stsat = [[] for i in range(maxsat)]

    with open(file,mode='r',errors='replace') as f_n:
        #
        # reading header of navigation file
        #
        while True:
            s_line = f_n.readline()
            # if 'RINEX VERSION' in s_line:
            #     ver = s_line[5:9]
            if 'END OF HEADER' in s_line:
                break
            if not s_line:
                break
        #
        # reading data of navigation file
        #
        krec = 0
        s_line = f_n.readline()
        while True:
            idsat[krec] = int(s_line[0:2])
            # iy = int(s_line[3:5])
            # imon = int(s_line[6:8])
            # iday = int(s_line[9:11])
            ih = int(s_line[12:14])
            imin = int(s_line[15:17])
            sec = float(s_line[18:22])
            tiempo[krec] = ih+imin/60.0+sec/3600.0
            s_line = f_n.readline()
            for irec in range(7):
                ii = irec*4
                if irec!=6:
                    oel[krec][ii] = float(s_line[3:3+15])*pow(10,int(s_line[19:19+3]))
                    oel[krec][ii+1] = float(s_line[3+19:3+19+15])*pow(10,int(s_line[19+19:19+19+3]))
                    oel[krec][ii+2] = float(s_line[3+19*2:3+19*2+15])*pow(10,int(s_line[19+19*2:19+19*2+3]))
                    oel[krec][ii+3] = float(s_line[3+19*3:3+19*3+15])*pow(10,int(s_line[19+19*3:19+19*3+3]))
                else:
                    oel[krec][ii] = float(s_line[3:3+15])*pow(10,int(s_line[19:19+3]))
                    oel[krec][ii+1] = 0.0
                    oel[krec][ii+2] = 0.0
                    oel[krec][ii+3] = 0.0
                s_line = f_n.readline()
            krec += 1
            if not s_line:
                break
        mrec = krec

    for krec in range(mrec):
        i1stsat[idsat[krec]].append(krec)
    return idsat,tiempo,oel,i1stsat


def load_navigation_file(file,rinux_ver):
    if rinux_ver=="3.02":
        return load_navi_file_v3_02(file)
    elif rinux_ver=="2.11":
        return load_navi_file_v2_11(file)
    else:
        print("Not supported")
        exit()

# G_stec,G_stec_bool,G_sat_bool,J_stec,J_stec_bool,J_sat_bool,p_o
def load_observ_file_v3_02(file):
    with open(file, "r") as o_f:
        lines = [s.rstrip() for s in o_f.readlines()]

    slines = [line.split() for line in lines]
    header_idx = 0
    # p_o is receiver ECEF XYZ
    p_o = np.zeros(3)
    G_l1num = -1
    G_l2num = -1
    J_l1num = -1
    J_l2num = -1

    for idx in range(len(lines)):
        if 'APPROX POSITION XYZ' in lines[idx]:
            p_o[0] = float(lines[idx][1:14])
            p_o[1] = float(lines[idx][15:28])
            p_o[2] = float(lines[idx][29:42])
        if 'OBS TYPES' in lines[idx]:
            gnss_type = lines[idx][0:1]
            if gnss_type == "G":  # GPS
                wave_sum = int(lines[idx][4:6])
                for i in range(wave_sum):
                    j = i % 13
                    k = i//13
                    wave = lines[idx+k][7+4*j:10+4*j]
                    if wave == "L1C":
                        G_l1num = i
                    elif wave == "L2W":
                        G_l2num = i
            if gnss_type == "J":
                wave_sum = int(lines[idx][4:6])
                for i in range(wave_sum):
                    j = i % 13
                    k = i//13
                    wave = lines[idx+k][7+4*j:10+4*j]
                    if wave == "L1C":
                        J_l1num = i
                    elif wave == "L2X":
                        J_l2num = i
        if 'END OF HEADER' in lines[idx]:
            header_idx = idx
            break
        if not lines[idx]:
            continue

    G_stec = np.zeros((maxsat, 2880))
    G_stec_bool = np.full((maxsat, 2880), False, dtype=bool)
    G_sat_bool = np.full(maxsat, False, dtype=bool)

    J_stec = np.zeros((maxsat, 2880))
    J_stec_bool = np.full((maxsat, 2880), False, dtype=bool)
    J_sat_bool = np.full(maxsat, False, dtype=bool)

    flag = False

    for idx in range(header_idx+1, len(slines)):
        # print(slines[idx])
        if "COMMENT" in lines[idx]:
            continue
        if not lines[idx]:
            continue
        if lines[idx][0:1] == ">":
            obs_year_str = lines[idx][4:6]
            if str(year2) != obs_year_str:
                flag = False
                continue
            else:
                ih = int(lines[idx][13:15])
                imin = int(lines[idx][16:18])
                sec = int(lines[idx][19:21])
                # time idx in observ
                iepoc = sec//30+2*imin+2*60*ih
                sat_sum = int(lines[idx][33:35])
                # print(ih, imin, sec, iepoc, sat_sum)
                # a = input()
                idx += 1
                flag = True
        if flag:
            gnss_type = lines[idx][0:1]
            isat = int(lines[idx][1:3])
            if gnss_type == "G":
                G_sat_bool[isat] = True
                wave_l1_bool = False
                wave_l2_bool = False
                str_l1 = lines[idx][3+16*G_l1num:17+16*G_l1num].strip()
                str_l2 = lines[idx][3+16*G_l2num:17+16*G_l2num].strip()
                if not str_l1 == "":
                    wave_l1 = float(str_l1)
                    wave_l1_bool = True
                if not str_l2 == "":
                    wave_l2 = float(str_l2)
                    wave_l2_bool = True
                if wave_l1_bool and wave_l2_bool:
                    G_stec[isat, iepoc] = 6.05 / (-0.65)*(vel_light*wave_l2/f_l2-vel_light*wave_l1/f_l1)
                    # print(wave_l1, wave_l2, G_stec[isat, iepoc])
                    G_stec_bool[isat, iepoc] = True
                    # a = input()
            if gnss_type == "J":
                J_sat_bool[isat] = True
                wave_l1_bool = False
                wave_l2_bool = False
                str_l1 = lines[idx][3+16*J_l1num:17+16*J_l1num].strip()
                str_l2 = lines[idx][3+16*J_l2num:17+16*J_l2num].strip()
                if not str_l1 == "":
                    wave_l1 = float(str_l1)
                    wave_l1_bool = True
                if not str_l2 == "":
                    wave_l2 = float(str_l2)
                    wave_l2_bool = True
                if wave_l1_bool and wave_l2_bool:
                    J_stec[isat, iepoc] = 6.05 / (-0.65)*(vel_light*wave_l2/f_l2-vel_light*wave_l1/f_l1)
                    # print(wave_l1, wave_l2, J_stec[isat, iepoc])
                    J_stec_bool[isat, iepoc] = True
                    # a = input()

    return G_stec,G_stec_bool,G_sat_bool,p_o # ,J_stec,J_stec_bool,J_sat_bool

# stec,stec_bool,sat_bool,p_o
def load_observ_file_v2_11(file):
    with open(file,mode='r',errors='replace') as f_o:
        #
        # reading header of observation file
        #
        p_o = np.zeros(3)
        l1num = -1
        l2num = -1
        # c1num = -1
        # p1num = -1
        # p2num = -1
        while True:
            s_line = f_o.readline()
            # if 'RINEX VERSION' in s_line:
            #     ver = s_line[5:9]
            if 'APPROX POSITION XYZ' in s_line:
                p_o[0] = float(s_line[1:14])
                p_o[1] = float(s_line[15:28])
                p_o[2] = float(s_line[29:42])
            if 'TYPES OF OBSERV' in s_line:
                wave_sum = int(s_line[4:6])
                for i in range(wave_sum):
                    j = i%9
                    wave = s_line[9+6*j:12+6*j]
                    if wave == " L1":
                        l1num = i
                    elif wave == " L2":
                        l2num = i
                    # elif wave == " C1":
                    #     c1num = i
                    # elif wave == " P1":
                    #     p1num = i
                    # elif wave == " P2":
                    #     p2num = i
                    if j == 8:
                        s_line = f_o.readline()
            if 'END OF HEADER' in s_line:
                break
            if not s_line:
                break
        # if p1num == -1 and c1num != -1:
        #     p1num = c1num

        #
        # reading data of observation file
        #
        stec = np.zeros((maxsat,2880))
        stec_bool = np.zeros((maxsat,2880),dtype=bool)
        sat_list = np.zeros(maxsat,dtype=np.int)
        sat_bool = np.full(maxsat,False,dtype=bool)
        while True:
            s_line = f_o.readline()
            while True:
                iy_str = s_line[1:3]
                if not s_line:
                    break
                if str(year2) != iy_str:
                    s_line = f_o.readline()
                else:
                    break
            if not s_line:
                break
            # imon = int(s_line[4:6])
            # iday = int(s_line[7:9])
            ih = int(s_line[10:12])
            imin = int(s_line[13:15])
            isec = int(s_line[16:18])
            # fsec = float(s_line[16:26])
            iepoc = isec//30+2*imin+2*60*ih
            # flag = int(s_line[28:29])
            # number of sat
            sat_sum = int(s_line[30:32])
            for isat in range(maxsat):
                sat_list[isat] = 0
            for i in range(sat_sum):
                j = i%12
                sat_num = int(s_line[33+3*j:35+3*j])
                # sat_list[i] is i th sat number
                sat_list[i] = sat_num
                sat_bool[sat_num] = True
                if j == 11 and sat_sum > i+1:
                    s_line = f_o.readline()
            for isat in range(maxsat):
                # sat_list[isat] is not empty
                if sat_list[isat] != 0:
                    sat_num = sat_list[isat]
                    wave_l1_bool = False
                    wave_l2_bool = False
                    for j in range(-(-wave_sum//5)):
                        s_line = f_o.readline()
                        if 5*j <= l1num and l1num < 5*(j+1):
                            k = l1num%5
                            str_l1 = s_line[1+16*k:14+16*k].strip()
                            if not str_l1 == "":
                                wave_l1 = float(str_l1)
                                wave_l1_bool = True
                        if 5*j <= l2num and l2num < 5*(j+1):
                            k = l2num%5
                            str_l2 = s_line[1+16*k:14+16*k].strip()
                            if not str_l2 == "":
                                wave_l2 = float(str_l2)
                                wave_l2_bool = True
                    if wave_l1_bool and wave_l2_bool:
                        stec[sat_num,iepoc] = 6.05/(-0.65)*(vel_light*wave_l2/f_l2-vel_light*wave_l1/f_l1)
                        stec_bool[sat_num,iepoc] = 1
                        # if abs(stec[isat,iepoc]-pre_stec[isat]) > threshold or iepoc == 0:
                        #     bias[isat,iepoc] = stec[isat,iepoc]
                        #     pre_bias[isat] = bias[isat,iepoc]
                        #     period_start_epoc[isat,period_sum[isat]] = iepoc
                        #     period_sum[isat] += 1
                        # else:
                        #     bias[isat,iepoc] = pre_bias[isat]
                        # period_end_epoc[isat,period_sum[isat]-1] = iepoc
                        # pre_stec[isat] = stec[isat,iepoc]
            if not s_line:
                break
    return stec,stec_bool,sat_bool,p_o

def load_observation_file(file: str,rinux_ver: str):
    if rinux_ver=="3.02":
        return load_observ_file_v3_02(file)
    elif rinux_ver=="2.11":
        return load_observ_file_v2_11(file)
    else:
        print("Not supported")
        exit()

def make_tomo_loadfile(path_nd,stec,stec_bool,sat_bool,i1stsat,p_o,isat):
        if sat_bool[isat]:
            with open(path_nd,mode='w') as f_nd:
                flag = 0
                for iepoc in range(2880):
                    if stec_bool[isat,iepoc]:
                        ttime = iepoc/120.0
                        # jrec = i1stsat[isat]
                        if len(i1stsat[isat])>1:
                            if flag == 0:
                                bias = stec[isat,iepoc]
                                bias = 0.0
                                pre_stec1 = 0.0
                                #
                                # writting header
                                #
                                f_nd.write("# PRN G{sat:02d}\n".format(sat=isat))
                                f_nd.write("# BIAS {b}\n".format(b=bias))
                                f_nd.write("# \n")
                                f_nd.write("# RINEX VER G_2.11\n")
                                f_nd.write("# FILE {s:>4}{d:03d}0.{y2:02d}o {s:>4}{d:03d}0.{y2:02d}n\n".format(s=sta,d=day,y2=year2))
                                f_nd.write("# \n")
                                f_nd.write("# RUN BY\n")
                                f_nd.write("# PROGRAM {p}\n".format(p="RNX2TOMO"))
                                f_nd.write("# UTCTIME {t}\n".format(t=datetime.datetime.now()))
                                f_nd.write("# \n")
                                f_nd.write("# WAVE FREQUENCY\n")
                                f_nd.write("# L1 {l1}\n".format(l1=f_l1))
                                f_nd.write("# L2 {l2}\n".format(l2=f_l2))
                                f_nd.write("# \n")
                                f_nd.write("# CYCLE SLIP THRESHOLD {th:8.6f}\n".format(th=threshold))
                                f_nd.write("# \n")
                                f_nd.write("# END OF HEADER\n")
                                flag = 1
                            min_d_t=10000.0
                            min_d_t_idx=-1
                            for i in range(len(i1stsat[isat])):
                                if abs(tiempo[i1stsat[isat][i]]-ttime)<min_d_t:
                                    min=abs(tiempo[i1stsat[isat][i]]-ttime)
                                    min_d_t_idx=i
                            time0 = tiempo[i1stsat[isat][min_d_t_idx]]
                            try:
                                xyz = gtxyz(ttime,time0,oel[i1stsat[isat][min_d_t_idx]])
                            except:
                                continue
                            stec1 = stec[isat,iepoc] - bias
                            if abs(stec1-pre_stec1)>threshold:
                                f_nd.write("# cycle slip detection\n")
                            f_nd.write('{t:.4f} {s:.4f} {x:.7f} {y:.7f} {z:.7f} {p_x:.7f} {p_y:.7f} {p_z:.7f}\n'.format(t=ttime,sat=isat,s=stec1,x=xyz[0],y=xyz[1],z=xyz[2],p_x=p_o[0],p_y=p_o[1],p_z=p_o[2]))
                            pre_stec1 = stec1



path_tec = "D:\\"

year4 = 2017
year2 = year4 % 100
days = [257,258,240]
# number of record

start=time.time()

for day in days:
    folder=path_tec+"1_tmp\\"+str(year4)+"\\"+str(day)+"/*"
    # print(folder)
    files=glob.glob(folder)
    stas=set()
    for file in files:
        stas.add(file[-12:-8])

    print(len(stas))
    idx=0
    for sta in stas:
        
        print("station:{s:4}, day:{d:03} is started. {t:.3f} {i}/{f}".format(s=sta,d=day,t=time.time()-start,i=idx,f=len(stas)))
        idx+=1
        # print(day,sta)

        #
        # readling observation file
        #
        idsat = np.zeros(maxrec,dtype=int)
        tiempo = np.zeros(maxrec)
        oel = np.zeros((maxrec,28))
        i1stsat = [[] for i in range(maxsat)]

        path_n = path_tec + "1_tmp/{y:04}/{d:03}/{s:4}{d:03}0.{y2:02}n".format(y=year4,y2=year2,s=sta,d=day)
        if not os.path.exists(path_n):
            continue
        idsat,tiempo,oel,i1stsat=load_navigation_file(path_n,"3.02")


        #
        # reading observation file
        #

        path_o = path_tec + "1_tmp/{y:04}/{d:03}/{s:4}{d:03}0.{y2:02}o".format(y=year4, y2=year2, s=sta, d=day)
        if not os.path.exists(path_o):
            continue
        G_stec,G_stec_bool,G_sat_bool,p_o=load_observation_file(path_o,"3.02")

        #
        # calculating satellite positions and outputting data
        #

        for isat in range(maxsat):
            if G_sat_bool[isat]:
                path_nw = "C:/data/{y:04}/{d:03}/G{sat:02}".format(y=year4,d=day,sat=isat)
                os.makedirs(path_nw,exist_ok=True)
                path_nd = path_nw + "/{s:4}.dat".format(s=sta)
                make_tomo_loadfile(path_nd,G_stec,G_stec_bool,G_sat_bool,i1stsat,p_o,isat)


        # print("station:{s:4}, day:{d:03} is finished. {t:.3f} {i}/{f}".format(s=sta,d=day,t=time.time()-start,i=idx,f=len(stas)))
