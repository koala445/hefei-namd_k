#!/usr/bin/python 
# -*- coding: utf-8 -*-

import os
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from mpl_toolkits import mplot3d 
import math
from scipy.optimize import curve_fit
from scipy import special

def main():
    plt.rcParams['font.family'] = 'Arial'  #'Times New Roman'
    #E = 1   #    initial energy elec: + ;hole: -  #0.6, 0.8, 1, 1.2, 1.4, 1.6,
    #delta_E = 0.00001 
    delta_E = 0.002
    Efermi = np.array([-4.528625, -5.2722035, -5.3783725, -5.416445 , -5.4291865  ,1 ,1 , # 0,1,2,3,4,5,6 
                       -4.566445 ,-5.414105, # 7 ,8
                       -5.4413085, -5.4174325, 1,  # 9 10 11
                       -5.352921325, -5.26030079,          #12,13
                       -5.481080755, -5.48377689, -5.48120359,-5.478049145   #14,15,16,17 
                       ])
    
    #Efermi =  -4.528625     # GR: 0 gy: 1 gdy:2 g3y:3 gty:4 g5y:5 g6y:6 ;;;
                             #gr-ABA :7; gdy-ABC :8
                              #a-gy:9 b-gy:10 12-gy:11
                             #4-hR3-gy :12   4-HL3-gy:13
                             # a-gdy:14,  a-g3y:15, a-g4y:16, a-g5y:17
    
    Efermi = Efermi[14]
  
   ## get_En_each_phonon()                  #test
    ##get_Couple_each_phonon ('EPPHTXT')    #test
    ##get_electronic_state_occupy_time()    # test
    
    #get_couple_Exyz('EPECTXT','EIGTXT','BASSEL')  #get complete state-state coupling data and sum to oneside #old version before 2024
    #get_couple_Exyz_new('EPELTXT','BASSEL')  #get complete state-state coupling data and sum to oneside
    #get_couple_E('EXYZ_coup.txt', E, delta_E)
    
    #plot_couple_E_imshow('EPELTXT', E, delta_E, figname='COUPLE_E_imshow.png') #'EPECTXT_E'
    ##plot_couple_tricontour('EPECTXT_E_3D.txt', Efermi,figname='COUPLE_E_tricontour.png')
    #plot_couple_tricontour('EXYZ_coup.txt', Efermi)   #as well
    
    #plot_couple_oneside_sum('coup_half-oneside_sum.txt','coup_oneside_sum.txt',Efermi)#'coup_diff_sum.txt','coup_same_sum.txt',Efermi, )
    #plot_couple_oneside_all('EXYZ_coup.txt','EXYZ_coup_diff.txt','EXYZ_coup_same.txt',Efermi)
    
    #plot_couple_contour('EPECTXT','EIGTXT', Efermi)  # well done! #old version before 2024
    #plot_couple_contour_new('EPELTXT','BASSEL', Efermi)  # well done!
    
    
    #plot_density_state('EIGTXT', Efermi)
    #get_occupation_band_index_time()      #useless
    
    ###########################  En-Time and lifetime or energy loss - initial level
    
    #watch_disturibution('SHPROP.40',Efermi)
    #plot_3D_Energy_Time('Energy-Time.txt')
    #fit_En_Time('Energy.txt',Efermi,figname='fit_Energy_Time.png')
    fit_En_Time_1_e('Energy.txt',Efermi,figname='fit_Energy_Time_1_e.png')
    #fit_OD_Time('occupy_OD.txt','Energy.txt',figname='fit_OD_Time.png')
    #fit_OD_Time_1_e('occupy_OD.txt','Energy.txt',figname='fit_OD_Time_1_e.png')
    #get_lifetime_loss('Energy-Time.txt')407
    #get_occu_elec_time('occu-elec-time.txt')
    
    
def get_couple_Exyz(file1,file2,file3):
    from collections import Counter  
# get EX EY E_coup data and index_xy
    EX=np.round(np.array(np.loadtxt(file2)),decimals=6)
    EY=np.round(np.array(np.loadtxt(file2)),decimals=6)
    EZ=np.round(np.array(np.loadtxt(file1)),decimals=9)
    index_x = np.loadtxt(file3, usecols=(1,),skiprows=1)
    index_y = np.loadtxt(file3, usecols=(1,),skiprows=1)
    
# save EXYZ_coup
    nb=int(np.sqrt(EZ.shape[0]))
    EXYZ_coup=np.empty(shape=(nb*nb,5))
    
    nd = Counter(index_x)[9]*Counter(index_x)[8] # n_diff = 2*nd
    EXYZ_coup_diff = np.empty(shape=(2*nd,3))
    EXYZ_coup_same = np.empty(shape=(nb*nb-2*nd,3))
    coup_diff_sum = np.zeros(shape=(nb,2))
    coup_same_sum = np.zeros(shape=(nb,2))
    
    n_all=0
    n_diff = 0
    n_same = 0
    for x in range(0,nb):
        coup_diff_sum[x,0] = EX[x]
        coup_same_sum[x,0] = EX[x]
        for y in range(0,nb):
            EXYZ_coup[n_all,0]=EX[x]
            EXYZ_coup[n_all,1]=EY[y]
            EXYZ_coup[n_all,2]=EZ[n_all]
            EXYZ_coup[n_all,3]=index_x[x]
            EXYZ_coup[n_all,4]=index_y[y]
            '''
            if  index_x[x] != index_y[y]: #judge diff band couple
                coup_diff_sum[x,1] += EZ[n_all]
                EXYZ_coup_diff[n_diff,0]=EX[x]
                EXYZ_coup_diff[n_diff,1]=EY[y]
                EXYZ_coup_diff[n_diff,2]=EZ[n_all]
                n_diff += 1
            else:
                coup_same_sum[x,1] += EZ[n_all]
                EXYZ_coup_same[n_same,0]=EX[x]
                EXYZ_coup_same[n_same,1]=EY[y]
                EXYZ_coup_same[n_same,2]=EZ[n_all]
                n_same += 1 
            '''
            n_all=n_all+1
    
    
    
    np.savetxt('EXYZ_coup.txt',EXYZ_coup,fmt='%1.6f %1.6f %1.9f %d %d ')
    #np.savetxt('EXYZ_coup_diff.txt',EXYZ_coup_diff,fmt='%1.6f %1.6f %1.9f')
    #np.savetxt('EXYZ_coup_same.txt',EXYZ_coup_same,fmt='%1.6f %1.6f %1.9f')
# save coup_oneside_sum.txt and coup_diff_sum.txt
    coup_each_sum = np.empty(shape=(nb,2))
    coup_half_each_sum = np.empty(shape=(nb,2))  # just upper or lower triangle coupling
    for i in range(0,nb):
        coup_each_sum[i,0] = EX[i]
        coup_each_sum[i,1] = np.sum(EZ[i*nb:(i+1)*nb])  
    for i in range(0,nb):
        coup_half_each_sum[i,0] = EX[i]
        each_coup = EZ[i*nb:(i+1)*nb]
        coup_half_each_sum[i,1] = np.sum(each_coup[0:i+1])
    
    np.savetxt('coup_oneside_sum.txt',coup_each_sum,fmt='%1.6f %1.9f')
    np.savetxt('coup_half-oneside_sum.txt',coup_half_each_sum,fmt='%1.6f %1.9f')
    #np.savetxt('coup_diff_sum.txt',coup_diff_sum,fmt='%1.6f %1.9f')
    #np.savetxt('coup_same_sum.txt',coup_same_sum,fmt='%1.6f %1.9f')
    #coup_diff = np.empty(shape=())
    #for j in range(0,nb*nb):
    #    if EXYZ_coup[i,3] != EXYZ_coup[i,4]:
    
    #print(i)
    #num_diff  = len(EXYZ_coup(EXYZ_coup[3] != EXYZ_coup[4]))
    #print(EXYZ_coup(EXYZ_coup[i][3] != EXYZ_coup[i][4]))
    #coup_diff = np.empty(shape=(num_diff,2))
    
# save coup_all and coup_everystate
    coup_all = np.sum(np.abs(EZ)) * 1000  #change unit to mev
    coup_everystate = np.average(np.abs(coup_each_sum[:,1])) * 1000  #change unit to mev 
    coup_max = np.max(np.abs(EZ)) * 1000 #change unit to mev
    with open('coup_analy.txt','w',encoding='utf-8') as f:
        f.write(str('coup_all :')+' ')
        f.write(str(coup_all)+'\n')
        f.write(str('coup_everystate :')+' ')
        f.write(str(coup_everystate)+'\n')
        f.write(str('coup_max :')+' ')
        f.write(str(coup_max)+'\n')
        f.close()
    print(len(EX),nb)

def get_couple_Exyz_new(file1,file2):
    from collections import Counter  
# get EX EY E_coup data and index_xy
    EX=np.round(np.array(np.loadtxt(file2, usecols=2,skiprows=1,dtype='float32' )),decimals=6) 
    EY=np.round(np.array(np.loadtxt(file2, usecols=2,skiprows=1,dtype='float32' )),decimals=6) 
    EZ=np.round(np.array(np.loadtxt(file1).flatten(),dtype='float32'),decimals=9)  


    index_x = np.loadtxt(file2, usecols=(1,),skiprows=1,dtype='float32') #index_which band
    index_y = np.loadtxt(file2, usecols=(1,),skiprows=1,dtype='float32')
    
# save EXYZ_coup
    nb=int(len(EX))
    print(nb,EZ)
    EXYZ_coup=np.empty(shape=(nb*nb,5))
    
    nd = Counter(index_x)[9]*Counter(index_x)[8] # n_diff = 2*nd
    EXYZ_coup_diff = np.empty(shape=(2*nd,3))
    EXYZ_coup_same = np.empty(shape=(nb*nb-2*nd,3))
    coup_diff_sum = np.zeros(shape=(nb,2))
    coup_same_sum = np.zeros(shape=(nb,2))
    
    n_all=0
    n_diff = 0
    n_same = 0
    for x in range(0,nb):
        coup_diff_sum[x,0] = EX[x]
        coup_same_sum[x,0] = EX[x]
        for y in range(0,nb):
            EXYZ_coup[n_all,0]=EX[x]
            EXYZ_coup[n_all,1]=EY[y]
            EXYZ_coup[n_all,2]=EZ[n_all]
            EXYZ_coup[n_all,3]=index_x[x]
            EXYZ_coup[n_all,4]=index_y[y]
            '''
            if  index_x[x] != index_y[y]: #judge diff band couple
                coup_diff_sum[x,1] += EZ[n_all]
                EXYZ_coup_diff[n_diff,0]=EX[x]
                EXYZ_coup_diff[n_diff,1]=EY[y]
                EXYZ_coup_diff[n_diff,2]=EZ[n_all]
                n_diff += 1
            else:
                coup_same_sum[x,1] += EZ[n_all]
                EXYZ_coup_same[n_same,0]=EX[x]
                EXYZ_coup_same[n_same,1]=EY[y]
                EXYZ_coup_same[n_same,2]=EZ[n_all]
                n_same += 1 
            '''
            n_all=n_all+1
    
    
    
    np.savetxt('EXYZ_coup.txt',EXYZ_coup,fmt='%1.6f %1.6f %1.9f %d %d ')
    #np.savetxt('EXYZ_coup_diff.txt',EXYZ_coup_diff,fmt='%1.6f %1.6f %1.9f')
    #np.savetxt('EXYZ_coup_same.txt',EXYZ_coup_same,fmt='%1.6f %1.6f %1.9f')
# save coup_oneside_sum.txt and coup_diff_sum.txt
    coup_each_sum = np.empty(shape=(nb,2))
    coup_half_each_sum = np.empty(shape=(nb,2))  # just upper or lower triangle coupling
    for i in range(0,nb):
        coup_each_sum[i,0] = EX[i]
        coup_each_sum[i,1] = np.sum(EZ[i*nb:(i+1)*nb])  
    for i in range(0,nb):
        coup_half_each_sum[i,0] = EX[i]
        each_coup = EZ[i*nb:(i+1)*nb]
        coup_half_each_sum[i,1] = np.sum(each_coup[0:i+1])
    
    np.savetxt('coup_oneside_sum.txt',coup_each_sum,fmt='%1.6f %1.9f')
    np.savetxt('coup_half-oneside_sum.txt',coup_half_each_sum,fmt='%1.6f %1.9f')
    #np.savetxt('coup_diff_sum.txt',coup_diff_sum,fmt='%1.6f %1.9f')
    #np.savetxt('coup_same_sum.txt',coup_same_sum,fmt='%1.6f %1.9f')
    #coup_diff = np.empty(shape=())
    #for j in range(0,nb*nb):
    #    if EXYZ_coup[i,3] != EXYZ_coup[i,4]:
    
    #print(i)
    #num_diff  = len(EXYZ_coup(EXYZ_coup[3] != EXYZ_coup[4]))
    #print(EXYZ_coup(EXYZ_coup[i][3] != EXYZ_coup[i][4]))
    #coup_diff = np.empty(shape=(num_diff,2))
    
# save coup_all and coup_everystate
    coup_all = np.sum(np.abs(EZ)) * 1000  #change unit to mev
    coup_everystate = np.average(np.abs(coup_each_sum[:,1])) * 1000  #change unit to mev 
    coup_max = np.max(np.abs(EZ)) * 1000 #change unit to mev
    with open('coup_analy.txt','w',encoding='utf-8') as f:
        f.write(str('coup_all :')+' ')
        f.write(str(coup_all)+'\n')
        f.write(str('coup_everystate :')+' ')
        f.write(str(coup_everystate)+'\n')
        f.write(str('coup_max :')+' ')
        f.write(str(coup_max)+'\n')
        f.close()
    print(len(EX),nb)
        

def get_couple_E(arr_name,E,delta_E):
# prepare  DataFrame
    
    EXYZ_coup = np.loadtxt(arr_name)
    Enn_coup = EXYZ_coup[:,0:2]
    EZ=  EXYZ_coup[:,2] 
    nb=int(np.sqrt(EZ.shape[0]))
    EZ_coup = EZ.reshape(nb*nb,1)
    #E = 1.5 #eV
    #delta_E = 0.05 #eV
    n = int(E/delta_E)
    Emin = EXYZ_coup[0,0] 
    
    for nk in range(0,n):
        Enk = nk*delta_E + Emin
        Enk1 = (nk+1)*delta_E + Emin
        E_average=(Enk+Enk1)/2
        Enn_coup[(Enn_coup >= Enk) & (Enn_coup < Enk1) ]=E_average
    EnnZ_coup=np.append(Enn_coup,EZ_coup,axis=1)   
    np.savetxt('EnnZ_coup',EnnZ_coup,fmt='%1.9f %1.9f %1.9f' )
    
#all kinds of dataframe
#2D data
    df=pd.DataFrame(EnnZ_coup,columns=['EX','EY','EZ'])
    #df.round({'EX': 6,'EY': 6, 'EZ':9})
    E_all = df.groupby(['EX','EY']).sum()
    E_mean = df.groupby(['EX','EY']).mean()
    df_E_mean = pd.DataFrame(E_mean).reset_index()
    E_max = df.groupby(['EX','EY']).max()

#arr_data
    arr_Emean=np.array(E_mean)
    arr_Emax=np.array(E_max)
    arr_Eall=np.array(E_all)
    arr_Emean=np.array(E_mean)
#save data\
    np.savetxt('EPECTXT_E',arr_Emean , fmt='%.10f')
    df_E_mean.to_csv('EPECTXT_E_3D.txt',sep='\t',index=False,float_format='%.9f')
    print( "E_mean-max:%f" % arr_Emean.max(),"E_mean-sum:%f" % arr_Emean.sum(),
          "E_max:%f" % arr_Emax.max(), "E_sum:%f" % EZ.sum(),
          )
    
    print("E_all-max:%f" % arr_Eall.max(),"E_all-sum:%f" % arr_Eall.sum(),
          )
   
def get_En_each_phonon ():
# select PHPROP.*
    filepath='./'
    filename_list = os.listdir(filepath)
    filelist=[]
    for filename in filename_list:
        if filename[0:6]=='PHPROP':
            filelist.append(filename)
    count = len(filelist)
    phonon_En = np.empty(shape=(count,2))
#  each phonon get how much energy from electron    
    regex = re.compile(r'\d+')        
    for file in filelist:
        n = [int(x) for x in regex.findall(file)][0]
        print(n)
        phonon_En[(n-1),0]= n
        phonon_En[(n-1),1]= np.sum(np.loadtxt(file,skiprows=24,usecols=(1,),dtype=float,))
       
    print(np.sum(phonon_En[:,1]))            
    np.savetxt('phonon_En.txt',phonon_En)
    print("\n phonon_En.txt has been saved.")

def get_Couple_each_phonon (file):
    
    couple_phonon = np.sum(np.loadtxt(file,dtype=float,),axis=0)
    index = np.arange(len(couple_phonon))+1
    Couple_phonon = np.vstack((index,couple_phonon)).T
    np.savetxt('Couple_phonon.txt',Couple_phonon)
    print("\n Couple_phonon.txt has been saved.")
    
    print(Couple_phonon)
    
def get_electronic_state_occupy_time():
    
    nu= np.loadtxt('INICON',usecols=(0,),dtype=int)
    print(nu,nu[0])
    SHPROP_first = np.loadtxt('SHPROP.'+str(nu), skiprows=24)
    time = SHPROP_first[:,0]
    
    for i in nu:
        SHPROP_ini = np.loadtxt('SHPROP.'+str(i), skiprows=24)
        index = np.where(SHPROP_first[0,:]==1)[0][1]     # [0][0] is initial time
        occupation_each = SHPROP_ini[:,index]
        if i == nu[0]:
            occupation_every = occupation_each
        else:
            occupation_every = np.vstack((occupation_every,occupation_each))
    occupation_mean = np.mean(occupation_every,axis=0)
    elec_state = np.vstack((time,occupation_mean)).T
    '''
    nu = int(np.loadtxt('INICON')[0,0])
    print(nu)
    SHPROP_first = np.loadtxt('SHPROP.'+str(nu), skiprows=24) 
    index = np.where(SHPROP_first[0,:]==1)[0][1]     # [0][0] is initial time 
    print(index)
    elec_state = SHPROP_first[:,(0,index)]
    '''
    np.savetxt('occupy_electronic_state.txt',elec_state)
    print('occupy_electronic_state.txt has been saved')
    
def get_occupation_band_index_time():
    
    nu_list = np.loadtxt('INICON')[:,0].astype(int)
    index_list = np.loadtxt('BASSEL', usecols=(1,),skiprows=1)
    index_8 = np.where(index_list==8)[0]
    index_9 = np.where(index_list==9)[0]
    time = len(np.loadtxt('SHPROP.' + str(nu_list[0])   , skiprows=24))
    
    all_occup_time = np.empty(shape=(len(nu_list),time,3))
    print(nu_list,index_list,index_8,time)
    for nu in nu_list:
        
        SHPROP_occupation = np.loadtxt('SHPROP.' + str(nu) , skiprows=24)[:,2:]
        #np.multiply(SHPROP_occupation,) 
        n=np.where(nu_list == nu)
       
        for i in range(len(SHPROP_occupation)):
            all_occup_time[n,i,0] = i+1
            all_occup_time[n,i,1] = np.sum(SHPROP_occupation[i,index_8])
            all_occup_time[n,i,2] = np.sum(SHPROP_occupation[i,index_9])
            
        
    aver_occup_time=np.average(all_occup_time,axis=0)
    print(aver_occup_time)
    np.savetxt('occupy_band_time.txt',aver_occup_time)
    
    
def  plot_density_state(file, Ef):
     
#get freq and scope to plot hist  qq
    data = np.round(np.loadtxt(file),decimals=6) - Ef
    num_bins= int((max(data)-min(data))/0.05)
    
    frenq,scope,patches= plt.hist(data, bins=num_bins, color=(0.378, 0.574, 0.781) ,rwidth=1)
#get frenq and mid-scope to plot line
    
    mid_scope = np.empty(shape=(len(scope)-1))      
    for i in range(len(mid_scope)):
        mid_value = (scope[i]+scope[i+1])/2
        mid_scope[i] = mid_value
#exclude 0 freq data and plot  
    index_nonzero = np.nonzero(frenq)
    plt.plot(mid_scope[index_nonzero],frenq[index_nonzero],color='red',linestyle='--')
    plt.grid(alpha=0.5,linestyle='-.')
    plt.xlabel('E-Ef (eV)',family='Times New Roman',fontsize=16,weight='bold') #fontweight='semibold')
    plt.ylabel('Frenquency of states',family='Times New Roman',fontsize=16,weight='bold') #fontweight='semibold')
    plt.xticks(fontsize=16,family='Times New Roman',weight='bold')
    plt.yticks(fontsize=16,family='Times New Roman',weight='bold')
    plt.tight_layout()
    plt.savefig('density_state.png', dpi=300)
    plt.show()

def plot_couple_E_imshow(file, E, delta_E, figname='COUPLE_E_imshow.png'):
    '''
    This function plots average couplings.

    Parameters:
    coup: ndarray, average coupling data in forms of coup[nb, nb]
    figname: string, output figure file name.
    '''
    
    coup = np.loadtxt(file)
    nb=int(np.sqrt(coup.shape[0]))
    coup.resize(nb,nb)
    
    fig = plt.figure()
    ax=plt.axes()
    figsize_x = 4.8
    figsize_y = 3.6 # in inches
    fig.set_size_inches(figsize_x, figsize_y)
# colormap
    cmap0 = plt.cm.jet
    cmap1 = mpl.colors.ListedColormap(['white'])
    color0=cmap0(np.linspace(0,1,400))
    color0_1=cmap0(np.linspace(0,1,80))
    color1=cmap1(np.linspace(0,1,10))
    newcolorsgroup=np.vstack((color0_1[15:29],color0_1[30:40],color0_1[41:48],color0[250:380]))
    newcmap=mpl.colors.ListedColormap(newcolorsgroup)
    #newcmap.set_under('white')
# plot imshow
    n = coup.shape[0]
    coup *= 1000.0 # change unit to meV
    Bmin = delta_E/2 ; Bmax = E+delta_E/2
    #cmin = 0.0; cmax = np.max(coup)
    cmin=round(np.min(coup),4);  cmax=round(np.max(coup),4)
    norm = mpl.colors.Normalize(cmin,cmax)
    plt.imshow(coup, cmap=newcmap, origin='lower', norm=norm,
        extent=(Bmin,Bmax,Bmin,Bmax), 
        interpolation='none')
    
    bounds= np.arange(cmin,cmax*1.01,(cmax-cmin)/5)
    cb=plt.colorbar(
        mpl.cm.ScalarMappable(cmap=newcmap,norm=norm),
        ax=ax,
        #extend='min',
        ticks=bounds,
        spacing='proportional',
        orientation='vertical',
        label='Coupling (meV)',
    )
    '''
    cb.ax.tick_params(axis='both',labelsize=10)
    plt.yticks(fontproperties='Times New Roman', fontsize=10,weight='bold') 
    plt.xticks(fontproperties='Times New Roman', fontsize=10,weight='bold')
    '''
    
    #cbar = plt.colorbar()
    # cbar.ax.set_title('   meV')
    #cbar.set_label('Coupling (meV)')
    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved."%figname)

def plot_couple_tricontour(file, E0,figname='COUPLE_E_tricontour.png'):

    x0=np.round(np.loadtxt(file,usecols=(0,)),decimals=6) - E0
    y0=np.round(np.loadtxt(file,usecols=(1,)),decimals=6) - E0
    z0=np.round(np.loadtxt(file,usecols=(2,)),decimals=9) * 1000 #change unit to mev
    
    
    fig = plt.figure()
    ax=plt.axes()
    figsize_x = 4.8
    figsize_y = 3.6 # in inches
    fig.set_size_inches(figsize_x, figsize_y)
    
# colormap
    cmap0 = plt.cm.jet
    cmap1 = mpl.colors.ListedColormap(['white'])
    color0=cmap0(np.linspace(0,1,400))
    color0_1=cmap0(np.linspace(0,1,80))
    color1=cmap1(np.linspace(0,1,10))
    newcolorsgroup=np.vstack((color0_1[15:29],color0_1[30:40],color0_1[41:48],color0[250:380]))
    newcmap=mpl.colors.ListedColormap(newcolorsgroup)
    
    cmin=round(np.min(z0),4); cmax=round(np.max(z0),4)
    
    norm = mpl.colors.Normalize(cmin,cmax)
    #norm=mpl.colors.LogNorm(cmin,cmax)
    bounds= np.arange(cmin,cmax*1.01,(cmax-cmin)/5)
    
    
    plt.tricontourf(x0,y0,z0,30,cmap=newcmap,norm=norm)
    #plt.scatter(x1, y1, s=1, color= 'grey' ,alpha=0.2)
    
    cb=plt.colorbar(
            mpl.cm.ScalarMappable(cmap=newcmap,norm=norm),
            ax=ax,
            #extend='min',
            ticks=bounds,
            spacing='proportional',
            orientation='vertical',
            label='Coupling (meV)',
        )
    
    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved."%figname)
    plt.show()

def plot_couple_contour(file1,file2, E0, figname='COUPLE_E_contour_arial.png'):
    x0=np.round(np.array(np.loadtxt(file2)),decimals=6) - E0
    y0=np.round(np.array(np.loadtxt(file2)),decimals=6) - E0
    z0=np.round(np.array(np.loadtxt(file1)),decimals=9) * 1000 #change unit to mev
   
    print(x0.shape)
    fig = plt.figure()
    plt.rc('font',family='arial',size=11)
    ax=plt.axes()
    ax.spines['bottom'].set_linewidth(1)                 #set frame
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['top'].set_linewidth(1)
    figsize_x = 4.8
    figsize_y = 3.6 # in inches
    fig.set_size_inches(figsize_x, figsize_y)
    
# colormap
    #cmap0 = plt.cm.jet
    #color0=cmap0(np.linspace(0,1,200))
    #color0_1=cmap0(np.linspace(0,1,80))
    #newcolorsgroup=np.vstack((color0_1[18:49],color0[125:180]))
    #newcmap=mpl.colors.ListedColormap(newcolorsgroup)    
    #newcmap = plt.cm.RdYlBu_r

######### map1_1
    '''   
    cmap0 = plt.cm.RdBu_r#RdYlBu_r##jet
    color0=cmap0(np.linspace(0,1,200))
    color0_1=cmap0(np.linspace(0,1,80))
    newcolorsgroup=np.vstack((color0_1[12:39],color0[101:175]))
    newcmap=mpl.colors.ListedColormap(newcolorsgroup)  
    '''
########### map1_2
    '''
    cmap0 = plt.cm.RdBu_r#RdYlBu_r##jet
    color0=cmap0(np.linspace(0,1,200))

    newcolorsgroup= color0  #color0[101:200]
    newcmap=mpl.colors.ListedColormap(newcolorsgroup)  
    '''
########### map1_2

    colorlist=[ (157/255,92/255,57/255),                         #black
        (254/255,147/255,118/255),      #red (242/255,106/255,17/255),
                (175/255,221/255,139/255),  # green (53/255,161/255,83/255)
               (160/255,201/255,229/255),(66/255,146/255,201/255)  ][::-1]
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("testCmap", colors=colorlist, N=512)

###########
    cmin=round(np.min(z0),2);  cmax=0.37 #cmax=0.54 for GR; 0.51 for GDY    #;
    norm = mpl.colors.Normalize(0,cmax)
       
# plot contourf    
    X,Y = np.meshgrid(x0,y0)
    print(X.shape, cmax,round(np.max(z0),2) )
    nb = len(x0)
    Z = np.resize(z0,(nb,nb))
    plt.contourf(X,Y,Z,30,cmap=newcmap,norm=norm)
    
    bounds= np.round(np.arange(0,cmax*1.0001,cmax/4),decimals=2) #bounds= np.linspace(cmin,cmax,5)
    #bounds = np.array(['0', '0.13','0.26','0.39', '0.52'])
    ##########colorbar    
    cb=plt.colorbar(
            mpl.cm.ScalarMappable(cmap=newcmap,norm=norm),
            ax=ax,
            #extend='min',
            ticks=bounds,
            spacing='proportional',
            orientation='vertical',
            label='Coupling (meV)',
        )
    font = {'family' : 'arial',#'Times New Roman',
        'color'  : 'black',
        'weight' : 'bold',
        'size'   : 14,
        }
    cb.ax.tick_params(labelsize=12)
    for i in range(len(bounds)):
        cb.ax.get_yticklabels()[i].set_fontweight("bold")
    #bounds = np.array(['0', '0.13','0.26','0.39', '0.52'])
    cb.ax.set_yticklabels(bounds )
    cb.set_label('Coupling (meV)',fontdict=font, labelpad = 6 ) #设置colorbar的标签字体及其大小
######### ticks labels   
    plt.xlabel('$\mathregular{E_{nk}}$ (eV)',family= 'arial', #'Times New Roman', 
               labelpad=2 , fontsize = 14, fontweight='bold')
    plt.ylabel('E$_{mk\'}$ (eV)',family= 'arial', #'Times New Roman',
               labelpad= 2, fontsize = 14, fontweight='bold')
    
    my_x = np.arange(0, 1.51, 0.25)      #for gdy    ##initial level     
    #my_x = np.linspace(-0.3, 0.3, 5)                  
    #my_xticks = np.array(['0.25', '0.5','0.75','1.0','1.25'])
    
    #my_x = np.arange(0.25, 1.5, 0.3)      #for gr    ##initial level                       
    #my_xticks = np.array(['0.25', '0.55','0.85','1.15','1.45'])
    my_xticks = np.array(['0','0.25', '0.5','0.75','1','1.25','1.5'])
    my_yticks = np.array(['0.00','0.25', '0.50','0.75','1.00','1.25','1.50'])
    xmin=np.min(x0); xmax=np.max(x0)
    ax.set_xlim(0,1.5)
    ax.set_ylim(-0.05,1.55)
    ax.xaxis.set_ticks(my_x)
    ax.xaxis.set_tick_params(pad=2.7, labelsize=13 )
    ax.yaxis.set_ticks(my_x)
    ax.yaxis.set_tick_params(pad=2, labelsize=13)
    ax.set_xticklabels(my_xticks, )
    ax.set_yticklabels(my_yticks, )
    for i in range(len(my_x)):
        ax.get_xticklabels()[i].set_fontweight("bold")
        ax.get_yticklabels()[i].set_fontweight("bold")
  
    
#plot line

    n1 = 0; n2 = 0
    k1 = 0.19 ;k2 = 0.27  # gr:0.195 #gdy 0.19 0.27 #a-gdy 0 0.05 #b-gdy 0.16 0.27 #gy  #g3y 
    #k1 = 0.19;k2 = 0.19 #gr
                             #a-gdy
    #k1= 0.16  ; k2 = 0.27       #b-gdy
    '''
    ####### make line inside data picture
    for i in x0:
        if i < np.max(x0)-k1:
            n1=n1+1
    for i in x0:
        if i < np.max(x0)-k2:
            n2=n2+1
    x1 = x0[0:n1] ; x2 = x0[0:n2]
    '''
    x1 = np.linspace(0,1.5,100); x2 = np.linspace(0,1.5,100)  ; #=my_x
    y1 = k1 + x1 ; y2 = k2 + x2
    plt.plot(x1,y1,color='black',linewidth=1.5, linestyle='--',dashes=(3,2),alpha=0.8)   
    plt.plot(x2,y2,color='black',linewidth=1.5, linestyle='--',dashes=(3,2),alpha=0.8)

######## plot diff_band points


#save picture  
    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved."%figname)
    plt.show()

def plot_couple_contour_new(file1,file2, E0, figname='COUPLE_E_contour_arial.png'):
    #x0=np.round(np.array(np.loadtxt(file2, usecols=2,skiprows=1 ),dtype=np.float32),decimals=6) - E0
    #y0=np.round(np.array(np.loadtxt(file2, usecols=2,skiprows=1 ),dtype=np.float32),decimals=6) - E0
    #z0=np.round(np.array(np.loadtxt(file1),dtype=np.float32),decimals=9) * 1000 #change unit to mev
    
    x0=np.array(np.loadtxt(file2, usecols=2,skiprows=1 ),dtype=np.float32) - E0
    y0=np.array(np.loadtxt(file2, usecols=2,skiprows=1 ),dtype=np.float32) - E0
    z0=np.array(np.loadtxt(file1),dtype=np.float32) * 1000 #change unit to mev
    
    print(x0.shape,z0.shape)
    fig = plt.figure()
    
    ax=plt.axes()
    ax.spines['bottom'].set_linewidth(1)                 #set frame
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['top'].set_linewidth(1)
    figsize_x = 4.8
    figsize_y = 3.6 # in inches
    fig.set_size_inches(figsize_x, figsize_y)
    
# colormap
    '''
    cmap0 = plt.cm.jet
    color0=cmap0(np.linspace(0,1,200))
    color0_1=cmap0(np.linspace(0,1,80))
    newcolorsgroup=np.vstack((color0_1[18:49],color0[125:180]))
    newcmap=mpl.colors.ListedColormap(newcolorsgroup)    
    '''
    colorlist=[ (157/255,92/255,57/255),                         #black
        (254/255,147/255,118/255),      #red (242/255,106/255,17/255),
                (175/255,221/255,139/255),  # green (53/255,161/255,83/255)
               (160/255,201/255,229/255),(66/255,146/255,201/255)  ][::-1]
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("testCmap", colors=colorlist, N=512)
    cmin=round(np.min(z0),4);  cmax= 0.03 #0.94   #round(np.max(z0),2) #cmax=0.54 #;
    norm = mpl.colors.Normalize(0,cmax)
       
# plot contourf    
    X,Y = np.meshgrid(x0,y0)
    print(X.shape, 'cmax:',round(np.max(z0),6) )
    nb = len(x0)
    #Z = np.resize(z0,(nb,nb))
    plt.contourf(X,Y,z0,30,cmap=newcmap,norm=norm)
    
    bounds=  np.round(np.arange(0,cmax*1.0001,cmax/4),decimals=2) #bounds= np.linspace(cmin,cmax,5)
    #bounds = np.array(['0', '0.13','0.26','0.39', '0.52'])
    ##########colorbar    
    cb=plt.colorbar(
            mpl.cm.ScalarMappable(cmap=newcmap,norm=norm),
            ax=ax,
            #extend='min',
            ticks=bounds,
            spacing='proportional',
            orientation='vertical',
            label='Coupling (meV)',
        )
    font = {
        'color'  : 'black',
        'weight' : 'bold',
        'size'   : 14,
        }
    cb.ax.tick_params(labelsize=12)
    for i in range(len(bounds)):
        cb.ax.get_yticklabels()[i].set_fontweight("bold")
    #bounds = np.array(['0', '0.13','0.26','0.39', '0.52'])
    cb.ax.set_yticklabels(bounds )
    cb.set_label('Coupling (meV)',fontdict=font, labelpad = 6 ) #设置colorbar的标签字体及其大小
######### ticks labels   
    plt.xlabel('$\mathregular{E_{nk}}$ (eV)', labelpad=2 , fontsize = 14, fontweight='bold')
    plt.ylabel('E$_{mk\'}$ (eV)', labelpad= 2, fontsize = 14, fontweight='bold')
    
    my_x = np.arange(0, 1.51, 0.25)      #for gdy    ##initial level     
    #my_x = np.linspace(-0.3, 0.3, 5)                  
    #my_xticks = np.array(['0.25', '0.5','0.75','1.0','1.25'])
    
    #my_x = np.arange(0.25, 1.5, 0.3)      #for gr    ##initial level                       
    #my_xticks = np.array(['0.25', '0.55','0.85','1.15','1.45'])
    my_xticks = np.array(['0','0.25', '0.5','0.75','1','1.25','1.5'])
    my_yticks = np.array(['0.00','0.25', '0.50','0.75','1.00','1.25','1.50']) 
    
    xmin=np.min(x0); xmax=np.max(x0)
    ax.set_xlim(0,1.5)
    ax.set_ylim(-0.05,1.55)
    ax.xaxis.set_ticks(my_x)
    ax.xaxis.set_tick_params(pad=2.7, labelsize=13 )
    ax.yaxis.set_ticks(my_x)
    ax.yaxis.set_tick_params(pad=2, labelsize=13)
    ax.set_xticklabels(my_xticks, )
    ax.set_yticklabels(my_yticks, )
    for i in range(len(my_x)):
        ax.get_xticklabels()[i].set_fontweight("bold")
        ax.get_yticklabels()[i].set_fontweight("bold")
  
    
#plot line

    n1 = 0; n2 = 0
    k1 = 0.14 ;k2 = 0.19 ; k3=0.27  # gr:0.195 #gdy 0.19 0.27 #a-gdy 0 0.05 #b-gdy 0.16 0.27 #gy  #g3y 
    #k1 = 0.19;k2 = 0.19 #gr
                             #a-gdy
    #k1= 0.16  ; k2 = 0.27       #b-gdy
    '''
    ####### make line inside data picture
    for i in x0:
        if i < np.max(x0)-k1:
            n1=n1+1
    for i in x0:
        if i < np.max(x0)-k2:
            n2=n2+1
    x1 = x0[0:n1] ; x2 = x0[0:n2]
    '''
    x1 = np.linspace(0,1.5,100); x2 = np.linspace(0,1.5,100)  ; x3 = np.linspace(0,1.5,100)#=my_x
    y1 = k1 + x1 ; y2 = k2 + x2 ; y3= k3+x3
    plt.plot(x1,y1,color='black',linewidth=1.5, linestyle='--',dashes=(3,2),alpha=0.8)   
    plt.plot(x2,y2,color='black',linewidth=1.5, linestyle='--',dashes=(3,2),alpha=0.8)
    #plt.plot(x3,y3,color='black',linewidth=1.5, linestyle='--',dashes=(3,2),alpha=0.8)

######## plot diff_band points


#save picture  
    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved."%figname)
    plt.show()

def plot_couple_oneside_sum(file1, file2, Ef):#file2,file3, Ef):
    fig = plt.figure()
    ax=plt.axes()
    figsize_x = 4.8
    figsize_y = 3.6 # in inches
    fig.set_size_inches(figsize_x, figsize_y)
    
    x0=np.loadtxt(file1,usecols=(0,),skiprows=1) -Ef #E-Efermi
    z0=np.loadtxt(file1,usecols=(1,),skiprows=1)*1000 # change unit to mev
    #z1=np.loadtxt(file2,usecols=(1,),skiprows=1)*1000
    #z2=np.loadtxt(file3,usecols=(1,),skiprows=1)*1000
    #plt.plot(x0,z0,linestyle='solid', linewidth=1,color='steelblue',)
    figname1='couple_half-oneside_sum.png'
    #plt.scatter(x0,z0,s=1,color='mediumblue') #
    plt.scatter(x0,z0,s=5,color='black') #
    plt.xlabel('E - Ef (eV)')
    #plt.xticks(np.arange(0,E*1.001,E/5))
    plt.ylabel('Coupling (meV)')
    plt.savefig(figname1, dpi=400)
    print("\n%s has been saved."%figname1)
    
    ###########################################################
    fig = plt.figure()
    ax=plt.axes()
    figsize_x = 4.8
    figsize_y = 3.6 # in inches
    fig.set_size_inches(figsize_x, figsize_y)
    figname2='couple_oneside_sum.png'
    
    x1=np.loadtxt(file2,usecols=(0,),skiprows=1) -Ef #E-Efermi
    z1=np.loadtxt(file2,usecols=(1,),skiprows=1)*1000 # change unit to mev
    figname2='couple_oneside_sum.png'
    #plt.scatter(x0,z0,s=1,color='mediumblue') #
    plt.scatter(x1,z1,s=5,color='black') #
    plt.xlabel('E - Ef (eV)')
    #plt.xticks(np.arange(0,E*1.001,E/5))
    plt.ylabel('Coupling (meV)')
    plt.savefig(figname2, dpi=400)
    print("\n%s has been saved."%figname2)
    plt.show()

def plot_couple_oneside_all(file1, file2, file3,  Ef,  ):
    fig = plt.figure()
    ax=plt.axes()
    figsize_x = 4.8
    figsize_y = 3.6 # in inches
    fig.set_size_inches(figsize_x, figsize_y)
    figname1='couple_oneside_all.png'
    x0=np.loadtxt(file1,usecols=(0,)) -Ef #E-Efermi
    z0=np.loadtxt(file1,usecols=(2,))*1000 # change unit to mev
    x1=np.loadtxt(file2,usecols=(0,)) -Ef #E-Efermi
    z1=np.loadtxt(file2,usecols=(2,))*1000 # change unit to mev
    x2=np.loadtxt(file3,usecols=(0,)) -Ef #E-Efermi
    z2=np.loadtxt(file3,usecols=(2,))*1000 # change unit to mev
    #plt.plot(x0,z0,linestyle='solid', linewidth=1,color='steelblue',)
    plt.scatter(x0,z0,s=5,color='black') #
    plt.xlabel('E - Ef (eV)')
    #plt.xticks(np.arange(0,E*1.001,E/5))
    plt.ylabel('Coupling (meV)')
    plt.savefig(figname1, dpi=400)
    print("\n%s has been saved."%figname1)
    
    fig = plt.figure()
    ax=plt.axes()
    figsize_x = 4.8
    figsize_y = 3.6 # in inches
    fig.set_size_inches(figsize_x, figsize_y)
    figname2='couple_diff_same.png'
    plt.scatter(x1,z1,s=5,color='red') #
    plt.scatter(x2,z2,s=5,color='blue',alpha=0.1) #
    #plt.yticks(np.arange(cmin,cmax*1.001,cmax/5))
    #plt.ylim(0,1)
    #plt.yscale('log')
    plt.xlabel('E - Ef (eV)')
    #plt.xticks(np.arange(0,E*1.001,E/5))
    plt.ylabel('Coupling (meV)')
    plt.savefig(figname2, dpi=400)
    print("\n%s has been saved."%figname2)
    plt.show()

def watch_disturibution(file1,Efermi, file2='EIGTXT'):
    data = np.loadtxt(file1)
    en = np.loadtxt(file2) - Efermi
    time = data[:,0]
    distur = data[:,2:]
    index = np.where(distur[0,:]==1)[0] 
    print(len(distur),index)
    
    fig = plt.figure()
    ax = plt.axes() 
    
    for i in range(300,1100,50):
        plt.plot(en, distur[i,:])
    
    plt.show()
    


def plot_3D_Energy_Time(file1, figname='3D_Energy_Time.png'):
    fig = plt.figure()
    ax  = plt.axes(projection='3d')
      
    data = np.loadtxt(file1)[:2000]              # select data range
    X = data[:,0] 
    num_en = np.arange(int(len(data[0,1:])/2)) + 1   #how many kinds of initial level
    Y = np.hstack((num_en,num_en))                   # initial level * materials number
    Z = data[:,1:]                                   # average energy in  initial level
    print(num_en,Y)
    colorgroup = ['red', 'blue']
    
    for i in range(len(Y)):
        x = X
        y = np.full(len(X),Y[i])
        z = Z[:, i]
        if i < len(Y)/2:            # select color from colorgroup
            j = 0
        else:
            j = 1 
        ax.plot(x, y, z, color = colorgroup[j] )
        
    my_x = np.arange(0, 2000, 500)          # time
    my_y = num_en                           
    my_z = np.arange(0, 1.6, 0.25)          #initial level
    my_xticks = np.arange(0,2.0,0.5)
    my_yticks = np.array(['0.5', '0.75','1.0','1.25','1.5'])
    
    ax.set_ylim(0.25,6)
    ax.set_xlim(0,2000)
    ax.set_zlim(0,1.6)
    ax.xaxis.set_ticks(my_x,)
    ax.xaxis.set_tick_params(rotation=45, pad=-6.5)
    ax.yaxis.set_ticks(my_y)
    ax.yaxis.set_tick_params(rotation=-10, pad=-3)
    ax.zaxis.set_ticks(my_z)
    ax.zaxis.set_tick_params(rotation=-13, pad=-1)
    
    ax.w_xaxis.set_pane_color((144/255, 238/255, 144/255, 1))
    ax.w_yaxis.set_pane_color((255/255, 250/255, 205/255, 1))
    ax.w_zaxis.set_pane_color((220/255, 220/255, 220/255, 1))
    
    ax.set_xticklabels(my_xticks, )
    ax.set_yticklabels(my_yticks, )
    
    for i in range(len(my_x)):
        ax.get_xticklabels()[i].set_fontweight("bold")
    for i in range(len(my_y)):
        ax.get_yticklabels()[i].set_fontweight("bold")
    for i in range(len(my_z)):
        ax.get_zticklabels()[i].set_fontweight("bold")
        
    plt.gca().set_box_aspect((10, 8, 8))     #x:y:z
    
    plt.xlabel('Time (ps)',family='Times New Roman',labelpad=-3) #fontweight='semibold')
    plt.ylabel('Initial level (eV)',family='Times New Roman',labelpad=-2,)    
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel('Energy (eV)', rotation = 85, labelpad=-1)
  
    plt.savefig(figname, dpi=600)
    print("\n%s has been saved."%figname)
    plt.show()
    
def fit_En_Time(file1, Efermi, figname='fit_Energy_Time.png'):
    data = np.loadtxt(file1)[:100000] 
    time = data[:,0]    # xdata
    en = data[:,1]-Efermi     # ydata
    xdata =time         #np.insert(time,0,0)
    ydata =en         #np.insert(en,0,en[0])
    
    fig = plt.figure()
    ax = plt.axes()
    
    plt.plot(xdata, ydata, 'o', color='b', markerfacecolor='w', label='Data', ms=5)
    
##############################################  
    index_special = np.where(ydata==np.max(ydata))[0]
    sigma=np.ones(len(xdata))
    sigma[index_special,] = 0.001 ###index_special+11]] = 0.001    
    param_bounds_5 = ([0,0,0,0,-np.inf],[np.inf,np.inf,np.inf,np.inf,np.inf])
    popt_5,pocv_5 = curve_fit(gauss_exponent, xdata, ydata, maxfev=4000, bounds = param_bounds_5,#
                              sigma=sigma)
    plt.plot(xdata, gauss_exponent(xdata, *popt_5), '-', lw=2, color = 'orange',
             label='gauss_exponent')
    print(#'fitting exponent_2 results:a1=%5.6f a2=%5.6f  b1=%5.6f b2=%5.6f  c=%5.6f ' % tuple
          'gauss_exponent:',(popt_5))  
        
##############################
    index_special = np.where(ydata==np.max(ydata))[0]
    sigma=np.ones(len(xdata))
    sigma[index_special,] = 0.001 ###index_special+11]] = 0.001     
    popt_2,pocv_2 = curve_fit(gauss_1, xdata, ydata, maxfev=14000,
                              sigma=sigma)
    plt.plot(xdata, gauss_1(xdata, *popt_2), 'g-', lw=2, 
             label='gauss_1')
    print(#'fitting exponent_2 results:a1=%5.6f a2=%5.6f  b1=%5.6f b2=%5.6f  c=%5.6f ' % tuple
          'gauss_1:',(popt_2))
    
    #print(
    #    (ydata[0]-ydata[np.where(xdata>=abs(popt_2[1]))[0][0]])/xdata[np.where(xdata>=abs(popt_2[1]))[0][0]]*1000, ) 
          
############################
    index_special = np.where(ydata==np.max(ydata))[0]
    sigma=np.ones(len(xdata))
    sigma[index_special,] = 0.001 ###index_special+11]] = 0.001     
    param_bounds_3 = ([0,0,0,0,-np.inf],[np.inf,np.inf,np.inf,np.inf,np.inf])    #([0,0,0,0,-np.inf],[np.inf,np.inf,np.inf,np.inf,np.inf])
    popt_3,pocv_3 = curve_fit(gauss_2, xdata, ydata, maxfev=14000, bounds = param_bounds_3,
                              sigma=sigma)
    plt.plot(xdata,  gauss_2(xdata, *popt_3), 'r-', lw=2, label='gauss_2',
             ) 
    #delta=%5.6f'
    print(#'fitting exponent_2 results:a1=%5.6f a2=%5.6f  b1=%5.6f b2=%5.6f  c=%5.6f ' % tuple
          'gauss_2:',(popt_3))
    
    
    print(
        'gauss_1-ELR:',(ydata[0]-ydata[np.where(xdata>=abs(popt_2[1]))[0][0]])/xdata[np.where(xdata>=abs(popt_2[1]))[0][0]]*1000,  
         'gauss_2-ELR1:',(ydata[0]-ydata[np.where(xdata>=abs(popt_3[2]))[0][0]])/xdata[np.where(xdata>=abs(popt_3[2]))[0][0]]*1000,
         'gauss_2-ELR2:', (ydata[0]-ydata[np.where(xdata>=abs(popt_3[3]))[0][0]])/xdata[np.where(xdata>=abs(popt_3[3]))[0][0]]*1000,
          )
###########
    ax.set_ylim(0,)
    plt.legend(loc='upper right', frameon=False, borderpad=0.3, handlelength=1.5, 
               prop={'weight' :'normal', 'family':'Times New Roman', 'size':'13'} )
    plt.xlabel('Time (fs)',family='Times New Roman', labelpad=2 , fontsize = 15, fontweight='normal')
    plt.ylabel('Energy (eV)',family='Times New Roman', labelpad= 2, fontsize = 15, fontweight='normal')
#save picture  
    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved."%figname)
    plt.show()

def fit_En_Time_1_e(file1, Efermi, figname='fit_Energy_Time_1_e.png'):
    data = np.loadtxt(file1)[:] 
    time = data[:,0]    # xdata
    en = data[:,1]-Efermi     # ydata
    xdata =time         
    ydata =en         
    
    fig = plt.figure()
    ax = plt.axes()
    plt.plot(xdata, ydata, 'o', color='b', markerfacecolor='w', label='Data', ms=5)
    
#############find 1/e*y    
    e_index = np.where(ydata >= (np.max(ydata)/math.e))[0][-1]
    e_y  = ydata[e_index]
    delta_y = ydata[0]-ydata[e_index]
    delta_x = xdata[e_index]
    delta_P = delta_y / delta_x *1000 ###unit ev/ps
    
    print('t=%d'%xdata[e_index], 'P=%f'%delta_P)
    
    plt.plot([0, delta_x], [e_y, e_y], c='b', linestyle='--')
    plt.plot([delta_x,delta_x], [0, e_y], c='b', linestyle='--')
    plt.show()

def fit_OD_Time_1_e(file1, file2,figname='fit_OD_Time_1_e.png'):
    data = np.loadtxt(file1)[:] 
    time = data[:,0]    # xdata
    OD = data[:,1]     # ydata
    xdata =time         
    ydata =OD       
    en_data = np.loadtxt(file2)[:] 
    en = en_data[:,1]
    
    fig = plt.figure()
    ax = plt.axes()
    plt.plot(xdata, ydata, 'o', color='b', markerfacecolor='w', label='Data', ms=5)
    
#############find 1/e*y    
    e_index = np.where(ydata >= (np.max(ydata)/math.e))[0][-1]
    e_y  = ydata[e_index]
    delta_y = en[0]-en[e_index]
    delta_x = xdata[e_index]
    delta_P = delta_y / delta_x *1000 ###unit ev/ps
    
    print('t=%d'%xdata[e_index], 'P=%f'%delta_P)
    
    plt.plot([0, delta_x], [e_y, e_y], c='b', linestyle='--')
    plt.plot([delta_x,delta_x], [0, e_y], c='b', linestyle='--')
    plt.show()

def fit_OD_Time(file1,file2,figname='fit_OD_Time.png'):
    data = np.loadtxt(file1)[:,:]
    time = data[:,0]    # xdata
    OD = data[:,1]     # ydata
    en_data = np.loadtxt(file2)[:] 
    en = en_data[:,1]
    
    xdata =  time  # np.insert(time,0,0)
    ydata =  OD  # np.insert(OD,0,OD[0]) / OD[0]  #norm
    
    xdata = time
    ydata = OD / OD[0]
    
    fig = plt.figure()
    ax = plt.axes()

    plt.plot(xdata, ydata, 'o', color='b', markerfacecolor='w', label='Data', ms=5)
##############################################  
    '''
    param_bounds_1 = ([0,0,0,0,0],[np.inf,np.inf,np.inf,np.inf,np.inf])
    popt_1,pocv_1 = curve_fit(gauss_exponent, xdata, ydata, maxfev=4000, bounds = param_bounds_1)
    plt.plot(xdata, gauss_exponent(xdata, *popt_1), 'y-', lw=4, 
             label='gauss_exponent')
    print(#'fitting exponent_2 results:a1=%5.6f a2=%5.6f  b1=%5.6f b2=%5.6f  c=%5.6f ' % tuple
          'gauss_exponent:',(popt_1))
    '''
    '''
##############################################  
    param_bounds_4 = ([0,0,0,0],[1,1,np.inf,np.inf])
    popt_4,pocv_4 = curve_fit(exponent_2, xdata, ydata, maxfev=4000, bounds = param_bounds_4)
    plt.plot(xdata, exponent_2(xdata, *popt_4), '-', lw=2, color = 'black',
             label='exponent_2')
    print(#'fitting exponent_2 results:a1=%5.6f a2=%5.6f  b1=%5.6f b2=%5.6f  c=%5.6f ' % tuple
          'exponent_2:',(popt_4))    

##############################################  
    param_bounds_5 = ([0,0,0,0,-np.inf],[1,1,np.inf,np.inf,np.inf])
    popt_5,pocv_5 = curve_fit(gauss_exponent, xdata, ydata, maxfev=4000, bounds = param_bounds_5)
    plt.plot(xdata, gauss_exponent(xdata, *popt_5), '-', lw=2, color = 'orange',
             label='gauss_exponent')
    print(#'fitting exponent_2 results:a1=%5.6f a2=%5.6f  b1=%5.6f b2=%5.6f  c=%5.6f ' % tuple
          'gauss_exponent:',(popt_5))  
    '''   
##############################################
    ##param_bounds_6 = ([0,0,0,0,0,0,-np.inf],[1,1,1,np.inf,np.inf,np.inf,np.inf])
    #popt_6,pocv_6 = curve_fit(gauss_exponent_2, xdata, ydata, maxfev=4000, bounds = param_bounds_6)
    #plt.plot(xdata, gauss_exponent_2(xdata, *popt_6), '-', lw=2, color = 'purple',
    #         label='gauss_exponent_2')
   # print(#'fitting exponent_2 results:a1=%5.6f a2=%5.6f  b1=%5.6f b2=%5.6f  c=%5.6f ' % tuple
    #      'gauss_exponent_2:',(popt_6))   


###############################   
    param_bounds_2 = ([0,0,-np.inf],[1,np.inf,np.inf])
    popt_2,pocv_2 = curve_fit(gauss_1, xdata, ydata, maxfev=4000,bounds = param_bounds_2)
    plt.plot(xdata, gauss_1(xdata, *popt_2), 'g-', lw=2, 
         label='gauss_1')
    print(#'fitting exponent_2 results:a1=%5.6f a2=%5.6f  b1=%5.6f b2=%5.6f  c=%5.6f ' % tuple
      'gauss_1:',(popt_2))


############################
    param_bounds_3 = ([0,0.5,0,0,0],[0.5,1,np.inf,np.inf,np.inf])
    popt_3,pocv_3 = curve_fit(gauss_2, xdata, ydata, maxfev=24000, bounds = param_bounds_3)
    plt.plot(xdata,  gauss_2(xdata, *popt_3), 'r-', lw=1, 
             label='gauss_2')
    #delta=%5.6f'
    print(#'fitting exponent_2 results:a1=%5.6f a2=%5.6f  b1=%5.6f b2=%5.6f  c=%5.6f ' % tuple
          'gauss_2:',(popt_3))
###########
    print(
        'gauss_1_P:',(en[0]-en[np.where(xdata>=abs(popt_2[1]))[0][0]])/xdata[np.where(xdata>=abs(popt_2[1]))[0][0]]*1000,  
          'gauss_2_P1:',(en[0]-en[np.where(xdata>=abs(popt_3[3]))[0][0]])/xdata[np.where(xdata>=abs(popt_3[3]))[0][0]]*1000,
    'gauss_2_P2:',(en[0]-en[np.where(xdata>=abs(popt_3[2]))[0][0]])/xdata[np.where(xdata>=abs(popt_3[2]))[0][0]]*1000)
    ax.set_ylim(0,)
    ax.set_xlim(0,)
    plt.legend(loc='upper right', frameon=False, borderpad=0.3, handlelength=1.5, 
               prop={'weight' :'normal', 'family':'Times New Roman', 'size':'13'} )
    plt.xlabel('Time (fs)',family='Times New Roman', labelpad=2 , fontsize = 15, fontweight='normal')
    plt.ylabel('OD ',family='Times New Roman', labelpad= 2, fontsize = 15, fontweight='normal')
#save picture  
    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved."%figname)
    plt.show()



def exponent_1(x, a1,  b1):
    return  a1 * np.exp(-x / b1) \
                #+5.3783725-5.162245

def exponent_2(x, a1, a2,  b1, b2):
         
    return  a1 * np.exp(-x / b1)+ (a2) * np.exp(-x / b2) 

def exponent_3(x, a1, a2, a3, b1, b2,b3, c ):      
    return  a1 * np.exp(-x / b1)+ a2 * np.exp(-x / b2)+ a3 * np.exp(-x / b3) +c

def exponent_4(x, a1, a2, a3,a4 , b1, b2,b3,b4, c ):     
    return  a1 * np.exp(-x / b1)+ a2 * np.exp(-x / b2)+ a3 * np.exp(-x / b3) +a4 * np.exp(-x / b4) +c

def sigmoid_exponent(x, a1, a2,a3, b1, b2,  c ):    
    return a1 * np.exp(-x / b1)+ a2 / (a3+np.exp( -x / b2 ) ) +c

def gauss_1(x,a1,b1,c):#c):
    
    return  a1* np.exp(-(x /  b1)**2 /2  )  +c  \
           #+5.3783725-5.162245

def gauss_2(x,a1,a2,b1,b2,c):
    
    return  a1 * np.exp(-(x /  b1)**2 /2 ) + a2 * np.exp(-(x /  b2)**2 /2 )+c  \
             #+5.3783725-5.162245  #for gdy ; 5.414105 - 5.35655516 #for gdy-ABC
        
def gauss_3(x,a1,a2,a3,b1,b2,b3,):
    return  a1 * np.exp(-(x /  b1)**2 ) +a2 * np.exp(-(x /  b2)**2 ) \
            +a3 * np.exp(-(x /  b3)**2 ) 

def gauss_exponent(x, a1, a2,   b1, b2 ,  c):
    return a1 * np.exp(-(x / b1)**2 )+a2 * np.exp(-x / b2)+c  \
           # +5.3783725-5.162245  #for gdy    5.414105 - 5.35655516 = 0.05754984 #for gdy-ABC

def gauss_exponent_2(x, a1, a2, a3,  b1, b2 , b3, c):
    return a1 * np.exp(-(x / b1)**2 )+ a2 * np.exp(-x / b2)+ (1-a1-a2) * np.exp(-x / b3)+c  \
    
def get_lifetime_loss(file1, figname='lifetime.txt'):
    
    data = np.loadtxt(file1)[:2000] 
    time = data[:,0]
    en = data[:,1:]
    num_all = np.arange(int(len(data[0,1:])))
    num_en = 5
    
   # def func_expon_decay(x, k, t ):
   #     return k*np.exp(-x/t)
       
    t = np.zeros(shape=(2,num_en), dtype=int)
    en_loss = np.zeros(shape=(2,num_en))
    k=0
    
    for i in num_all:
        en_ini = en[0,i]
        
        if i < num_en:
            en_cb = 5.3783725-5.162245
            y = (en_ini - en_cb)*np.exp(-1) + en_cb
            j = 0
            t_index = np.where(en[:,i] < y )[0][0]
            t[j][i]= t_index + 1
            en_loss[j,i] = np.around((en_ini -en[t_index,i])/t[j,i] *1000, 4)
        else:
            en_cb = 0
            y = (en_ini - en_cb)*np.exp(-1)
            j = 1
            t_index = np.where(en[:,i] < y )[0][0]
            t[j][k]= t_index + 1
            en_loss[j,k] = np.around((en_ini -en[t_index,i])/t[j,k] *1000, 4)
            k= k+1
    
    
    seq = np.arange(5)+1
    t = np.insert(t.T/1000, 0, seq, axis=1)
    
    lifetime_loss = np.hstack((t,en_loss.T))
    
    #lifetime_loss = np.insert(lifetime_loss, 0, seq, axis=1)
    
    np.savetxt('lifetime_loss.txt',lifetime_loss)
    
    ############## plot 
    x = lifetime_loss[:,0]
    y = lifetime_loss[:,1:]
    fig = plt.figure()
    ax1 = fig.subplots()
    ax2 = ax1.twinx()
    
    line1=ax1.plot(x,y[:,0],'s-',color='r',label='GDY')
    line2=ax1.plot(x,y[:,1],'o-',color='b',label='GR')
    line3=ax2.plot(x,y[:,2],'s--',color='r',label='3')
    line4=ax2.plot(x,y[:,3],'o--',color='b',label='4')
    
    ax1.set_xlabel('E$_{ini}$ (eV)',family='Times New Roman',labelpad=-1)
    ax1.set_ylabel('Lifetime (ps)',family='Times New Roman')
    ax2.set_ylabel('Energy loss rate (eV/ps)',family='Times New Roman')    
    
    my_x = np.arange(1, 6, 1)          # time
    my_y = np.arange(0, 2.2, 0.5)
    my_y2 = np.arange(0, 4, 1)        
    my_xticks = np.array(['0.5', '0.75','1.0','1.25','1.5'])
    #my_yticks = np.array(['0.5', '0.75','1.0','1.25','1.5'])
    ax1.set_xlim(0.5,5.5)
    ax1.set_ylim(0,2.2)
    ax2.set_ylim(0,4)
    
    ax1.xaxis.set_ticks(my_x)
    ax1.xaxis.set_tick_params(pad=1)
    ax1.yaxis.set_ticks(my_y)
    ax1.yaxis.set_tick_params(pad=2)
    ax2.yaxis.set_ticks(my_y2)
    ax2.yaxis.set_tick_params(pad=2)
    ax1.set_xticklabels(my_xticks, )
    #ax1.set_yticklabels(my_yticks, )
    #ax2.set_yticklabels(my_yticks, )
    plt.gca().set_box_aspect(0.8)     #x:y
    
    legend_font = {
    'family': 'Times New Roman',  # 字体
    'style': 'normal',
    #'size': font_size,  # 字号
    'weight': "normal",  # 是否加粗，不加粗
}
    
    
    ax1.legend(
    ##[line1,('-'),line2,('-')],
    ['GDY','GR'],
    bbox_to_anchor=(0.5, 0.8),
    loc='lower center',  # 图例的底部中央位置在图像上部居中
    frameon=False,  # 不显示图例框线
    prop=legend_font
)
    ax2.legend(
    #,
    ['t','ELR'],
    bbox_to_anchor=(0.5, 0.65),
    loc='lower center',  # 图例的底部中央位置在图像上部居中
    frameon=False,  # 不显示图例框线
    prop=legend_font
)
    
    
    plt.show()

def get_occu_elec_time(file1, figname='occu-elec-time.txt'):    
    data = np.loadtxt(file1)[:2000] 
    time = data[:,0]
    en = data[:,1:]
    num_all = np.arange(int(len(data[0,1:])))
    num_en = 5
    
    t = np.zeros(shape=(2,num_en), dtype=int)
    
    k=0
    
    for i in num_all:
        if i < num_en:
            j = 0
            t_index = np.where(en[:,i] < 0.01 )[0][0]
            t[j][i]= t_index + 1
            
        else:
            j = 1
            t_index = np.where(en[:,i] < 0.01 )[0][0]
            t[j][k]= t_index + 1
            k= k+1
        print(t_index)
    seq = np.arange(5)+1
    t = np.insert(t.T/1000, 0, seq, axis=1)
    np.savetxt('occu_time.txt',t)
    
    
    
if __name__=='__main__':
    main()

