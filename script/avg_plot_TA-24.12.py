# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 21:32:46 2023

@author: daifulong
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec

from mpl_toolkits import mplot3d 
import math
from scipy import special
from scipy.optimize import curve_fit
from brokenaxes import brokenaxes

def main():
    plt.rcParams['font.family'] = 'Arial' #'Times New Roman'
    np.set_printoptions(precision=5)
    
    matrial = np.array([ 2, 3, ])  #2 for GR; 3 for GDY
    func_select = 3
    #avg_ref_temp()
    #avg_temp()
    #plot_TA_all('avg_ref_temp.txt')
    plot_TA_2_1('avg_ref_temp.txt')
    
    #fitting_TA('5259_OD_time_norm_zero.txt',func_select)   #change file name
    #fitting_TA('4700_OD_time_norm_zero.txt',func_select)   #change file name
    #fitting_TA('5297_OD_time_norm_zero.txt',func_select)   #change file name
    
def avg_ref_temp(prefix = 'Ref_temp', figname = 'avg_ref_temp.txt'):
    filepath = './'
    filename_list = os.listdir(filepath)
    filelist = []
    for filename in filename_list:
        if filename[0:8] == prefix:
            filelist.append(filename)

    tem_shape = np.loadtxt(filelist[0] ).shape
    nu = len(filelist)
    avg_ref_temp = np.zeros((nu,tem_shape[0],tem_shape[1])) #larger than last required data
    n=0
    for file in filelist:
        sample = np.loadtxt(file )
        avg_ref_temp[n,:] = sample
        n=n+1
        if n==11:             #abandon bad data, 16 for gdy;
            print("select the top %d temps"%n)
            break      
    avg_ref_temp = np.sum(avg_ref_temp , 0) / n      
    np.savetxt('avg_ref_temp.txt', avg_ref_temp)  #mean for required data
    print("\n%s has been saved."%figname)
    
def avg_temp(prefix = 'temp', figname = 'avg_temp.txt'):
    filepath = './'
    filename_list = os.listdir(filepath)
    filelist = []
    for filename in filename_list:
        if filename[0:4] == prefix:
            filelist.append(filename)

    tem_shape = np.loadtxt(filelist[0] ).shape
    nu = len(filelist)
    avg_temp = np.empty((nu,tem_shape[0],tem_shape[1]))
    n=0
    for file in filelist:
        sample = np.loadtxt(file )
        avg_temp[n,:] = sample
        n=n+1
        if n==11:             #abandon bad data, 16 for gdy;
            print("select the top %d temps"%n)
            break
    avg_temp = np.sum(avg_temp , 0) / n    
    np.savetxt('avg_temp.txt', avg_temp)
    print("\n%s has been saved."%figname)
    
def plot_TA_all(file, figname='TA_all.png'):
    all_data = np.loadtxt(file)          
    freq_range = all_data[1:,0]
    freq_intensity = all_data[1:,1]
    time = all_data[0, 2:]
    data = all_data[1:,2:]
    
    #OD_max_index = np.max(data)
    #time_max = time[OD_max_index] 
    
    n_start = 18 ; n_end = 130     #change relaxation time
    x0 = freq_range
    y0 = time[n_start:n_end]
    z0 = data[:,n_start:n_end] *1000  # mΔOD
    
    index_all = []
    for i in range(len(x0)):
        index_all.insert(-1,np.where(z0[i] == np.max(z0[i]))[0])
    
    index_0 = np.min(index_all)   #make y0 zero 
    n_0 = y0[index_0]
    print(index_0,'n0=',n_0)
    y0 = y0-n_0
    
    
    fig = plt.figure()
    plt.rc('font',family='Times New Roman',size=11)
    ax=plt.axes()
    ax.spines['bottom'].set_linewidth(1)                 #set frame
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['top'].set_linewidth(1)
    figsize_x = 4.8
    figsize_y = 3.6 # in inches
    fig.set_size_inches(figsize_x, figsize_y)
# colormap
    cmap0 = plt.cm.jet
    color0=cmap0(np.linspace(0,1,200))
    color0_1=cmap0(np.linspace(0,1,80))
    newcolorsgroup=np.vstack((color0_1[18:49],color0[125:180]))
    newcmap=mpl.colors.ListedColormap(newcolorsgroup)    
    cmin=round(np.min(z0),2);  cmax=round(np.max(z0),2)
    norm = mpl.colors.Normalize(0,cmax)
       
# plot contourf    
    X,Y = np.meshgrid(x0,y0)
    print(x0.shape, X.shape, z0.shape, np.min(z0) )
    
    plt.contourf(X,Y,z0.T,50,cmap=newcmap,norm=norm)
    
    #bounds= np.round(np.arange(0,cmax,(cmax-0)/4),decimals=2)
    bounds = np.linspace(0,cmax,5)
    #bounds = np.array(['0', '0.13','0.26','0.39', '0.52'])
    ##########colorbar    
    cb=plt.colorbar(
            mpl.cm.ScalarMappable(cmap=newcmap,norm=norm),
            ax=ax,
            #extend='min',
            ticks=bounds,
            spacing='proportional',
            orientation='vertical',
            
        )
    font = {'family' : 'Times New Roman',
        'color'  : 'black',
        'weight' : 'bold',
        'size'   : 14,
        }
    cb.ax.tick_params(labelsize=12)
    for i in range(len(bounds)):
        cb.ax.get_yticklabels()[i].set_fontweight("bold")
    #bounds = np.array(['0', '0.13','0.26','0.39', '0.52'])
    cb.ax.set_yticklabels(bounds )
    cb.set_label('m$\Delta$OD',fontdict=font, labelpad = 6 ) #设置colorbar的标签字体及其大小    
 ######### ticks labels   
    plt.xlabel('Wavelength (nm)',family='Times New Roman', labelpad=1 , fontsize = 15, fontweight='bold')
    plt.ylabel('Time (ps)',family='Times New Roman', labelpad= 2, fontsize = 15, fontweight='bold')
    
    
    xmin=np.min(x0); xmax=np.max(x0)
    ymin=np.min(y0); ymax=np.max(y0) 
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(-0.6,12)
    my_x = np.arange(4800,5300, 100)      #
    my_y = np.arange(0,14,2)                   
    #my_yticks = np.arange(n_0,ymax,1)
    plt.minorticks_on() # 显示副刻度线
    plt.tick_params(top=False,bottom=True,left=True,right=False) 
    plt.rcParams['xtick.direction'] = 'in'#将x周的刻度线方向设置向内
    ax.xaxis.set_minor_locator(MultipleLocator(50)) # 设置 X 轴上的副刻度线之间的间隔为0.4
    ax.yaxis.set_minor_locator(MultipleLocator(1)) # 设置 Y 轴上的副刻度线之间的间隔为0.4
    
    ax.xaxis.set_ticks(my_x)
    ax.xaxis.set_tick_params(pad=3, labelsize=13 )
    ax.yaxis.set_ticks(my_y)
    ax.yaxis.set_tick_params(pad=2, labelsize=13)
    #ax.set_xticklabels(my_xticks, )
    #ax.set_yticklabels(my_yticks )
    for i in range(len(my_x)):
        ax.get_xticklabels()[i].set_fontweight("bold")
    for i in range(len(my_y)):
        ax.get_yticklabels()[i].set_fontweight("bold")   
 #save picture  
    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved."%figname)
    plt.show()   

###############################################################  get max freq-intensity
    fig = plt.figure()
    ax = plt.axes()
    plt.plot(freq_range, freq_intensity)
    fig = plt.figure()
    ax = plt.axes()
    for i in range(len(freq_range)):
        plt.scatter(freq_range[i], np.max(z0[i,:])/1000)
    
    index = np.where(freq_intensity == np.max(freq_intensity))
    index_max = np.where(data == np.max(data))
    wave_max = freq_range[index_max[0]]
    print('wave_max=%f'%wave_max , index_max, 'z0_max=%f'%np.max(z0), 'freq_sum=%f'%freq_range[index])

######### get OD_time.txt
    data_max = data[index_max[0],:].flatten()
    OD_max = data_max  *1000
    
    OD_time = np.vstack((time,OD_max)).T
    np.savetxt('%d_OD_time.txt'%wave_max, OD_time)
    fig = plt.figure()
    ax=plt.axes()
    plt.plot(time, OD_max)
    
########## get OD_time_zero.txt
    time_zero = time - time[index_max[1]]        #### make time 0
    OD_time_zero = np.vstack((time_zero,OD_max)).T
    np.savetxt('%d_OD_time_zero.txt'%wave_max, OD_time_zero)
    fig = plt.figure()
    ax=plt.axes()
    plt.plot(time_zero, OD_max)
    #ax.set_xlim(-4,10)  

########## get OD_time_norm.txt
    OD_max_norm = OD_max  / np.max(z0)
    
    OD_time_norm_zero = np.vstack((time_zero,OD_max_norm)).T
    np.savetxt('%d_OD_time_norm_zero.txt'%wave_max, OD_time_norm_zero)
    fig = plt.figure()
    ax=plt.axes()
    plt.plot(time_zero, OD_max_norm)
    #ax.set_xlim(-4,10) 

######################## get 5297_OD_time.txt or last data
    data_last = data[-1,:].flatten()
    index_last = np.where(data_last == np.max(data_last))[0]
    time_last = time - time[index_last]    
    
    data_last_norm= data_last/np.max(data_last)
    OD_last_time = np.vstack((time_last,data_last)).T
    OD_last_time_norm = np.vstack((time_last,data_last_norm)).T
    np.savetxt('5297_OD_time_zero.txt',OD_last_time)
    np.savetxt('5297_OD_time_norm_zero.txt',OD_last_time_norm)
    fig = plt.figure()
    ax=plt.axes()
    plt.plot(time_last, data_last)
    print('5297_OD_time.txt is saved')

######################## get 5259_OD_time.txt
    data_select = data[-5,:].flatten()*1000
    index_select = np.where(data_select == np.max(data_select))[0]
    time_select =time - time[index_select]
    
    data_select_norm = data_select/np.max(data_select)
    OD_select_time = np.vstack((time_select,data_select)).T
    OD_select_time_norm = np.vstack((time_select,data_select_norm)).T
    np.savetxt('5259_OD_time_zero.txt',OD_select_time)
    np.savetxt('5259_OD_time_norm_zero.txt',OD_select_time_norm)
    fig = plt.figure()
    ax=plt.axes()
    plt.plot(time_select, data_select)
    print('5259_OD_time_zero.txt is saved')

##################### get 4700_OD_time_norm_zero.txt
    data_select_4700 = data[0,:].flatten()*1000
    index_select_4700 = np.where(data_select_4700 == np.max(data_select_4700))[0]
    time_select_4700 =time - time[index_select_4700]
    
    data_select_norm_4700 = data_select_4700/np.max(data_select_4700)
    OD_select_time_4700 = np.vstack((time_select_4700,data_select_4700)).T
    OD_select_time_norm_4700 = np.vstack((time_select_4700,data_select_norm_4700)).T
    np.savetxt('4700_OD_time_zero.txt',OD_select_time_4700)
    np.savetxt('4700_OD_time_norm_zero.txt',OD_select_time_norm_4700)
    fig = plt.figure()
    ax=plt.axes()
    plt.plot(time_select_4700, data_select_4700)
    print('4700_OD_time_zero.txt is saved')  

 


def plot_TA_2_1(file, figname='TA_2_1_arial.svg'):
    
    all_data = np.loadtxt(file)          
    freq_range = all_data[1:,0]
    freq_intensity = all_data[1:,1]
    time = all_data[0, 2:]
    data = all_data[1:,2:]
    
    OD_max_index = np.where(data== np.max(data))
    OD_max_index_time = np.where(data== np.max(data))[1]
    #print(OD_max_index )
    time_max = time[OD_max_index_time]
    freq_max = freq_range[OD_max_index[0]]
    
    #n_start = 18 ; n_end = 130     #change relaxation time
    x0 = freq_range                #wavelength
    y0 = time - time_max   #y0 = time[n_start:n_end]   # change time of max_OD to 0
    z0 = data*1000      #z0 = data[:,n_start:n_end] *1000  # mΔOD
    #print(time_max,'OD_max:%f'%(y0[OD_max_index_time]), z0[0,OD_max_index_time],'freq_max:%f'%(freq_max))
    #index_all = []
    #for i in range(len(x0)):
    #    index_all.insert(-1,np.where(z0[i] == np.max(z0[i]))[0])
    
    #index_0 = np.min(index_all)   #make y0 zero 
    #n_0 = y0[index_0]
    #print(index_0,'n0=',n_0)
    #y0 = y0-n_0
 ########################################## plt.figure
    fig = plt.figure()
    #ax1 = fig.add_subplot(221, )
    #ax1 = plt.subplot2grid((3,3),(0,0),colspan=3,rowspan=1)
    #ax2 = fig.add_subplot(223)
    #ax3 = fig.add_subplot(224)
    ax1 = plt.subplot2grid((17,11),(0,0),colspan=10,rowspan=11)  ###!!! 网格 ，起始位置(行，列)， 占用几列，几行 
    ax2 = plt.subplot2grid((17,11),(11,0),colspan=10,rowspan=6)
    ax3 = plt.subplot2grid((17,11),(0,10),colspan=1,rowspan=10) #color bar
    
    #fig,(ax1,ax2)=plt.subplots(2,1,sharex=True) #共享X轴
    #fig,ax=plt.subplots(2,1)#,sharex=True)
    
    #ax1 = ax[0]
    #ax2 = ax[1]
    #ax3 = ax[0,1]
    #ax2 = ax1.twinx()   #共用一个图
    #gs = gridspec.GridSpec(2, 1,  height_ratios=[3,2]   ) #width_ratios=[1,1],
    #ax1 = fig.add_subplot(gs[0,0])
    #ax2 = fig.add_subplot(gs[1,0])
    fig.subplots_adjust(hspace=0,wspace=1,left=0.14, bottom=0.14, right=0.86, top=0.9)
    plt.rc('font',family='Arial',size=11)
    
    ax1.spines['bottom'].set_linewidth(1)                 #set frame
    ax1.spines['left'].set_linewidth(1)
    ax1.spines['right'].set_linewidth(1)
    ax1.spines['top'].set_linewidth(1)
    ax2.spines['bottom'].set_linewidth(1)                 #set frame
    ax2.spines['left'].set_linewidth(1)
    ax2.spines['right'].set_linewidth(1)
    ax2.spines['top'].set_linewidth(1)
    figsize_x = 6.2 #6.2   #5.6 #for 2024.7.20
    figsize_y = 5.2 #7.6   #4.6 # in inches
    fig.set_size_inches(figsize_x, figsize_y)
# colormap
    #cmap0 = plt.cm.jet
    #color0=cmap0(np.linspace(0,1,200))
    #color0_1=cmap0(np.linspace(0,1,80))
    #newcolorsgroup=np.vstack((color0_1[18:49],color0[125:180]))
    #newcmap=mpl.colors.ListedColormap(newcolorsgroup)    
    
    newcmap= plt.cm.gist_heat_r  ### hot ;hot_r
    cmin=round(np.min(z0),2);  cmax=round(np.max(z0),2)
    norm = mpl.colors.Normalize(cmin,cmax)
       
################### plot contourf for ax1
    X,Y = np.meshgrid(x0,y0)
    print(x0.shape, X.shape, z0.shape, np.min(z0) )
    
    ax1.contourf(X,Y,z0.T,50,cmap=newcmap,norm=norm,color='b')
       

###################  plot ΔOD_WAVE curve for ax2

########## 6 selected moments of TA_signals 
    time_target = np.array([-0.1, 0, 0.1, 1, 10, 100])   #!!!
    time_index = []
    #print(np.where(y0 < time_target[0])[0][-1],range(len(time_target)))
    my_color = ['#b4403e','#256ea2','#ea841e','#399335','#603c87','#9d5c39' ]  ##[(180/255,64/255,62/255)]
    for i in range(len(time_target)):
        #print(time_target[i])
        time_index.append( np.where(y0 <= time_target[i])[0][-1])
    
    #print(time_index)
    for j in range(len(time_index)):
        label_format = '{:.2f}'.format(time_target[j]).rstrip("0").rstrip(".")
        ax2.plot(x0,z0[:,time_index[j]],label=label_format+' ps',color=my_color[j] )#label='%.2f ps'%(time_target[j]))
    
    ax2.legend(
        loc = (0.05,0.65),  ##'upper left',
        ncol = 3,            # 每行几个图例
        frameon = False,      #去掉图例边框
        prop = {'size':13,
                'weight':'bold',
                'family': 'Arial',
            },
        labelspacing = 0.2,
        columnspacing= 0.5,
        handletextpad=0.5, #图例标记和标签之间的距离
        )
 
    
 ######### ticks label for ax1 ax2 
    #ax1.set_xticks([])
    #ax1.xaxis.set_visible(False)
    ax1.set_ylabel('Time (ps)',family='Arial', labelpad= 2, fontsize = 17, fontweight='bold')
    ax2.set_xlabel('Wavelength (nm)',family='Arial', labelpad=1 , fontsize = 17, fontweight='bold')
    ax2.set_ylabel('m$\Delta$OD',family='Arial', labelpad=1 , fontsize = 17, fontweight='bold')
    
    xmin=np.min(x0); xmax=np.max(x0)
    ymin=np.min(y0); ymax=np.max(y0) 
    
    ax1.set_ylim(-0.6,12)      #!!!
    ax1.set_xlim(xmin,xmax)
    ax2.set_xlim(xmin,xmax)
    ax2.set_ylim(-2,45)          ########(-0.2,2.2)       ##!!! (-2,45) for GR (-0.2,3.5) for GDY
    my_x = np.arange(4800,5300, 100)      #
    my_y1 = np.arange(0,14,2)
    my_y2 = np.arange(0,45,10) #np.arange(0,3.5,1) for GDY #np.linspace(0,45,10)  for GR                   
    #my_yticks = np.arange(n_0,ymax,1)
    plt.minorticks_on() # 显示副刻度线
    plt.tick_params(top=False,bottom=True,left=True,right=False) 
    plt.rcParams['xtick.direction'] = 'in'#将x周的刻度线方向设置向内
    ax1.xaxis.set_minor_locator(MultipleLocator(50))
    ax2.xaxis.set_minor_locator(MultipleLocator(50)) # 设置 2-X 轴上的副刻度线之间的间隔为50
    ax1.yaxis.set_minor_locator(MultipleLocator(1)) # 设置 1-Y 轴上的副刻度线之间的间隔为1
    ax2.yaxis.set_minor_locator(MultipleLocator(5))# 设置 2-Y 轴上的副刻度线之间的间隔为5 for GR# 0.5 for GDY
       
    ax1.yaxis.set_ticks(my_y1)
    ax1.yaxis.set_tick_params(pad=3, labelsize=16)
    ax2.xaxis.set_ticks(my_x)
    ax2.xaxis.set_tick_params(pad=3, labelsize=16 )
    ax2.yaxis.set_ticks(my_y2)
    ax2.yaxis.set_tick_params(pad=3, labelsize=16)
    #ax.set_xticklabels(my_xticks, )
    #ax.set_yticklabels(my_yticks )
    for i in range(len(my_x)):
        ax2.get_xticklabels()[i].set_fontweight("bold")
    for i in range(len(my_y1)):
        ax1.get_yticklabels()[i].set_fontweight("bold") 
    for i in range(len(my_y2)):
        ax2.get_yticklabels()[i].set_fontweight("bold")    

    
########## ticks label for ax3   of colorbar
    #bounds= np.round(np.arange(0,cmax,(cmax-0)/4),decimals=2)
    bounds = np.round(np.linspace(cmin,cmax,2),decimals=2)
    #bounds = np.array(['0', '0.13','0.26','0.39', '0.52'])
    ##########colorbar    
    cb=plt.colorbar(
            mpl.cm.ScalarMappable(cmap=newcmap,norm=norm),
            #ax = ax3,
            cax=ax3,
            #extend='min',
            ticks=bounds,
            spacing='proportional',
            orientation= 'vertical' #'horizontal',#'vertical'
            
        )
    font = {'family' : 'Arial',
        'color'  : 'black',
        'weight' : 'bold',
        'size'   : 14,
        }
    cb.ax.tick_params(labelsize=14)
    for i in range(len(bounds)):
        cb.ax.get_yticklabels()[i].set_fontweight("bold")
    #bounds = np.array(['0', '0.13','0.26','0.39', '0.52'])
    cb.ax.set_yticklabels(bounds, )
    cb.set_label('m$\Delta$OD',fontdict=font, labelpad = 8, ) #设置colorbar的标签字体及其大小    
    
        
 #save picture  
    #plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved."%figname)
    plt.show() 
    
    
    
    
def fitting_TA(file,func_select,figname='fitting_TA_arial.svg'):
    x0=np.loadtxt(file,usecols=0)
    y0=np.loadtxt(file,usecols=1)
    
    index_max= np.where(y0 == np.max(y0))[0]
    n_0 = x0[index_max]              # change zero point considering u on the func_select
    print(index_max,n_0)
    
    n_s = 16
    n_e = 133    #  time range 
    xdata = x0[n_s:n_e] - n_0       # change time range   20-245 -15.96
    ydata = y0[n_s:n_e]          # same
    
    #bax = brokenaxes(xlims=((-1,10),(40,55)), wspace=0.2)          ###### breakpoint
    fig = plt.figure()  
    ax = plt.axes()
    ax.spines['bottom'].set_linewidth(1)                 #set frame
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['top'].set_linewidth(1)
    figsize_x = 4.8
    figsize_y = 3.6 # in inches
    fig.set_size_inches(figsize_x, figsize_y)
    
    plt.plot(xdata, ydata, 'o', color='b', markerfacecolor='w', label='Experimental Data', ms=5)
    #bax.plot(xdata, ydata, 'o', color='b', markerfacecolor='w', label='Experimental Data', ms=5)
################    exponent_3 * gauss
    if func_select == 3:
        index_special = np.where(ydata==np.max(ydata))[0]
        sigma=np.ones(len(xdata))
        sigma[index_special] = 1  ###index_special+11]] = 0.001 
        #print(ydata[index_special-1])
        #param_bounds_3 = ([-1,-1,-1,0,0,0,-np.inf],[1,1,1,np.inf,np.inf,np.inf,np.inf])
        popt_3,pcov_3 = curve_fit(func_erf_3, xdata, ydata, maxfev=54000, # bounds = param_bounds_3, 
                       sigma=sigma)
        plt.plot(xdata, func_erf_3(xdata, *popt_3), 'r-', lw=2, 
                 label='Fitting Curve')
        #bax.plot(xdata, func_erf_3(xdata, *popt_3), 'r-', lw=2, 
        #         label='Fitting Curve')
        #delta=%5.6f'
        perr = np.sqrt(np.diag(pcov_3))
        print('fitting exponent_3 results:a1=%5.6f a2=%5.6f a3=%5.6f b1=%5.6f b2=%5.6f b3=%5.6f u=%5.6f ' % tuple(popt_3),)
             #'perr=',perr)
        print('perr=',perr)
        fit_data = np.vstack((xdata,func_erf_3(xdata, *popt_3))).T
        np.savetxt('fit_data.txt', fit_data)
        ini_data = np.vstack((xdata,ydata)).T
        np.savetxt('ini_data.txt',ini_data )
        np.savetxt('fitting exponent_3 results.txt',popt_3 )
    if func_select == 2:
        index_special = np.where(ydata==np.max(ydata))[0]
        sigma=np.ones(len(xdata))
        sigma[index_special] = 1 
        #print(ydata[index_special-1])
        #param_bounds_2 = ([0,0,0,0,-np.inf],[1,1,np.inf,np.inf,np.inf])
        popt_2,pcov_2 = curve_fit(func_erf_2, xdata, ydata, maxfev=54000,###bounds=param_bounds_2,
                          sigma=sigma        )
        plt.plot(xdata, func_erf_2(xdata, *popt_2), 'r-', lw=2, 
             label='Fitting Curve' )
        #bax.plot(xdata, func_erf_2(xdata, *popt_2), 'r-', 
        #     label='Fitting Curve' )
        #delta=%5.6f'
        print(pcov_2)
        perr = np.sqrt(np.diag(pcov_2))
        print('fitting exponent_2 results:a1=%5.6f a2=%5.6f  b1=%5.6f b2=%5.6f u=%5.6f ' % tuple(popt_2),
              'perr=',perr)
        fit_data = np.vstack((xdata,func_erf_2(xdata, *popt_2))).T
        np.savetxt('fit_data.txt', fit_data)
        ini_data = np.vstack((xdata,ydata)).T
        np.savetxt('ini_data.txt',ini_data )
        np.savetxt('fitting exponent_3 results.txt',popt_2 )
        
################ ticks labels  
 
    plt.xlabel('Time (ps)',family='Times New Roman', labelpad=2 , fontsize = 15, fontweight='bold')
    plt.ylabel('Normalized $\Delta$OD',family='Times New Roman', labelpad= 2, fontsize = 15, fontweight='bold')
    #ax.set_xlim(-2,12)
    print('start_time=%f,end_time=%f'%(xdata[0],xdata[-1]))
    plt.savefig(figname, dpi=400)
'''    
    my_x = np.arange(0, 52, 10)                          
    my_xticks = np.array(['0', '10','20','30','40','50'])
    my_y = np.arange(0, 1.1, 0.2)
    my_yticks = np.array(['0', '0.2','0.4','0.6','0.8','1.0'])
    
    xmin=xdata[0];xmax=xdata[-1]
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(-0.1,1.1)
    ax.xaxis.set_ticks(my_x)
    ax.xaxis.set_tick_params(pad=3, labelsize=13 )
    ax.yaxis.set_ticks(my_y)
    ax.yaxis.set_tick_params(pad=2, labelsize=13)
    ax.set_xticklabels(my_xticks, )
    ax.set_yticklabels(my_yticks, )
    for i in range(len(my_x)):
        ax.get_xticklabels()[i].set_fontweight("bold")
    for i in range(len(my_y)):
        ax.get_yticklabels()[i].set_fontweight("bold")    
    
    plt.legend(loc='upper right', frameon=False, borderpad=0.3, handlelength=1.5, 
               prop={'weight' :'bold', 'family':'Times New Roman', 'size':'13'}  )
    plt.text(30,0.65, s='pump = 620nm\nprobe = 5259nm',family='Times New Roman',
             fontsize=13,fontweight='bold')
    #np.savetxt('opt_TA.txt',)
  
#save picture  
    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved."%figname)
    plt.show()
'''
   
def func_erf_3(x, a1, a2, a3, b1, b2, b3, u ):
    #u= -0.02                   #location
    delta_0 = 0.1500   #FWHM
    delta = delta_0/(2*math.sqrt(2*math.log(2)))
    #a1=4.841526 ;a2=0.622949;b1=0.013966; b2=0.953623
    #b3=100
    #print(delta)
    #gauss = 1/(delta*np.sqrt(2)*math.pi)*np.exp(-(x-u)/(2*delta**2))
    #exponent_3 = a1 * np.exp(-x / b1)+ a2 * np.exp(-x / b2)+ a3 * np.exp(-x / b3)     
    #return a1*np.exp(-x/b1)*np.exp((u+delta**2/b1/2)/b1)*(1+special.erf((x-(u+delta**2/b1))/(math.sqrt(2)*delta))) \
    #       + a2*np.exp(-x/b2)*np.exp((u+delta**2/b2/2)/b2)*(1+special.erf((x-(u+delta**2/b2))/(math.sqrt(2)*delta)))\
    #        +  a3*np.exp(-x/b3)*np.exp((u+delta**2/b3/2)/b3)*(1+special.erf((x-(u+delta**2/b3))/(math.sqrt(2)*delta))) 
   
    return a1/2*np.exp(-x/b1)*np.exp((u+delta**2/b1/2)/b1)*(1+special.erf((x-(u+delta**2/b1))/(math.sqrt(2)*delta))) \
           + a2/2*np.exp(-x/b2)*np.exp((u+delta**2/b2/2)/b2)*(1+special.erf((x-(u+delta**2/b2))/(math.sqrt(2)*delta)))\
            +  a3/2*np.exp(-x/b3)*np.exp((u+delta**2/b3/2)/b3)*(1+special.erf((x-(u+delta**2/b3))/(math.sqrt(2)*delta))) 

def func_erf_2(x, a1, a2, b1, b2, u ):
    #u= -0.02                   #location
    delta_0 = 0.080     #FWHM
    delta = delta_0/(2*math.sqrt(2*math.log(2)))
    
    return a1/2*np.exp(-x/b1)*np.exp((u+delta**2/b1/2)/b1)*(1+special.erf((x-(u+delta**2/b1))/(math.sqrt(2)*delta))) \
           + a2/2*np.exp(-x/b2)*np.exp((u+delta**2/b2/2)/b2)*(1+special.erf((x-(u+delta**2/b2))/(math.sqrt(2)*delta)))



def func_erf_4(x, a1, a2, a3, a4, b1, b2, b3, b4, u  ):
     #u= -0.02                   #location
     delta_0 = 0.120     #FWHM
     delta = delta_0/(2*math.sqrt(2*math.log(2)))
     #print(delta)
     #gauss = 1/(delta*np.sqrt(2)*math.pi)*np.exp(-(x-u)/(2*delta**2))
     #exponent_3 = a1 * np.exp(-x / b1)+ a2 * np.exp(-x / b2)+ a3 * np.exp(-x / b3)     
     return a1*np.exp(-x/b1)*np.exp((u+delta**2/b1/2)/b1)*(1+special.erf((x-(u+delta**2/b1))/(math.sqrt(2)*delta))) \
            + a2*np.exp(-x/b2)*np.exp((u+delta**2/b2/2)/b2)*(1+special.erf((x-(u+delta**2/b2))/(math.sqrt(2)*delta)))\
            + a3*np.exp(-x/b3)*np.exp((u+delta**2/b3/2)/b3)*(1+special.erf((x-(u+delta**2/b3))/(math.sqrt(2)*delta)))\
            + a4*np.exp(-x/b4)*np.exp((u+delta**2/b4/2)/b4)*(1+special.erf((x-(u+delta**2/b4))/(math.sqrt(2)*delta))) 
   

if __name__=='__main__':
    main()    