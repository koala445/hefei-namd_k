# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 00:02:17 2024

@author: daifulong
"""

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
from scipy import interpolate

def main():
    plt.rc('font',family= 'arial')   #'Times New Roman',size=11)
    potsize = 6.33/40  # =y*x # 0.004 (y)* x (x=10)*n = 6.33um   # k =10000 multiplication
    R = -0.5  #V
    plot_PT_dR_image('PT_dR.txt',potsize)
    plot_PT_dR_line_image('PT_dR_line.txt',R)
    
    plot_PT_R_image('PT_R.txt',potsize)
    
    plot_PT_dR_R_image('PT_dR.txt','PT_R.txt',potsize)
    plot_PT_dR_R_line_image('PT_dR_R_line.txt')
    
def plot_PT_dR_image(file1,potsize,figname='PT_dR.svg'):
       
    z0 = np.array(np.loadtxt(file1),dtype=np.float32)/10000*1000000  # uV   
    z0 = z0[18:34,14:30].T                   #!!!change here to modify
    nb = len(z0)
    x0 = y0 = np.linspace(0,nb-1,nb) *potsize
    #z0 = interpolate.interp2d(x0,y0,z0,kind='cubic')
    
    print(x0,z0.shape)
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
    #cmap0 = plt.cm.jet
    #color0=cmap0(np.linspace(0,1,200))
    #color0_1=cmap0(np.linspace(0,1,80))
    #newcolorsgroup=np.vstack((color0_1[18:49],color0[125:180]))
    #newcmap=mpl.colors.ListedColormap(newcolorsgroup)        
    #newcmap = plt.cm.jet
    
    ###map1-2
    colorlist=[ (157/255,92/255,57/255),                         #black
        (254/255,147/255,118/255),      #red (242/255,106/255,17/255),
                (175/255,221/255,139/255),  # green (53/255,161/255,83/255)
               #(160/255,201/255,229/255), #white blue
               (66/255,146/255,201/255),
               (37/255,110/255,162/255) ][::-1]
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("testCmap", colors=colorlist, N=512)
    cmin=np.min(z0);  cmax=round(np.max(z0),1)    #round(np.max(z0),2) #cmax=0.54 #;
    norm = mpl.colors.Normalize(0,cmax)  
       
# plot contourf    
    X,Y = np.meshgrid(x0,y0)
    print(X.shape, cmax,round(np.max(z0),2) )
    plt.contourf(X,Y,z0,30,cmap=newcmap,norm=norm,linewidths=0)  
    plt.contourf(X,Y,z0,30,cmap=newcmap,norm=norm,linewidths=0)   # plot twice to rase contour
    #bounds= np.round(np.arange(0,cmax*1.4,cmax/4),decimals=2) #bounds= np.linspace(cmin,cmax,5)
    bounds= np.round(np.linspace(0,cmax,5),decimals=1)
    ##########colorbar    
    cb=plt.colorbar(
            mpl.cm.ScalarMappable(cmap=newcmap,norm=norm),
            ax=ax,
            #extend='min',
            ticks=bounds,
            spacing='proportional',
            orientation='vertical',
            label='deltaR (uV)',
        )
    font = {'family' : 'arial',
        'color'  : 'black',
        'weight' : 'bold',
        'size'   : 14,
        }
    cb.ax.tick_params(labelsize=12)
    cb.formatter.set_powerlimits((2, 0))
    for i in range(len(bounds)):
        cb.ax.get_yticklabels()[i].set_fontweight("bold")
    cb.ax.set_yticklabels(bounds )
    cb.set_label('ΔR (\u03BCV)',fontdict=font, labelpad = 6 ) #设置colorbar的标签字体及其大小
######### ticks labels   
    plt.xlabel('X (\u03BCm)',family='arial', labelpad=2 , fontsize = 14, fontweight='bold')
    plt.ylabel('Y (\u03BCm)',family='arial', labelpad= 2, fontsize = 14, fontweight='bold')
    
    #my_x = np.arange(0, 1.51, 0.25)      #for gdy    ##initial level     
    my_x = np.round(np.linspace(0, (nb-1)*potsize, 5),decimals=2)                  
    #my_xticks = np.array(['0.25', '0.5','0.75','1.0','1.25'])
    
    #my_x = np.arange(0.25, 1.5, 0.3)      #for gr    ##initial level                       

    my_xticks = np.array(['0','0.59', '1.19','1.78','2.37'])
    my_yticks = np.array(['0','0.59', '1.19','1.78','2.37'])
    
    xmin=np.min(x0); xmax=np.max(x0)
    ax.set_xlim(0,(nb-1)*potsize)
    ax.set_ylim(0,(nb-1)*potsize)
    #ax.set_xlim(0,1.5)
    #ax.set_ylim(-0.05,1.55)
    ax.xaxis.set_ticks(my_x)
    ax.xaxis.set_tick_params(pad=2.7, labelsize=13 )
    ax.yaxis.set_ticks(my_x)
    ax.yaxis.set_tick_params(pad=2, labelsize=13)
    ax.set_xticklabels(my_xticks, )
    ax.set_yticklabels(my_yticks, )
    for i in range(len(my_x)):
        ax.get_xticklabels()[i].set_fontweight("bold")
        ax.get_yticklabels()[i].set_fontweight("bold")

###########plot line
    nx = 0 ; ny = 6     ####!!!!!!  change here   !!!!!!
    if nx==0 and ny != 0:
        z1 = z0[ny,nx:]
    if ny ==0 and nx != 0: 
        z1 = z0[ny:,nx]
    
    
    if nx == 0 and ny != 0:
        x1 = x0 ; y1 = ([ny * potsize]*nb)   # ---
        plt.plot(x1,y1,color='black',linewidth=1.5, linestyle='--',dashes=(3,2),alpha=0.8)   

    if ny == 0 and nx != 0:
        x1 = ([nx * potsize]*nb) ; y1 = y0  # |||||
        plt.plot(x1,y1,color='black',linewidth=1.5, linestyle='--',dashes=(3,2),alpha=0.8)   
    print(y1,z1)
#save PT_image picture  
    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved."%figname)
    plt.show()
    
####### plot PT_line_picture
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
    figname_line = 'PT_dR_line.png'
    if nx == 0 and ny != 0:
        x_line = x1
        plt.plot(x_line,z1,color='black',linewidth=1.5, linestyle='-',dashes=(3,2),alpha=0.8)       
    if ny == 0 and nx != 0:
        x_line = y1
        plt.plot(x_line,z1,color='black',linewidth=1.5, linestyle='-',dashes=(3,2),alpha=0.8)       
#save PT_line_image picture  
    PT_line = np.vstack((x_line, z1))
    np.savetxt('PT_dR_line.txt',PT_line)
    plt.tight_layout()
    plt.savefig(figname_line, dpi=400)
    print("\n%s has been saved."%figname_line)
    plt.show()
    
def plot_PT_dR_line_image(file1,R):
    data_all = np.array(np.loadtxt(file1),dtype=np.float32)  
    #z0 = z0[5:20,4:19]
    n_s=0
    n_e=60
    x = data_all[0,n_s:n_e]
    x = x - x[0] 
    y = data_all[1,n_s:n_e]/R
####### plot PT_line_picture
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
    figname_line ='PT_dR_line_select.png'    
    plt.plot(x,y,color='black',linewidth=1.5, linestyle='-',dashes=(3,2),alpha=0.8)           
    
#save PT_line_R_image picture  
    PT_dR_line_select = np.vstack((x, y))
    np.savetxt('PT_dR_line_select.txt',PT_dR_line_select)
    plt.tight_layout()
    plt.savefig(figname_line, dpi=400)
    print("\n PT_dR_line_select has been saved.")
    plt.show()    

def plot_PT_R_image(file1,potsize):
           
    z0 = np.array(np.loadtxt(file1),dtype=np.float32)  # V   
    z0 = z0[18:34,14:30].T                  #!!!!
    z0 = np.round(z0,decimals=2)
    nb = len(z0)
    x0 = y0 = np.linspace(0,nb-1,nb) *potsize
    
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
    #cmap0 = plt.cm.jet
    #color0=cmap0(np.linspace(0,1,200))
    #color0_1=cmap0(np.linspace(0,1,80))
    #newcolorsgroup=np.vstack((color0_1[18:49],color0[125:180]))
    #newcmap=mpl.colors.ListedColormap(newcolorsgroup)        
    #newcmap = plt.cm.jet
    
    ####map1-2
    colorlist=[ (157/255,92/255,57/255),                         #black
        (254/255,147/255,118/255),      #red (242/255,106/255,17/255),
                (175/255,221/255,139/255),  # green (53/255,161/255,83/255)
               #(160/255,201/255,229/255), #white blue
               (66/255,146/255,201/255),
               (37/255,110/255,162/255) ][::-1]
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("testCmap", colors=colorlist, N=512)
    cmin=np.min(z0);cmax=np.max(z0)    #round(np.max(z0),2) #cmax=0.54 #;
    #cmin=0    ####!!! change here
    
    norm = mpl.colors.Normalize(cmin,cmax)
       
# plot contourf    
    X,Y = np.meshgrid(x0,y0)
    print( cmin, cmax, )
    plt.contourf(X,Y,z0,30,cmap=newcmap,norm=norm)
    plt.contourf(X,Y,z0,30,cmap=newcmap,norm=norm)
    #bounds= np.round(np.arange(0,cmax,cmax/5),decimals=1) #bounds= np.linspace(cmin,cmax,5)
    bounds= np.linspace(cmin,cmax,5)
    print(bounds)
    ##########colorbar    
    cb=plt.colorbar(
            mpl.cm.ScalarMappable(cmap=newcmap,norm=norm),
            ax=ax,
            #extend='min',
            ticks=bounds,
            spacing='proportional',
            orientation='vertical',
            label='R (V)',
        )
    font = {'family' : 'arial',
        'color'  : 'black',
        'weight' : 'bold',
        'size'   : 14,
        }
    cb.ax.tick_params(labelsize=12)
    for i in range(len(bounds)):
        cb.ax.get_yticklabels()[i].set_fontweight("bold")
    #cb.ax.set_yticklabels(bounds)
    cb.set_label('R (V)',fontdict=font, labelpad = 6 ) #设置colorbar的标签字体及其大小
######### ticks labels   
    plt.xlabel('X (\u03BCm)',family='arial', labelpad=2 , fontsize = 14, fontweight='bold')
    plt.ylabel('Y (\u03BCm)',family='arial', labelpad= 2, fontsize = 14, fontweight='bold')
    
    #my_x = np.arange(0, 1.51, 0.25)      #for gdy    ##initial level     
    my_x = np.round(np.linspace(0, (nb-1)*potsize, 5),decimals=2)                  
    #my_xticks = np.array(['0.25', '0.5','0.75','1.0','1.25'])
    
    #my_x = np.arange(0.25, 1.5, 0.3)      #for gr    ##initial level                       
    my_xticks = np.array(['0','0.59', '1.19','1.78','2.37'])
    my_yticks = np.array(['0','0.59', '1.19','1.78','2.37'])
    
    xmin=np.min(x0); xmax=np.max(x0)
    ax.set_xlim(0,(nb-1)*potsize)
    ax.set_ylim(0,(nb-1)*potsize)
    #ax.set_xlim(0,1.5)
    #ax.set_ylim(-0.05,1.55)
    ax.xaxis.set_ticks(my_x)
    ax.xaxis.set_tick_params(pad=2.7, labelsize=13 )
    ax.yaxis.set_ticks(my_x)
    ax.yaxis.set_tick_params(pad=2, labelsize=13)
    ax.set_xticklabels(my_xticks, )
    ax.set_yticklabels(my_yticks, )
    for i in range(len(my_x)):
        ax.get_xticklabels()[i].set_fontweight("bold")
        ax.get_yticklabels()[i].set_fontweight("bold")

    figname = 'PT_R.svg'
#save PT_image picture  
    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved."%figname)
    plt.show()


    
def plot_PT_dR_R_image(file_dR,file_R,potsize,figname='PT_dR_R.svg'):
    z1 = np.array(np.loadtxt(file_dR),dtype=np.float32)   # uV    
    z2 = np.array(np.loadtxt(file_R),dtype=np.float32)  # V    
    z0 = z1/z2
    z0 = z0[18:34,14:30].T                       #!!!change here to modify
    nb = len(z0)
    x0 = y0 = np.linspace(0,nb-1,nb) *potsize
    
    print(x0)
    fig = plt.figure()
    
    ax=plt.axes()
    ax.spines['bottom'].set_linewidth(1)                 #set frame
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['top'].set_linewidth(1)
    figsize_x = 4.8
    figsize_y = 3.8 # in inches
    fig.set_size_inches(figsize_x, figsize_y)
    
# colormap
    #cmap0 = plt.cm.jet
    #color0=cmap0(np.linspace(0,1,200))
    #color0_1=cmap0(np.linspace(0,1,80))
    #newcolorsgroup=np.vstack((color0_1[18:49],color0[125:180]))
    #newcmap=mpl.colors.ListedColormap(newcolorsgroup)        
    #newcmap = plt.cm.jet
    ####### map1-2
    colorlist=[ (157/255,92/255,57/255),                         #black
        (254/255,147/255,118/255),      #red (242/255,106/255,17/255),
                (175/255,221/255,139/255),  # green (53/255,161/255,83/255)
               #(160/255,201/255,229/255), #white blue
               (66/255,146/255,201/255),
               (37/255,110/255,162/255) ][::-1]
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("testCmap", colors=colorlist, N=512)
    cmin=round(np.min(z0),1); cmax =round(np.max(z0),1) #2 # cmax=round(np.max(z0),1)    #round(np.max(z0),2) #cmax=0.54 #;
    norm = mpl.colors.Normalize(0,2)
       
# plot contourf    
    X,Y = np.meshgrid(x0,y0)
    print(X.shape, cmax,round(np.max(z0),2) )
    plt.contourf(X,Y,z0,30,cmap=newcmap,norm=norm) 
    plt.contourf(X,Y,z0,30,cmap=newcmap,norm=norm)
    #bounds= np.round(np.arange(0,cmax*1.1,cmax/4),decimals=2) #bounds= np.linspace(cmin,cmax,5)
    bounds= np.round(np.linspace(0,2,5),decimals=2) #
    ##########colorbar    
    cb=plt.colorbar(
            mpl.cm.ScalarMappable(cmap=newcmap,norm=norm),
            ax=ax,
            #extend='min',
            ticks=bounds,
            #spacing='proportional',
            orientation='vertical',##orientation='vertical',horizontal
            label='ΔR/R (uV)',
        )
    font = {'family' : 'arial',
        'color'  : 'black',
        'weight' : 'bold',
        'size'   : 14,
        }
    cb.ax.tick_params(labelsize=12)
    for i in range(len(bounds)):
        cb.ax.get_yticklabels()[i].set_fontweight("bold")  ##cb.ax.get_xticklabels()[i].set_fontweight("bold")
    cb.ax.set_yticklabels(bounds ) #cb.ax.set_xticklabels(bounds )
    cb.set_label('ΔR/R',fontdict=font, labelpad = 6 ) #设置colorbar的标签字体及其大小
######### ticks labels   
    plt.xlabel('X (\u03BCm)',family='arial', labelpad=2 , fontsize = 14, fontweight='bold')
    plt.ylabel('Y (\u03BCm)',family='arial', labelpad= 2, fontsize = 14, fontweight='bold')
    
    #my_x = np.arange(0, 1.51, 0.25)      #for gdy    ##initial level     
    my_x = np.round(np.linspace(0, (nb-1)*potsize, 5),decimals=2)                  
    #my_xticks = np.array(['0.25', '0.5','0.75','1.0','1.25'])
    
    #my_x = np.arange(0.25, 1.5, 0.3)      #for gr    ##initial level                       
    my_xticks = np.array(['0','0.59', '1.19','1.78','2.37'])
    my_yticks = np.array(['0','0.59', '1.19','1.78','2.37'])
    
    xmin=np.min(x0); xmax=np.max(x0)
    ax.set_xlim(0,(nb-1)*potsize)
    ax.set_ylim(0,(nb-1)*potsize)
    #ax.set_xlim(0,1.5)
    #ax.set_ylim(-0.05,1.55)
    ax.xaxis.set_ticks(my_x)
    ax.xaxis.set_tick_params(pad=2.7, labelsize=13 )
    ax.yaxis.set_ticks(my_x)
    ax.yaxis.set_tick_params(pad=2, labelsize=13)
    ax.set_xticklabels(my_xticks, )
    ax.set_yticklabels(my_yticks, )
    for i in range(len(my_x)):
        ax.get_xticklabels()[i].set_fontweight("bold")
        ax.get_yticklabels()[i].set_fontweight("bold")
    
    #plt.annotate('×10$^{-4}$',xy=(6.5,6.5))
    plt.text(0.79, 0.908, '$× 10^{-4}$', fontsize=9,weight='bold',family='arial',
             transform=plt.gcf().transFigure, color='black')
    #ax.set_xticks([]), ax.set_yticks([])
    #ax.xaxis.set_visible(False), ax.yaxis.set_visible(False)
###########plot line
    nx =0 ; ny = 6     ####!!!!!!  change here   !!!!!!  21 for GR3-5 8 for GDY5
    if nx==0 and ny != 0:
        z1 = z0[ny,nx:]
    if ny ==0 and nx != 0: 
        z1 = z0[ny:,nx]
    
    
    if nx == 0 and ny != 0:
        x1 = x0 ; y1 = ([ny * potsize]*nb)   # ---
        plt.plot(x1,y1,color='black',linewidth=1.5, linestyle='--',dashes=(3,2),alpha=0.8)   

    if ny == 0 and nx != 0:
        x1 = ([nx * potsize]*nb) ; y1 = y0  # |||||
        plt.plot(x1,y1,color='black',linewidth=1.5, linestyle='--',dashes=(3,2),alpha=0.8)   
    print(y1,z1)
#save PT_image picture  
    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved."%figname)
    plt.show()
    
####### plot PT_line_picture
    fig = plt.figure()
    
    ax=plt.axes()
    ax.spines['bottom'].set_linewidth(1)                 #set frame
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['top'].set_linewidth(1)
    figsize_x = 4.8
    figsize_y = 3.6 # in inches
    fig.set_size_inches(figsize_x, figsize_y)
    figname_line = 'PT_line_dR_R.png'
    if nx == 0 and ny != 0:
        x_line = x1
        plt.plot(x_line,z1,color='black',linewidth=1.5, linestyle='-',dashes=(3,2),alpha=0.8)       
    if ny == 0 and nx != 0:
        x_line = y1
        plt.plot(x_line,z1,color='black',linewidth=1.5, linestyle='-',dashes=(3,2),alpha=0.8)       
#save PT_line_image picture  
    PT_line = np.vstack((x_line, z1))
    np.savetxt('PT_dR_R_line.txt',PT_line)
    plt.tight_layout()
    plt.savefig(figname_line, dpi=400)
    print("\n%s has been saved."%figname_line)
    plt.show()

def plot_PT_dR_R_line_image(file1):
    data_all = np.array(np.loadtxt(file1),dtype=np.float32)  
    n_s=0
    n_e=100
    x0 = data_all[0,n_s:n_e]
    x0 = x0 - x0[0] 
    y0 = data_all[1,n_s:n_e]

####### interpolate
    xp = np.linspace(x0[0],x0[-1],100)
    yp = interpolate.interp1d(x0,y0,kind='quadratic')     #ks = ['zero', 'slinear', 'quadratic', 'cubic']
    print(yp)
####### plot PT_line_picture
    fig = plt.figure()
    
    ax=plt.axes()
    ax.spines['bottom'].set_linewidth(1)                 #set frame
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['top'].set_linewidth(1)
    figsize_x = 4.8
    figsize_y = 3.6 # in inches
    fig.set_size_inches(figsize_x, figsize_y)
    figname_line = 'PT_line_dR_R_select.png'    
    plt.plot(x0,y0,color='black',linewidth=1.5, linestyle='-',dashes=(3,2),alpha=0.8,label='origin')
    plt.plot(xp,yp(xp),label='inter')           
    
#save PT_line_R_image picture  
    PT_dR_R_line = np.vstack((x0, y0))
    PT_dR_R_interpolate_line = np.vstack((xp, yp(xp)))
    np.savetxt('PT_dR_R_interpolate_line_select.txt',PT_dR_R_interpolate_line.T)
    np.savetxt('PT_dR_R_line_select.txt',PT_dR_R_line.T)
    plt.tight_layout()
    plt.savefig(figname_line, dpi=400)
    print("\n PT_line_R has been saved.")
    plt.show()        
if __name__=='__main__':
    main()