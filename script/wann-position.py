# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 09:22:04 2022

@author: daifulong
"""

import re
import numpy as np


struc_name= 'b-gdy-POSCAR.vasp' #change for different struc
file_name='wann-struc-b.txt'



with open(file_name,'r') as f:
    count= len(open(file_name,'r').readlines())
    wann_position = np.empty(shape=(count,3))
    n=0
    for line in f.readlines():
        string=re.findall(r'[(](.*?)[)]',line)
    
        string = string[0].split(',')
        
        x=np.float64(string[0])
        y=np.float64(string[1])
        z=np.float64(string[2])
    
        print(x,y,z)
        wann_position[n,0]=x
        wann_position[n,1]=y
        wann_position[n,2]=z
        n=n+1
        
np.savetxt('wann_position_b.txt',wann_position)
        


    
    
    
    
