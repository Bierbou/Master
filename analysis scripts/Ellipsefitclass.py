# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 13:41:02 2015

@author: ga56pan
"""


import math 
import numpy as np
from numpy.linalg import eig, inv


class Ellipsefit(object):
    def __init__(self):
        self.data = []
    def rect(self,r, theta):
        """theta in degrees

        returns tuple; (float, float); (x,y)
        """
        x = r * math.cos(math.radians(theta))
        y = r * math.sin(math.radians(theta))
        return x,y
    
    def fit_ellipse(self,x_input,y_input):
        x = x_input[:,np.newaxis]
        y = y_input[:,np.newaxis]
        D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
        S = np.dot(D.T,D)
        C = np.zeros([6,6])
        C[0,2] = C[2,0] = 2; C[1,1] = -1
        E, V =  eig(np.dot(inv(S), C))
        n = np.argmax(np.abs(E))
        a = V[:,n]
        return a
        
        
    def ellipse_center(self,a):
        a,b,c,d,f = a[0], a[1]/2, a[2], a[3]/2, a[4]/2 
        num = b*b-a*c
        x0=(c*d-b*f)/num
        y0=(a*f-b*d)/num
        return np.array([x0,y0])

   
    def ellipse_semi_axis_length(self, a ):
        a,b,c,d,f,g = a[0], a[1]/2, a[2], a[3]/2, a[4]/2, a[5]
        up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
        down1=(b**2-a*c)*(np.sqrt((a-c)**2 +4*b**2)-(a+c))
        down2=(b**2-a*c)*(-np.sqrt((a-c)**2 +4*b**2)-(a+c))
        res1=np.sqrt(up/down1)
        res2=np.sqrt(up/down2)
        return np.array([res1, res2])


    def ellipse_angle_of_rotation(self, a ):
        a,b,c = a[0], a[1]/2, a[2]
        return 0.5*np.arctan(2*b/(a-c))


    def ellipse_angle_of_rotation2(self, a ):
        a,b,c = a[0],a[1]/2, a[2],
        if b == 0 and a<c:
            return 0            
        if b == 0 and a>c:
            return np.pi/2
        else:
            if a < c:
                return np.arctan(2*b/(a-c))/2
            else:
                return np.pi/2 + np.arctan(2*b/(a-c))/2



















































