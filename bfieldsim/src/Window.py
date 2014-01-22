# -*- coding: utf-8 -*-
"""
Created on Sun May 19 15:35:49 2013

TODO:
add interpolation/compression features
separate out the specific window definition

@author: Will
"""

#from collections import namedtuple
import numpy as np

wmin = -50 #um
wmax = 50 #um
#wmin = -13 * 512
#wmax = 13 * 512
wnum = 2048

wmin = wmin * 1e-6 #m
wmax = wmax * 1e-6 #m



# Window = namedtuple('Window', ['min','max','num_cells', 'window'])

# window = Window(wmin, wmax, wnum, np.linspace(wmin, wmax, wnum))

class Point():
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Window():
   def __init__(self, xmin, xmax, num_cells):
       self.xmin = xmin
       self.xmax = xmax
       self.num_cells = num_cells
       self.window = np.linspace(xmin, xmax, num_cells)
       self.cell_size = (self.xmax - self.xmin) / self.num_cells
	   
window = Window(wmin, wmax, wnum)

class Window_2D():
    def __init__(self, x_0, y_0, x_f, y_f, x_num, y_num):
        self.xmin = x_0
        self.ymin = y_0
        self.xmax = x_f
        self.ymax = y_f
        self.x_num = x_num
        self.y_num = y_num
        self.num_cells = x_num*y_num
        self.window = np.array([Point(x,y) for x in np.linspace(self.xmin, self.xmax, self.x_num)
                                  for y in np.linspace(self.ymin, self.ymax, self.y_num)])
        getxf = lambda pt: pt.x
        getyf = lambda pt: pt.y
        getxv = np.vectorize(getxf)
        getyv = np.vectorize(getyf)
        self.X = getxv(self.window)
        self.Y = getyv(self.window)
        self.window = self.window.reshape((y_num, x_num))
#        print self.window
#        print 'window shape'
#        print self.window.shape
#        print y_num
#        print x_num
        self.X = self.X.reshape((y_num, x_num))
#        print self.X
        self.Y = self.Y.reshape((y_num, x_num))