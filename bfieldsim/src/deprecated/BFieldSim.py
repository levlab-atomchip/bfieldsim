# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 22:36:44 2013

Magnetic Trap Simulator, Python edition
This is a Python port of the most up-to-date simulator I had written in 
MatLab as of 7/4/10. It simulates
the combined fields from the macrowires and microwires.

@author: Will
"""



## Trap parameters
# The script needs:
#       -the bias field (B_ybias) in Gauss
#       -the central wire length (l) in microns
#       -the trap height (h) in microns
#       -the wire current (I) in amps
#
# Formulae are provided for deriving the necessary bias field to realize a
# particular trap height. Trap length must be provided, plus any two of
# (B_ybias, I, and h).


# Axes are defined so that z is up (out of the chip plane), x is along the
# central axis of the z trap ( = ODT axis) and y follows the right-hand
# rule.

#The z ground plane is the top of the MACOR surface.

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from AtomChip import *
import numpy as np
import acwires
import ThinWireSpecs as WireSpecs
import matplotlib.pyplot as plt
from math import pi, sqrt

clear = "\n"*100

#Simulation Mode
fieldplot = 1
plotxz = 0

## Plot parameters
plot_on=1
#plot freq (1) or G/cm^2 (0)?
#plot gradient of x (1) or y (0)?
freq=1 
xplot=0
# all distance in units of meters
# resolution=a/4000 # resolution of the plots, meters
# plotleft=-a_ax/50 # boundary of plots, meters
# plotright=a_ax/50
# plottop = a/50
# plotbottom = -a/50

#resolution=a/5 # resolution of the plots, meters
#plotleft=-a_ax # boundary of plots, meters
#plotright=a_ax
#plottop = a
#plotbottom = -a

resolution=0.0001 # resolution of the plots, meters
plotleft=-0.005 # boundary of plots, meters
plotright=0.005
plottop = 0.0025
plotbottom = -0.0025

#resolution=0.1 # resolution of the plots, meters
#plotleft=-5 # boundary of plots, meters
#plotright=5
#plottop = 2
#plotbottom = -2
# 
# resolution=a/400 # resolution of the plots, meters
# plotleft=-a_ax/8 # boundary of plots, meters
# plotright=a_ax/8
# plottop = a/8
# plotbottom = -a/8
# 
# resolution=a/2000 # resolution of the plots, meters
# plotleft= -l_med/4 # boundary of plots, meters
# plotright= l_med/4
# plottop = a/40
# plotbottom = -a/40

# resolution=a/100000 # resolution of the plots, meters
# plotleft=-l_med/500 # boundary of plots, meters
# plotright=l_med/500
# plottop = a/500
# plotbottom = -a/500

# For traps along y-axis
# plotleft=-10*a # boundary of plots, microns
# plotright=10*a
# plottop = 3*a
# plotbottom = -a

## Find trap

if plotxz == 1:
    n = 100 #how many points to check
#    z_range = np.zeros(2,n)
    z = 1e-3
    z_range = np.linspace(z*.8,z*1.5,n) #meters
else:
    n = 100
#    z_range = np.zeros(2,n)
    z_range = np.linspace(z*.9,z*1.1,n) #meters

z_spacing = z_range[1]-z_range[0] #meters


#for ii = 1:n
B_tot_trap = np.array([])
B_bias = np.array((B_xbias, B_ybias, B_zbias))
#print B_bias
x_trap = 0
y_trap = 0

for ii in xrange(n):
    tot_field = np.array((0.0,0.0,0.0))
    for wire in WireSpecs.allwires:
        if wire.current != 0:
#            print wire.bfieldcalc(x_trap, y_trap, z_range[ii])
            this_field = wire.bfieldcalc(x_trap, y_trap, z_range[ii])
#            print this_field
            tot_field = tot_field + this_field
#            print wire.name
#            print wire.current
#            print 1e6*tot_field
#    print tot_field
    tot_field_norm = np.linalg.norm(tot_field + B_bias)
#    print tot_field_norm
    B_tot_trap = np.append(B_tot_trap,tot_field_norm)
#    z_range(2,ii) = B_tot_center #in Tesla



#correct for gravity
#gravity_equivalent_field = (m_Rb87 * g)/(mu_B) * z_range
#B_tot_trap = B_tot_trap - gravity_equivalent_field
#print B_tot_trap
#print B_tot_trap*1e4
min_B_tot = np.min(B_tot_trap)
min_index = np.argmin(B_tot_trap)

trap_height = z_range[min_index]
z_trap = trap_height
#if fieldplot == 0:
print 'Trap Height is %2.0f'%1e6*trap_height
#trap_height = h*1e-6
#trap_height = 3e-3

plt.plot(z_range*1e3, B_tot_trap*1e4)
plt.xlabel('Z axis (mm)') #Standard axis labelling
plt.ylabel('Effective |B|, (G)')
plt.show()

## Field calculations

if fieldplot == 1:
    
    # permute coordinates to create slices in other planes.
    
    if plotxz == 1:
        
        
        x, z = np.meshgrid(np.arange(plotleft, plotright, resolution), z_range)
        B_tot = np.zeros(x.shape)
        for coords in np.ndenumerate(x):
            tot_field = np.array((0.0,0.0,0.0))
            for wire in WireSpecs.allwires:
                this_field = wire.bfieldcalc(x[coords[0]], 
                                             y_trap, 
                                             z[coords[0]])
#                print wire.name
#                print this_field*1e4
                tot_field += this_field
#                print tot_field
#            print coords[1]
#            print z[coords[0]]
#            print tot_field * 1e4
            tot_field_norm = np.linalg.norm(tot_field + B_bias)
            B_tot[coords[0]] = tot_field_norm
#            B_tot[coords[0]] = (tot_field_norm 
#                                - (m_Rb87 * g)/(mu_B) * z[coords][1])B_tot[coords[0] = tot_field_norm
#        B_tot_grav_corr = B_tot - (m_Rb87 * g)/(mu_B) * z
        print B_tot
        
        fig = plt.figure()
        ax = fig.gca(projection='3d')

#        ax.plot_trisurf(x*1e3,z*1e3,B_tot*1e4, cmap=cm.jet, linewidth=0.2)
        ax.plot_surface(x*1e3,z*1e3,B_tot*1e4, rstride=8, cstride=8, alpha=0.3)
        plt.xlabel('X axis (mm)')
        plt.ylabel('Z axis (mm)') #standard axis labelling
        ax.set_zlabel('B field (G)')
        plt.show()
#        meshc(x*1e3,z*1e3,B_tot*1e4)
#        
#        plt.show()
        
    else:
        x, y =np.meshgrid(np.arange(plotleft, plotright, resolution), 
                          np.arange(plotbottom, plottop, resolution))
#        print x.shape
        B_tot = np.zeros(x.shape)
        npoints = B_tot.size
        npoints_complete = 0
        for coords in np.ndenumerate(x):    
#            print coords
#            print x[coords][1]
#            print coords[1]
#            print y[coords][1]
#            print y[coords[0]]
            tot_field = np.array((0.0,0.0,0.0))
            i = 1
            num_wires = len(WireSpecs.allwires)
            for wire in WireSpecs.allwires:
#                print clear
#                print "%d of %d points complete\n"%(npoints_complete, npoints)
#                print "%2.3f percent complete"%(100.0*npoints_complete/npoints)
#                print "Wire %d of %d"%(i, num_wires)
                i += 1
                this_field = wire.bfieldcalc(x[coords[0]], 
                                              y[coords[0]], 
                                                z_trap)
                tot_field += this_field
#            print tot_field
            tot_field_norm = np.linalg.norm(tot_field + B_bias)
            B_tot[coords[0]] = tot_field_norm
#            B_tot[coords[0]] = (tot_field_norm 
#                                - (m_Rb87 * g)/(mu_B) * z_trap) 
                                #gravity correction
            npoints_complete += 1
#        print B_tot
        
        
        
        #find the trap
        min_field = np.min(B_tot)
        min_ind = np.unravel_index(B_tot.argmin(), B_tot.shape)
        x_ind = min_ind[0]
        y_ind = min_ind[1]
        x_trap = x[x_ind, 0]
        y_trap = y[0, y_ind]
        
        GradBx,GradBy = np.gradient(B_tot,resolution, resolution)
        GGradBx, GGradByx = np.gradient(GradBx,resolution, resolution)
        GGradBxy,GGradBy = np.gradient(GradBy,resolution, resolution)
        ## Extract and print frequency
        # The first part attempts to extract the transverse and longitudinal
        # frequencies by fitting to a paraboloid and extracting the principal values
        traploc = [x_ind, y_ind]
#        print traploc
#        print GGradBy
        a = abs(GGradBy[(traploc[0],traploc[1])])
        b = abs(.5*(GGradByx[(traploc[0],traploc[1])] 
                + GGradBxy[(traploc[0],traploc[1])]))
        c = abs(GGradBx[(traploc[0],traploc[1])])
#        print a
#        print b
#        print c
        Principal_Matrix = np.array([[a, b], [b, c]])
        values, vectors = np.linalg.eig(Principal_Matrix)
        freqs = []
        for val in values:
            freqs.append((1.2754)*sqrt(abs(val))) #Hz
        f_transverse = max(freqs)
        f_longitudinal = min(freqs)
        omega_transverse = 2*pi*f_transverse
        omega_longitudinal = 2*pi*f_longitudinal
        
        # Try to extract f_transverse by fitting a parabola to z_range
        grad_z = np.gradient(B_tot_trap, z_spacing)
        ggrad_z = np.gradient(grad_z, z_spacing)
        ddBdzz = ggrad_z[min_index]
        f_z = 1.2754*sqrt(abs(ddBdzz))
        omega_z = 2*pi*f_z
        
        B_ybias = B_ybias*1e4 #Gauss
        B_tot_center = B_tot[(traploc[0],traploc[1])]*1e4 #Gauss
        trap_height = trap_height*1e6 #microns
        
        ## Other Trap Parameters
        ODT_temp = 1e-6 #Kelvin, a guess
        f_rad_ODT = 40 #Hz, extracted from Lin
        f_ax_ODT = .3 #Hz, extracted from Lin
        
        AC_trap_temp = (((f_z**2 * f_longitudinal)
                        /(f_rad_ODT**2 * f_ax_ODT))**(1/3)
                        * ODT_temp)
#        cloud_length = (2*1e6*(k_B * AC_trap_temp 
#                        / (m_Rb87 * omega_longitudinal**2))**.5) #microns
        cloud_width = (2*1e6*(k_B * AC_trap_temp 
                        / (m_Rb87 * omega_z**2))**.5) #microns
        
        ## Plotting
        if plot_on==1:
#            meshc(x*1e3,y*1e3,B_tot*1e4)
#            plt.xlabel('X axis (mm)')
#            plt.ylabel('Y axis (mm)') #standard axis labelling
#            plt.zlabel('B field (G)')
#            plt.show()
            
            fig = plt.figure()
            ax = fig.gca(projection='3d')

#            ax.plot_trisurf(x*1e3,y*1e3,B_tot*1e4, cmap=cm.jet, linewidth=0.2)
            ax.plot_surface(x*1e3,y*1e3,B_tot*1e4, 
                            rstride=8, 
                            cstride=8, 
                            alpha=0.3)
            plt.xlabel('X axis (mm)')
            plt.ylabel('Y axis (mm)') #standard axis labelling
#            plt.zlabel('B field (G)')
            ax.set_zlabel('B field (G)')
            plt.show()
            
            #     figure(14)
            #     plot(z_range(1,:)*1e3, ggrad_z)
            #     xlabel('Z axis (\mum)') #Standard axis labelling
            #     ylabel('2nd Derivative of Field Strength')
