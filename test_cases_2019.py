import os

import sys

cur_path = os.path.abspath(__file__)

cur_path = os.path.dirname(cur_path)

sys.path.append(os.path.join(cur_path, ".."))

import scipy.ndimage.interpolation as sni

import h5py

import pickle

from sklearn.gaussian_process import GaussianProcessRegressor

from mpl_toolkits.mplot3d.art3d import Poly3DCollection


import mayavi.mlab as mlab

from tvtk.api import tvtk, write_data

import matplotlib.pyplot as plt

import pydem

from pydem import *

from math import exp

def desirability(y):

    weights = np.asarray([5,3,3,3], dtype=float)

    weights /= np.sum(weights)

    mins = [0.395, 2.15, 140, 0.048]

    maxes = [1.28, 3.28, 165, 0.107]

    minimize = [False, True, True, True]

    # minimize = [True, False, False, False]

    desirability = np.ones(y.shape[0])

    for i in range(len(weights)):

        if minimize[i]:

            temp = (maxes[i] - y[:,i])/(maxes[i] - mins[i])

        else:

            temp = (y[:,i] - mins[i])/(maxes[i] - mins[i])

        temp[temp < 0] = 0

        temp[temp > 1] = 1

        temp **= weights[i]

        desirability *= temp

    # desirability = np.power(desirability, 1.0/np.sum(weights))

    return desirability

    

def responses(all_points):

    vars = sympy.symbols('TON TOFF IP SV WF WT')

    TON = vars[0]

    TOFF= vars[1]

    IP  = vars[2]

    SV  = vars[3]

    WF  = vars[4]

    WT  = vars[5]

    VAR = [TON,TOFF,IP,SV,WF,WT]

    

    f = [0] * 4

    fmax = [0] * 4

    fmin = [0] * 4



    f = [0]*4

    # Machining Rate (MR)

    f1 = 11.06243 + (0.021958 * TON) - (0.4023 * TOFF) - (0.039917 * IP) - (0.045079 * SV) + (3.06872 * 10.0**-3.0 * TOFF**2.0) + (1.79134 * 10.0**-5.0 * IP**2.0) + (2.30469 * 10.0**-4.0 * TON * IP) + (1.96875 * 10.0**-4.0 * TOFF * IP) + (7.39583 * 10.0**-4.0 * TOFF * SV)

    f[0] = SymbolicFunction(f1, VAR)

    # Surface Roughness (SR)

    f2 = 10.87640 - (0.038333 * TON) - (0.1284 * TOFF) - (0.072771 * IP) - (0.008667 * SV) + (0.1247 * WF) - (0.000365432 * WT) + (9.72222 * 10.0**-4.0 * TOFF**2.0) - (8.75000 * 10.0**-3.0 * WF**2.0)  + (1.91358 * 10.0**-7.0 * WT**2.0) + (6.48438 * 10.0**-4.0 * TON * IP)

    f[1] = SymbolicFunction(f2, VAR)

    # Dimensional Deviation (DD)

    f3 = 150.52609+ (1.43750  * TON) - (4.0869 * TOFF) - (0.23551  * IP) - (0.46806  * SV) - (1.8670 * WF) - (0.015278    * WT) + (0.034827 * TOFF**2.0) + (8.23864 * 10.0**-4.0 * IP**2.0) + (0.12542 * WF**2.0) + (2.77778 * 10.0**-4.0 * SV * WT)

    f[2] = SymbolicFunction(f3, VAR)

    # Wire Wear Ratio (WWR)

    # old one

    # f4 = 0.0 + (0.049740 * TON) - (0.00104365 * 10.0**-3.0 * TOFF) + (2.06250 * 10.0**-4.0 * IP) + (8.05208 * 10.0**-4.0 * SV) + (6.06636 * 10.0**-3.0 * WT) + (0.080469 * 10.0**-4.0 * TON**2.0) - (3.16358 * 10.0**-8.0 * WT**2.0) - (1.34375 * 10.0**-3.0 * TON * SV) 

    # Updated to reflect correct equation (found via matlab)

    f4 = 2.2831    - (0.041917 * TON) - (0.0018264 * TOFF) + (0.00020625 * IP) + (0.0073896 * SV) + (6.0664 * 10.0**-5.0 * WT) + (2.0117 * 10.0**-4.0 * TON**2.0) - (3.1636 * 10.0**-8.0 * WT**2.0) - (6.7187 * 10.0**-5.0 * TON * SV)

    f[3] = SymbolicFunction(f4, VAR)

    objectives = []

    for f_temp in f:

        objectives.append(f_temp.get_y(all_points).flatten())

    objectives = np.transpose(np.asarray(objectives))

    return objectives

    

def run_tests():

    scale0 = False

    scale1 = True

    scale2 = False

    scale3 = False
    

    # switch to compare the different error metrics

    em_chosen = hdemi



    total_start_time = time.time()



# Test plotting functions

# Working Correctly!



    # max_len = 10

    # values = np.arange(max_len)

    # names = ["blah", "crap", "sadfkad", "aef", "iou", "wow", "qwefa"]

    # for i in range(5,6):

        # test_points = np.zeros((max_len, i))

        # test_bound = np.zeros((max_len, i))

        # for j in range(i):

            # test_points[:,j] = np.arange(max_len)

            # test_bound[:,j] = np.arange(max_len) + max_len

        # plot(test_points, values=values, bnd=test_bound, names=names)

        # plot(test_points, bnd=test_bound, names=names)

        # plot(test_points, values=values, names=names)

    # plot_combinations(test_points, values, test_bound, names)

    # exit()

    

# Scale 0
# Lvl z
# Working!

    if scale0:

        print('Starting Scale 0')

        ferriteFraction=(np.arange(0.1, 1.0, 0.05))

        ferriteSize=(np.arange(1.0,30.0,1.0))

        pearliteSpacing=(np.arange(0.05,0.35,0.02))

        print(ferriteSize)

        print(ferriteFraction)

        print(pearliteSpacing)

        ie1=0.01

        ie2=0.01

        ie3=0.01
        
        Mn = 0.7

        vars = sympy.symbols('FERRITEFRACTION FERRITESIZE PEARLITESPACING')

        FERRITEFRACTION = vars[0]

        FERRITESIZE = vars[1]

        PEARLITESPACING = vars[2]

        VAR=[FERRITEFRACTION,FERRITESIZE,PEARLITESPACING]

        #f    = -0.853260 + 0.0248455 * PEARLITESPACING + 0.000808578 * PEARLITESPACING * FERRITEFRACTION + 0.000391126 * PEARLITESPACING * FERRITESIZE % thesis values

        #fmax = -0.728569 + 0.0242107 * PEARLITESPACING + 0.000835134 * PEARLITESPACING * FERRITEFRACTION + 0.000393725 * PEARLITESPACING * FERRITESIZE

        #fmin = -0.977951 + 0.0254804 * PEARLITESPACING + 0.000782021 * PEARLITESPACING * FERRITEFRACTION + 0.000388527 * PEARLITESPACING * FERRITESIZE
        
        #8 variables for Mechanical Properties
        #For now: hold chemical composition + Ferrite Transformation T constant
        #[Mn] = 0.11, [N] = 0.007, [C] = 0.18, [Si] = 0.36, [P] = 0.019, [Cu] = 0.08, Er = 0, Tmf = 700
        
        Si = 0.36
        
        N = 0.007
        
        C = 0.18
        
        P = 0.019
        
        Cu = 0.08
        
        Tmf = 700

        f = [0] * 3

        fmax = [0] * 3

        fmin = [0] * 3
        
        #Variables: FERRITEFRACTION, FERRITESIZE, PEARLITESPACING
        
        #Yield Strength
        
        f[1] = SymbolicFunction(63*Si + 425*(N**0.5) + (FERRITEFRACTION**(1/3))*(35 + 58*Mn) + 17*((0.001*FERRITESIZE)**-0.5) + (1 - FERRITEFRACTION**(1/3))*(179 + 3.9*PEARLITESPACING**-0.5), VAR) # Mean Response

        fmax[1] = SymbolicFunction(62.6 + 26.1*Mn + 60.2*Si + 759*P + 212.9*Cu + 3286*N + 19.7*((0.001*FERRITESIZE)**-0.5), VAR)

        fmin[1] = SymbolicFunction(FERRITEFRACTION*(77.7 + 59.9*Mn + 9.1*((0.001*FERRITESIZE)**-0.5)) + 478*(N**0.5) + 1200*P + (1 - FERRITEFRACTION)*(145.5 + 3.5*PEARLITESPACING**-0.5), VAR)

        #Tensile Strength
        
        f_0 = FERRITEFRACTION*(20 + 2440*(N**0.5) + 18.5*((0.001*FERRITESIZE)**-0.5)) + 750*(1 - FERRITEFRACTION) + 3*(1 - (FERRITEFRACTION**0.5))*(PEARLITESPACING**-0.5) + 92.5*Si
        
        f[0] = SymbolicFunction(f_0, VAR)
        
        fmax[0] = SymbolicFunction(f_0 * 1.05, VAR)
        
        fmin[0] = SymbolicFunction(f_0 * 0.95, VAR)
        
        #Hardness
        
        f_2 = FERRITEFRACTION*(361 - 0.357*Tmf + 50*Si) + 175*(1 - FERRITEFRACTION)
        
        f[2] = SymbolicFunction(f_2, VAR)
        
        fmax[2] = SymbolicFunction(f_2 * 1.05, VAR)
        
        fmin[2] = SymbolicFunction(f_2 * 0.95, VAR)
        
        fa = [fmax]

        fa.append(fmin)

        fa.append(f)

        # names = ['$HD_{EMI}:Level$ $4_{I=1.5 \\mathrm{MPa-ms}} | Regr_4$',  # Chart Title

        #         '$T^{\circ} (\mathrm{MPa})$',             # xlabel

        #         '$E_{dis} (\mathrm{kJ}/\mathrm{m}^2)$',    # ylabel

        #         '$t_{panel} (\mathrm{mm})$']             # zlabel

        names = ['$Level$ $4$',  # Chart Title

                          '$X_{f} (\mathrm{})$',             # xlabel

                          '$D_{a} (\mathrm{Micro m})$',    # ylabel

                          '$S_{o} (\mathrm{Micro m})$']             # zlabel

        xs = [ferriteFraction, ferriteSize, pearliteSpacing]

        dxs = [ie1, ie2, ie3]

        bound = [[0,0,0], [750,330,170]]
        
        #bound = [[130],[sys.maxsize]

        bound = PrismaticBoundary(bound)

        start_time = time.time()
        
        feas_values, bound_scale0 = idem(fa, xs, dxs, bound, em=em_chosen) #The order of columns of bound_scale0 may be incorrect

        print(time.time()-start_time)
        
        with open("scale_0_bounds.p", "wb") as f:

            pickle.dump(bound_scale0, f)
            
        with open("feasible_0", "w") as f:
            
            for item in feas_values:
                f.write("%s\n" % item)
            

        ax = plot(bound_scale0.feasible_points, bnd=bound_scale0.boundary_points, names=names)
        
        

        # ax.azim = 29

        ax.azim = -119

        ax.elev = 22

        plt.savefig('Scale_0_Python', dpi=plt.gcf().dpi)

        plt.show()

        plt.close()



# Scale 1 

# Not Working!



    if scale1:

        print('Starting Scale 1')

        Mn        = (np.arange(0.3,1.5,0.05))

        CR        = (np.arange(5.0,110.0,1))

        D         = (np.arange(10.0,120.0,2.0))
                
        Si = 0.36
        
        N = 0.007
        
        C = 0.18
        
        P = 0.019
        
        Cu = 0.08
        
        Tmf = 700
        
        Er = 0

        ie1=0.01

        ie2=0.01

        ie3=0.01

        ge=0.02

        vars = sympy.symbols('MN COOLRATE GRAINSIZE')

        MN = vars[0]

        COOLRATE = vars[1]

        GRAINSIZE = vars[2]

        VAR=[MN,COOLRATE,GRAINSIZE]

        f = [0] * 3

        fmax = [0] * 3

        fmin = [0] * 3


        #ferriteFraction----------------------------------------------------------------------

        f_1 = 1 - (0.206 - 0.117*MN - 0.0005*COOLRATE - 0.00113*GRAINSIZE + 0.248*C + 0.00032*MN*COOLRATE + 0.000086*MN*GRAINSIZE + 0.9539*MN*C - 4.259*(10**-6)*COOLRATE*GRAINSIZE + 0.00726*COOLRATE*C + 0.0023*GRAINSIZE*C - 0.0305*(MN**2) - 0.0000056*(COOLRATE**2) + 4.859*(10**-6)*(GRAINSIZE**2) + 0.79*(C**2))
        
        f[1] =    SymbolicFunction(f_1, VAR)

        fmax[1] = SymbolicFunction(f_1*(1 + ge), VAR)

        fmin[1] = SymbolicFunction(f_1*(1 - ge), VAR)
        

        #ferriteSize---------------------------------------------------------------------------
        
        #Ceq = (C + MN)/6

        f_0 = (1 - 0.45*(Er**0.5))*((-0.4 + 6.37*(C + MN)/6) + (24.2 - 59*(C + MN)/6)*(COOLRATE**-0.5) + 22.0*(1 - ((2.718**(-0.015))**GRAINSIZE)))
        
        f[0] = SymbolicFunction(f_0, VAR)

        fmax[0] = SymbolicFunction((1.0 + ge) * f_0, VAR)

        fmin[0] = SymbolicFunction((1.0 - ge) * f_0, VAR)
        

        #pearliteSpacing-----------------------------------------------------------------------
        
        f_2 = 0.1307 + 1.027*C - 1.993*(C**2) - 0.1108*MN + 0.0305*(COOLRATE**-0.52)
        
        f[2] = SymbolicFunction(f_2, VAR)

        fmax[2] = SymbolicFunction(f_2*(1.0 + ge), VAR)

        fmin[2] = SymbolicFunction(f_2*(1.0 - ge), VAR)
        
        #---------------------------------------------------------------------------------------

        fa = [fmax]

        fa.append(fmin)

        fa.append(f)

        names = ['$Level$ $3',

                '$Mn (%)$',

                '$Cool Rate (K/min)$',

                '$Grain Size (Micro m)$']

        xs = [Mn, CR, D]

        dxs = [ie1, ie2, ie3]

        start_time = time.time()

        with open("scale_0_bounds.p", "rb") as f:

            bound_scale0 = pickle.load(f)
            
        with open("file_0.txt", "w") as f:
            
            for item in bound_scale0.boundary_points:
                f.write("%s\n" % item)

        # bound_scale0.set_ignore_bounds(np.asarray([[0,0,0],[sys.maxsize, sys.maxsize, sys.maxsize]]))

        with open("file_01.txt", "w") as f:
            
            for item in bound_scale0.boundary_points:
                f.write("%s\n" % item)
                
       
        bound_MWP = [[13.0, 0.30, 0.05], [29.0, 0.95, 0.33]]
        bound_MWP = PrismaticBoundary(bound_MWP)

        #LW the order of the columns of bound_scale0 might be incorrect
        
        feas_values, bound_scale1 = idem(fa, xs, dxs, bound_MWP, em=em_chosen) #LW come back

        plot_combinations(bound_scale1.feasible_points, feas_values, bound_scale1.boundary_points, names)

        print(time.time()-start_time)

        with open("scale_1_bounds.p", "wb") as f:

            pickle.dump(bound_scale1, f)
            
        
        
        ax = plot(bound_scale1.feasible_points, bnd=bound_scale1.boundary_points, names=names)
        
        plt.show()

        # ax.azim = 29

        ax.azim = -119

        ax.elev = 22

        plt.savefig('Scale_1_Python', dpi=plt.gcf().dpi)
        
        plt.close()


# Scale 2 

# Working Correctly!



    bound_types = ['Convex', 'Concave']

    ignore_concavities = [True, False]

 

    if scale2:

        print('Starting Scale 2')

        wcmr=(np.arange(0.15,0.45,0.05))

        vp=(np.arange(0.01,0.51,0.05))

        rm=(np.arange(0.1,30.1,5))

        ft=8

        ie1=0.05

        ie2=0.05

        ie3=0.1

        fe=0.2

        vars = sympy.symbols('WCMR VP RM')

        WCMR = vars[0]

        VP = vars[1]

        RM = vars[2]

        VAR=[WCMR,VP,RM]

        f=[SymbolicFunction(0.177 * (99.3 / WCMR * (1 - VP) / (RM**.5)) ** 0.74, VAR)]

        fmax=[SymbolicFunction(0.216 * (99.3 / WCMR * (1 - VP) / (RM**.5)) ** 0.74, VAR)]

        fmin=[SymbolicFunction(0.144 * (99.3 / WCMR * (1 - VP) / (RM**.5)) ** 0.74, VAR)]

        

        fa = [fmax]

        fa.append(fmin)

        fa.append(f)

        names = ['$Level$ $2_{\\itf_t = %d \\mathrm{MPa}}$' % ft,

            '$w/cm$', '$V_p$', '$r_{pore} (\\mathrm{nm})$']

        xs = [wcmr, vp, rm]

        dxs = [ie1, ie2, ie3]

        bound = [[ft], [sys.maxsize]]

        bound = PrismaticBoundary(bound)

        start_time = time.time()
        

        for (bound_type, ig_conc) in zip(bound_types, ignore_concavities):

            feas_values, bound_scale2 = idem(fa, xs, dxs, bound, em=em_chosen, ignore_concavity=ig_conc)

            with open("scale_2_bounds_%s.p" % bound_type, "wb") as f:

                pickle.dump(bound_scale2, f)

        # plot_combinations(bound_scale2.feasible_points, feas_values, bound_scale2.boundary_points, names)

        feas = bound_scale2.feasible_points

        bnd = bound_scale2.boundary_points

        check_point = np.asarray([[0.3, 0.1, 3]])

        print("Check error point to see if inner before exclude: %s" % str(bound_scale2.is_inner(check_point)))

        if hasattr(bound_scale2, 'make_exclude_points') and callable(bound_scale2.make_exclude_points):

            bound_scale2.make_exclude_points(check_point)

        print("Check error point to see if inner after exclude: %s" % str(bound_scale2.is_inner(check_point)))

        check_point = np.asarray([[0.2, 0.2, 0.1]])

        print("Is other inner point still inner: %s" % str(bound_scale2.is_inner(check_point)))

        print(time.time()-start_time)

        

        ax = plot(bound_scale2.feasible_points, bnd=bound_scale2.boundary_points, names=names)

        check_points = np.asarray([[0.3, 0.1, 3], [0.3, 0.1, 2]])

        ax.scatter(check_points[:,0], check_points[:,1], check_points[:,2], c='r', s=40, alpha=1, marker='s')

        # ax.set_yticklabels([])

        ax.set_ylim([0,.5])

        ax.set_yticks([0,.25,.5])

        ax.set_xticks(np.arange(.1,.51,.1))

        ax.set_zticks(np.arange(0,9,2))

        # plt.locator_params(axis='y',nbins=2)

        ax.azim = -78

        ax.elev = 7.5

        plt.savefig('Scale_2_Python', dpi=plt.gcf().dpi)

        # plt.show()

        plt.close()



        names.append('Distance to boundary')

        xs=(np.arange(0.15,0.45,0.05))

        ys=(np.arange(0.01,0.51,0.01))

        zs=(np.arange(0.1,8,.3))

        points = np.zeros((0,3))

        dists = np.zeros((0,2,3))

        for x in xs:

            for y in ys:

                for z in zs:

                    temp_point = np.asarray([[x,y,z]])

                    if bound_scale2.is_inner(temp_point)[0]:

                        points = np.concatenate((points, temp_point), axis=0)

                        if x == .2 and y == .06 and z == .4:

                            print('blah')

                        dists = np.concatenate((dists, bound_scale2.bound_dist(temp_point)[np.newaxis,...]), axis=0)

        for i in range(2):

            for d in range(3):



                pydem.plot(points, values=dists[:,i,d], names=names)

                plt.savefig('scale2_dists_%d_%d' % (i,d))

        

# Scale 3

    if scale3:

        print( "**********Starting Scale 3**********" )



        vcem=(np.arange(0.1,0.3,0.02))

        vsf=(np.arange(0.03,0.08,0.01))

        wcmr=(np.arange(0.15,0.32,0.02))

        ft=8.0

        tcure=90.0

        ie1=0.05

        ie2=0.05

        ie3=0.05

        fe=0.0

        ge=0.1

        he=0.15

        vars = sympy.symbols('VCEM VSF WCMR')

        VCEM = vars[0]

        VSF = vars[1]

        WCMR = vars[2]

        VAR=[VCEM,VSF,WCMR]

        rhoCem=3150.0

        rhoWater=1000.0

        rhoSF=2200.0

        rhoAgg=2700.0

        costCem=0.081

        costWater=4.0

        costSF=0.88

        costAgg=0.013

        f = [0] * 3

        fmax = [0] * 3

        fmin = [0] * 3

        f[0] = SymbolicFunction(WCMR, VAR)

        fmax[0] = SymbolicFunction((1.0 + fe) * WCMR, VAR)

        fmin[0] = SymbolicFunction((1.0 - fe) * WCMR, VAR)

        f1 = -0.00398 - 0.000167 * tcure + 0.201 * VCEM + 0.193 * VSF + 0.298 * WCMR - 0.000761 * tcure * VCEM - 0.000735 * tcure * VSF - 0.00198 * tcure * WCMR - 0.315 * VCEM * VCEM - 0.558 * VCEM * VSF + 1.37 * VCEM * WCMR + 1.08 * VSF * WCMR - 0.165 * WCMR * WCMR

        f[1] = SymbolicFunction(f1, VAR)

        fmax[1] = SymbolicFunction((1 + ge) * f1, VAR)

        fmin[1] = SymbolicFunction((1 - ge) * f1, VAR)

        f2 = 70.9 - 0.76 * tcure - 71.5 * VCEM - 91.1 * VSF - 16.8 * WCMR + 1.1 * tcure * VCEM + 1.33 * tcure * VSF + 0.307 * tcure * WCMR - 31.9 * VCEM * VCEM - 309 * VCEM * VSF - 65.2 * VCEM * WCMR - 64.5 * VSF * VSF - 78.3 * VSF * WCMR

        f[2] = SymbolicFunction(f2, VAR)

        fmax[2] = SymbolicFunction((1 + he) * f2, VAR)

        fmin[2] = SymbolicFunction((1 - he) * f2, VAR)

        fa = [fmax]

        fa.append(fmin)

        fa.append(f)

        #names = ['$HD_{EMI}: Level$ $1_{T_{cure} = %d ^\circ\mathrm{C} | f_t = %d \mathrm{MPa}}$' % (tcure,ft),

        #        '$V_{cem} (\mathrm{m}^{3})$',

        #        '$V_{sf} (\mathrm{m}^{3})$',

        #        '$w/cm$']

        names = ['$Level$ $1_{T_{cure} = %d ^\circ\mathrm{C} | f_t = %d \mathrm{MPa}}$' % (tcure, ft),

                 '$V_{cem}$',

                 '$V_{sf}$',

                 '$w/cm$']

        xs = [vcem, vsf, wcmr]

        dxs = [ie1, ie2, ie3]

        collections_3d = []

        mesh_colors = [(1,0,0), (0,0,1)]

        for i, bound_type in enumerate(bound_types):

            with open("scale_2_bounds_%s.p" % bound_type, "rb") as f:

                bound_scale2 = pickle.load(f)

            start_time = time.time()

            feas_values, bound_scale3 = idem(fa, xs, dxs, bound_scale2, em=em_chosen, ignore_concavity=True, ignore_boundary=False)

            print(time.time()-start_time)

            # ax = plot(bound_scale3.feasible_points, bnd=bound_scale3.boundary_points, names=names)

            ax = plot(bound_scale3.feasible_points, bound_scale3.boundary_points, names=names)

            plt.locator_params(axis='y',nbins=5)

            ax.set_xlim([.1,.32])

            # print(ax.get_ylim())

            # print(ax.get_zlim())

            ax.set_ylim([.02,.08])

            ax.set_zlim([.14,.32])

            ax.set_zticks([.15, .2, .25, .3])

            plt.gcf().tight_layout(pad=2)

            plt.savefig('Scale_3_Python_%s' % bound_type, dpi=plt.gcf().dpi)

            # plt.show()

            plt.close()



            problem_pt = [.028, .05, .17]



            if hasattr(bound_scale3, 'make_exclude_points') and callable(bound_scale3.make_exclude_points):

                triangles = bound_scale3.get_external_faces()

                mesh = tvtk.PolyData(points=bound_scale3.bound.points, polys=triangles)

                write_data(mesh, 'bound_surface' + bound_type)



            # x = bound_scale3.bound.points[:,0]

            # y = bound_scale3.bound.points[:,1]

            # z = bound_scale3.bound.points[:,2]



            # mlab.triangular_mesh(x, y, z, triangles, color=mesh_colors[i], opacity=0.5, representation='surface')

            # mlab.triangular_mesh(x, y, z, triangles, color=(0,0,0), opacity=1, representation='wireframe')





        # mlab.axes(extent=[.1,.32,.02,.08,.14,.32],ranges=[.1,.32,.02,.08,.14,.32], line_width=3.0, color=(0,0,0))

        # mlab.xlabel('V_cem')

        # mlab.ylabel('V_sf')

        # mlab.zlabel('w/cm')

        # mlab.show()

        # mlab.savefig('Scale_3_Python_mayavi.pdf', size=[3000,3000], magnification='auto')

    print('Total Execution time: ' + str(time.time() - total_start_time))





if __name__ == "__main__":

    run_tests()