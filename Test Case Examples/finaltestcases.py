import os

import sys

cur_path = os.path.abspath(__file__)

cur_path = os.path.dirname(cur_path)

sys.path.append(os.path.join(cur_path, ".."))

import scipy.ndimage.interpolation as sni

import h5py


import pickle

from sklearn.gaussian_process import GaussianProcessRegressor

from mpl_toolkits.mplot3d.art3d import Poly3DCollection #LW This is not compatible with Python 3

# import mayavi as mlab #LW I could not get the vtk to install so I can't install mayavi.

# from tvtk.api import tvtk, write_data

import numpy as np

import matplotlib.pyplot as plt

import pydem

from pydem import *

from simpy import *




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

    desirability = np.power(desirability, 1.0/np.sum(weights))

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

    scale0 = False #LW Working

    scale1 = False #LW Working

    scale2 = True #LW Working

    scale3 = False #LW Working

    wedm = False #LW Working

    wedm_optimal = False #LW Workingq

    test_numeric = False #LW Working.

    test_gaussian = False
 
 
    test_cpfem = False #Need file SinglePt_data_all_Ti64B.hdf5

    test_holes = False #LW Working. Returns a figure, but algorithm in find_boundary had to be adjusted. 
    
    test_multi_volumes = False #LW Working. There are lines with type error problems that have been commented out. 

    test_projected_dists = False #LW Working.
    
    test_expand_plotting = False #LW Working without line 3446 in pydem.py

    test_sparse_input = False #LW Gaussian working but not idem

    test_idce = False #LW Working

    test_dynamic_bound = False #LW Working

    test_concavity_theorem = False #LW Working
    

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

# Working Correctly!

    if scale0:

        print('Starting Scale 0')

        strength=(np.arange(10.0,22.0,2.0))

        dissipatedenergy=(np.arange(20.0,105.0,5.0)) 

        thickness=(np.arange(39.0,69.0,6.0))

        print(dissipatedenergy)

        print(strength)

        print(thickness)

        ie1=0.1

        ie2=0.1

        ie3=0.1

        fe=0.1

        vars = sympy.symbols('STRENGTH DISSIPATEDENERGY THICKNESS')

        STRENGTH = vars[0]

        DISSIPATEDENERGY = vars[1]

        THICKNESS = vars[2]

        VAR=[STRENGTH,DISSIPATEDENERGY,THICKNESS]

        #f    = -0.853260 + 0.0248455 * THICKNESS + 0.000808578 * THICKNESS * STRENGTH + 0.000391126 * THICKNESS * DISSIPATEDENERGY % thesis values

        #fmax = -0.728569 + 0.0242107 * THICKNESS + 0.000835134 * THICKNESS * STRENGTH + 0.000393725 * THICKNESS * DISSIPATEDENERGY

        #fmin = -0.977951 + 0.0254804 * THICKNESS + 0.000782021 * THICKNESS * STRENGTH + 0.000388527 * THICKNESS * DISSIPATEDENERGY

        f = -0.857 + 0.0262 * THICKNESS + 0.000651 * THICKNESS * STRENGTH + 0.000422 * THICKNESS * DISSIPATEDENERGY # IMMI paper values


        fmax = 1.1*f

        fmin = 0.9*f

        fa = [[SymbolicFunction(fmax, VAR)]]

        fa.append([SymbolicFunction(fmin, VAR)])

        fa.append([SymbolicFunction(f, VAR)])

        # names = ['$HD_{EMI}:Level$ $4_{I=1.5 \\mathrm{MPa-ms}} | Regr_4$',  # Chart Title

        #         '$T^{\circ} (\mathrm{MPa})$',             # xlabel

        #         '$E_{dis} (\mathrm{kJ}/\mathrm{m}^2)$',    # ylabel

        #         '$t_{panel} (\mathrm{mm})$']             # zlabel

        names = ['$Level$ $4$',  # Chart Title

                          '$T^{\circ} (\mathrm{MPa})$',             # xlabel

                          '$E_{dis} (\mathrm{MPa-mm})$',    # ylabel

                          '$t_{panel} (\mathrm{mm})$']             # zlabel

        xs = [strength, dissipatedenergy, thickness]

        dxs = [ie1, ie2, ie3]

        bound = [[1.5], [sys.maxsize]] 

        bound = PrismaticBoundary(bound)

        start_time = time.time()

        feas_values, bound_scale0 = idem(fa, xs, dxs, bound, em=em_chosen)
        
        print((time.time()-start_time))

        with open("scale_0_bounds.p", "wb") as f:

            pickle.dump(bound_scale0, f)
        print('shape',bound_scale0.feasible_points.shape)
        
    
        ax = plot(bound_scale0.feasible_points, bnd=bound_scale0.boundary_points, names=names)
        
        
        # ax.azim = 29

        ax.azim = -119

        ax.elev = 22
        

        plt.savefig('Scale_0_Python', dpi=plt.gcf().dpi)

        plt.show()

        plt.close()



# Scale 1 

# Working Correctly!



    if scale1:

        print('Starting Scale 1')

        ft        = (np.arange(5.00,9.50,0.5))

        pitch     = (np.arange(6.00,39.0,3.0))

        thickness = (np.arange(39.0,69.0,6.0))

        vf=0.019

        ie1=0.1

        ie2=0.1

        ie3=0.1

        fe=0.1

        ge=0.2

        vars = sympy.symbols('FT PITCH THICKNESS')

        FT = vars[0]

        PITCH = vars[1]

        THICKNESS = vars[2]

        VAR=[FT,PITCH,THICKNESS]

        f = [0] * 3

        fmax = [0] * 3

        fmin = [0] * 3

        f[1] =    SymbolicFunction(0.166 + (4320.0 * vf) - (62.4 * vf * PITCH), VAR)
        
        fmax[1] = SymbolicFunction((1.0 + fe) * (0.166 + (4320.0 * vf) - (62.4 * vf * PITCH)), VAR)
         
        fmin[1] = SymbolicFunction((1.0 - fe) * (0.166 + (4320.0 * vf) - (62.4 * vf * PITCH)), VAR)

       

        f_0 = 3.47 + 0.569*(1-vf)*FT + 4380*(vf-0.00866)*PITCH**-0.937

        f[0] = SymbolicFunction(f_0, VAR)
        

        fmax[0] = SymbolicFunction((1.0 + ge) * f_0, VAR)
       
        fmin[0] = SymbolicFunction((1.0 - ge) * f_0, VAR)
       


        # if vf >= 0.0075:

        #     f_0 = (0.089 * (1.0 - vf) * FT) + (1300.0 * (vf - 0.0075) * (PITCH ** - 0.3)) # thesis function

        #     f_0 = 0.166 + 4320*vf

        #     f[0] = SymbolicFunction(f_0, VAR)

        #     fmax[0] = SymbolicFunction((1.0 + ge) * f_0, VAR)

        #     fmin[0] = SymbolicFunction((1.0 - ge) * f_0, VAR)

        # else:

        #     f[0] = SymbolicFunction((0.85 * (1.0 - vf) * FT) + (1300.0 * (0.0) * (PITCH ** - 0.3)), VAR)

        #     fmax[0] = SymbolicFunction((1.0 + ge) * ((0.85 * (1.0 - vf) * FT) + (1300.0 * (0.0) * (PITCH ** - 0.3))), VAR)

        #     fmin[0] = SymbolicFunction((1.0 - ge) * ((0.85 * (1.0 - vf) * FT) + (1300.0 * (0.0) * (PITCH ** - 0.3))), VAR)

        f[2] = SymbolicFunction(THICKNESS, VAR)
        
        fmax[2] = SymbolicFunction(THICKNESS, VAR)
        
        fmin[2] = SymbolicFunction(THICKNESS, VAR)
        
        
        
       
        fa = [fmax]
        
        fa.append(fmin)

        fa.append(f)
        


        names = ['$Level$ $3 V_f= %d $' % vf,

                '${\\itf_t} (MPa)$',

                '$Pitch (mm)$',

                '$t_{panel} (mm)$']

        xs = [ft, pitch, thickness]

        dxs = [ie1, ie2, ie3]

        start_time = time.time()

        with open("scale_0_bounds.p", "rb") as f:

            bound_scale0 = pickle.load(f)
            
        with open("file_0.txt", "w") as f:
            
            for item in bound_scale0.boundary_points:
                 f.write("%s\n" % item)

        bound_scale0.set_ignore_bounds(np.asarray([[0,0,sys.maxsize],[sys.maxsize, sys.maxsize, sys.maxsize]]))

        with open("file_01.txt", "w") as f:
            
            for item in bound_scale0.boundary_points:
                f.write("%s\n" % item)
        
        test_x = np.asarray([[6.5, 6, 45]])

        test_x = np.asarray([[7, 6, 45]])

        test_x = np.asarray([[5, 6, 51]])

        dx = test_x * dxs
        
        
        
        # hdemi(fa, test_x, dx, bound_scale0) #LW not functioning quite correctly
        feas_values, bound_scale1 = idem(fa, xs, dxs, bound_scale0, em=em_chosen)

        plot_combinations(bound_scale1.feasible_points, feas_values, bound_scale1.boundary_points, names)

        print((bound_scale1.boundary_points))
        
        print((time.time()-start_time))

        with open("scale_1_bounds.p", "wb") as f:

            pickle.dump(bound_scale1, f)
        
        ax = plot(bound_scale1.feasible_points, bnd=bound_scale1.boundary_points, names=names)

        # ax.azim = 29

        ax.azim = -119

        ax.elev = 22

        plt.savefig('Scale_1_Python', dpi=plt.gcf().dpi)
       
        plt.show()
       
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

                pickle.dump(bound_scale2, f) #LW changed from cPickle to pickle

        plot_combinations(bound_scale2.feasible_points, feas_values, bound_scale2.boundary_points, names)

        feas = bound_scale2.feasible_points


        bnd = bound_scale2.boundary_points

        check_point = np.asarray([[0.3, 0.1, 3]])

        print(("Check error point to see if inner before exclude: %s" % str(bound_scale2.is_inner(check_point))))

        if hasattr(bound_scale2, 'make_exclude_points') and callable(bound_scale2.make_exclude_points):

            bound_scale2.make_exclude_points(check_point)

        print(("Check error point to see if inner after exclude: %s" % str(bound_scale2.is_inner(check_point))))

        check_point = np.asarray([[0.2, 0.2, 0.1]])

        print(("Is other inner point still inner: %s" % str(bound_scale2.is_inner(check_point))))

        print((time.time()-start_time))

     

        ax = plot(bound_scale2.feasible_points, bnd=bound_scale2.boundary_points, names=names)

        check_points = np.asarray([[0.3, 0.1, 3], [0.3, 0.1, 2]])

        ax.scatter(check_points[:,0], check_points[:,1], check_points[:,2], c='r', s=40, alpha=1, marker='s')

        # ax.set_yticklabels([])

        ax.set_ylim([0,.5])

        ax.set_yticks([0,.25,.5])

        ax.set_xticks(np.arange(.1,.51,.1))

        ax.set_zticks(np.arange(0,9,2))

        plt.locator_params(axis='y',nbins=2)

        ax.azim = -78

        ax.elev = 7.5

        plt.savefig('Scale_2_Python', dpi=plt.gcf().dpi)

        plt.show()

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
        
        plt.show()
        plt.close()

        

# Scale 3

    if scale3:

        print ("**********Starting Scale 3**********")



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

            print((time.time()-start_time))

            ax = plot(bound_scale3.feasible_points, bound_scale3.boundary_points, names=names)

            plt.locator_params(axis='y',nbins=5)

            ax.set_xlim([.1,.32])

            # print(ax.get_ylim())

            # print(ax.get_zlim())

            ax.set_ylim([.02,.08])

            ax.set_zlim([.14,.32])

            ax.set_zticks([.15, .2, .25, .3])

            plt.gcf().tight_layout(pad=2)

            plt.savefig('Scale_3_Python_%s' % bound_type, dpi=plt.gcf().dpi) #LW The savefig function causes remapping error. 

            plt.show()

            plt.close()



            problem_pt = [.028, .05, .17]



            if hasattr(bound_scale3, 'make_exclude_points') and callable(bound_scale3.make_exclude_points):

                triangles = bound_scale3.get_external_faces()

                mesh = tvtk.PolyData(points=bound_scale3.bound.points, polys=triangles)

                write_data(mesh, 'bound_surface' + bound_type)



            x = bound_scale3.bound.points[:,0]

            y = bound_scale3.bound.points[:,1]

            z = bound_scale3.bound.points[:,2]



            mlab.triangular_mesh(x, y, z, triangles, color=mesh_colors[i], opacity=0.5, representation='surface')

            mlab.triangular_mesh(x, y, z, triangles, color=(0,0,0), opacity=1, representation='wireframe')





        mlab.axes(extent=[.1,.32,.02,.08,.14,.32],ranges=[.1,.32,.02,.08,.14,.32], line_width=3.0, color=(0,0,0))

        mlab.xlabel('V_cem')

        mlab.ylabel('V_sf')

        mlab.zlabel('w/cm')

        mlab.show()

        mlab.savefig('Scale_3_Python_mayavi.pdf', size=[3000,3000], magnification='auto')



    if wedm:

        print("**********Starting WEDM Optimization**********")

        # 2 values per variable

        # Ran in 141 seconds

    #    ton   = (np.arange(112.0, 128.0, 8.0))     # Pulse on Time (112 to 120)

    #    toff  = (np.arange(44.0, 68.0, 12.0))      # Pulse off time (44 to 56)

    #    ip    = (np.arange(120.0, 280.0, 80.0))    # Peak Current (120 to 200)

    #    sv    = (np.arange(40.0, 80.0, 20.0))      # Spark Gap Voltage (40 to 60)

    #    wf    = (np.arange(4.0, 16.0, 6.0))        # Wire Feed (4 to 10)

    #    wt    = (np.arange(500.0, 2300.0, 900.0))  # Wire Tension (500 to 1400)

        # 3 values per variable

        ton   = (np.arange(112.0, 122.0, 2.0))     # Pulse on Time (112 to 120)

        toff  = (np.arange(44.0, 60.0, 4.0))       # Pulse off time (44 to 56)

        ip    = (np.arange(120.0, 220.0, 20.0))    # Peak Current (120 to 200)

        sv    = (np.arange(40.0, 65.0, 5.0))      # Spark Gap Voltage (40 to 60)

        wf    = (np.arange(4.0, 12.0, 2.0))        # Wire Feed (4 to 10)

        wt    = (np.arange(500.0, 1700.0, 300.0))  # Wire Tension (500 to 1400)

        # print(ton)

        # print(toff)

        # print(ip)

        # print(sv)

        # print(wf)

        # print(wt)

        # ie1=0.05

        # ie2=0.05

        # ie3=0.05

        # ie4=0.05

        # ie5=0.05

        # ie6=0.05

        # fe1=0.025

        # fe2=0.025

        # fe3=0.025

        # fe4=0.025

        ie1=0.025

        ie2=0.025

        ie3=0.025

        ie4=0.025

        ie5=0.025

        ie6=0.025

        fe1=0.025

        fe2=0.025

        fe3=0.025

        fe4=0.025

        vars = sympy.symbols('TON TOFF IP SV WF WT')

        TON = vars[0]

        TOFF= vars[1]

        IP  = vars[2]

        SV  = vars[3]

        WF  = vars[4]

        WT  = vars[5]

        VAR = [TON,TOFF,IP,SV,WF,WT]

        

    #    MR = 11.06243 + (0.021958 * TON) - (0.4023 * TOFF) - (0.039917 * IP) - (0.045079 * SV) + (3.06872 * 10.0**-3.0 * TOFF**2.0) + (1.79134 * 10.0**-5.0 * IP**2.0) + (2.30469 * 10.0**-4.0 * TON * IP) + (1.96875 * 10.0**-4.0 * TOFF * IP) + (7.39583 * 10.0**-4.0 * TOFF * SV)

    #    SR = 10.87640 - (0.038333 * TON) - (0.1284 * TOFF) - (0.072771 * IP) - (0.008667 * SV) + (0.1247 * WF) - (0.000365432 * WT) + (9.72222 * 10.0**-4.0 * TOFF**2.0) - (8.75000 * 10.0**-3.0 * WF**2.0)  + (1.91358 * 10.0**-7.0 * WT**2.0) + (6.48438 * 10.0**-4.0 * TON * IP)

    #    DD = 150.52609+ (1.43750  * TON) - (4.0869 * TOFF) - (0.23551  * IP) - (0.46806  * SV) - (1.8670 * WF) - (0.015278    * WT) + (0.034827 * TOFF**2.0) + (8.23864 * 10.0**-4.0 * IP**2.0) + (0.12542 * WF**2.0) + (2.77778 * 10.0**-4.0 * SV * WT)

    #    WWR = 0.0     + (0.049740 * TON) - (0.00104365 * 10.0**-3.0 * TOFF) + (2.06250 * 10.0**-4.0 * IP) + (8.05208 * 10.0**-4.0 * SV) + (6.06636 * 10.0**-3.0 * WT) + (0.080469 * TON**2.0) - (3.16358 * 10.0**-8.0 * WT**2.0) - (1.34375 * 10.0**-3.0 * TON * SV)       

    #    f    = [(MR * SR * DD * WWR)**(1.0/4.0)]

    #    fmax = [(1+fe) * (MR * SR * DD * WWR)**(1.0/4.0)]

    #    fmin = [(1-fe) * (MR * SR * DD * WWR)**(1.0/4.0)]



        f = [0] * 4

        fmax = [0] * 4

        fmin = [0] * 4



        # Machining Rate (MR)

        f1 = 11.06243 + (0.021958 * TON) - (0.4023 * TOFF) - (0.039917 * IP) - (0.045079 * SV) + (3.06872 * 10.0**-3.0 * TOFF**2.0) + (1.79134 * 10.0**-5.0 * IP**2.0) + (2.30469 * 10.0**-4.0 * TON * IP) + (1.96875 * 10.0**-4.0 * TOFF * IP) + (7.39583 * 10.0**-4.0 * TOFF * SV)

        f[0] = SymbolicFunction(f1, VAR)

        fmax[0] = SymbolicFunction((1.0 + fe1) * f1, VAR)

        fmin[0] = SymbolicFunction((1.0 - fe1) * f1, VAR)



        # Surface Roughness (SR)

        f2 = 10.87640 - (0.038333 * TON) - (0.1284 * TOFF) - (0.072771 * IP) - (0.008667 * SV) + (0.1247 * WF) - (0.000365432 * WT) + (9.72222 * 10.0**-4.0 * TOFF**2.0) - (8.75000 * 10.0**-3.0 * WF**2.0)  + (1.91358 * 10.0**-7.0 * WT**2.0) + (6.48438 * 10.0**-4.0 * TON * IP)

        f[1] = SymbolicFunction(f2, VAR)

        fmax[1] = SymbolicFunction((1.0 + fe2) * f2, VAR)

        fmin[1] = SymbolicFunction((1.0 - fe2) * f2, VAR)

        # Dimensional Deviation (DD)

        f3 = 150.52609+ (1.43750  * TON) - (4.0869 * TOFF) - (0.23551  * IP) - (0.46806  * SV) - (1.8670 * WF) - (0.015278    * WT) + (0.034827 * TOFF**2.0) + (8.23864 * 10.0**-4.0 * IP**2.0) + (0.12542 * WF**2.0) + (2.77778 * 10.0**-4.0 * SV * WT)

        f[2] = SymbolicFunction(f3, VAR)

        fmax[2] = SymbolicFunction((1.0 + fe3) * f3, VAR)

        fmin[2] = SymbolicFunction((1.0 - fe3) * f3, VAR)

        # Wire Wear Ratio (WWR)

        # Updated to reflect correct equation (found via matlab)

    #    f[3] = 0.0      + (0.049740 * TON) - (0.00104365 * 10.0**-3.0 * TOFF) + (2.06250 * 10.0**-4.0 * IP) + (8.05208 * 10.0**-4.0 * SV) + (6.06636 * 10.0**-3.0 * WT) + (0.080469 * 10.0**-4.0 * TON**2.0) - (3.16358 * 10.0**-8.0 * WT**2.0) - (1.34375 * 10.0**-3.0 * TON * SV) 

    #    fmax[3] = (1.0 + fe4) * (0.0      + (0.049740 * TON) - (0.00104365 * 10.0**-3.0 * TOFF) + (2.06250 * 10.0**-4.0 * IP) + (8.05208 * 10.0**-4.0 * SV) + (6.06636 * 10.0**-3.0 * WT) + (0.080469 * 10.0**-4.0 * TON**2.0) - (3.16358 * 10.0**-8.0 * WT**2.0) - (1.34375 * 10.0**-3.0 * TON * SV))

    #    fmin[3] = (1.0 - fe4) * (0.0      + (0.049740 * TON) - (0.00104365 * 10.0**-3.0 * TOFF) + (2.06250 * 10.0**-4.0 * IP) + (8.05208 * 10.0**-4.0 * SV) + (6.06636 * 10.0**-3.0 * WT) + (0.080469 * 10.0**-4.0 * TON**2.0) - (3.16358 * 10.0**-8.0 * WT**2.0) - (1.34375 * 10.0**-3.0 * TON * SV)) 

        f4 = 2.2831    - (0.041917 * TON) - (0.0018264 * TOFF) + (0.00020625 * IP) + (0.0073896 * SV) + (6.0664 * 10.0**-5.0 * WT) + (2.0117 * 10.0**-4.0 * TON**2.0) - (3.1636 * 10.0**-8.0 * WT**2.0) - (6.7187 * 10.0**-5.0 * TON * SV)

        f[3] = SymbolicFunction(f4, VAR)

        fmax[3] = SymbolicFunction((1.0 + fe4) * f4, VAR)

        fmin[3] = SymbolicFunction((1.0 - fe4) * f4, VAR)

        fa = [fmax]

        fa.append(fmin)

        fa.append(f)

        names = ['$WEDM$',          # Chart Title

                '$T_{on} (\mu \mathrm{s})$',   # 0-dir

                '$T_{off} (\mu \mathrm{s})$',  # 1-dir

                '$Ip (\mathrm{A})$',          # 2-dir

                '$SV (\mathrm{V})$',          # 3-dir

                '$WF (\mathrm{m}/\mathrm{min})$',      # 4-dir

                '$WT (\mathrm{g})$'           # 5-dir

                ]            

        xs = [ton, toff, ip, sv, wf, wt]

        dxs = [ie1, ie2, ie3, ie4, ie5, ie6]

        bound = [[0.395, 2.15, 140.0, 0.048], [1.28, 3.28, 165.0, 0.107]]

    #    bound = [[0.395, 2.15, 140.0], [1.28, 3.28, 165.0]]

        bound = PrismaticBoundary(bound)

        start_time = time.time()

        feas_values, bound_WEDM = idem(fa, xs, dxs, bound, em=em_chosen, ignore_concavity=True, ignore_boundary=True)



        temp_out = "wedm_bounds.p"

        if not os.path.isfile(temp_out): #LW Edit Not sure of the function of this line, but it was causing the bound_WEDM file to be blank.

            with open(temp_out, "wb") as f:

                pickle.dump(bound_WEDM, f)#LW changed from cPickle to pickle
                
        # axis_ranges = [[112, 115, 118, 121], [44, 50, 55, 60], [120, 150, 180, 210]]

        axis_ranges = xs[:3]

        # plot_expand_dims_subplot(bound_WEDM.feasible_points, bnd=bound_WEDM.boundary_points, names=names, axis_ranges=axis_ranges, discrete_dims=True)

        

        optimal_kumar_40 = np.asarray([[112, 44, 120.02, 46.51, 9.55, 1400]])

        optimal_kumar_43 = np.asarray([[112, 44.06, 177.12, 40, 4, 968]])

        

        feas = bound_WEDM.feasible_points

        

        filter_40 = np.ones(feas.shape[0], dtype=bool)

        filter_40 = np.logical_and(filter_40, feas[:,0]==112)

        filter_40 = np.logical_and(filter_40, feas[:,1]==44)

        filter_40 = np.logical_and(filter_40, feas[:,2]==120)

            

        feas_40 = feas[filter_40, 3:]

        hdemi_40 = feas_values[filter_40]

        print((bound_WEDM.is_inner(optimal_kumar_40)))
        
        ax = plot(feas_40, names=([names[0]] + names[4:]), values=hdemi_40)

        ax.scatter(optimal_kumar_40[:,3], optimal_kumar_40[:,4], optimal_kumar_40[:,5], c='g', marker='D', alpha=1, s=60)

        ax.set_xticks(np.arange(38, 56, 4))

        ax.set_yticks(np.arange(3, 12, 2))

        plt.tight_layout(pad=2.5)

        plt.savefig('wedm_opt_inner', dpi=plt.gcf().dpi) #LW Savefig function causes remapping error.

        plt.show()

        plt.close()

        

        filter_36 = np.ones(feas.shape[0], dtype=bool)

        filter_36 = np.logical_and(filter_36, feas[:,0]==112)

        filter_36 = np.logical_and(filter_36, feas[:,1]==44)

        filter_36 = np.logical_and(filter_36, feas[:,2]==180)

            

        feas_36 = feas[filter_36, 3:]

        hdemi_36 = feas_values[filter_36]

        print((bound_WEDM.is_inner(optimal_kumar_43)))

        ax = plot(feas_36, names=([names[0]] + names[4:]), values=hdemi_36)

        ax.scatter(optimal_kumar_43[:,3], optimal_kumar_43[:,4], optimal_kumar_43[:,5], c='r', marker='s', alpha=1, s=60)

        ax.set_xticks(np.arange(38, 56, 4))

        ax.set_yticks(np.arange(3, 12, 2))

        plt.tight_layout(pad=2.5)

        plt.savefig('wedm_opt_outer', dpi=plt.gcf().dpi) #LW savefig function causes remapping errors. 

        plt.show()

        plt.close()



        print((time.time()-start_time))

    

    if wedm_optimal:

        ton   = np.linspace(112.0, 122.0, 10)   # Pulse on Time (112 to 120)

        toff  = np.linspace(44.0, 56.0, 10)     # Pulse off time (44 to 56)

        ip    = np.linspace(120.0, 200.0, 10)   # Peak Current (120 to 200)

        sv    = np.linspace(40.0, 60.0, 10)     # Spark Gap Voltage (40 to 60)

        wf    = np.linspace(4.0, 10.0, 10)      # Wire Feed (4 to 10)

        wt    = np.linspace(500.0, 1400.0, 10)  # Wire Tension (500 to 1400)

        

        all_inputs = (ton, toff, ip, sv, wf, wt)

        n = np.prod([x.size for x in all_inputs])

        all_points = np.array(np.meshgrid(*all_inputs)).T.reshape(-1,len(all_inputs))

        print((all_points.shape))

        print("Starting load boundary")
        

        with open("wedm_bounds.p", "rb") as f: #LW Edit changed 'r' to 'rb' for bytes compatibility.

            bound_WEDM = pickle.load(f) #LW EDIT changed cPickle to pickle

        print("Finished load boundary")

        start_time = time.time()

        inners = bound_WEDM.is_inner(all_points)

        print(("Finding inner points took %f" % (time.time()-start_time)))

        start_time = time.time()

        all_points = all_points[inners,:]

        print(("Trimming list took %f" % (time.time()-start_time)))

        objectives = responses(all_points)

        start_time = time.time()

        

        desirabilityvar = desirability(objectives) #LW EDIT Changed the variable because it references itself, but it takes like three hours to run, so i'll come back

        print(("Evaluating all points took %f" % (time.time()-start_time)))

        with open("wedm_values.p", "wb") as f:

            pickle.dump([all_points, desirabilityvar], f) #LW EDIT change cPickle to pickle

        print((np.max(desirability)))

        

        

    if test_numeric:

        f = lambda x: 1-np.linalg.norm(x)

        f = [[NumericFunction(f)]]

        xs = [np.linspace(0,5,10), np.linspace(0,5,10)]

        dxs = [0.01, 0.01]

        bound = [[0],[10]]

        bound = PrismaticBoundary(bound)

        feas_values, bound_numeric = idem(f, xs, dxs, bound)

        plot_combinations(bound_numeric.feasible_points, feas_values, bound_numeric.boundary_points)
        
        
    if test_gaussian:

        #gp = GaussianProcessRegressor(theta0=.1, thetaL=.001, thetaU=1.0) #LW this needs an updated equivalent. 
        #gp = GaussianProcessOptimizer(obj_func, 0.1, .001)
        gp = GaussianProcessRegressor()
        
        in_points = 3

        xs = np.linspace(0,10,in_points)

        ys = np.linspace(0,10,in_points)

        xs, ys = np.meshgrid(xs, ys)

        xs = np.reshape(xs, (xs.size,1))

        ys = np.reshape(ys, (ys.size,1))

        xs = np.concatenate((xs,ys), axis=1)

        response = np.zeros((len(xs)))

        for i in range(len(response)):

            response[i] = 1 - np.linalg.norm(xs[i,:])

        gp.fit(xs, response)

        in_x = xs

        in_y = response

        plot_points = 20

        xs = np.linspace(0,10,plot_points)

        ys = np.linspace(0,10,plot_points)

        xs, ys = np.meshgrid(xs, ys)

        xs = np.reshape(xs, (xs.size,1))

        ys = np.reshape(ys, (ys.size,1))

        xs = np.concatenate((xs,ys), axis=1)

        response, MSE = gp.predict(xs,True)
        
        sigma = np.sqrt(MSE)

        upper = response + 2*sigma

        lower = response - 2*sigma

        h = plt.figure()

        ax = h.add_subplot(111, projection='3d')

        ax.scatter(xs[:,0], xs[:,1], response, c='k')

        ax.scatter(xs[:,0], xs[:,1], upper, c='r')

        ax.scatter(xs[:,0], xs[:,1], lower, c='b')

        real_response = np.zeros((len(xs)))

        for i in range(len(response)):

            real_response[i] = 1 - np.linalg.norm(xs[i,:])

        ax.scatter(xs[:,0], xs[:,1], real_response, c='g')

        plt.show()

        

        f = lambda x: gp.predict(x)[0]

        f_upper = lambda x: kriging_upper(gp, x)

        f_lower = lambda x: kriging_lower(gp, x)

        #f = [[NumericFunction(f_upper)], [NumericFunction(f_lower)], [NumericFunction(f)]]

        f = [[NumericFunction(f)]]

        xs = [np.linspace(0,5,10), np.linspace(0,5,10)]

        dxs = [0.01, 0.01]

        bound = [[0],[10]]

        bound = PrismaticBoundary(bound)
   
        feas_values, bound_numeric = idem(f, xs, dxs, bound) #LW This isn't quite working. 
        # plot_combinations(bound_numeric.feasible_points, feas_values, bound_numeric.boundary_points)
        
        
    if test_sparse_input:

        nums = [7, 18, 6]

        mask = np.zeros(tuple(nums), dtype=bool)

        mask[:,9:12,:] = True

        test_ignore_volume = PrismaticBoundary(np.asarray([[5,50,30], [30,75,80]]))

        strength = np.linspace(10.0,22.0,nums[0])

        dissipatedenergy = np.linspace(20.0,105.0,nums[1])

        thickness = np.linspace(39.0,69.0,nums[2])

        print(dissipatedenergy)

        print(strength)

        print(thickness)

        ie1=0.1

        ie2=0.1

        ie3=0.1

        fe=0.1

        vars = sympy.symbols('STRENGTH DISSIPATEDENERGY THICKNESS')

        STRENGTH = vars[0]

        DISSIPATEDENERGY = vars[1]

        THICKNESS = vars[2]

        VAR=[STRENGTH,DISSIPATEDENERGY,THICKNESS]

        f    = -0.853260 + 0.0248455 * THICKNESS + 0.000808578 * THICKNESS * STRENGTH + 0.000391126 * THICKNESS * DISSIPATEDENERGY

        fmax = -0.728569 + 0.0242107 * THICKNESS + 0.000835134 * THICKNESS * STRENGTH + 0.000393725 * THICKNESS * DISSIPATEDENERGY

        fmin = -0.977951 + 0.0254804 * THICKNESS + 0.000782021 * THICKNESS * STRENGTH + 0.000388527 * THICKNESS * DISSIPATEDENERGY

        fa = [[SymbolicFunction(fmax, VAR)]]

        fa.append([SymbolicFunction(fmin, VAR)])

        fa.append([SymbolicFunction(f, VAR)])

        names = ['$HD_{EMI}:Level$ $0_{I=1.5 \\mathrm{MPa-ms}} | Regr_4$',  # Chart Title

                '$\sigma_t (\mathrm{MPa})$',             # xlabel

                '$E_{dissipated} (\mathrm{kJ}/\mathrm{m}^2)$',    # ylabel

                '$Thickness (\mathrm{mm})$']             # zlabel

        xs = [strength, dissipatedenergy, thickness]

        dxs = [ie1, ie2, ie3]

        bound = [[1.5], [sys.maxsize]]

        bound = PrismaticBoundary(bound)

        start_time = time.time()

        feas_values, bound_scale0 = idem(fa, xs, dxs, bound, em=em_chosen, ignore_region=test_ignore_volume)

        print(time.time()-start_time)

        ax = plot(bound_scale0.feasible_points, bnd=bound_scale0.boundary_points, names=names)

        # ax.azim = 29

        ax.azim = -119

        ax.elev = 22

        check_points = np.asarray([[15, 70, 65], [18, 70, 65], [18, 70, 60]])

        print(bound_scale0.is_inner(check_points))

        plt.show()

        locations = []

        for x in np.linspace(10,22,20):

            for y in np.linspace(20,105,40):

                for z in np.linspace(39,69,20):

                    locations.append(np.atleast_2d([x,y,z]))

        locations = np.concatenate(locations, axis=0)

        valid = bound_scale0.is_inner(locations)

        valid_locs = locations[valid]

        print(valid_locs.shape)

        plot(valid_locs, names=names)

        plt.show()

      
    


        

    if test_cpfem:



        print("**********Starting CPFEM Optimization**********")

        file_singlept = h5py.File('SinglePt_data_all_Ti64B.hdf5', 'r')

        data_cb = np.asarray(file_singlept['data_cb'])

        data_cs = np.asarray(file_singlept['data_cs'])

        data_cu = np.asarray(file_singlept['data_cu'])

        data_mu = np.asarray(file_singlept['data_mu'])

        fip_cb = np.log(data_cb[:,3])

        fip_cs = np.log(data_cs[:,3])

        fip_cu = np.log(data_cu[:,3])

        mod = data_mu[:,3]

        yld = data_mu[:,4]

        # phi1_fit = data_mu[:,0]

        # phi0_fit = data_mu[:,1]

        # phi2_fit = data_mu[:,2]

        # xs_fit = np.transpose(np.asarray((phi1_fit, phi0_fit, phi2_fit)))

        phi1_fit = (np.linspace(0.0, 360.0, 25)) * np.pi / 180.0

        phi0_fit = (np.linspace(0.0, 180.0, 13)) * np.pi / 180.0

        phi2_fit = (np.linspace(0.0, 360.0, 25)) * np.pi / 180.0

        xs_fit = (phi1_fit, phi0_fit, phi2_fit)

        dxs = [0.001, 0.001, 0.001]



        # Approximately 60 seconds before adding the interaction terms. Not a very good fit.

        #phi1_fit = data_mu[:,0]

        #phi0_fit = data_mu[:,1]

        #phi2_fit = data_mu[:,2]

        #xs_fit = np.transpose(np.asarray((phi1_fit, phi0_fit, phi2_fit)))

        #phi1 = (np.linspace(0.0, 360.0, 10)) * np.pi / 180.0

        #phi0 = (np.linspace(0.0, 180.0, 10)) * np.pi / 180.0

        #phi2 = (np.linspace(0.0, 360.0, 10)) * np.pi / 180.0

        #xs = [phi1, phi0, phi2]

        #ranges = (np.max(xs_fit,axis=0) - np.min(xs_fit, axis=0))/2.0

        #basis = make_fourier(xs_fit, ranges)

        #model1 = np.linalg.lstsq(basis, fip_cb)[0]

        #model2 = np.linalg.lstsq(basis, fip_cs)[0]

        #model3 = np.linalg.lstsq(basis, fip_cu)[0]

        #model4 = np.linalg.lstsq(basis, mod)[0]

        #model5 = np.linalg.lstsq(basis, yld)[0]

        #f = [0] * 5

        #f[0] = lambda x: eval_fourier(np.atleast_2d(x),ranges,model1)

        #f[1] = lambda x: eval_fourier(np.atleast_2d(x),ranges,model2)

        #f[2] = lambda x: eval_fourier(np.atleast_2d(x),ranges,model3)

        #f[3] = lambda x: eval_fourier(np.atleast_2d(x),ranges,model4)

        #f[4] = lambda x: eval_fourier(np.atleast_2d(x),ranges,model5)



        # Approximately 400 seconds. A good fit

        #phi1 = (np.linspace(0.0, 360.0, 10)) * np.pi / 180.0

        #phi0 = (np.linspace(0.0, 180.0, 10)) * np.pi / 180.0

        #phi2 = (np.linspace(0.0, 360.0, 10)) * np.pi / 180.0

        #xs = [phi1, phi0, phi2]

        #f = [0] * 5

        #f[0] = lambda x: si.interpn(xs_fit, np.reshape(fip_cb, (25,13,25)), x, bounds_error=False, fill_value=np.mean(fip_cb))[0]

        #f[1] = lambda x: si.interpn(xs_fit, np.reshape(fip_cs, (25,13,25)), x, bounds_error=False, fill_value=np.mean(fip_cs))[0]

        #f[2] = lambda x: si.interpn(xs_fit, np.reshape(fip_cu, (25,13,25)), x, bounds_error=False, fill_value=np.mean(fip_cu))[0]

        #f[3] = lambda x: si.interpn(xs_fit, np.reshape(mod, (25,13,25)), x, bounds_error=False, fill_value=np.mean(mod))[0]

        #f[4] = lambda x: si.interpn(xs_fit, np.reshape(yld, (25,13,25)), x, bounds_error=False, fill_value=np.mean(yld))[0]



        # Approximately 400 seconds. A good fit. Will allow you to wrap around the data

        phi1 = (np.linspace(0.0, 360.0, 10)) / 360.0 * 24.0

        phi0 = (np.linspace(0.0, 180.0, 10)) / 180.0 * 12.0

        phi2 = (np.linspace(0.0, 360.0, 10)) / 360.0 * 24.0

        xs = [phi1, phi0, phi2]

        f = [0] * 5

        

        sample_points = (phi1_fit, phi0_fit, phi2_fit)

        sample_sizes = list(map(len, sample_points))

        

        f[0] = lambda x: sni.map_coordinates(np.reshape(fip_cb, sample_sizes), np.transpose(x.reshape(1,len(x))), mode='wrap')[0]

        f[1] = lambda x: sni.map_coordinates(np.reshape(fip_cs, sample_sizes), np.transpose(x.reshape(1,len(x))), mode='wrap')[0]

        f[2] = lambda x: sni.map_coordinates(np.reshape(fip_cu, sample_sizes), np.transpose(x.reshape(1,len(x))), mode='wrap')[0]

        f[3] = lambda x: sni.map_coordinates(np.reshape(mod, sample_sizes), np.transpose(x.reshape(1,len(x))), mode='wrap')[0]

        f[4] = lambda x: sni.map_coordinates(np.reshape(yld, sample_sizes), np.transpose(x.reshape(1,len(x))), mode='wrap')[0]

        

        n = np.prod([x.size for x in sample_points])

        sample_points = np.array(np.meshgrid(*sample_points, indexing='ij')).T.reshape(-1,len(sample_points))

        

        plot_vals = [fip_cs, fip_cu, fip_cb, yld, mod]

        plot_names = ['FIP_s', 'FIP_u', 'FIP_b', 'y', 'mod']

        plot_labels = ['$\mathrm{ln}(FIP_s)$', '$\mathrm{ln}(FIP_u)$', '$\mathrm{ln}(FIP_b)$', '$\sigma_y$', '$E$']

        old_pad = pydem.pad_dist

        old_font = pydem.f_size

        pydem.pad_dist = -5

        pydem.f_size = 36

        for (name, vals, label) in zip(plot_names, plot_vals, plot_labels):

            names = ['', '', '', '', '']

            # names = [r'', # Chart Title

            # r'$\phi_{1} (rad)$', # 0-dir

            # r'$\Phi (rad)$', # 1-dir

            # r'$\phi_{2} (rad)$', # 2-dir

            # label

            # ]

            ax = plot(sample_points, values=vals, names=names, disp_cbar=False)
   
            fig = plt.gcf()

            # labels = [ax.title, ax.xaxis.label, ax.yaxis.label, ax.zaxis.label, fig.axes[1].xaxis.label, fig.axes[1].yaxis.label]

            # for label in labels:

                # label.set_fontsize(24)

            # fig.axes[1].tick_params(labelsize=24) # set colorbar label sizes

            ax.axes.xaxis.set_ticklabels([])

            ax.axes.yaxis.set_ticklabels([])

            ax.axes.zaxis.set_ticklabels([])

            plt.tight_layout(pad=1.5)

            plt.savefig('Ti64_' + name + "_no_labels", dpi=fig.dpi)

        

        # Approximately 300 seconds. A good fit

        #phi1 = (np.linspace(0.0, 360.0, 10)) * np.pi / 180.0

        #phi0 = (np.linspace(0.0, 180.0, 10)) * np.pi / 180.0

        #phi2 = (np.linspace(0.0, 360.0, 10)) * np.pi / 180.0

        #xs = [phi1, phi0, phi2]

        #f = [0] * 5

        #model1 = si.RegularGridInterpolator(xs_fit, np.reshape(fip_cb, (25,13,25)), bounds_error=False, fill_value=np.mean(fip_cb))

        #model2 = si.RegularGridInterpolator(xs_fit, np.reshape(fip_cs, (25,13,25)), bounds_error=False, fill_value=np.mean(fip_cs))

        #model3 = si.RegularGridInterpolator(xs_fit, np.reshape(fip_cu, (25,13,25)), bounds_error=False, fill_value=np.mean(fip_cu))

        #model4 = si.RegularGridInterpolator(xs_fit, np.reshape(mod, (25,13,25)), bounds_error=False, fill_value=np.mean(mod))

        #model5 = si.RegularGridInterpolator(xs_fit, np.reshape(yld, (25,13,25)), bounds_error=False, fill_value=np.mean(yld))

        #f[0] = lambda x: model1(x)[0]

        #f[1] = lambda x: model2(x)[0]

        #f[2] = lambda x: model3(x)[0]

        #f[3] = lambda x: model4(x)[0]

        #f[4] = lambda x: model5(x)[0]



        # Approximately ?? seconds. Never initialized for some reason. It takes an extremely long amount of time to fit each model

        # phi1_fit = data_mu[:,0]

        # phi0_fit = data_mu[:,1]

        # phi2_fit = data_mu[:,2]

        # xs_fit = np.asarray((phi1_fit, phi0_fit, phi2_fit))

        # model1 = Pipeline([('poly', PolynomialFeatures(degree=2)),('linear', LinearRegression(fit_intercept=False))])

        # model2 = Pipeline([('poly', PolynomialFeatures(degree=2)),('linear', LinearRegression(fit_intercept=False))])

        # model3 = Pipeline([('poly', PolynomialFeatures(degree=2)),('linear', LinearRegression(fit_intercept=False))])

        # model4 = Pipeline([('poly', PolynomialFeatures(degree=2)),('linear', LinearRegression(fit_intercept=False))])

        # model5 = Pipeline([('poly', PolynomialFeatures(degree=2)),('linear', LinearRegression(fit_intercept=False))])

        # model1 = model1.fit(xs_fit, np.log(fip_cb))

        # model2 = model2.fit(xs_fit, np.log(fip_cs))

        # model3 = model3.fit(xs_fit, np.log(fip_cu))

        # model4 = model4.fit(xs_fit, np.log(mod))

        # model5 = model5.fit(xs_fit, np.log(yld))

        #

        # f = [0] * 5

        # f[0] = lambda x: model1.predict([x])[0]

        # f[1] = lambda x: model2.predict([x])[0]

        # f[2] = lambda x: model3.predict([x])[0]

        # f[3] = lambda x: model4.predict([x])[0]

        # f[4] = lambda x: model5.predict([x])[0]

        f = [[NumericFunction(f[0]), NumericFunction(f[1]), NumericFunction(f[2]), NumericFunction(f[3]), NumericFunction(f[4])]]

        # bound = [[-sys.maxsize-1, -sys.maxsize-1, -sys.maxsize-1, 130000.0, 1000.0], [0.0, 0.0, 0.0, sys.maxsize, sys.maxsize]]

        # bound = [[-10.0, -10.0, -10.0, 130000.0, 1000.0], [0.0, 0.0, 0.0, sys.maxsize, sys.maxsize]]

        bound = [[-sys.maxsize, -sys.maxsize, -sys.maxsize, 130000.0, 1000.0], [-9.0, -12.0, -12.0, sys.maxsize, sys.maxsize]]

        bound = PrismaticBoundary(bound)

        names = [r'$_{Ti64\beta}$', # Chart Title

        r'$\phi_{1} (rad)$', # 0-dir

        r'$\Phi (rad)$', # 1-dir

        r'$\phi_{2} (rad)$', # 2-dir

        ]

        feas_values, bound_numeric = idem(f, xs, dxs, bound, ignore_concavity=True, ignore_boundary=True)

        plot_combinations(bound_numeric.feasible_points, bnd=bound_numeric.boundary_points, names=names)

        ax = plt.gca()

        ax.azim = -34.96

        ax.elev = 30.0

        # labels = [ax.title, ax.xaxis.label, ax.yaxis.label, ax.zaxis.label]

        # for label in labels:

            # label.set_fontsize(24)

        ax.axes.xaxis.set_ticklabels([])

        ax.axes.yaxis.set_ticklabels([])

        ax.axes.zaxis.set_ticklabels([])



        plt.tight_layout(pad=1.5)

        plt.savefig('Ti64_feasible_revised', dpi=plt.gcf().dpi)

        plt.show()

        plt.close()

        

        pydem.pad_dist = old_pad

        pydem.f_size = old_font

        

    if test_holes:

        rad_lower = 5

        rad_upper = 10

        rad_int = 5

        deg_int = 11

        

        # num sphere dimensions, only 2 or 3 for now

        n = 2

        

        points = np.zeros((rad_int*deg_int**2,n))

        rads = np.linspace(rad_lower, rad_upper, rad_int)

        thetas = np.linspace(0, np.pi*2, deg_int)

        phis = np.linspace(-np.pi/2, np.pi/2, deg_int)

        for i in range(rad_int):

            for j in range(deg_int):

                r = rads[i]

                t = thetas[j]

                if n == 3:

                    for k in range(deg_int):

                        

                        p = phis[k]

                        points[i*deg_int**2 + j*deg_int + k] = [np.cos(t)*np.cos(p)*r, np.sin(t)*np.cos(p)*r, np.sin(p)*r]

                else:

                    points[i*deg_int + j] = [np.cos(t)*r, np.sin(t)*r]

        print((points.shape))

        temp = np.zeros((0,n))

        # bound_hole = ConcaveBoundary(points, temp)

        bound_hole = ConcaveBoundary(points, temp)

        if isinstance(bound_hole, ConcaveBoundary):

            centroids = bound_hole.simplex_centroids()

            excludes = np.linalg.norm(centroids, axis=1)

            excludes = excludes < rad_lower

            excludes = np.where(excludes)[0].flatten()

            bound_hole.exclude_simplices(excludes)

        points_test = np.zeros((20**n,n))

        for i in range(20):

            for j in range(20):

                if n == 3:

                    for k in range(20):

                        points_test[i*400+j*20+k] = [i-10, j-10, k-10]

                else:

                    points_test[i*20+j] = [i, j]

        interior_points = bound_hole.is_inner(points_test)

        interior = np.zeros((0,n))

        bound_dist = np.zeros((0,2,n))

        

        test_point = np.asarray([5,4])

        bound_hole.bound_dist(np.atleast_2d(test_point))

        # exit()

        

        for i in range(len(interior_points)):

            if interior_points[i]:

                temp_interior = points_test[i]


                interior = np.concatenate((interior, temp_interior[np.newaxis]), axis=0)

                
                temp = bound_hole.bound_dist(np.atleast_2d(temp_interior)) #LW These three lines cause typeerror's.
                temp = bound_hole.bound_dist(np.atleast_2d(temp_interior))
                temp = temp[np.newaxis,:,:]

                bound_dist = np.concatenate((bound_dist, temp), axis=0)

        base_names = ["", "", ""]

        for i in range(n):

            for j in range(2):

                val_name = "Distance to boundary in " + ("negative" if j==0 else "positive") + (" x" if i==0 else " y") + " direction"

                temp_names = base_names + [val_name]

                plot(interior, values=bound_dist[:,j,i], names=temp_names)

                plt.show()

                

    if test_projected_dists:

        # test various boolean operations with a MultiBoundary

        # test with a square from (0,0) to (1,1) and circle radius one centered at (0,0)

        rads = 1000

        feasible = np.asarray([[0.5*np.cos(i*2.0*np.pi/rads), 0.5*np.sin(i*2.0*np.pi/rads)] for i in range(rads)])

        feasible_2 = np.asarray([[np.cos(i*2.0*np.pi/rads), np.sin(i*2.0*np.pi/rads)] for i in range(rads)])

        feasible = np.concatenate((feasible, feasible_2), axis=0)

        boundary = np.asarray([[2*np.cos(i*2.0*np.pi/rads), 2*np.sin(i*2.0*np.pi/rads)] for i in range(rads)])

        boundary_1 = ConcaveBoundary(feasible, boundary)

        

        centroids = boundary_1.simplex_centroids()

        remove_centroids = []

        for i, c in enumerate(centroids):

            dist = np.linalg.norm(c)

            if dist > 0.5 and dist < 1:

                remove_centroids.append(i)

        boundary_1.exclude_simplices(remove_centroids)

        

        points_test = np.zeros((0,2))

        test_values = np.zeros((0,2,2))

        for j in np.linspace(-2,2,50):

            for k in np.linspace(-2,2,50):

                temp = np.asarray([[j,k]])

                if boundary_1.is_inner(temp)[0]:

                    points_test = np.concatenate((points_test, temp), axis=0)

                    output_range = np.concatenate((temp - .02, temp + .02), axis=0)

                    temp = boundary_1.projected_point_distance(output_range)

                    temp = np.reshape(temp, (1,2,2))

                    test_values = np.concatenate((test_values, temp), axis=0)

        for j in range(2):

            for k in range(2):

                plot(points_test, values=test_values[:,j,k])

                plt.show()

    

                

    if test_multi_volumes:

        # test various boolean operations with a MultiBoundary

        # test with a square from (0,0) to (1,1) and circle radius one centered at (0,0)

        bound = np.asarray([[0,0], [1.5,1.5]])

        boundary_2 = PrismaticBoundary(bound)

        print((boundary_2.is_inner(np.asarray([[.5,.5], [.5,2], [2,2]]))))

        rads = 100

        rad_1 = 0.25

        rad_2 = 0.5

        rad_3 = 1

        feasible = np.asarray([[rad_1*np.cos(i*2.0*np.pi/rads), rad_1*np.sin(i*2.0*np.pi/rads)] for i in range(rads)])

        feasible_2 = np.asarray([[rad_2*np.cos(i*2.0*np.pi/rads), rad_2*np.sin(i*2.0*np.pi/rads)] for i in range(rads)])

        feasible = np.concatenate((feasible, feasible_2), axis=0)

        boundary = np.asarray([[rad_3*np.cos(i*2.0*np.pi/rads), rad_3*np.sin(i*2.0*np.pi/rads)] for i in range(rads)])

        boundary_1 = ConcaveBoundary(feasible, boundary)

        

        # centroids = boundary_1.simplex_centroids()

        # remove_centroids = []

        # for i, c in enumerate(centroids):

            # dist = np.linalg.norm(c)

            # if dist > 0.5 and dist < 1:

                # remove_centroids.append(i)

        # boundary_1.exclude_simplices(remove_centroids)





        # test_points = np.zeros((0,2))

        # for j in np.linspace(-2,2,100):

        #     for k in np.linspace(-2,2,100):

        #         temp = np.asarray([[j,k]])

        #         if boundary_1.is_inner(temp)[0]:

        #             test_points = np.concatenate((test_points, temp), axis=0)

        # plot(test_points)

        # plt.show()





        # new_bounds = boundary_1.split_disconnected_volumes()

        # print(len(new_bounds))

        # for i in range(len(new_bounds)):

            # temp_bound = new_bounds[i]

            # test_points = np.zeros((0,2))

            # for j in np.linspace(-2,2,100):

                # for k in np.linspace(-2,2,100):

                    # temp = np.asarray([[j,k]])

                    # if temp_bound.is_inner(temp)[0]:

                        # test_points = np.concatenate((test_points, temp), axis=0)

            # plot(test_points)

            # plt.show()



        multi = MultiBoundary(boundary_1, boundary_2, 3)

        points_test = []

        for j in np.linspace(-2,2,100):

            for k in np.linspace(-2,2,100):

                temp = np.asarray([[j,k]])

                if multi.bound_2.is_inner(temp)[0]:

                    points_test.append(temp)

        points_test = np.concatenate(points_test, axis=0)

        plot(points_test)

        plt.show()

        plt.close()





        merge_types = [0,1,2]

        names = ['', 'x', 'y', 'Distance to boundary decreasing in x']

        fig_names = ['MultiBoundary_Intersect', 'MultiBoundary_Union', 'MultiBoundary_Difference']

        for i in merge_types:

            multi = MultiBoundary(boundary_1, boundary_2, i+1)

            points_test = []

            test_values = []

            for j in np.linspace(-2,2,50):

                for k in np.linspace(-2,2,50):

                    temp = np.asarray([[j,k]])

                    if multi.is_inner(temp)[0]:

                        points_test.append(temp)
#               
#                         temp = multi.bound_dist(temp, dim=0, dire=0) #LW These lines were causing a type error, "only integer scalar arrays can be converted to a scalar index"
# 
#                         temp = np.reshape(temp, (1,2,2))

                        test_values.append(temp)

            test_values = np.concatenate(test_values, axis=0)

            points_test = np.concatenate(points_test, axis=0)

            for j in range(1):

                for k in range(1):

                    ax = plot(points_test, values=test_values[:,j], names=names) #LW edit This was originally [:,j,k], but there were too many indices for the array.
                    #This may need to be changed back when other functions are fixed. 

                    ax.set_xlim([-2,2])

                    ax.set_ylim([-2,2])

                    fig = plt.gcf()

                    fig.axes[1].yaxis.label.set_fontsize(12)

                    plt.savefig(fig_names[i],dpi=fig.dpi)

                    plt.show()

                    plt.close()

               

    if test_expand_plotting:

        dim_range = 10

        num_dims = 4

        feas = np.meshgrid(*((np.arange(dim_range), )*4))

        new_feas = np.zeros((dim_range**num_dims,0))
        
        for i in range(len(feas)):

            new_feas = np.concatenate((new_feas, np.reshape(feas[i], (feas[i].size, 1))), axis=1)

        print((new_feas.shape))

        names = ["blah", "x1", "x2", "x3", "x4", "x5", "x6"]

        axis_ranges = [[0,5,10], [0,5,10], [0,10]]

        plot_expand_dims_subplot(new_feas, values=None, bnd=None, names=names, axis_ranges=axis_ranges)

        

    if test_sparse_input:

        nums = [7, 18, 6]

        mask = np.zeros(tuple(nums), dtype=bool)

        mask[:,9:12,:] = True

        test_ignore_volume = PrismaticBoundary(np.asarray([[5,50,30], [30,75,80]]))

        strength = np.linspace(10.0,22.0,nums[0])

        dissipatedenergy = np.linspace(20.0,105.0,nums[1])

        thickness = np.linspace(39.0,69.0,nums[2])

        print(dissipatedenergy)

        print(strength)

        print(thickness)

        ie1=0.1

        ie2=0.1

        ie3=0.1

        fe=0.1

        vars = sympy.symbols('STRENGTH DISSIPATEDENERGY THICKNESS')

        STRENGTH = vars[0]

        DISSIPATEDENERGY = vars[1]

        THICKNESS = vars[2]

        VAR=[STRENGTH,DISSIPATEDENERGY,THICKNESS]

        f    = -0.853260 + 0.0248455 * THICKNESS + 0.000808578 * THICKNESS * STRENGTH + 0.000391126 * THICKNESS * DISSIPATEDENERGY

        fmax = -0.728569 + 0.0242107 * THICKNESS + 0.000835134 * THICKNESS * STRENGTH + 0.000393725 * THICKNESS * DISSIPATEDENERGY

        fmin = -0.977951 + 0.0254804 * THICKNESS + 0.000782021 * THICKNESS * STRENGTH + 0.000388527 * THICKNESS * DISSIPATEDENERGY

        fa = [[SymbolicFunction(fmax, VAR)]]

        fa.append([SymbolicFunction(fmin, VAR)])

        fa.append([SymbolicFunction(f, VAR)])

        names = ['$HD_{EMI}:Level$ $0_{I=1.5 \\mathrm{MPa-ms}} | Regr_4$',  # Chart Title

                '$\sigma_t (\mathrm{MPa})$',             # xlabel

                '$E_{dissipated} (\mathrm{kJ}/\mathrm{m}^2)$',    # ylabel

                '$Thickness (\mathrm{mm})$']             # zlabel

        xs = [strength, dissipatedenergy, thickness]

        dxs = [ie1, ie2, ie3]

        bound = [[1.5], [sys.maxsize]] #LWEDIT maxint has been changed to maxsize in Python 3

        bound = PrismaticBoundary(bound)

        start_time = time.time()
 
        feas_values, bound_scale0 = idem(fa, xs, dxs, bound, em=em_chosen, ignore_region=test_ignore_volume)

        print((time.time()-start_time))

        ax = plot(bound_scale0.feasible_points, bnd=bound_scale0.boundary_points, names=names)

        # ax.azim = 29

        ax.azim = -119

        ax.elev = 22

        check_points = np.asarray([[15, 70, 65], [18, 70, 65], [18, 70, 60]])

        print((bound_scale0.is_inner(check_points)))

        plt.show()

        locations = []

        for x in np.linspace(10,22,20):

            for y in np.linspace(20,105,40):

                for z in np.linspace(39,69,20):

                    locations.append(np.atleast_2d([x,y,z]))

        locations = np.concatenate(locations, axis=0)

        valid = bound_scale0.is_inner(locations)

        valid_locs = locations[valid]

        print((valid_locs.shape))

        plot(valid_locs, names=names)

        plt.show() 



    if test_idce:

        wcmr = (np.arange(0.15, 0.45, 0.05))

        vp = (np.arange(0.01, 0.51, 0.05))

        rm = (np.arange(0.1, 30.1, 5))

        ft = 8

        ie1 = 0.05

        ie2 = 0.05

        ie3 = 0.1

        fe = 0.2

        vars = sympy.symbols('WCMR VP RM')

        WCMR = vars[0]

        VP = vars[1]

        RM = vars[2]

        VAR = [WCMR, VP, RM]

        f = [SymbolicFunction(0.177 * (99.3 / WCMR * (1 - VP) / (RM ** .5)) ** 0.74, VAR)]

        fmax = [SymbolicFunction(0.216 * (99.3 / WCMR * (1 - VP) / (RM ** .5)) ** 0.74, VAR)]

        fmin = [SymbolicFunction(0.144 * (99.3 / WCMR * (1 - VP) / (RM ** .5)) ** 0.74, VAR)]



        fa = [fmax]

        fa.append(fmin)

        fa.append(f)

        names = ['$Level$ $2_{\\itf_t = %d \\mathrm{MPa}}$' % ft,

                 '$w/cm$', '$V_p$', '$r_{pore} (\\mathrm{nm})$']

        xs = [wcmr, vp, rm]

        dxs = [ie1, ie2, ie3]

        bound = [[ft], [sys.maxsize]]

        bound = PrismaticBoundary(bound)

        temp_type = ['Concave']

        temp_concav = [False]



        fun_var = [[0, 1, 2]]



        for (bound_type, ig_conc) in zip(temp_type, temp_concav):

            feas_values, bound_scale2 = idem(fa, xs, dxs, bound, em=em_chosen, ignore_concavity=ig_conc,

                                             fun_var=fun_var)

        plot_combinations(bound_scale2.feasible_points, feas_values, bound_scale2.boundary_points, names)

        plt.show()



    if test_dynamic_bound:

        # build 1D, 2D and use UHPC as 3D test functions

        # 1D is special case because cannot build simplices

        # plot the distance functions and interior functions for known values



        # 2D

        inner_diam = 1

        outer_diam = 4

        search = outer_diam + 1

        search_inc = 0.25

        increment = 2 * search_inc

        searches = np.arange(-search, search + .001, search_inc)

        grid = np.arange(-outer_diam, outer_diam + .001, increment)

        x, y = np.meshgrid(grid, grid)

        feas = np.concatenate((np.atleast_2d(x.flatten()), np.atleast_2d(y.flatten())), axis=0)

        feas = np.transpose(feas)

        n = feas.shape[1]

        temp = np.linalg.norm(feas, axis=1)

        feas = feas[np.logical_and(temp >= inner_diam, temp <= outer_diam)]



        bnd = np.zeros((0, n))

        xs = np.unique(feas[:, 0])

        for x in xs[1:-1]:

            y1 = (outer_diam ** 2 - x ** 2) ** 0.5

            y2 = -y1

            bnd = np.concatenate((bnd, np.asarray([[x, y1], [x, y2]])))

        bnd = np.concatenate((bnd, bnd[:, [1, 0]]), axis=0)

        # plt.scatter(feas[:, 0], feas[:, 1])

        # plt.scatter(bnd[:, 0], bnd[:, 1], marker='*')

        # plt.show()



        feass = []

        dists = []

        bound = pydem.DynamicBoundary(feas, bnd, [grid, grid])

        for x in searches:

            for y in searches:

                point = np.asarray([[x,y]])

                if bound.is_inner(point)[0]:

                    if x == -3.5 and y == 0:

                        print('')

                    feass.append(point)

                    dists.append(bound.bound_dist(point)[1,1])

        feass = np.concatenate(feass, axis=0)



        plt.close('all')



        plt.figure(1)

        plt.scatter(feass[:, 0], feass[:, 1], c=dists, s=100)

        # plt.scatter(bnd[:, 0], bnd[:, 1], marker='*')

        plt.show()



    if test_concavity_theorem:

        points = np.zeros((8,3))

        points[:4,:] = np.asarray([[0,0,0], [0,1,0], [1,0,0], [1,1,0]])

        points[4:,:] = points[:4,:]

        for i in range(10000):

            points[4:,2] = np.random.rand(4)

            temp = ConvexHull(points)

            if temp.vertices.shape[0] != points.shape[0]:

                print('failed!')

        print(points)



        points = np.zeros((6,3))

        points[:3, :] = np.asarray([[0, 0, 0], [0, 1, 0], [1, 0, 0]])

        points[3:, :] = points[:3, :]

        for i in range(10000):

            points[3:, 2] = np.random.rand(3)

            temp = ConvexHull(points)

            if temp.vertices.shape[0] != points.shape[0]:

                print('failed!')

        print(points)



    print(('Total Execution time: ' + str(time.time() - total_start_time)))





if __name__ == "__main__":

    run_tests()