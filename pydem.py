import numpy as np #LW Is the current version installed?

import time

from scipy.spatial import ConvexHull #LW is the current version installed?

from scipy.spatial import Delaunay

import sympy

from scipy.spatial.qhull import QhullError

from scipy import interpolate


try:

    import matplotlib.pyplot as plt

    from mpl_toolkits.mplot3d import Axes3D

    from matplotlib.text import Text

    enable_plotting = True

except ImportError:

    enable_plotting = False

import numdifftools as nd  

import sys

import itertools

import scipy.interpolate as si

from matplotlib.colors import LinearSegmentedColormap

import collections



hdemi_thresh = 1

hdemi_tol = 0.005

ht_p = hdemi_thresh + hdemi_tol

ht_m = hdemi_thresh - hdemi_tol

cdict = {'red':   [(0.0,  0.0, 0.0),

                   (1.0,  157/255.0, 157/255.0)],



         'green': [(0.0,  0.0, 40/255.0),

                   (1.0,  123/255.0, 123/255.0)],



         'blue':  [(0.0,  0.0, 82/255.0),

                   (1.0,  72/255.0, 72/255.0)]}

custom_cmap = LinearSegmentedColormap('BlueGold', cdict)

pad_dist = 13

f_size = 16

base_fig_size = (6,4)

fig_dpi = 300





class Boundary(object):



    def __init__(self, *args, **kwargs):

        self.holes_enabled = False

        self.ignore_bounds = None

  #LW Why is there no return or raise here?  

    def set_ignore_bounds(self, bounds):

        """

        Allow user to set arbitrarily large distances in different directions



        Boundary searches will not occur in directions with non-zero values

        Shape is (2, n_dimensions). First row is for searches decreasing in the ith dimension for [0,i-1].



        Parameters

        ----------

        bounds : array_like

            values to return in each ignored direction (for non-zero values)



        Returns

        -------



        """

        bounds = np.array(bounds)

        self.ignore_bounds = bounds

     #LW Why is there no return or raise here?   

    def set_holes_enabled_value(self, holes_enabled):

        """

        Set boolean value for holes enabled in search.



        Enabling this feature will check all simplices along search path to ensure none are missing. Significantly slows code.



        Parameters

        ----------

        holes_enabled : bool

            If True assume holes exist within the feasible space



        Returns

        -------



        """

        self.holes_enabled = holes_enabled



    def is_inner(self, points):

        """



        Parameters

        ----------

        points : array_like

            List of points in n-space to check for location within this Boundary



        Returns

        -------

        inner_values : array

            List of values for each input point. True for those interior of Boundary, False otherwise



        """

        raise NotImplementedError()

        return [False]*len(points)

        

    def bound_dist(self, point, **kwargs):

        """



        Parameters

        ----------

        point : array_like

            Single point in the n-space to compute distances to Boundary.

        kwargs



        Returns

        -------

        dists : ndarray

            (2,n) array of distances. [0,i] is decreasing in i+1th dimension. [1,i] is decreasing in i+1th dimension



        """

        raise NotImplementedError()

        dists = np.zeros((2,point.shape[1])) - 1

        return dists



    def extents(self):

        """



        Returns

        -------

        extents : ndarray

            (2,n) array of the minimum and maximum extents of this Boundary in each of the n dimensions



        """

        raise NotImplementedError()



    def points(self):

        """



        Returns

        -------

        points : ndarray

            (m,n) array of m points associated with this Boundary which at the very least will enclose the convex hull.

        """

        raise NotImplementedError()



    def find_boundary(self, path):

        """

        Use bisection method to determine a point between two endpoints which lies on the boundary.

        Number of bisection iterations is currently fixed at 10.



        Parameters

        ----------

        path : array_like

            Two (or more) points which are the start and end points of the path to search along. Shape (n_points, n_dims)



        Returns

        -------

        bound_point : approximate location of boundary. If no intermediate boundary, returns the starting location



        """
        
        k_max = 10

        inners = self.is_inner(path)
        
        outsides = np.where(np.asarray(inners) == False)
 
        if len(outsides) < 1:

            print("No path points are outside the boundary, no intermediate boundary")

            return np.atleast_2d(path[0])

        if len(outsides) == len(path):

            print("No path points are inside the boundary, no intermediate boundary")

            return np.atleast_2d(path[0])

        if not inners[0]:

            path = path[::-1]

            inners = inners[::-1]
           
            outsides = np.where(np.asarray(inners) == False)

        if not inners[0]:

            print("No endpoint is inside the boundary, unsure of direction to proceed")

            return np.atleast_2d(path[0])

        

        u_lower = 0.0
        

        # use linear interpolation along path to make a little easier
    
        # tck, u = interpolate.splprep(np.transpose(path[:outsides[0]+1]), k=min([path.shape[0]-1,3]))
        
        
        tck, u = si.splprep(np.transpose(path), k=min([path.shape[0]-1,3]))
        
        #LW EDIT colon was causing type error with the indices. The +1 was out of range of the index.
        # K was not greater than path[outsides[0]] because this was only 1 and k = 1. It's more likely there were not outside points,
        # but this does return a figure. 
        

        # in binary search guarantee to be searching true->false for boundary
        u_upper = 1.0
        u_m = (u_lower + u_upper)/2.0

        eval_point = np.atleast_2d(np.asarray(si.splev(u_m, tck)))

        k = 0

        total = 0

        while k < k_max:

            if "simplex_neighbors" in dir(self) and self.holes_enabled:

                iters = 0

                k_sub = 0

                temp_low = np.atleast_2d(np.asarray(si.splev(u_lower, tck)))

                while True:

                    if iters >= k_max or self.simplex_neighbors(eval_point, temp_low):

                        break

                    u_m = (u_lower + u_m)/2.0

                    eval_point = np.atleast_2d(np.asarray(si.splev(u_m, tck)))

                    iters += 1

                    k_sub = 1

                k -= k_sub

            if self.is_inner(eval_point)[0]:

                u_lower = u_m

            else:

                u_upper = u_m

            u_m = (u_lower + u_upper)/2.0

            eval_point = np.atleast_2d(np.asarray(si.splev(u_m, tck)))

            if total > 100:

                # print(path)

                # print(u_lower)

                # print(u_m)

                # print(u_upper)

                # print("I was looping a lot in find_boundary")

                break

            k += 1

            total += 1

        return eval_point



    def invert_boundary(self, new_extents):

        """

        Construct a new Boundary which is the inverse of the current Boundary, i.e. is_inner will be the inverse of previous evaluations.



        Parameters

        ----------

        new_extents : array_like

            (2, n_dims) defines the hyperrectangular region which the new Boundary object will be valid over.



        Returns

        -------

        inv_bound : ConcaveBoundary

            The union of this new Boundary and the input bound will contain the hyperrectangular region of new_extents



        """

        new_points = self.points()

        temp = PrismaticBoundary(new_extents)

        new_points = np.concatenate((new_points, temp.points()), axis=0)

        inv_bound = ConcaveBoundary(new_points, np.zeros((0, new_points.shape[1])))

        simplex_centroids = inv_bound.simplex_centroids()

        inner_simplex = self.is_inner(simplex_centroids)

        exclude_simplices = np.arange(simplex_centroids.shape[0])

        exclude_simplices = exclude_simplices[inner_simplex]

        inv_bound.exclude_simplices(exclude_simplices)

        return inv_bound

        

        

class MultiBoundary(Boundary):

    """

    Support all Boundary operations using logical operations on multiple Boundary objects.

    Currently supported operations are union, intersection, and difference.

    May be used recursively to define any number of regions.

    """

    INTERSECT_TYPE = 1

    UNION_TYPE = 2

    DIFFERENCE_TYPE = 3

    

    def __init__(self, bound_1, bound_2, bool_type, *args, **kwargs):

        super(MultiBoundary, self).__init__(*args, **kwargs)

        self.bool_type = bool_type

        self.bound_1 = bound_1

        self.bound_2 = bound_2

        if self.bool_type == MultiBoundary.DIFFERENCE_TYPE:

            new_extents = self.extents()

            self.bound_2 = self.bound_2.invert_boundary(new_extents)



    def is_inner(self, points):

        inners_1 = self.bound_1.is_inner(points)

        inners_2 = self.bound_2.is_inner(points)

        if self.bool_type == MultiBoundary.INTERSECT_TYPE:

            return np.logical_and(inners_1, inners_2)

        elif self.bool_type == MultiBoundary.UNION_TYPE:

            return np.logical_or(inners_1, inners_2)

        elif self.bool_type == MultiBoundary.DIFFERENCE_TYPE:

            return np.logical_and(inners_1, inners_2)

        return [False]*len(points)



    def extents(self):

        ex_1 = self.bound_1.extents()

        ex_2 = self.bound_2.extents()

        mask = ex_1[1] < ex_2[1]

        ex_1[1][mask] = ex_2[1][mask]

        mask = ex_1[0] > ex_2[0]

        ex_1[0][mask] = ex_2[0][mask]

        return ex_1



    def points(self):

        return np.concatentate([self.bound_1.points(), self.bound_2.points()], axis=0)

    

    def bound_dist(self, point, dim=None, dire=None):

        dists_1 = self.bound_1.bound_dist(point, dim=dim, dire=dire)

        dists_2 = self.bound_2.bound_dist(point, dim=dim, dire=dire)

        

        if self.bool_type == MultiBoundary.INTERSECT_TYPE:

            mask = dists_1 > dists_2

            dists_1[mask] = dists_2[mask]

        elif self.bool_type == MultiBoundary.UNION_TYPE:

            n = point.shape[1]

            dists_1 = np.zeros(dists_1.shape)

            if dim is not None:

                space = [dim]

            else:

                space = range(n)

            if dire is not None:

                directions = [dire]

            else:

                directions = range(2)

            for d in space:

                for i in directions:

                    new_point = np.copy(point)

                    while True:

                        if self.bound_1.is_inner(new_point):

                            dists = self.bound_1.bound_dist(new_point, dim=d, dire=i)

                        elif self.bound_2.is_inner(new_point):

                            dists = self.bound_2.bound_dist(new_point, dim=d, dire=i)

                        else:

                            break

                        dists = dists[i]

                        if i == 0:

                            dists = -dists

                        new_point[0,d] += dists[d]

                    dists_1[i, d] = np.max(np.abs(new_point-point))

        elif self.bool_type == MultiBoundary.DIFFERENCE_TYPE:

            mask = dists_1 > dists_2

            dists_1[mask] = dists_2[mask]



        return dists_1





class PrismaticBoundary(Boundary):

    """

    Describe hyperrectangular regions in space for all Boundary functions.

    Is only defined between a minimum and maximum value along each dimension.

    """

    def __init__(self, *args, **kwargs):

        super(PrismaticBoundary, self).__init__(*args, **kwargs)

        new_bound = np.array(args[0])

        if new_bound.shape[0] != 2:

            raise ValueError("Prismatic bounds must be array-like (2,dims) to store upper and lower bounds")

        if np.any(new_bound[0] >= new_bound[1]):

            raise ValueError("Lower bound must be less than upper bound for all dimensions")

        self.bound = new_bound



    def is_inner(self, points):

        if len(points.shape) == 1:

            points = np.reshape(points, (1,len(points)))

        if points.shape[1] != self.bound.shape[1] :

            raise ValueError("dimension of points (%d) must match dimension of bounds (%d)" % (points.shape[1], self.bound.shape[1]))

        temp = np.logical_and(points > self.bound[0,:], points < self.bound[1,:])

        temp = np.all(temp, axis=1)

        return temp



    def extents(self):

        return np.copy(self.bound)



    def points(self):

        lims = self.extents()

        lims = np.transpose(lims)

        bound_vertex = np.meshgrid(*tuple(lims))

        bound_vertex = [np.atleast_2d(coords.flatten()) for coords in bound_vertex]

        bound_vertex = np.concatenate(bound_vertex, axis=0)

        bound_vertex = np.transpose(bound_vertex)

        return bound_vertex



    def bound_dist(self, point, **kwargs):

        dists = np.concatenate((point-self.bound[0,:], self.bound[1,:]-point), axis=0)

        return dists



            

class ConcaveBoundary(Boundary):

    """

    Boundary object with support for several additional features.

    Uses Delaunay Triagnulation and to reduce the higher dimensional space into simplices.

    Each simplex may then be removed from the hull using exclude_simplices to better describe the underlying hypervolume.

    """

    

    def __init__(self, feasible_points, boundary_points, *args, **kwargs):

        super(ConcaveBoundary, self).__init__(*args, **kwargs)

        feasible_points = np.copy(feasible_points)

        boundary_points = np.copy(boundary_points)

        bound = np.concatenate((feasible_points, boundary_points), axis=0)

        min_pct = np.min(np.max(bound, axis=0)-np.min(bound, axis=0))*.005

        # Qbb Qc Qz Qx 

        bound = Delaunay(bound, qhull_options="QJ%f" % min_pct)

        print("Finished Delaunay")

        # bound = Delaunay(bound)

        self.bound = bound

        self.simplex_is_inner = []

        self.init_inner_list()

        self.exclude_points = None

        self.boundary_points = boundary_points

        self.feasible_points = feasible_points



    def extents(self):

        extents = np.zeros((2,self.feasible_points.shape[1]))

        extents[0] = np.min(self.bound.points, axis=0)

        extents[1] = np.max(self.bound.points, axis=0)

        return extents



    def points(self):

        return self.bound.points



    def init_inner_list(self):

        self.simplex_is_inner = np.ones(self.bound.simplices.shape[0] + 1, dtype=bool)

        self.simplex_is_inner[-1] = False

    

    def simplex_centroids(self):

        nodes = self.bound.points

        elements = self.bound.simplices

        centroids = np.average(nodes[elements], axis=1)

        return centroids

            

    def exclude_simplices(self, simplex_list):

        self.simplex_is_inner[simplex_list] = False

        # for s in simplex_list:

        #     self.simplex_is_inner[s] = True

                    

    def make_exclude_points(self, ex_points):

        if self.exclude_points is None:

            self.exclude_points = ex_points

            self.init_inner_list()

        else:

            self.exclude_points = np.concatenate((self.exclude_points, ex_points), axis=0)

        locs = self.bound.find_simplex(ex_points)

        self.simplex_is_inner[locs] = False



    def is_inner(self, points, valid_simplices=None):

        simplices = self.bound.find_simplex(points)

        inners = self.simplex_is_inner[simplices]

        return inners

        

    def projected_point_distance(self, output_range, nominal_output=None):

        if nominal_output is None:

            nominal_output = np.mean(output_range, axis=0)

        n = output_range.shape[1]

        dists = np.zeros((2,n))

        all_points = self.boundary_points

        simplex_points = self.bound.simplices[np.logical_not(self.simplex_is_inner[:-1])].flatten()

        simplex_points = np.unique(simplex_points)

        hole_bound_points = self.bound.points[simplex_points]

        all_points = np.concatenate((all_points, hole_bound_points), axis=0)

        for d in range(n):

            valid_band_points = np.ones(all_points.shape[0], dtype=bool)

            for d2 in range(n-1):

                d2 = (d + d2 + 1) % n

                valid_band_points = np.logical_and(valid_band_points, all_points[:,d2] >= output_range[0,d2])

                valid_band_points = np.logical_and(valid_band_points, all_points[:,d2] <= output_range[1,d2])

            band_points = all_points[valid_band_points]

            band_dim = band_points[:,d]

            temp_upper = band_dim[band_dim > nominal_output[d]]

            if not len(temp_upper):

                dists[1, d] = -1

            else:

                dists[1, d] = np.min(temp_upper) - nominal_output[d]

            temp_lower = band_dim[band_dim < nominal_output[d]]

            if not len(temp_lower):

                dists[0, d] = -1

            else:

                dists[0, d] = nominal_output[d] - np.max(temp_lower)

            

        if self.ignore_bounds is not None:

            dists[np.where(self.ignore_bounds)] = sys.maxsize

            

        return dists

        

    def simplex_neighbors(self, point_1, point_2):

        simp_1 = self.bound.find_simplex(point_1)

        simp_2 = self.bound.find_simplex(point_2)

        if simp_1 == simp_2:

            return True

        n_list = self.bound.neighbors

        if simp_1 != -1:

            n_1 = n_list[simp_1]

        else:

            n_1 = []

        if simp_2 != -1:

            n_2 = n_list[simp_2]

        else:

            n_2 = []

        are_neighbors = simp_1 in n_2 or simp_2 in n_1

        return are_neighbors

                            

    def bound_dist(self, point, dim=None, dire=None):
       
        n = point.shape[1]

        extra = 0.01

        dim_min = np.min(self.bound.points, axis=0)

        dim_max = np.max(self.bound.points, axis=0)

        dim_range = dim_max - dim_min

        dim_max += dim_range*extra

        dim_min -= dim_range*extra

        dists = np.zeros((2,n))

        if not self.is_inner(point)[0]:

            dists[0,0] = -1

            return dists

        if dim is not None:

            space = [dim]

        else:

            space = range(n)

        if dire is not None:

            directions = [dire]

        else:

            directions = range(2)

        for d in space:

            d_new = [dim_min[d], dim_max[d]]

            for i in directions:

                if self.ignore_bounds is None or not self.ignore_bounds[i, d]:

                    x_2 = np.copy(point)

                    x_2[0,d] = d_new[i]

                    temp_points = np.concatenate((point, x_2), axis=0)
                    
                    x_bound = self.find_boundary(temp_points)

                    dists[i,d] = np.linalg.norm(point - x_bound)

                else:

                    dists[i,d] = sys.maxsize

        return np.abs(dists)

        

    def get_connected_volumes(self):

        neigh_list = np.copy(self.bound.neighbors)

        out_simplices = np.where(np.logical_not(self.simplex_is_inner[:-1]))[0]

        for s in out_simplices:

            neigh_list[s,:] = -1

            neigh_list[neigh_list==s] = -1

        print(len(neigh_list))

        trees = {}

        simplex_tree = [i for i in range(len(neigh_list))]

        for i in range(len(neigh_list)):

            trees[i] = {}

            trees[i][i] = True

        for i in range(len(neigh_list)):

            for j in range(len(neigh_list[i])):

                n = neigh_list[i,j]

                if n != -1 and n not in trees[simplex_tree[i]]:

                    neigh_tree = trees[simplex_tree[n]].keys()

                    for k in neigh_tree:

                        temp_tree = trees[simplex_tree[k]]

                        temp_tree.pop(k, None)

                        if len(temp_tree) == 0:

                            trees.pop(simplex_tree[k], None)

                        simplex_tree[k] = simplex_tree[i]

                        trees[simplex_tree[i]][k] = True



        for s in out_simplices:

            trees.pop(int(s), None)

        return trees

        

    def split_disconnected_volumes(self):

        """



        :return:

        """

        trees = self.get_connected_volumes()

        if len(trees) == 1:

            return [self]

        points = self.bound.points

        new_bounds = []

        print("*************************************************")

        print("%d disconnected volumes detected" % len(trees))

        for key in trees:

            simplices = trees[key].keys()

            temp_points = points[np.unique(self.bound.simplices[simplices])]

            print(temp_points.shape)

            print("Boundary contains %d points from the initial hull" % len(temp_points))

            new_bound = ConcaveBoundary(temp_points, np.zeros((0,temp_points.shape[1])))

            centroids = new_bound.simplex_centroids()

            valid_centroids = self.is_inner(centroids, trees[key])

            invalid_simplices = np.where(np.logical_not(valid_centroids))[0]

            new_bound.exclude_simplices(invalid_simplices)

            new_bounds.append(new_bound)

        print("*************************************************")

        return new_bounds



    def get_external_faces(self):

        ex_faces = np.zeros((0,self.bound.points.shape[1]))

        # can speed this up by vectorizing, but lazy for now

        neighbors = self.bound.neighbors

        simplices = self.bound.simplices

        for i in range(self.bound.simplices.shape[0]):

            inner_simplex = self.simplex_is_inner[i]

            temp_simplex = simplices[i]

            for j in range(neighbors.shape[1]):

                n = neighbors[i,j]

                already_explored = n < i and n > -1

                not_included = not self.simplex_is_inner[n]

                if inner_simplex and not already_explored and not_included:

                    ex_faces = np.concatenate((ex_faces, [np.concatenate((temp_simplex[:j], temp_simplex[(j+1):]))]), axis=0)

        return ex_faces





class DynamicBoundary(Boundary):

    def __init__(self, feasible_points, boundary_points, xs, *args, **kwargs):

        super(DynamicBoundary, self).__init__(*args, **kwargs)

        feasible_points = np.copy(feasible_points)

        boundary_points = np.copy(boundary_points)

        self.boundary_points = boundary_points

        self.feasible_points = feasible_points



        self.xs = [np.array(xs[i]) for i in range(len(xs))]

        n = feasible_points.shape[1]

        type_list = [feasible_points, boundary_points]

        grid_map = {}



        for k in range(len(type_list)):

            temp_points = type_list[k]

            for i in range(temp_points.shape[0]):

                temp_point = temp_points[i]

                indices = [np.where(np.logical_and(temp_point[d] >= xs[d][0:-1], temp_point[d] <= xs[d][1:]))[0] for d in range(n)]

                indices = np.meshgrid(*tuple(indices))

                indices = [np.atleast_2d(coords.flatten()) for coords in indices]

                indices = np.concatenate(indices, axis=0)

                indices = np.transpose(indices)

                for index in indices:

                    index = tuple(index)

                    if index not in grid_map:

                        grid_map[index] = [[], []]

                    grid_map[index][k].append(i)

        self.grid_map = grid_map



    def is_inner(self, points):

        inners = []

        n = points.shape[1]

        for point in points:

            indices = [np.where(np.logical_and(point[d] >= self.xs[d][0:-1], point[d] <= self.xs[d][1:]))[0] for d in

                       range(n)]

            indices = np.meshgrid(*tuple(indices))

            indices = [np.atleast_2d(coords.flatten()) for coords in indices]

            indices = np.concatenate(indices, axis=0)

            indices = np.transpose(indices)

            is_in = False

            for index in indices: # this loop should only every execute once. Either points is inside single grid location, else in multiple and each should return true

                # verify previous statement

                index = tuple(index)

                if self.is_full_grid_ind(index):

                    is_in = True

                    break

                elif index in self.grid_map:

                    grid_points = self.grid_map[index]

                    grid_points = np.concatenate((self.feasible_points[grid_points[0]], self.boundary_points[grid_points[1]]) , axis=0)

                    if grid_points.shape[0] >= n + 1:

                        try:

                            temp_hull = Delaunay(grid_points)

                            loc = temp_hull.find_simplex(point)

                            is_in = loc != -1

                            if is_in:

                                break

                        except QhullError:

                            is_in = False

            inners.append(is_in)

        return inners



    def is_gridded(self, points, ignore_dim=[]):

        points = np.asarray(points)

        grids = np.ones(points.shape[0], dtype=bool)

        search_dims = set(range(points.shape[1])) - set(ignore_dim)

        for d in search_dims:

            grids = np.logical_and(grids, np.any(points[:,d]==self.xs[d],axis=1))

        return grids



    def is_full_grid_ind(self, grid_index):

        if grid_index in self.grid_map:

            if len(self.grid_map[grid_index][0]) == 2**self.feasible_points.shape[1]:

                return True

            return False

        return False



    def bound_int(self, point, ray, grid_index):

        """

        Compute boundary intersection with the search point and ray. Assumes that we are in a grid index with overall

        bounndary, but will not error otherwise.

        Parameters

        ----------

        point

        ray

        grid_index



        Returns

        -------



        """

        n = point.size

        if grid_index not in self.grid_map:

            return np.zeros(point.shape)

        vol_points = self.feasible_points[self.grid_map[grid_index][0], :]



        if self.boundary_points.size:

            vol_points = np.concatenate((vol_points, self.boundary_points[self.grid_map[grid_index][1], :]), axis=0)



        # degenerate cases to allow function to work when query point is on a grid face

        remain_dims = list(range(n)) #LWEDIT change range object to list
        search_dim = np.where(ray)[0][0]

        lost_dims = []

        for d in (list(range(n-1, search_dim, -1)) + list(range(search_dim-1, -1, -1))): #LWEIDT Python 3 does return ranges lists. List function added.

            is_dim = np.any(point[d] == self.xs[d])

            if is_dim:

                remain_dims.pop(d)

                lost_dims.append(d)



        if len(remain_dims) == 1:

            grid_dims = list(range(search_dim)) + list(range(search_dim+1,n)) #LWEDIT Changes range object to list object.

            dist = 0

            pt = np.copy(point)

            for temp in vol_points:

                temp_dist = (temp[search_dim] - point[search_dim])*ray[search_dim]

                if np.all(temp[grid_dims] == point[grid_dims]) and temp_dist > dist:

                    dist = temp_dist

                    pt = temp

            return pt



        retain_points = []

        for i in range(vol_points.shape[0]):

            temp_point = vol_points[i]

            if np.all(np.abs(temp_point[lost_dims] - point[lost_dims]) < 2**-30):

                retain_points.append(i)

        vol_points = vol_points[:,remain_dims]

        vol_points = vol_points[retain_points]

        orig_point = np.copy(point)

        point = orig_point[remain_dims]

        ray = ray[remain_dims]

        orig_search_dim = search_dim

        search_dim = np.where(ray)[0][0]

        n = vol_points.shape[1]



        if vol_points.shape[0] >= n + 1:

            try:

                temp_hull = ConvexHull(vol_points)

                temp_delaunay = Delaunay(vol_points)

                dist = 0

                pt = np.copy(point)

                for facet in temp_hull.equations:

                    if abs(np.dot(ray, facet[:-1])) > 2**-30:

                        temp_int = point + ray*(-1*facet[-1]-np.dot(point, facet[:-1]))/np.dot(ray,facet[:-1])

                        temp_simplex = temp_delaunay.find_simplex(temp_int)

                        temp_dist = (temp_int[search_dim] - point[search_dim])*ray[search_dim]

                        if temp_simplex > -1 and temp_dist > dist:

                            dist = temp_dist

                            pt = temp_int

                orig_point[remain_dims] = pt

                return orig_point

            except QhullError:

                raise RuntimeWarning('Convex hull not well formed for bound_int search')

        ray_int = orig_point

        dire = np.sign(ray[search_dim])

        if dire > 0:

            remains = vol_points[:,search_dim] >= point[search_dim]

        else:

            remains = vol_points[:,search_dim] <= point[search_dim]

        vol_points = vol_points[remains]

        ray_int[orig_search_dim] = dire * np.amin(vol_points[:,search_dim] * dire)

        return ray_int



    def bound_dist(self, point, dim=None, dire=None):

        """

           Note that we can get a very very good idea of absolute minimum over the area desired by checking all exterior corners

           and the minimum distance of the points contained within the region desired (HDEMI projection)

        """

        point = np.copy(point)

        point = np.atleast_2d(point)

        n = point.shape[1]

        dists = np.zeros((2,n))

        if not self.is_inner(point)[0]:

            dists[0,0] = -1

            return dists

        point = point.flatten()

        ignore_search = np.zeros((2,n), dtype=bool)

        if dim is not None:

            space = [dim]

        else:

            space = range(n)

        if dire is not None:

            directions = [dire]

        else:

            directions = range(2)



        indices = [np.where(np.logical_and(point[d] >= self.xs[d][0:-1], point[d] <= self.xs[d][1:]))[0] for d in

                   range(n)]

        start_index = [indices[d][0] for d in range(n)]



        for d in space:

            for i in directions:

                if self.ignore_bounds is None or not  self.ignore_bounds[i,d]:

                    grid_index = start_index[:]

                    ray = np.zeros((n,))

                    ray[d] = 2*i - 1

                    bound_pt = np.copy(point)

                    explore_bound = True

                    while explore_bound:

                        if grid_index[d] < 0 or grid_index[d] >= (len(self.xs[d]) - 1):

                            explore_bound = False

                            break

                        elif self.is_full_grid_ind(tuple(grid_index)):

                            bound_pt[d] = self.xs[d][grid_index[d] + i]

                        else:

                            bound_pt = self.bound_int(point, ray, tuple(grid_index))

                            if abs(self.xs[d][grid_index[d] + i] - bound_pt[d]) > 2**-30:

                                explore_bound = False

                                break

                        grid_index[d] += 2*i - 1

                    dists[i,d] = np.linalg.norm(bound_pt - point)

                else:

                    dists[i,d] = sys.maxsize

        return np.abs(dists)



class AnonymousFunction:

    def get_y(self, x):

        """



        Parameters

        ----------

        x : array_like

            point (or list of points) to evaluate based on nominal input configurations



        Returns

        -------

        y : array_like

            evaluates objective function to perform mapping for the associated dimension in output space

            process must be repeated over multiple AnonymousFunction objects to accomplish full mapping



        """

        raise NotImplementedError()

        return 0



    def get_dy(self, x, dx):

        """



        Parameters

        ----------

        x : array_like

            Nominal input configuration to evaluate

        dx : array_like

            Symmetric variation along each axis based on uncertainty



        Returns

        -------

        dy : float

            Maximum potential variation based on the dot product of the absolute value of the gradient and the absolute value of the dx vector



        """

        raise NotImplementedError()

        return 0





class NumericFunction(AnonymousFunction):



    """

    This class is an implementation of AnonymousFunction which uses numdifftools to numerically differentiate objective functions.

    """

    def __init__(self, f):

        """



        Parameters

        ----------

        f : function

            Single objective function which accepts a list of values as an input

        """

        self.function = f

        self.d_function = nd.Gradient(f)

    

    def get_y(self, x):

        """

        Parameters

        ----------

        x : array_like

            point (or list of points) to evaluate based on nominal input configurations

            While the generic AnonymousFunction accepts a list of points indexed by the first axis, we only evaluate the first point in that list currently



        Returns

        -------

        y : array_like

            output of the objective function (mapping) for the assigned output dimension



        """

        return self.function(x[0])

    

    def get_dy(self, x, dx):

        """

        Parameters

        ----------

        x : array_like

            Nominal input configuration to evaluate

            Note that we must maintain consistency with a (1,n) point being passed in but gradient evaluation will fail

            if numdifftools sees a 2d point being passed in so we pass in the first row

        dx : array_like

            Symmetric variation along each axis based on uncertainty



        Returns

        -------

        dy : float

            Maximum potential variation based on the dot product of the absolute value of the gradient and the absolute value of the dx vector

        """

        partials = np.abs(self.d_function(x[0]))

        dx = np.abs(dx)

        return np.dot(partials.flatten(), dx.flatten())

    

    

class SymbolicFunction(AnonymousFunction):

    """

    This class is an implementation of AnonymousFunction which uses sympy to evaluate and differentiate symbolic objective functions.

    """

    def __init__(self, f, symbols):

        """



        Parameters

        ----------

        f : sympy expression

            Construct using symbol objects, e.g. f = -0.853260 + 0.0248455 * THICKNESS + 0.000808578 * THICKNESS * STRENGTH + 0.000391126 * THICKNESS * DISSIPATEDENERGY

        symbols : list of sympy symbols

            These can be constructed using the sympy.symbols function, e.g. vars = sympy.symbols('STRENGTH DISSIPATEDENERGY THICKNESS')

        """

        dims = len(symbols)
        
        

        self.function = sympy.lambdify(symbols, f, modules='numpy')

        self.derivatives = [0]*dims

        for i in range(dims):

            d_f = sympy.diff(f,symbols[i])

            self.derivatives[i] = sympy.lambdify(symbols, d_f, modules='numpy')
            
        
        
        

        

    def get_y(self, x):
    # def __repr__(self,x):

        """



        Parameters

        ----------

        x : array_like

            point (or list of points) to evaluate based on nominal input configurations



        Returns

        -------

        y : array_like

            evaluates objective function to perform mapping for the associated dimension in output space

            process must be repeated over multiple AnonymousFunction objects to accomplish full mapping



        """
       
        
        x = np.atleast_2d(x)

        temp_inp = np.split(x, x.shape[1], axis=1)
        
       

        return self.function(*temp_inp) #Lw added return

    

    def get_dy(self, x, dx):

        """



        Parameters

        ----------

        x : array_like

            Nominal input configuration to evaluate

        dx : array_like

            Symmetric variation along each axis based on uncertainty



        Returns

        -------

        dy : float

            Maximum potential variation based on the dot product of the absolute value of the gradient and the absolute value of the dx vector



        """

        d_f = 0

        dx = dx.flatten().tolist()

        for d in range(len(dx)):

            d_f += np.abs(self.derivatives[d](*x[0]) * dx[d])

        return d_f



def n_dim_norm(points):

    """

    Create normal vector for n-dimensional

    Parameters

    ----------

    points



    Returns

    -------



    """

    points = np.asarray(points)

    n = points.shape[1]

    norm = np.zeros((1,n))

    if n != points.shape[0]:

        raise RuntimeWarning('Invalid number of points provided for normal function')

        return norm

    vecs = points[1:,:] - points[0,:]

    for d in range(n):

        norm[d] = (-1)**(d)*np.linalg.det(vecs[range(d) + range(d+1,n)])

    return norm

    

def vor(y, dy, bound):

    """

    Valid Output Region

    Possible replacement error margin for the Hyperdimensional error margin index. Evaluates the outer extents to see if all of the output region (expressed as a hyperrectangle) lies within the feasible space, and is thus a robust space.

    Could currently fail due to gaps in the feasible space which are contained within the extents of the output region. Can be addressed efficiently by means of KD-tree or similar.



    Parameters

    ----------

    y : ndarray

        Bounding outputs to evaluate, (n_functions, n_dims)

        Nominal output should always be last.

    dy : ndarray

        Variation along each output dimension (n_functions, n_dims)

        Nominal variation should always be last

    bound : Boundary

        Feasible region representation



    Returns

    -------

    error_margin : float

        Uses same convention as HDemi for compatibility of logic. -1 if not in feasible region, 1.5 otherwise



    """

    n = y.shape[1]

    fs = y.shape[0]

    o_lim = np.zeros((n,2))



    for d in range(n):

        poss_max = []

        poss_min = []

        for i in range(fs):

            poss_max.append(y[i,d] + dy[i,d])

            poss_min.append(y[i,d] - dy[i,d])

        o_lim[d,:] = np.asarray([np.min(poss_min),np.max(poss_max)])

    

    # construct all outermost extents + the center of output range

    eval_points = np.meshgrid(*tuple(o_lim))

    eval_points = [np.atleast_2d(coords.flatten()) for coords in eval_points]

    eval_points = np.concatenate(eval_points, axis=0)

    eval_points = np.transpose(eval_points)

    temp = np.zeros((1,n))

    temp[0] = y[-1]

    eval_points = np.concatenate((eval_points, temp), axis=0)

    

    valid_outputs = bound.is_inner(eval_points)

    

    if np.all(valid_outputs):

        return 1.5

    return -1

        

        

def miv(y, dy, bound):

    """

    Maximum Independent Variation

    Possible replacement error margin for the Hyperdimensional error margin index. Evaluates the mid-points of the outer extents to see if all of the output region (expressed as a hyperrectangle) lies within the feasible space, and is thus a robust space.

    Could currently fail due to gaps in the feasible space which are contained within the extents of the output region. Can be addressed efficiently by means of KD-tree or similar.



    Parameters

    ----------

    y : ndarray

        Bounding outputs to evaluate, (n_functions, n_dims)

        Nominal output should always be last.

    dy : ndarray

        Variation along each output dimension (n_functions, n_dims)

        Nominal variation should always be last

    bound : Boundary

        Feasible region representation



    Returns

    -------

    error_margin : float

        Uses same convention as HDemi for compatibility of logic. -1 if not in feasible region, 1.5 otherwise



    """

    n = y.shape[1]

    fs = y.shape[0]

    eval_points = np.zeros((2*n+1,n))

    for d in range(n):

        poss_max = []

        poss_min = []

        for i in range(fs):

            poss_max.append(y[i,d] + dy[i,d])

            poss_min.append(y[i,d] - dy[i,d])

        eval_points[:,d] = y[-1,d]

        eval_points[d*2+1,d] = np.min(poss_min)

        eval_points[d*2+2,d] = np.max(poss_max)

    valid_outputs = bound.is_inner(eval_points)

    if np.all(valid_outputs):

        return 1.5

    return -1





def hdemi(y, dy, bound, projection=False):

    """

    Compute the Hyperdimensional error margin index. See "An Inductive Design Exploration Method for Robust Multiscale Materials Design" Choi, McDowell, Allen, Rosen, Mistree.

    Requires the distance computations from Boundary objects which may, or may not, be subject to significant errors if holes exist, etc.

    Distance computations are also relatively expensive, so recommend exploring mlv or vor alternatives.



    Parameters

    ----------

    y : ndarray

        Bounding outputs to evaluate, (n_functions, n_dims)

        Nominal output should always be last.

    dy : ndarray

        Variation along each output dimension (n_functions, n_dims)

        Nominal variation should always be last

    bound : Boundary

        Feasible region representation



    Returns

    -------

    error_margin : float

        Return -1 if value is within the feasible region, otherwise compute the full hyperdimensional error margin index which is the minimum ratio of the output range distance to the distance to the boundary along that dimension.

        HD_{EMI} = \begin{cases} min_{i} \left[ \frac{\left\Vert \bar{y}-b_i \right\Vert}{\Delta y_{i}} \right], & \mbox{for } \bar{y} \in \Omega \\ -1, & \mbox{for } \bar{y} \notin \Omega \end{cases}



    """
   
    n = y.shape[1]

    fs = y.shape[0]

    o_dist = np.zeros((2,n))
 
    
    for d in range(n):

        poss_max = []

        poss_min = []

        for i in range(fs):

            poss_max.append(y[i,d] + dy[i,d]) #(0,0), (1,0), (2,0), (0,1), (1,1)...

            poss_min.append(y[i,d] - dy[i,d])

        o_dist[:,d] = np.asarray([np.min(poss_min),np.max(poss_max)])

    if projection and "projected_point_distance" in dir(bound):

        dists = bound.projected_point_distance(o_dist, nominal_output=y[None,-1,:])

    else:

        dists = bound.bound_dist(y[None,-1,:])

    o_dist[0,:] -= y[-1,:]

    o_dist[1,:] -= y[-1,:]

    

    div_0_mask = o_dist == 0

    o_dist[div_0_mask] = dists[div_0_mask]/10000.0

    if np.all(dists>=0):

        return np.min(np.abs(dists/o_dist))

    return -1



def get_y(f, x, x_err, fun_var=None, y_db=None, xs=None):

    n_bnd = len(f)

    n_out = len(f[0])



    eval_fun_var = y_db is not None and fun_var is not None and xs is not None



    y = np.zeros((len(f), n_out))

    dy = np.zeros(y.shape)

    y2 = np.zeros(y.shape)

    dy2 = np.zeros(y.shape)

    for d in range(n_out):



        failed = False

        if eval_fun_var:

            temp = [np.where(x[0,i] == xs[i])[0] for i in fun_var[d]]

            ind_db = [temp[i][0] if temp[i].size>0 else -1 for i in range(len(temp))]

            if np.any([ind_db[i]<0 for i in range(len(ind_db))]):

                failed = True

            else:

                ind_db = tuple([slice(None), slice(None)] + ind_db)

                temp = y_db[d][ind_db]

                y[:, d] = temp[:, 0]

                dy[:, d] = temp[:, 1]

        if not eval_fun_var or failed:

            dxs = x_err * x[0]

            for j in range(n_bnd):

                anon_f = f[j][d]

                y[j, d] = anon_f.get_y(x)

                dy[j, d] = abs(anon_f.get_dy(x, dxs))



    return y, dy



def find_boundary(f, x_1, x_2, x_err, bound, f_0=[], em=hdemi, fun_var=None, y_db=None, xs=None):

    """

    Find an approximation of the boundary representing the feasible region. Uses bisection search between a known interior point and exterior point.



    Parameters

    ----------

    f : array_like

        Mapping functions (number of bounding functions with nominal value last, number of output variables)

    x_1 : ndarray

        Robust nominal input to search between

    x_2 : ndarray

        Non-robust nominal input to search between

    x_err : ndarray

        Relative uncertainty in each dimension

    bound : Boundary

        Feasible region in the output space

    f_0 : list, optional

        Previously computed evaluations of f array for both x_1 and x_2 to reduce total computation slightly

    em : function, optional

        Error margin function which may be evaluated em(f, x, dx, bound) and return expected values < 1 for non-robust solutions and values >= 1 for robust solutions.

    fun_var : list of lists, optional

        List of which variables are used in each mapping functions. Interior lists will be integer valued.

        Used to minimize the number of function evaluations performed. If provided, each unique combination of inputs will only be evaluated once.

    y_db : array_like, optional

        List of ndarray. Unique combinations of gridded function values for each function. Must be provided with fun_var

    xs : list of lists, optional

        Discrete values to sample in each dimension. Unless ingore_region is specified, the full factorial combination space of these values will be explored.



    Returns

    -------

    x_m : ndarray

       Approximation of the boundary point based on 10 bisection iterations or error margin values in range (1-tolerance, 1+tolerance).



    """

    k_max = 7

    if f_0:

        h_1 = f_0[0]

        h_2 = f_0[1]

    else:

        y, dy = get_y(f, x_1, x_err, fun_var, y_db, xs)

        h_1 = em(y, dy, bound)

        y, dy = get_y(f, x_2, x_err, fun_var, y_db, xs)

        h_2 = em(y, dy, bound)

    # error check outer ranges just in case...

    if h_1 < ht_p and h_1 > ht_m:

        return x_1, None

    elif h_2 < ht_p and h_2 > ht_m:

        return x_2, None

    x_m = (x_1 + x_2)/2.0

    for k in range(k_max):

        y, dy = get_y(f, x_m, x_err, fun_var, y_db, xs)

        h_m = em(y, dy, bound)

        if h_m < ht_p and h_m > ht_m:

            break

        elif (h_1 > ht_p and h_m > ht_p) or (h_1 < ht_m and h_m < ht_m):

            x_1 = x_m

            h_1 = h_m

        else:

            x_2 = x_m

            h_2 = h_m

        x_m = (x_1 + x_2)/2.0

    if h_1 < ht_m:

        exclude_point = x_1

    else:

        exclude_point = x_2

    return x_m, exclude_point





def idem(f, xs, x_err, objective_bound, ignore_concavity=False, ignore_boundary=False, em=vor, ignore_region=None, fun_var=None):

    """

    Inductive Design Exploration Method. See "An Inductive Design Exploration Method for Robust Multiscale Materials Design" Choi, McDowell, Allen, Rosen, Mistree.

    Perform all critical steps to result in a feasible region boundary for the current design level based on constraints from the previous level, mapping functions, and uncertainties.

    Initially screen a discretized input space based on xs which will satisfy the robust design threshold.

    Several additional checks may be completed such as finding the boundary between the discrete feasible points and the infeasible region along each dimension,

    potential concavities in the feasible space, and ignoring potentially feasible combinations due to additional constraints.



    Parameters

    ----------

    f : array_like

        Mapping functions (number of bounding functions with nominal value last, number of output variables)

    xs : list of lists

        Discrete values to sample in each dimension. Unless ingore_region is specified, the full factorial combination space of these values will be explored.

    x_err : list

        Relative uncertainty for each input dimension

    objective_bound : Boundary

        Describes the region which satisfies design criteria for the previous level (current output region)

    ignore_concavity : bool, optional

        Controls the logic flow for evaluating concave regions. Compute the feasibility of the centroid of each simplex in the Delaunay Triangulation.

        If centroid is feasible, simplex remains in the feasible space, else it is removes and is_inner will evaluate to False within this simplex.

        Defaults to False (evaluate concavity).

    ignore_boundary : bool, optional

        Controls the logic flow for evaluating boundary points (HDemi ~ 1). If False, will execute explore_boundary and probe for boundary locations between feasible and infeasible configurations.

        Defaults to False (evaluate boundaries).

    em : function, optional

        Error margin function which may be evaluated em(f, x, dx, bound) and return expected values < 1 for non-robust solutions and values >= 1 for robust solutions.

    ignore_region : Boundary, optional

        Initial configurations within this region will be ignored for feasibility computation.

        Once simplices are constructed for ConcaveBoundary of the current feasible space, centroids are examined to ensure they do not lie in this ignore_region.

    fun_var : list of lists, optional

        List of which variables are used in each mapping functions. Interior lists will be integer valued.

        Used to minimize the number of function evaluations performed. If provided, each unique combination of inputs will only be evaluated once.



    Returns

    -------

    feas_values : ndarray

        (n_feas_points) shape array containing the error margin values for the evaluated feasible points

    bound : ConcaveBoundary

        ConcaveBoundary made up of the feasible and boundary points (if found). Concavities may or may not be removed based on selected

    """



    n = len(xs)

    n_bnd = len(f)

    n_out = len(f[0])



    eval_fun_var = fun_var is not None and len(fun_var) == n_out

    y_db = None

    if eval_fun_var:



        y_db = []

        for i in range(n_out):

            sub_ind = [len(xs[j]) for j in fun_var[i]]

            db_shape = [n_bnd, 2] + sub_ind

            temp_db = np.zeros(tuple(db_shape))

            sub_num = np.prod(sub_ind)

            x = np.zeros((1, n))

            for j in range(n_bnd):

                anon_f = f[j][i]

                for k in range(sub_num):

                    index = np.unravel_index(k, sub_ind)

                    for ii in range(len(fun_var[i])):

                        inp_d = fun_var[i][ii]

                        x[0,inp_d] = xs[inp_d][index[ii]]

                    dxs = x_err * x

                    y = anon_f.get_y(x)

                    dy = abs(anon_f.get_dy(x, dxs))

                    temp_db[(j, 0) + index] = y

                    temp_db[(j, 1) + index] = dy

            y_db.append(temp_db)



    feasible = find_feasible(f, xs, x_err, objective_bound, em=em, ignore_region=ignore_region, fun_var=fun_var, y_db=y_db)

    feas = feasible.points

    feas_values = feasible.robustness



    if ignore_boundary:

        print('Ignoring boundary computations')

        bnd = np.zeros((0,n))

    else:

        bound_vals = explore_boundary(f, xs, x_err, objective_bound, feasible.explored_robustness, em=em, ignore_region=ignore_region, fun_var=fun_var, y_db=y_db)

        bnd = bound_vals.points



    # if len(feas) > 0 and n > 1:

    #     start_time = time.time()

    #     boundary = ConcaveBoundary(feas, bnd)

    #     if ignore_concavity:

    #         print('Not evaluating concavity')

    #     else:

    #         fix_concavity(f, x_err, objective_bound, boundary, em=em, ignore_region=ignore_region)

    #     print('Finished constructing feasible hull, elapsed time=%02.02f seconds' % (time.time()-start_time))

    #     return feas_values, boundary

    # elif len(feas) > 0:

    #     bound = np.zeros((2,n))

    #     temp = np.concatenate((feas, bnd), axis=0)

    #     bound[0,:] = np.min(temp, axis=0)

    #     bound[1,:] = np.max(temp, axis=0)

    #     bound = PrismaticBoundary(bound)

    #     bound.feasible_points = feas

    #     bound.boundary_points = bnd

    #     return feas_values, bound

    if len(feas) > 0:

        start_time = time.time()

        boundary = DynamicBoundary(feas, bnd, xs)

        print('Finished constructing feasible hull, elapsed time=%02.02f seconds' % (time.time()-start_time))

        return feas_values, boundary

    else:

        return feas_values, None

        



def find_feasible(f, xs, x_err, objective_bound, em=vor, ignore_region=None, fun_var=None, y_db=None):

    """

    Perform the initial step of IDEM which requires the evaluation of the discretized region using an error metric.

    Parameters

    ----------

    f : array_like

        Mapping functions (number of bounding functions with nominal value last, number of output variables)

    xs : list of lists

        Discrete values to sample in each dimension. Unless ingore_region is specified, the full factorial combination space of these values will be explored.

    x_err : list

        Relative uncertainty for each input dimension

    objective_bound : Boundary

        Describes the region which satisfies design criteria for the previous level (current output region)

    em : function

        Error margin function which may be evaluated em(f, x, dx, bound) and return expected values < 1 for non-robust solutions and values >= 1 for robust solutions.

    ignore_region : Boundary, optional

        Initial configurations within this region will be ignored for feasibility computation.

    fun_var : list of lists, optional

        List of which variables are used in each mapping functions. Interior lists will be integer valued.

        Used to minimize the number of function evaluations performed. If provided, each unique combination of inputs will only be evaluated once.

    y_db : array_like

        List of ndarray. Unique combinations of gridded function values for each function. Must be provided with fun_var



    Returns

    -------

    points : ndarray

        Array of points in n-dimensional space (n_points, n_dims) for those configurations found to be robust

    robustness : ndarray

        Array of robustness values for all robust solutions based on the error margin computed.

    explored_robustness : ndarray

        Total array of all explored points and their computed error margin values.

    """

    start_time = time.time()

    feas_dims = tuple(map(len, xs))

    n = len(xs)

    robust = np.zeros(np.prod(feas_dims))

    x_err = np.array(x_err) #this turns the list x_err into an array. 

    x_err = np.reshape(x_err, (1,n)) #LW added the .ma. to update command, didn't seem to make a difference though



    n_bnd = len(f)

    n_out = len(f[0])



    # build database of function values and variation for combinations of inputs

    # (n_bound_funs, 2 (y, dy), n_inp_1, n_inp_2, ...)



    # evaluate feasible region

    for i in range(len(robust)):

        index = np.unravel_index(i, feas_dims)



        x = np.zeros((1,n))

        for j in range(n):

            x[0,j] = (xs[j][index[j]])

        if ignore_region is None or not ignore_region.is_inner(x):

            y, dy = get_y(f, x, x_err, fun_var, y_db, xs)



            temp = em(y, dy, objective_bound)

            robust[i] = temp

        else:

            robust[i] = -1

        if len(robust) > 10 and i%(len(robust)/10)==0: #LW Is this supposed to be float or integer division? 

            print('Feasible: %02.0f%% complete, elapsed time=%02.02f seconds' % (i/float(len(robust))*100, time.time()-start_time))

    print('Finished finding feasible points, elapsed time=%02.02f seconds' % (time.time()-start_time))        

    feas_mask = robust > ht_p

    print('%10d Feasible points found' % feas_mask.sum())

    indices = np.arange(len(robust))

    feas_index = indices[feas_mask]

#     feas_index = np.unravel_index(feas, feas_dims) #LW error: UnboundLocal Error

    total_indices = tuple(xs)

    locs = np.meshgrid(*total_indices,indexing='ij')

    new_locs = np.zeros((locs[0].size, len(locs)))

    for i in range(len(locs)):

        new_locs[:,i] = locs[i].flatten()

    feas_robust = robust[feas_mask]

    feas = new_locs[feas_index,:]



    Feasible = collections.namedtuple('Feasible', 'points robustness explored_robustness')

    f = Feasible(feas, feas_robust, robust)

    return f





def explore_boundary(f, xs, x_err, objective_bound, x_robust, em=vor, ignore_region=None, fun_var=None, y_db=None):

    """

    Examine robust configurations with non-robust neighbors in all dimensions.

    If any such neighbor pair is found, the  boundaries are explored using find_boundary to solve for the point where the error margin approximates the true boundary within some tolerance.

    This function can require significant additional computational effort depending on the initial gridding density and the dimensionality.



    Parameters

    ----------

    f : array_like

        Mapping functions (number of bounding functions with nominal value last, number of output variables)

    xs : list of lists

        Discrete values to sample in each dimension. Unless ingore_region is specified, the full factorial combination space of these values will be explored.

    x_err : list

        Relative uncertainty for each input dimension

    objective_bound : Boundary

        Describes the region which satisfies design criteria for the previous level (current output region)

    em : function

        Error margin function which may be evaluated em(f, x, dx, bound) and return expected values < 1 for non-robust solutions and values >= 1 for robust solutions.

    ignore_region : Boundary, optional

        Do not explore boundary between infeasible_region feasible points and robust configurations.

    fun_var : list of lists, optional

        List of which variables are used in each mapping functions. Interior lists will be integer valued.

        Used to minimize the number of function evaluations performed. If provided, each unique combination of inputs will only be evaluated once.

    y_db : array_like

        List of ndarray. Unique combinations of gridded function values for each function. Must be provided with fun_var





    Returns

    -------

    points : ndarray

        (n_bound_points, n_dim) Array of the explored boundary points, all within tolerance of the true boundary location based on error metric.



    """

    feas_dims = tuple(map(len, xs))

    start_time = time.time()

    n = len(xs)

    bnd = []

    for k in range(len(x_robust)):

        ind_1 = np.unravel_index(k, feas_dims)

        x_1 = np.reshape([xs[i][ind_1[i]] for i in range(n)], (1,n))
        f_1 = x_robust[k]
        
        
        if f_1 <= ht_p and f_1 > ht_m:

            bnd.append(x_1)
        

        elif f_1 > ht_p:

            for i in range(n):

                dim_diffs = np.zeros((n,2))

                dim_diffs[i,:] = np.asarray([[-1,1]])

                for j in range(2):

                    ind_2 = tuple(np.asarray(ind_1)+dim_diffs[:,j])

                    ind_2 = tuple(map(int, ind_2)) #LW needs list or tupble before it. map object won't work for if statement

                    if np.min(ind_2) > -1 and ind_2[i] < feas_dims[i]:

                        flat_ind = np.ravel_multi_index(ind_2, feas_dims)

                        f_2 = x_robust[flat_ind]

                        x_2 = np.asarray([xs[ii][ind_2[ii]] for ii in range(n)])

                        if f_2 < ht_m and (ignore_region is None or not ignore_region.is_inner(x_2)):

                            found_bound, outside_point = find_boundary(f, x_1, x_2, x_err, objective_bound, f_0=[f_1, f_2], em=em, fun_var=fun_var, y_db=y_db, xs=xs)

                            found_bound = np.reshape(found_bound, (1,n))

                            bnd.append(found_bound)
                            
                            

        if len(x_robust) > 10 and k%(len(x_robust)/10)==0:

            print('Boundary: %02.0f%% complete, time=%02.02f' % (k/float(len(x_robust))*100, time.time()-start_time))

    print('Finished finding boundary points, elapsed time=%02.02f seconds' % (time.time()-start_time))

    print('%6.00f Boundary points found' % len(bnd))
    

    bnd = np.concatenate(bnd, axis=0)



    Bound = collections.namedtuple("Bound", "points") # idea to later add the objective function values here

    b = Bound(bnd)

    return b





def fix_concavity(f, x_err, objective_bound, boundary, em=vor, ignore_region=None, plot_concave=False):

    """

    Examine all simplex centroids in the Delaunay Triangulation to determine if they lie within feasible regions.

    Obviously still limited in terms of finding infeasible, concave regions due to fineness of meshing related to initial discretization.

    Larger performance hit in higher dimensional space due to the larger number of simplices constructed.



    Parameters

    ----------

    f : array_like

        Mapping functions (number of bounding functions with nominal value last, number of output variables)

    x_err : list

        Relative uncertainty for each input dimension

    objective_bound : Boundary

        Describes the region which satisfies design criteria for the previous level (current output region)

    boundary : ConcaveBoundary

        Triangulation to remove simplexes from using the aforementioned criteria

    em : function

        Error margin function which may be evaluated em(f, x, dx, bound) and return expected values < 1 for non-robust solutions and values >= 1 for robust solutions.

    ignore_region : Boundary, optional

        Perform additional filtering if this region is supplied to ensure that is_inner for the constructed boundary does not overlap with the ignored region.

    plot_concave : bool, optional

        Plot the resultant space using a gridded view to ensure proper construction of space. Default False



    Returns

    -------



    """

    start_time = time.time()

    centroids = boundary.simplex_centroids()

    print('%d simplex centroids to evaluate for concavity feasibility' % len(centroids))

    remove_simplex = []



    n_bnd = len(f)

    n_out = len(f[0])



    for i, c in enumerate(centroids):

        y, dy = get_y(f, [c], x_err)

        if em(y, dy, objective_bound) < ht_m or (ignore_region is not None and ignore_region.is_inner(c)):

            remove_simplex.append(i)

        if len(centroids) > 10 and i%(len(centroids)/10)==0:

            print('Concavity Check: %02.0f%% complete, time=%02.02f' % (i/float(len(centroids))*100, time.time()-start_time))

    print('Will remove %d simplices due to concavity' % len(remove_simplex))

    boundary.exclude_simplices(remove_simplex)

    if plot_concave and enable_plotting:

        points = centroids[remove_simplex]

        plot(points)

        plt.show()



        

def get_sub_dim_range(feas, axis_ranges, range_index, values=None, bnd=None):

    """

    Get points within ranges along each dimension (hyperrectangular region) for use in plotting functions.

    Note the hyperrectangle is only defined in the first n dimensions, where n is the length of axis_ranges.

    Dimensions not specified are assumed to be valid on :math:'\\left[-\\infty, \\infty\\right]'



    Parameters

    ----------

    feas : ndarray

        (n_points, n_dim) array of feasible points to plot

    axis_ranges : array_like

        list of lists describing the ranges possible along each axis. Ranges are expressed by adjacent indices.

        i.e. feas[1,0] is valid if the value is in range [axis_ranges[0][0], axis_ranges[0][1])

    range_index : list

        list of integer indices into the axis_ranges array for the valid combination of axis_values to use.

    values : list, optional

        robustness (or other) values associated with the feasible points.

    bnd : array_like, optional





    Returns

    -------

    feas : ndarray

        indexed values of the feas array input meeting criteria

    values : ndarray

        indexed values of the values array input matching the feasible points which satisfied criteria

        None if no values were supplied

    bnd : ndarray

        indexed values of the bnd array input meeting criteria



    """

    feas_mask = np.ones(feas.shape[0])

    exp_dims = len(axis_ranges)

    for j in range(exp_dims):

        feas_mask = np.logical_and(feas_mask, feas[:,j] >= axis_ranges[j][range_index[j]])

        feas_mask = np.logical_and(feas_mask, feas[:,j] < axis_ranges[j][range_index[j]+1])

    feas = feas[feas_mask,exp_dims:]

    if values is not None:

        values = values[feas_mask]

    else:

        values = None

    if bnd is not None:

        bnd_mask = np.ones(bnd.shape[0])

        for j in range(exp_dims):

            bnd_mask = np.logical_and(bnd_mask, bnd[:,j] >= axis_ranges[j][range_index[j]])

            bnd_mask = np.logical_and(bnd_mask, bnd[:,j] < axis_ranges[j][range_index[j]+1])

        bnd = bnd[bnd_mask, exp_dims:]

    else:

        bnd = None

    return feas, values, bnd



    

def get_sub_dim_discrete(feas, axis_values, value_index, values=None, bnd=None):

    """

    Get points with discrete values (not ranges) to plot.

    Select points whose location in the first n dimensions is exactly that of the criteria described by axis_values and value_index.



    Parameters

    ----------

    feas : ndarray

        (n_points, n_dim) array of feasible points to plot

    axis_values : array_like

        list of lists describing the discrete values possible along each axis.

        The length of the first dimension of axis_values is the number of dimensions to select for exact values (flatten into separate plots) while plotting.

        values are compared exactly, use get_sub_dim_range for values which may fall inside a range.

    value_index : list

        list of integer indices into the axis_values array for the valid combination of axis_values to use.

        All combinations of feas[:,i] == axis_values[i][value_index[i]] must be true for points to be selected.

        Must be of same length or longer than axis_values.

    values : list, optional

        robustness (or other) values associated with the feasible points.

    bnd : array_like, optional





    Returns

    -------

    feas : ndarray

        indexed values of the feas array input meeting criteria

    values : ndarray

        indexed values of the values array input matching the feasible points which satisfied criteria

        None if no values were supplied

    bnd : ndarray

        indexed values of the bnd array input meeting criteria



    """

    feas_mask = np.ones(feas.shape[0])

    exp_dims = len(axis_values)

    for j in range(exp_dims):

        feas_mask = np.logical_and(feas_mask, feas[:,j] == axis_values[j][value_index[j]])

    feas = feas[feas_mask,exp_dims:]

    if values is not None:

        values = np.asarray(values)

        values = values[feas_mask]

    else:

        values = None

    if bnd is not None:

        bnd_mask = np.ones(bnd.shape[0])

        for j in range(exp_dims):

            bnd_mask = np.logical_and(bnd_mask, bnd[:,j] == axis_values[j][value_index[j]])

        bnd = np.asarray(bnd)

        bnd = bnd[bnd_mask, exp_dims:]

    else:

        bnd = None

    return feas, values, bnd

    

    

def plot_expand_dims_subplot(feas, values=None, bnd=None, names=[], axis_ranges=[], discrete_dims=False):

    """

    Best method so far to plot higher dimensional spaces. Expanded dimensions are those who are not plotted on the individual subplot, but rather expanded over a sequence of subplots.

    That is, a 2D grid of subplots may be used to plot the first two dimensions of a series of points, with each subplot being used to plot the remainder of dimensions.

    Currently support expanding 3 dimensions, that is create multiple grid images that may be flipped through in sequence to provide 3-dimensional plot expansion.

    Always saves plots since they typically are of such resolution as to not display nicely on the user's screen.



    Parameters

    ----------

    feas : ndarray

        (n_points, n_dim) array of feasible points to plot

    values : list, optional

        List of values associated with the feasible points. Colorbar will be created to show variation of these values.

    bnd : ndarray, optional

        (n_bound_points, n_dim) array of boundary points. Will always be plotted in black to differentiate from feasible points.

    names : list, optional

        List of various plotting names that must be supplied. In order; plot title, 1-axis, 2-axis (if exists), 3-axis (if exists), 4-axis (if exists), 5-axis (if exists), ... n-axis (if exists), and values (if exists).

    axis_ranges : list of lists, optional

        First list iterates over axes to be expanded in plotting. Internal lists provide values to be used for plotting these axis.

        e.g. [[0,1], [2,4,6]] will expand the first two  and create a 2x3 grid of subplots. If discrete_dims is False, feas and bnd points between 0 and 1 on the first axis and between 2 and 4 on the second axis on subplot (2,3,1)

        Default to [], which shortcircuits to normal plotting.

    discrete_dims : bool, optional

        If True, each value in axis_ranges is treated as a discrete selection from feas and bnd to plot, not as ranges. Not recommened to use if boundary points are supplied since these are rarely found along discrete values. Defaults to False.

    Returns

    -------



    """

    exp_dims = len(axis_ranges)

    if exp_dims == 0:

        plot(feas,values,bnd,names)

        plt.show()

    if exp_dims > 3:

        raise ValueError("Currently only support expansion of three dimensions or less")
    
    if discrete_dims:

        dim_num_subplots = np.asarray(list(map(len, axis_ranges))) #LW edit fixes length error in line 3320.

    else:

        dim_num_subplots = np.asarray([len(axis_ranges[i])-1 for i in range(exp_dims)])
    
    y_subplots = dim_num_subplots[-2] if len(dim_num_subplots) > 1 else 1

    x_subplots = dim_num_subplots[-1]

    x_size = base_fig_size[0]*x_subplots

    y_size = base_fig_size[1]*y_subplots

    prev_fig_num = 0

    num_plots = np.prod(dim_num_subplots)

    

    if discrete_dims:

        col_names = [names[exp_dims] + " $ = %.3e$" %

            (axis_ranges[-1][j]) for j in range(dim_num_subplots[-1])]

    else:

        col_names = [names[exp_dims] + " $ \in [%.3e,%.3e)$" %

            (axis_ranges[-1][j], axis_ranges[-1][j+1]) for j in range(dim_num_subplots[-1])]

    if exp_dims > 2:

        if discrete_dims:

            fig_names = [names[1] + " $ = %.3e$" % 

                (axis_ranges[0][j]) for j in range(dim_num_subplots[0])]

        else:

            fig_names = [names[1] + " $ \in [%.3e,%.3e)$" % 

                (axis_ranges[0][j], axis_ranges[0][j+1])  for j in range(dim_num_subplots[0]-1)]

    else:

        fig_names = ['']

    if exp_dims > 1:

        if discrete_dims:

            row_names = [names[exp_dims-1] + " $ = %.3e$" %

                (axis_ranges[-2][j]) for j in range(dim_num_subplots[-2])]

        else:

            row_names = [names[exp_dims-1] + " $\in [%.3e,%.3e)$" %

                (axis_ranges[-2][j], axis_ranges[-2][j+1]) for j in range(dim_num_subplots[-2])]

    else:

        row_names = ['']

    

    for i in range(num_plots):

        range_index = np.unravel_index(i, dim_num_subplots)

        fig_num = 0 if len(axis_ranges) < 3 else range_index[0]

        sub_num = i - fig_num*np.prod(dim_num_subplots[1:]) + 1

        subplot = (y_subplots, x_subplots, sub_num)

        if discrete_dims:

            new_feas, new_values, new_bnd = get_sub_dim_discrete(feas, axis_ranges, range_index, values=values, bnd=bnd)

        else:

            new_feas, new_values, new_bnd = get_sub_dim_range(feas, axis_ranges, range_index, values=values, bnd=bnd)

        if sub_num == 1:

            fig = plt.figure(fig_num, figsize=(x_size, y_size), dpi=fig_dpi)

        new_names = [names[0]]

        new_names.extend(names[(1+exp_dims):])

        ax = plot(new_feas, new_values, new_bnd, new_names, fig, subplot=subplot)

        sub_num -= 1

        if sub_num % x_subplots == 0:

            ax.annotate(row_names[sub_num//x_subplots], xy=(0,.5), xytext=(-6*pad_dist, 0), #LW EDIT Single division sign changed to double for float/integer division.

                xycoords='axes fraction', textcoords='offset points',

                ha='center', va='center', rotation=90, fontsize=1.5*f_size)

        if sub_num / x_subplots == 0:

            ax.annotate(col_names[sub_num], xy=(.5,1), xytext=(0,pad_dist*2),

                xycoords='axes fraction', textcoords='offset points',

                ha='center', va='center', fontsize=1.5*f_size)

            

        if prev_fig_num != fig_num or i == (num_plots-1):

            prev_fig = plt.figure(prev_fig_num)

            # add figure title

            if exp_dims == 3:

               
                prev_fig.text(.5, .93, fig_names[prev_fig_num], transform=prev_fig.transFigure, fontsize=2*f_size, ha='center', va='bottom')

            prev_fig.tight_layout()

            prev_fig.subplots_adjust(left=.1, top=.85)

            prev_fig.savefig("expanded_dims_fig_%d" % prev_fig_num, dpi=fig.dpi)

            plt.close(prev_fig)

#             prev_fig_num = fig_num #LW This line causes test_expand_plotting to not function. 



def plot_combinations(feas, values=None, bnd=None, names=[], max_dim=3, save_fig=False):

    """

    Plot various combinations of input variables. Take points in higher dimensions space, restrict them to max_dim number of variables in several lower dimension plots.

    Additional dimensions will be flattened into the lower dimensional plot.

    e.g. 5 dimensional points plotted in a maximum 3D plot will produce 10 unique plots displaying the combinations of 3 variables on a 3D plot to display the 5D points provided.

    Lower dimensions plots are made using the plot function.



    Parameters

    ----------

    feas : ndarray

        (n_points, n_dim) array of feasible points to plot

    values : list, optional

        List of values associated with the feasible points. Colorbar will be created to show variation of these values.

    bnd : ndarray, optional

        (n_bound_points, n_dim) array of boundary points. Will always be plotted in black to differentiate from feasible points.

    names : list, optional

        List of various plotting names that must be supplied. In order; plot title, 1-axis, 2-axis (if exists), 3-axis (if exists), 4-axis (if exists), 5-axis (if exists), and values (if exists).

        All values which are not necessary collapse the order, e.g. a labelled 2-axis plot with values would be ordered [1-axis label, 2-axis label, value label (will be assigned to Colorbar).

    max_dim : int, optional

        Controls the number of dimensions that will be shown on each plot

        default is a 3D plot

    save_fig : bool, optional

        Control whether to save the figures or plot to the first n figures where n is the number of plot combinations

        Saved figures will be named according to standard conventions for this module

        default False

    Returns

    -------



    """

    """ Plot all possible combinations of max_dim number of dimensions.

        The remaining dimensions are projected onto the lower dimensional subspace """

        

    n = feas.shape[1]

    if len(bnd) > 0 and bnd.shape[1] != n:

        raise ValueError("boundary and feasible region must exist in same n-dimensional space")

    

    if n <= max_dim:

        plot(feas, values, bnd, names, save_fig=save_fig)

    else:

        dims = range(n)

        comb = itertools.combinations(dims, max_dim)

        for i, c in enumerate(comb):

            print(c)

            print(i)

            feas_sub = feas[:,c]

            bnd_sub = bnd

            if len(bnd) > 0:

                bnd_sub = bnd[:,c]

            if len(names) > n:

                names_sub = [0]*(max_dim+1)

                names_sub[0] = names[0]

                for j in range(max_dim):

                    names_sub[j+1] = names[c[j]+1]

            else:

                names_sub = []

            fig = plt.figure(i, figsize=base_fig_size, dpi=fig_dpi)

            plot(feas_sub, values, bnd_sub, names_sub, fig=fig, save_fig=save_fig)





def size_scale(values, s_min, s_max):

    """



    Parameters

    ----------

    values : ndarray

        values to be displayed using the size of the scatter points

    s_min : float

        minimum value this set of values should be compared to

    s_max : float

        maximum value this set of values should be compared to



    Returns

    -------

    sizes : ndarray

        arbitrary scaling of values which should be appropriate for linearly scaled data and be visually distinct



    """

    return 30 + 200*(values-s_min)/(s_max-s_min)





def plot(feas, values=None, bnd=None, names=None, fig=None, save_fig=False, subplot=(1,1,1), disp_cbar=True):

    """

    Base plotting function for IDEM related visulaizations. Incorporates plotting of feasible points and boundary points.

    Feasible points may be colored based on objective values or robustness measures, for example.



    Plotting will automatically adjust from 2D to 3D plots. For points in n-dimensional space where n < 3, all plots are 2D. Colors are assigned based on values, if provided.

    Plotting in 3D will automatically adjust to display up to 5 dimensions, where the 4th is shown as size, and the 5th as color. This is not advised for gridded spaces with n > 3 since points will overlap in 3D coordinates.



    Parameters

    ----------

    feas : ndarray

        (n_points, n_dim) array of feasible points to plot. n_dim <= 5

    values : list, optional

        List of values associated with the feasible points. Colorbar will be created to show variation of these values.

    bnd : ndarray, optional

        (n_bound_points, n_dim) array of boundary points. Will always be plotted in black to differentiate from feasible points.

    names : list, optional

        List of various plotting names that must be supplied. In order; plot title, 1-axis, 2-axis (if exists), 3-axis (if exists), 4-axis (if exists), 5-axis (if exists), and values (if exists).

        All values which are not necessary collapse the order, e.g. a labelled 2-axis plot with values would be ordered [1-axis label, 2-axis label, value label (will be assigned to Colorbar).

    fig : matplotlib.figure or int, optional

        Preexisting figure instance or figure number to plot on. If not supplied, a new figure will be created.

    save_fig : bool, optional

        Save the figure when done. Default to False.

    subplot : location, optional

        See matplotlib documentation for valid subplot location indices. Default to plotting over the entire figure, subplot (1,1,1)

    disp_cbar : bool, optional

        Determine if the colobar will be displayed (only used if point values ar provided). Default to True.

    Returns

    -------

    ax : matplotlib.axes

        Axis that was used to draw the current plot for further manipulation

    """

    """ names in order (title, axis 1, ..., axis n, value) """

    norm_size = 30

    

    fig_num = 1

    if type(fig) == int:

        fig_num = fig

        fig = None
        



    if fig is None:

        fig = plt.figure(fig_num, figsize=base_fig_size, dpi=fig_dpi)

        fig.clf()
        


    if not enable_plotting:

        print("plotting has been disabled, please ensure matplotlib is installed correctly")

        return

    if len(feas.shape) < 2:

        raise ValueError("Feasible points must be in array_like with second dimension = spatial dimension")

    if feas.shape[0] == 0:

        print("No data to plot")

        ax = fig.add_subplot(*subplot, facecolor='white') #LW changes "axisbg" to facecolor

        ax.get_xaxis().set_visible(False)

        ax.get_yaxis().set_visible(False)

        ax.patch.set_visible(False)

        ax.axis('off')

        return ax


    values = np.array([]) if values is None else np.asarray(values)

    bnd = np.zeros((0,feas.shape[1])) if bnd is None else np.asarray(bnd)

    names = [] if names is None else names



    n = feas.shape[1] 

    # cm = plt.cm.get_cmap('spring')

    # cm = custom_cmap

    cm = plt.cm.get_cmap('viridis')

    cbar = None

    



    if bnd.shape[0] > 0 and bnd.shape[1] != n:

        raise ValueError("boundary and feasible region must exist in same n-dimensional space")



    if n < 2:

        feas = np.reshape(feas, (len(feas), 1))

        feas = np.concatenate((feas, np.ones(feas.shape)), axis=1)

        bnd = np.reshape(bnd, (len(bnd), 1))

        bnd = np.concatenate((bnd, np.ones(bnd.shape)), axis=1)

    if n < 3:

        ax = fig.add_subplot(*subplot, facecolor='white') #LW changed "axisbg" to facecolor

        if len(values) == feas.shape[0]:

            colors = values

            if len(names) > 3:

                clabel = names[3]

            else:

                clabel = '$HD_{EMI}$'

        else:

            colors = 'b'

        im = ax.scatter(feas[:,0], feas[:,1], c=colors, cmap=cm, edgecolors='none')

        if colors != 'b' and disp_cbar: #LWedit changed is not to !=

            cbar = fig.colorbar(im, ax=ax, pad=0.15)

            cbar.set_label(clabel)

            im.set_clim(np.min(values), np.max(values))

        if len(bnd) > 0:
         
            ax.scatter(bnd[:,0], bnd[:,1], marker='^', c='k')

        if len(names) > 0:

            ax.set_title(names[0], fontsize=f_size*1.25)

        if len(names) > 1:

            ax.set_xlabel(names[1], fontsize=f_size, labelpad=pad_dist)

        if len(names) > 2:

            ax.set_ylabel(names[2], fontsize=f_size, labelpad=pad_dist)

        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):

            item.set_fontsize(f_size)

            item.set_color("k")

        

    elif n < 6:

        ax = fig.add_subplot(*subplot, projection='3d', facecolor='white') #LW changes "axisbg" to facecolor

        #plt.hold('on') #LW 'hold' has been removed from Python 3. The default is on. 

        if len(values) > 0:

            c_min = np.min(values)

            c_max = np.max(values)

        else:

            c_min = 0

            c_max = 1

        if n < 5 and len(values) == feas.shape[0]:

            colors = values

            if len(names) > (n+1):

                clabel = names[n+1]

            else:

                clabel = '$HD_{EMI}$'

        elif n == 5:

            # get min and max values for color dimension

            c_min = np.min(feas[:,4])

            c_max = np.max(feas[:,4])
            

            if len(bnd) > 0:

                c_min = min(c_min, np.min(bnd[:,4]))

                c_max = max(c_min, np.max(bnd[:,4]))

            colors = feas[:,4]

            clabel = names[5]

        else:

            colors = 'b'

        if n > 3:

            s_min = np.min(feas[:,3])

            s_max = np.max(feas[:,3])

            if len(bnd) > 0:

                s_min = min(s_min, np.min(bnd[:,3]))

                s_max = max(s_min, np.max(bnd[:,3]))

            sizes = size_scale(feas[:,3], s_min, s_max)

        else:

            sizes = norm_size
        
        # im = ax.scatter(feas[:,0], feas[:,1], feas[:,2], s=sizes, c=colors, vmin=c_min, vmax=c_max, edgecolor='none', cmap=cm, alpha=1)
        im = ax.scatter(feas[:,0], feas[:,1], feas[:,2], s=sizes) #LWEDIT
        
        if colors != 'b' and disp_cbar: #LWedit changed "is not" to !=

            cbar = fig.colorbar(im, ax=ax, pad=.15)

            cbar.set_label(clabel)

            im.set_clim(c_min, c_max)

        if len(bnd) > 0:

            if n == 5:

                colors = bnd[:,4]

            else:

                colors = 'k'

            if n > 3:

                sizes = size_scale(bnd[:,3], s_min, s_max)

            else:

                sizes = norm_size
            
            # ax.scatter(bnd[:,0], bnd[:,1], bnd[:,2], marker='^', s=sizes, c=colors, vmin=c_min, vmax=c_max, edgecolor='none', cmap=cm, alpha=1)
            ax.scatter(bnd[:,0], bnd[:,1], bnd[:,2],s=sizes) #LWEDIT
        if len(names) > 0:

            ax.set_title(names[0], fontsize=f_size*1.25)

        if len(names) > 1:

            ax.set_xlabel(names[1], fontsize=f_size, labelpad=pad_dist)

        if len(names) > 2:

            ax.set_ylabel(names[2], fontsize=f_size, labelpad=pad_dist)

        if len(names) > 3:

            ax.set_zlabel(names[3], fontsize=f_size, labelpad=pad_dist)

#         plt.hold('off') #LW, the default is "on", but hold no longer functions. Do we need an equivalent coding for this? 

    else:

        raise ValueError("Plots of greater than 5 dimensions are not supported")

    

    to_change_labels = ax.get_xticklabels() + ax.get_yticklabels()

    if n > 2:

        to_change_labels += ax.get_zticklabels()

    if cbar is not None:

        to_change_labels += cbar.ax.get_yticklabels()

        for temp_text in cbar.ax.findobj(match=Text, include_self=False):

            to_change_labels.append(temp_text)

    for item in (to_change_labels):

        item.set_fontsize(f_size)

        item.set_color("k")

    ax.tick_params(labelsize=f_size*0.75)

    fig.tight_layout()

    if save_fig:
        
        plt.savefig('feasible_plot_fig_%d' % fig.number, dpi=fig.dpi)

    return ax





def kriging_upper(gp, x):

    y, mse = gp.predict(x, eval_MSE=True)

    vals = y + 2*np.sqrt(mse)

    return vals[0]





def kriging_lower(gp, x):

    y, mse = gp.predict(x, eval_MSE=True)

    vals = y - 2*np.sqrt(mse)

    return vals[0]





def make_fourier(points, ranges, order=5):

    constants = np.ones((points.shape[0], 1))

    temp = np.copy(points)

    for i in range(temp.shape[1]):

        temp[:, i] /= ranges[i]

    basis = np.concatenate((constants, points, np.cos(np.pi * temp), np.sin(np.pi * temp)), axis=1)

    #    basis = np.concatenate((constants, np.cos(np.pi*temp), np.sin(np.pi*temp)), axis=1)

    for i in range(order - 1):

        basis = np.concatenate((basis, np.cos((i + 2) * np.pi * temp), np.sin((i + 2) * np.pi * temp)), axis=1)

    return basis





def eval_fourier(points, ranges, coeff):

    return np.dot(make_fourier(points, ranges), coeff)





def mesh_2_list(x, y):

    temp = np.column_stack([x.flatten(), y.flatten()])

    print(temp.shape)

    return temp