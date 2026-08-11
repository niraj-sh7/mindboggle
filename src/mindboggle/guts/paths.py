#!/usr/bin/env python
"""
Curves, tracks, skeletons connecting surface mesh vertices.

Authors:
    - Arno Klein, 2012-2016  (arno@mindboggle.info)  http://binarybottle.com
    - Yrjo Hame, 2012  (yrjo.hame@gmail.com)

Copyright 2016,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def connect_points_erosion(
    S,
    neighbor_lists,
    outer_anchors,
    inner_anchors=[],
    values=[],
    erode_ratio=0.1,
    erode_min_size=10,
    save_steps=[],
    save_vtk="",
    background_value=-1,
    verbose=False,
):
    """
    >>> pass
"""
    import numpy as np

    from mindboggle.guts.mesh import extract_edge, find_endpoints, topo_test
    from mindboggle.guts.segment import segment_regions

    # Make sure arguments are numpy arrays:
    if not isinstance(S, np.ndarray):
        S = np.array(S)
    erode_by_value = False
    if len(values) and 0 <= erode_ratio < 1:
        erode_by_value = True
        if not isinstance(values, np.ndarray):
            values = np.array(values)

    keep = outer_anchors[:]
    keep.extend(inner_anchors)
    remove_endpoints = True

    if save_steps:
        from mindboggle.mio.vtks import rewrite_scalars

        S0 = S.copy()

    # ------------------------------------------------------------------------
    # Iteratively remove simple points:
    # ------------------------------------------------------------------------
    if verbose:
        print(f"  Remove up to {erode_ratio} of edge vertices per iteration")
    complex = []
    count = -1
    exist_simple = True
    while exist_simple:
        exist_simple = False
        if verbose or save_steps:
            count += 1

        # --------------------------------------------------------------------
        # Only consider updating vertices that are on the edge of the
        # region and are not among the indices to keep or known simple points:
        # --------------------------------------------------------------------
        indices = np.where(background_value != S)[0].tolist()
        edge = extract_edge(indices, neighbor_lists)
        if edge:
            edge = np.array(list(set(edge).difference(complex)))
            len_edge = np.shape(edge)[0]
            if len_edge:
                # ------------------------------------------------------------
                # Segment edge vertices into separate connected groups:
                # ------------------------------------------------------------
                edge_segs = segment_regions(
                    edge,
                    neighbor_lists,
                    1,
                    [],
                    False,
                    False,
                    [],
                    [],
                    [],
                    "",
                    background_value,
                    verbose,
                )

                edge_seg_numbers = [
                    x for x in np.unique(edge_segs) if x != background_value
                ]
                if verbose:
                    len_numbers = len(edge_seg_numbers)
                    if len_numbers > 1:
                        print(
                            f"    {count}: {len_edge} edge points in {len_numbers} segments"
                        )
                    else:
                        print(f"    {count}: {len_edge} edge points")
                first_seg = True
                for edge_seg_number in edge_seg_numbers:
                    edge_seg = np.where(edge_segs == edge_seg_number)[0]
                    edge_seg = np.array(list(set(edge_seg).difference(keep)))
                    len_edge_seg = np.shape(edge_seg)[0]
                    if len_edge_seg:
                        # ----------------------------------------------------
                        # Remove topologically simple points
                        # in order of lowest to highest values:
                        # ----------------------------------------------------
                        ntests = len_edge_seg
                        if erode_by_value and ntests > erode_min_size:
                            Isort = np.argsort(values[edge_seg])
                            edge_seg = edge_seg[Isort]
                            if erode_ratio > 0:
                                ntests = int(len_edge_seg * erode_ratio) + 1

                        for index in edge_seg[0:ntests]:
                            # Test to see if each index is a simple point:
                            simple, d = topo_test(index, S, neighbor_lists)

                            # If a simple point, remove and run again:
                            # (Note: Must remove at each iteration)
                            if simple:
                                S[index] = background_value
                                exist_simple = True
                            # Else store to exclude in future:
                            else:
                                complex.append(index)

                        # If no simple points, test all of the indices:
                        if not exist_simple and erode_by_value:
                            if verbose:
                                print("    No simple points")
                            for index in edge_seg[ntests::]:
                                simple, d = topo_test(index, S, neighbor_lists)
                                # If a simple point, remove and run again:
                                if simple:
                                    S[index] = background_value
                                    exist_simple = True
                                # Else store to exclude in future:
                                else:
                                    complex.append(index)

                        # Save incremental VTK files for debugging:
                        if count in save_steps and first_seg:
                            IDs = background_value * np.ones(len(values))
                            IDs[indices] = values[indices]
                            rewrite_scalars(
                                save_vtk,
                                "edge" + str(count) + ".vtk",
                                IDs,
                                "edges",
                                S0,
                                background_value,
                            )
                        first_seg = False

                # ------------------------------------------------------------
                # Remove branches by iteratively removing endpoints:
                # ------------------------------------------------------------
                if remove_endpoints:
                    indices = np.where(background_value != S)[0].tolist()
                    endpts = True
                    while endpts:
                        endpts = find_endpoints(indices, neighbor_lists)
                        if endpts:
                            endpts = [x for x in endpts if x not in outer_anchors]
                            if endpts:
                                S[endpts] = background_value
                                indices = list(set(indices).difference(endpts))

    skeleton = indices

    return skeleton


def connect_points_hmmf(
    indices_points,
    indices,
    L,
    neighbor_lists,
    wN_max=1.0,
    do_erode=True,
    background_value=-1,
    verbose=False,
):
    """
    Connect mesh vertices with a skeleton of 1-vertex-thick curves using HMMF.

    The goal of this algorithm is to assign each vertex a locally optimal
    Hidden Markov Measure Field (HMMF) value and to connect vertices according
    to a cost function that penalizes vertices that do not have high likelihood
    values and have HMMF values different than their neighbors.

    We initialize the HMMF values with likelihood values normalized to the
    interval (0.5, 1.0] (to guarantee correct topology) and take those values
    that are greater than the likelihood threshold (1 for each anchor point).

    We iteratively update each HMMF value if it is near the likelihood
    threshold such that a H_step makes it cross the threshold,
    and the vertex is a "simple point" (its addition/removal alters topology).

    Parameters for computing the cost and cost gradients:

        ``wL``: weight influence of likelihood on the cost function
        ``wN``: weight influence of neighbors on the cost function
        ``H_step``: the amount that the HMMF values are incremented

    Parameters
    ----------
    indices_points : list of integers
        indices of vertices to connect (should contain > 1)
    indices : list of integers
        indices of vertices through which to connect points
    L : numpy array of floats
        likelihood values for all vertices in mesh
    neighbor_lists : list of lists of integers
        indices to neighboring vertices for each vertex in mesh
    wN_max : float
        maximum neighborhood weight (trust prior more for smoother fundi)
    do_erode : bool
        erode to create skeleton?
    background_value : integer
        background value
    verbose : bool
        print statements?

    Returns
    -------
    skeleton : list of integers
        indices to vertices connecting the points

    Examples
    --------
    >>> # Connect vertices according to (usually likelihood) values in a fold:
    >>> import numpy as np
    >>> from mindboggle.guts.paths import connect_points_hmmf
    >>> from mindboggle.guts.paths import find_outer_endpoints
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.guts.compute import median_abs_dev
    >>> from mindboggle.guts.mesh import find_neighbors_from_file
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> url1 = urls['left_mean_curvature']
    >>> url2 = urls['left_travel_depth']
    >>> url3 = urls['left_folds']
    >>> curv_file = fetch_data(url1, '', '.vtk')
    >>> depth_file = fetch_data(url2, '', '.vtk')
    >>> folds_file = fetch_data(url3, '', '.vtk')
    >>> curvs, name = read_scalars(curv_file, True, True)
    >>> depths, name = read_scalars(depth_file, True, True)
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> L = curvs * depths
    >>> print(np.array_str(L[0:5], precision=5, suppress_small=True))
    [-0.11778 -0.35642 -0.80759 -0.25654 -0.04411]
    >>> neighbor_lists = find_neighbors_from_file(curv_file)
    >>> background_value = -1
    >>> # Limit number of folds to speed up the test:
    >>> limit_folds = True
    >>> if limit_folds:
    ...     fold_numbers = [4] #[4, 6]
    ...     indices = [i for i,x in enumerate(folds) if x in fold_numbers]
    ...     i0 = [i for i,x in enumerate(folds) if x not in fold_numbers]
    ...     folds[i0] = background_value
    ... else:
    ...     indices = range(len(L))
    >>> # Outer anchors:
    >>> min_separation = 10
    >>> verbose = False
    >>> indices_points, tracks = find_outer_endpoints(indices, neighbor_lists,
    ...                             L, depths, min_separation,
    ...                             background_value, verbose)
    >>> wN_max = 2.0
    >>> do_erode = True
    >>> skeleton = connect_points_hmmf(indices_points, indices, L,
    ...     neighbor_lists, wN_max, do_erode, background_value, verbose)
    >>> skeleton[0:10]
    [50324, 51535, 52692, 52698, 52699, 52706, 52707, 52708, 52717, 52725]

    Write out vtk file and view (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> from mindboggle.mio.vtks import rewrite_scalars # doctest: +SKIP
    >>> folds_copy = np.copy(folds) # doctest: +SKIP
    >>> folds[skeleton] = 100 # doctest: +SKIP
    >>> folds[indices_points] = 120 # doctest: +SKIP
    >>> rewrite_scalars(depth_file, 'connect_points_hmmf.vtk',
    ...                 folds, 'skeleton', folds, -1) # doctest: +SKIP
    >>> plot_surfaces('connect_points_hmmf.vtk') # doctest: +SKIP

    """
    import numpy as np

    from mindboggle.guts.mesh import topo_test
    from mindboggle.guts.paths import connect_points_erosion

    # Make sure argument is a numpy array
    if not isinstance(L, np.ndarray):
        L = np.array(L)

    # ------------------------------------------------------------------------
    # Parameters:
    # ------------------------------------------------------------------------
    # Cost and cost gradient parameters:
    wN_min = 0.0  # minimum neighborhood weight
    # wN_max = 2.0  # maximum neighborhood weight (trust prior more for smoother fundi)
    H_step = 0.1  # step down HMMF value

    # Parameters to speed up optimization and for termination of the algorithm:
    grad_min = 0.1  # minimum gradient factor
    grad_max = 1.0  # maximum gradient factor
    slope_exp = 2
    rate_factor = 0.9
    min_cost_change = 0.0001  # minimum change in the sum of costs
    n_tries_no_change = 3  # number of loops without sufficient change
    min_count = 50  # min. iterations (to overcome initial increasing costs)
    max_count = 300  # maximum number of iterations (in case no convergence)

    # Miscellaneous parameters:
    print_interval = 10

    def compute_costs(
        likelihoods, hmmfs, hmmfs_neighbors, numbers_of_neighbors, wN, Z=[]
    ):
        """
        Cost function for penalizing unlikely fundus curve vertices.

        This cost function penalizes vertices with low fundus likelihood values,
        and whose Hidden Markov Measure Field (HMMF) values differ from
        their neighbors:

        cost = hmmf * (1.1 - likelihood) +
               wN * sum(abs(hmmf - hmmf_neighbors)) / size_neighbors

        term 1 promotes high likelihood values
        term 2 promotes smoothness of the HMMF values

        Note: 1.1 is used instead of 1 to ensure that there is a cost
              for all points, even for those with likelihoods close to 1.

        Parameters
        ----------
        likelihoods : numpy array of floats
            likelihood values in interval [0,1]
        hmmf : numpy array of floats
            HMMF values
        hmmf_neighbors : numpy array of floats
            HMMF values of neighboring vertices for each vertex
        numbers_of_neighbors : numpy array of integers
            number of neighbors for each vertex
        wN : float
            weight influence of neighbors on cost (term 2)
        Z : numpy array of integers in {0,1} (same shape as hmmf_neighbors)
            to remove zero-padded neighborhood elements after hmmf subtraction

        Returns
        -------
        costs : numpy array of floats
            cost values

        """
        import numpy as np

        if all(numbers_of_neighbors):
            # Subtract each HMMF value from its neighbors:
            diff = abs(hmmfs - hmmfs_neighbors)

            # Since the above neighbors may have been padded,
            # remove the extra differences by multiplying them times zero:
            diff = diff * Z

            # Compute the cost for each vertex:
            costs = (
                hmmfs * (1.1 - likelihoods)
                + wN * np.sum(diff, axis=0) / numbers_of_neighbors
            )
        else:
            raise OSError("No HMMF neighbors to compute cost.")

        return costs

    # ------------------------------------------------------------------------
    # Initialize all Hidden Markov Measure Field (HMMF) values with
    # likelihood values (except 0) normalized to the interval (0.5, 1.0]
    # (to guarantee correct topology). Assign a 1 for each anchor point.
    # This influences surrounding vertex neighborhoods.
    # Note: 0.5 is the class boundary threshold for the HMMF values.
    H = np.zeros(len(L))
    H_new = (L + 1.000001) / 2
    H_new[L == 0.0] = 0
    H_new[H_new > 1.0] = 1
    H[H_new > 0.5] = H_new[H_new > 0.5]
    H[indices_points] = 1
    H_new = H.copy()
    H_tests = H.copy()

    # Find the HMMF values for the neighbors of each vertex:
    N = neighbor_lists
    N_sizes = np.array([len(x) for x in N])
    max_num_neighbors = max(N_sizes[indices])
    N_array = np.zeros((max_num_neighbors, len(L)))
    for index in indices:
        N_array[0 : N_sizes[index], index] = N[index]
    N_array_shape = np.shape(N_array)
    N_flat = np.ravel(N_array)
    N_flat_list = N_flat.tolist()
    H_N = np.reshape(H[[int(x) for x in N_flat_list]], N_array_shape)
    ind_flat = [i for i, x in enumerate(N_flat_list) if x > 0]
    len_flat = len(N_flat_list)

    # A zero in N calls H[0], so remove zero-padded neighborhood elements:
    Z = np.zeros((max_num_neighbors, len(L)))
    for index in indices:
        Z[0 : N_sizes[index], index] = 1

    # Assign cost values to each vertex (for indices):
    C = np.zeros(len(L))
    C[indices] = compute_costs(
        L[indices], H[indices], H_N[:, indices], N_sizes[indices], wN_max, Z[:, indices]
    )
    npoints = len(indices)

    # Loop until count reaches max_count or until end_flag equals zero
    # (end_flag allows the loop to continue a few times even if no change):
    count = 0
    end_flag = 0
    wN = wN_max
    gradient_factor = grad_min
    while end_flag < n_tries_no_change and count < max_count:
        # Select indices with a positive HMMF value:
        V = [indices[i] for i, x in enumerate(H[indices]) if x > 0.0]

        # Update neighborhood H values:
        # H_N = np.reshape(H[N_flat_list], N_array_shape)
        H_N = np.zeros(len_flat)
        H_N[ind_flat] = H[[int(x) for x in N_flat[ind_flat]]]
        H_N = np.reshape(H_N, N_array_shape)

        # Compute the cost gradient for the HMMF values:
        H_decr = H - H_step
        H_decr[H_decr < 0] = 0.0
        C_decr = compute_costs(L[V], H_decr[V], H_N[:, V], N_sizes[V], wN, Z[:, V])
        H_tests[V] = H[V] - gradient_factor * (C[V] - C_decr)
        H_tests[H_tests < 0] = 0.0
        H_tests[H_tests > 1] = 1.0

        # For each index:
        for index in V:
            # Do not update anchor point costs:
            if index not in indices_points:
                # Update a vertex HMMF value if it is away from the threshold:
                update = True
                # Or if it crosses the threshold and is a topologically
                # "simple point" (0.5 not considered part of the fundus):
                if H[index] > 0.5 >= H_tests[index]:
                    update, n_in = topo_test(index, H_new, N)
                elif H[index] <= 0.5 < H_tests[index]:
                    update, n_in = topo_test(index, 1 - H_new, N)
                if update:
                    H_new[index] = H_tests[index]

        # Update the cost values:
        C[V] = compute_costs(L[V], H_new[V], H_N[:, V], N_sizes[V], wN, Z[:, V])

        # Sum the cost values across all vertices and tally the number
        # of HMMF values greater than the threshold.
        # After iteration 1, compare current and previous values.
        # If the values are similar, increment end_flag:
        costs = sum(C[V].tolist())
        npoints_thr = len([x for x in H[V].tolist() if x > 0.5])

        # Terminate the loop if there are insufficient changes:
        if count > 0:
            delta_cost = (costs_previous - costs) / npoints
            delta_points = npoints_thr_previous - npoints_thr
            if delta_points == 0:
                if delta_cost < min_cost_change and count > min_count:
                    end_flag += 1
            else:
                end_flag = 0

            # Display information every n_mod iterations:
            if verbose and not np.mod(count, print_interval):
                print(
                    f"      Iteration {count}: {delta_points} crossing threshold "
                    f"(wN={wN:0.3f}, grad={gradient_factor:0.3f}, cost={delta_cost:0.3f})"
                )

            # Increment the gradient factor and decrement the neighborhood factor
            # so that spacing is close in early iterations and far apart later:
            factor = (count / np.round(rate_factor * max_count)) ** slope_exp
            if gradient_factor < grad_max:
                gradient_factor = factor * (grad_max - grad_min) + grad_min
            if wN > wN_min:
                wN = wN_max - factor * (wN_max - wN_min)

        # Reset for next iteration:
        costs_previous = costs
        npoints_thr_previous = npoints_thr
        H = H_new

        count += 1

    if verbose:
        print("      Updated hidden Markov measure field (HMMF) values")

    # Skeletonize:
    if do_erode:
        # Threshold the resulting array:
        S = background_value * np.ones(len(L))
        S[indices] = H[indices]
        S[S > 0.5] = 1.0
        S[S <= 0.5] = background_value

        skeleton = connect_points_erosion(
            S,
            neighbor_lists=N,
            outer_anchors=indices_points,
            inner_anchors=[],
            values=H,
            erode_ratio=0.5,
            erode_min_size=10,
            save_steps=[],
            save_vtk="",
            background_value=background_value,
            verbose=verbose,
        )
        if verbose:
            npoints_thr = len([x for x in S if x != background_value])
            print(
                f"      Removed {int(npoints_thr - len(skeleton))} points to create one-vertex-thin "
                "skeletons"
            )
    else:
        # Threshold the resulting array:
        S = np.zeros(len(L))
        S[indices] = H[indices]
        S[S > 0.5] = 1.0
        S[S <= 0.5] = 0.0
        skeleton = [i for i, x in enumerate(S.tolist()) if x == 1]
        if verbose:
            print(
                f"      Removed {int(sum(S.tolist()) - len(skeleton))} points to create one-vertex-thin "
                "skeletons"
            )

    return skeleton


def smooth_skeletons(
    skeletons,
    bounds,
    vtk_file,
    likelihoods,
    wN_max=1.0,
    do_erode=True,
    save_file=False,
    output_file="",
    background_value=-1,
    verbose=False,
):
    """
    Smooth skeleton by dilation followed by connect_points_hmmf().

    Steps ::
        1. Segment skeleton into separate sets of connected vertices.
        2. For each skeleton segment, extract endpoints.
        3. Dilate skeleton segment.
        4. Connect endpoints through dilated segment by connect_points_hmmf().
        5. Store smoothed output from #4.

    Parameters
    ----------
    skeletons : list of integers
        skeleton number for each vertex
    bounds : list of integers
        region number for each vertex; constrains smoothed skeletons
    vtk_file : string
        file from which to extract neighboring vertices for each vertex
    likelihoods : list of integers
        fundus likelihood value for each vertex
    wN_max : float
        maximum neighborhood weight (trust prior more for smoother skeletons)
    do_erode : bool
        erode skeleton?
    save_file : bool
        save output VTK file?
    output_file : string
        output VTK file
    background_value : integer or float
        background value
    verbose : bool
        print statements?

    Returns
    -------
    skeletons : list of integers
        skeleton numbers for all vertices
    n_skeletons :  integer
        number of skeletons
    skeletons_file : string (if save_file)
        name of output VTK file with skeleton numbers

    Examples
    --------
    >>> # Smooth skeleton to extract fundus from one or more folds:
    >>> import numpy as np
    >>> from mindboggle.mio.plots import plot_surfaces
    >>> from mindboggle.guts.paths import smooth_skeletons
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.guts.compute import median_abs_dev
    >>> from mindboggle.guts.mesh import find_neighbors_from_file
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> curv_file = fetch_data(urls['left_mean_curvature'], '', '.vtk')
    >>> depth_file = fetch_data(urls['left_travel_depth'], '', '.vtk')
    >>> folds_file = fetch_data(urls['left_folds'], '', '.vtk')
    >>> fundus_file = fetch_data(urls['left_fundus_per_fold'], '', '.vtk')
    >>> curvs, name = read_scalars(curv_file, True, True)
    >>> depths, name = read_scalars(depth_file, True, True)
    >>> vtk_file = curv_file
    >>> likelihoods = depths * curvs
    >>> [np.float("{0:.{1}f}".format(x, 5)) for x in likelihoods[0:5]]
    [-0.11778, -0.35642, -0.80759, -0.25654, -0.04411]
    >>> bounds, name = read_scalars(folds_file, True, True)
    >>> skeletons, name = read_scalars(fundus_file, True, True)
    >>> background_value = -1
    >>> # Limit number of folds to speed up the test:
    >>> limit_folds = True
    >>> if limit_folds:
    ...     fold_numbers = [7] #[4, 6]
    ...     i0 = [i for i,x in enumerate(bounds) if x not in fold_numbers]
    ...     bounds[i0] = background_value
    ...     skeletons[i0] = background_value
    >>> wN_max = 1.0
    >>> do_erode = True
    >>> save_file = True
    >>> output_file = 'smooth_skeletons.vtk'
    >>> verbose = False
    >>> smoothed_skeletons, n_skeletons, skel_file = smooth_skeletons(skeletons,
    ...     bounds, vtk_file, likelihoods, wN_max, do_erode, save_file,
    ...     output_file, background_value, verbose)
    >>> np.where(np.array(smoothed_skeletons)!=-1)[0][0:8]
    array([112572, 113453, 113454, 113469, 114312, 114313, 114325, 115087])

    Write out vtk file and view (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> from mindboggle.mio.vtks import rewrite_scalars # doctest: +SKIP
    >>> iskels = [i for i,x in enumerate(smoothed_skeletons)
    ...           if x != background_value] # doctest: +SKIP
    >>> bounds[iskels] = 100 # doctest: +SKIP
    >>> rewrite_scalars(depth_file, 'smooth_skeletons_no_background.vtk',
    ...                 bounds, 'skeleton', bounds, -1) # doctest: +SKIP
    >>> plot_surfaces('smooth_skeletons.vtk') # doctest: +SKIP

    """

    import os
    from time import time

    import numpy as np

    from mindboggle.guts.mesh import dilate, find_endpoints, find_neighbors_from_file
    from mindboggle.guts.paths import connect_points_hmmf
    from mindboggle.guts.segment import segment_regions
    from mindboggle.mio.vtks import rewrite_scalars

    t0 = time()

    neighbor_lists = find_neighbors_from_file(vtk_file)
    indices = np.where(bounds != background_value)[0]
    npoints = len(bounds)

    # ------------------------------------------------------------------------
    # Loop through skeletons:
    # ------------------------------------------------------------------------
    unique_IDs = [x for x in np.unique(skeletons) if x != background_value]
    n_skeletons = len(unique_IDs)
    if n_skeletons == 1:
        sdum = ""
    else:
        sdum = "s"
    if verbose:
        print(f"Smooth {n_skeletons} skeleton{sdum}...")
    Z = background_value * np.ones(npoints)
    smoothed_skeletons = Z.copy()
    for ID in unique_IDs:
        skeleton = [i for i, x in enumerate(skeletons) if x == ID]
        if verbose:
            print(f"  Skeleton {int(ID)}:")

        # --------------------------------------------------------------------
        # Segment skeleton vertices into separate connected groups:
        # --------------------------------------------------------------------
        skel_segs = segment_regions(
            skeleton,
            neighbor_lists,
            1,
            [],
            False,
            False,
            [],
            [],
            [],
            "",
            background_value,
            verbose,
        )

        skel_seg_numbers = [x for x in np.unique(skel_segs) if x != background_value]
        len_numbers = len(skel_seg_numbers)
        if verbose and len_numbers > 1:
            print(f"    {len_numbers} segments")
        for skel_seg_number in skel_seg_numbers:
            skel_seg = np.where(skel_segs == skel_seg_number)[0].tolist()

            # ----------------------------------------------------------------
            # Find endpoints:
            # ----------------------------------------------------------------
            endpoints = find_endpoints(skel_seg, neighbor_lists)

            # ----------------------------------------------------------------
            # Dilate the skeleton within the bounds:
            # ----------------------------------------------------------------
            nedges = 2
            if verbose:
                print("    Dilate skeleton within bounds...")
            dilated = dilate(skel_seg, nedges, neighbor_lists)
            dilated = list(set(dilated).intersection(indices))
            if dilated:
                # ------------------------------------------------------------
                # Set undilated likelihoods to background to keep neighbors:
                # ------------------------------------------------------------
                L = Z.copy()
                L[dilated] = likelihoods[dilated]

                # ------------------------------------------------------------
                # Smoothly re-skeletonize the dilated skeleton:
                # ------------------------------------------------------------
                if verbose:
                    print("    Smoothly re-skeletonize dilated skeleton...")
                new_skeleton = connect_points_hmmf(
                    endpoints,
                    dilated,
                    L,
                    neighbor_lists,
                    wN_max,
                    do_erode,
                    background_value,
                    verbose,
                )

                ## Plot overlap of dilated and pre-/post-smoothed skeleton:
                # from mindboggle.mio.plots import plot_surfaces
                # D = background_value * np.ones(npoints)
                # D[dilated]=1; D[skel_seg]=2; D[new_skeleton]=3
                # rewrite_scalars(vtk_file, 'test.vtk', D, 'D', bounds)
                # plot_surfaces('test.vtk')

                # ------------------------------------------------------------
                # Store skeleton:
                # ------------------------------------------------------------
                smoothed_skeletons[new_skeleton] = ID
    if verbose:
        print(f"  ...Smoothed {n_skeletons} skeleton{sdum} ({time() - t0:.2f} seconds)")

    # ------------------------------------------------------------------------
    # Return skeletons, number of skeletons, and file name:
    # ------------------------------------------------------------------------
    smoothed_skeletons = smoothed_skeletons.tolist()

    if save_file:
        if output_file:
            skeletons_file = output_file
        else:
            skeletons_file = os.path.join(os.getcwd(), "smooth_skeletons.vtk")
            rewrite_scalars(
                vtk_file,
                skeletons_file,
                smoothed_skeletons,
                "smoothed_skeletons",
                [],
                background_value,
            )
    else:
        skeletons_file = None

    return smoothed_skeletons, n_skeletons, skeletons_file


def track_segments(seed, segments, neighbor_lists, values, sink, background_value=-1):
    """
    >>> pass
"""
    import numpy as np

    if not sink:
        import sys

        sys.exit("Missing sink vertices.")

    track = [seed]
    for isegment, segment in enumerate(segments):
        # Find the seed's neighborhood N in the segment:
        N = neighbor_lists[seed]
        N = [x for x in N if values[x] != background_value]
        segment_set = set(segment)
        N_segment = list(segment_set.intersection(frozenset(N)))
        if N:
            # Add the neighborhood vertex with the maximum value to the track:
            if N_segment:
                seed = N_segment[np.argmax(values[N_segment])]
                track.append(seed)

                # If the track has run into the region's border, return the track:
                # if list(set(N_segment).intersection(sink)):
                if seed in sink:
                    return track

            # If there is no neighbor in the new segment,
            # back up to the previous segment:
            elif isegment > 0:
                bridge = []
                max_bridge = 0
                N_previous = list(set(N).intersection(segments[isegment - 1]))
                for Np in N_previous:
                    N_next = list(segment_set.intersection(neighbor_lists[Np]))
                    if N_next:
                        if np.max(values[N_next]) > max_bridge:
                            seed = N_next[np.argmax(values[N_next])]
                            bridge = [Np, seed]
                            max_bridge = np.max(values[N_next])
                if bridge:
                    track.extend(bridge)

                    # If the track has run into the region's border, return the track:
                    if seed in sink:
                        return track

        # If there is no neighborhood for the seed, return None:
        else:
            return None

    # If the track remains empty or does not reach the border, return the track:
    return None


def find_outer_endpoints(
    indices,
    neighbor_lists,
    values,
    values_seeding,
    min_separation=10,
    background_value=-1,
    verbose=False,
):
    """
    >>> pass
"""
    import numpy as np

    from mindboggle.guts.mesh import find_neighborhood
    from mindboggle.guts.paths import track_segments
    from mindboggle.guts.segment import extract_borders, segment_rings

    # ------------------------------------------------------------------------
    # Settings:
    # ------------------------------------------------------------------------
    # For finding outward tracks:
    do_threshold = True
    remove_fraction = 0.5  # Remove fraction of the surface (see below)
    do_filter_tracks = True

    # Initialize R, T, S, V:
    R = indices[:]
    T = []
    S = np.array(values_seeding)
    V = np.array(values)

    # ------------------------------------------------------------------------
    # Extract region boundary:
    # ------------------------------------------------------------------------
    B = np.ones(len(V))
    B[indices] = 2
    borders, foo1, foo2 = extract_borders(list(range(len(B))), B, neighbor_lists)

    # ------------------------------------------------------------------------
    # Initialize seeds with vertices at the median-depth boundary:
    # ------------------------------------------------------------------------
    if do_threshold:
        thresholdS = np.median(S[indices])  # + np.std(S[indices])
        indices_high = [x for x in indices if S[x] >= thresholdS]
        # Make sure threshold is within the maximum values of the boundary:
        if list(frozenset(indices_high).intersection(borders)):
            do_threshold = False
        else:
            if verbose:
                print(f"  Initialize seeds at {thresholdS:.2f} (median of fold depth)")
    # ------------------------------------------------------------------------
    # Or initialize seeds with vertices at the shrunken region boundary:
    # ------------------------------------------------------------------------
    if not do_threshold:
        thresholdS = remove_fraction * np.max(S[indices])
        if verbose:
            print(f"  Initialize seeds at {1 - remove_fraction:.2f} of fold depth")

    # Extract threshold boundary vertices as seeds:
    indices_high = [x for x in indices if S[x] >= thresholdS]
    B = np.ones(len(S))
    B[indices_high] = 2
    seeds, foo1, foo2 = extract_borders(list(range(len(S))), B, neighbor_lists)

    # ------------------------------------------------------------------------
    # Segment the mesh from the seeds iteratively toward the boundary:
    # ------------------------------------------------------------------------
    R = list(frozenset(R).difference(indices_high))
    R = list(frozenset(R).difference(seeds))
    segments = segment_rings(R, seeds, neighbor_lists, 1, background_value)

    # Run tracks from the seeds through the segments toward the boundary:
    if verbose:
        print(
            f"    Track through {len(segments)} concentric segments ({len(R)} vertices) "
            f"from threshold {thresholdS:0.2f}"
        )
    for seed in seeds:
        track = track_segments(
            seed, segments, neighbor_lists, V, borders, background_value
        )
        if track:
            T.append(track)

    # ------------------------------------------------------------------------
    # Filter the tracks in two ways:
    # 1. Keep tracks that have a high median track value.
    # 2. Filter tracks so that there are no two track endpoint vertices
    # within a given number of edges between each other.  We select the track
    # with higher median value when their endpoints are close.
    # ------------------------------------------------------------------------
    if do_filter_tracks and T:
        if verbose:
            print(f"    Filter {len(T)} tracks")

        # Compute median track values:
        Tvalues = [np.median(V[x]) for x in T]

        # Keep tracks with a high median track value:
        # background = np.median(V[indices])
        # background = np.median(Tvalues) + np.std(Tvalues)
        background = np.median(V[R]) + np.std(V[R])
        Ihigh = [i for i, x in enumerate(T) if Tvalues[i] > background]
        T = [T[i] for i in Ihigh]
        Tvalues = [Tvalues[i] for i in Ihigh]

        # Gather endpoint vertex indices:
        E = [x[-1] for x in T]

        # Loop through endpoints:
        E2 = []
        T2 = []
        while E:
            # Find endpoints close to the first endpoint:
            near = find_neighborhood(neighbor_lists, [E[0]], min_separation)
            Isame = [i for i, x in enumerate(E) if x == E[0]]
            Inear = [i for i, x in enumerate(E) if x in near]
            if Inear or len(Isame) > 1:
                Inear.extend(Isame)

                # Select endpoint with the maximum median track value:
                Imax = Inear[np.argmax([Tvalues[i] for i in Inear])]
                E2.append(E[Imax])
                T2.append(T[Imax])

                # Keep data for points that are not nearby for the next loop:
                Ikeep = [i for i in range(len(E)) if i not in Inear]
                E = [E[i] for i in Ikeep]
                T = [T[i] for i in Ikeep]
                Tvalues = [Tvalues[i] for i in Ikeep]

            # If nothing nearby, simply keep the endpoint:
            else:
                E2.append(E[0])
                T2.append(T[0])
                del E[0]
                del T[0]
                del Tvalues[0]

        endpoints = E2
        endtracks = T2

        if len(T2) == 1:
            s = "track"
        else:
            s = "tracks"
        if verbose:
            print(f"    Retain {len(T2)} {s}")

    else:
        # Gather endpoint vertex indices:
        endpoints = [x[-1] for x in T]
        endtracks = T

    return endpoints, endtracks


def find_max_values(points, values, min_separation=10, thr=0.5):
    """
    >>> pass
"""
    from operator import itemgetter

    import numpy as np

    # Make sure arguments are numpy arrays:
    if not isinstance(points, np.ndarray):
        points = np.array(points)

    # Sort values and find indices for values above the threshold:
    L_table = [[i, x] for i, x in enumerate(values)]
    L_table_sort = np.transpose(sorted(L_table, key=itemgetter(1)))[:, ::-1]
    IL = [int(L_table_sort[0, i]) for i, x in enumerate(L_table_sort[1, :]) if x > thr]

    # Initialize special points list with the index of the maximum value,
    # remove this value, and loop through the remaining high values:
    if IL:
        highest = [IL.pop(0)]
        for imax in IL:
            # Determine if there are any special points
            # near to the current maximum value vertex:
            i = 0
            found = 0
            while i < len(highest) and found == 0:
                # Compute Euclidean distance between points:
                D = np.linalg.norm(points[highest[i]] - points[imax])

                # If distance less than threshold, consider the point found:
                if min_separation > D:
                    found = 1
                ## Compute directional distance between points if close:
                # elif D < max_distance:
                #    dirV = np.dot(points[indices_special[i]] - points[imax],
                #               min_directions[highest[i]])
                #    # If distance less than threshold, consider point found:
                #    if np.linalg.norm(dirV) < min_separation:
                #        found = 1

                i += 1

            # If there are no nearby special points,
            # assign the maximum value vertex as a special point:
            if not found:
                highest.append(imax)
    else:
        highest = []

    return highest


# def track_values(seed, indices, neighbor_lists, values, sink=[]):
#     """
#     Build a track from a seed vertex along increasing vertex values of a mesh.
#
#     This function builds a track from an initial seed vertex through
#     increasing value vertices of a surface mesh, optionally terminating
#     at any of a set of sink vertices.
#
#     NOTE :: Untested!
#
#     Parameters
#     ----------
#     seed : integer
#         index to initial seed vertex from which to grow a track
#     indices : list of integers
#         indices of vertices through which to connect points
#     neighbor_lists : list of lists of integers
#         indices to neighboring vertices for each vertex
#     values : numpy array of floats
#         values for all vertices that help to guide a track
#     sink : list of integers
#         indices for vertices that end a track (optional)
#
#     Returns
#     -------
#     track : list of integers
#         indices of ordered vertices for a single track
#
#     Examples
#     --------
#     >>> # Track from deepest point in a fold to its boundary:
#     >>> import numpy as np
#     >>> from mindboggle.mio.vtks import read_scalars
#     >>> from mindboggle.guts.mesh import find_neighbors_from_file
#     >>> from mindboggle.guts.paths import track_values
#     >>> from mindboggle.mio.fetch_data import prep_tests
#     >>> urls, fetch_data = prep_tests()
#     >>> url1 = urls['left_travel_depth']
#     >>> url2 = urls['left_folds']
#     >>> depth_file = fetch_data(url1, '', '.vtk')
#     >>> folds_file = fetch_data(url2, '', '.vtk')
#     >>> depths, name = read_scalars(depth_file, True, True)
#     >>> folds, name = read_scalars(folds_file, True, True)
#     >>> # Limit number of folds to speed up the test:
#     >>> limit_folds = True
#     >>> if limit_folds:
#     ...     background_value = -1
#     ...     fold_numbers = [4] #[4, 6]
#     ...     indices = [i for i,x in enumerate(folds) if x in fold_numbers]
#     ...     i0 = [i for i,x in enumerate(folds) if x not in fold_numbers]
#     ...     folds[i0] = background_value
#     >>> neighbor_lists = find_neighbors_from_file(depth_file)
#     >>> seed = indices[np.argmax(depths[indices])]
#     >>> seed
#     65804
#     >>> values = depths
#     >>> [float("{0:.{1}f}".format(x, 5)) values[0:5]]
#     [0.02026, 0.06009, 0.12859, 0.04564, 0.00774]
#     >>> sink = []
#     >>> track = track_values(seed, indices, neighbor_lists, values, sink)
#     >>> track[0:10]
#     [65804, 64495, 65805, 67118, 64494, 64495, 67098, 65789]
#
#     View track in fold on surface (skip test):
#
#     >>> from mindboggle.mio.vtks import rewrite_scalars
#     >>> from mindboggle.mio.plots import plot_surfaces
#     >>> folds_copy = np.copy(folds)
#     >>> folds_copy[track] = 10
#     >>> folds_copy[seed] = 15
#     >>> rewrite_scalars(depth_file, 'track.vtk', folds_copy, 'track', folds)
#     >>> plot_surfaces('track.vtk')
#
#     """
#     import numpy as np
#
#     track = [seed]
#     value = values[seed]
#     prev_neighbors = track
#     prev_value = 0.0
#     while value >= prev_value:
#         prev_value = value
#
#         # Find the seed's neighborhood N within indices:
#         N = neighbor_lists[seed]
#         N = list(frozenset(N).intersection(indices))
#         N = [x for x in N if x not in prev_neighbors]
#         if N:
#
#             # Add the neighborhood vertex with maximum value to the track:
#             seed = N[np.argmax(values[N])]
#             track.append(seed)
#             value = values[seed]
#             prev_neighbors += N
#
#             # If the track has run into the sink vertices, return track:
#             #if list(frozenset(N_segment).intersection(sink)):
#             if sink and seed in sink:
#                 return track
#
#         # # If there is no neighborhood for the seed, widen the neighborhood:
#         # else:
#         #     # Find the neighborhood N's neighborhood N2 within indices:
#         #     N2 = []
#         #     for n in N:
#         #         N2.append(neighbor_lists[n])
#         #     N2 = list(frozenset(N2).intersection(indices))
#         #     N2 = [x for x in N2 if x not in past_seeds]
#         #     if N2:
#         #
#         #         # Add the neighborhood vertex with maximum value to the track:
#         #         seed = N2[np.argmax(values[N2])]
#         #         track.append(seed)
#         #         value = values[seed]
#         #         prev_neighbors += N2
#         #
#         #         # If the track has run into the sink vertices, return track:
#         #         if sink and seed in sink:
#         #             return track
#         #
#         #     # Widen the neighborhood once more:
#         #     else:
#         #         # Find the neighborhood N's neighborhood N2 within indices:
#         #         N3 = []
#         #         for n in N2:
#         #             N3.append(neighbor_lists[n])
#         #         N3 = list(frozenset(N3).intersection(indices))
#         #         N3 = [x for x in N3 if x not in past_seeds]
#         #         if N3:
#         #
#         #             # Add the neighborhood vertex with maximum value to the track:
#         #             seed = N3[np.argmax(values[N3])]
#         #             track.append(seed)
#         #             value = values[seed]
#         #             prev_neighbors += N3
#         #
#         #             # If the track has run into the sink vertices, return track:
#         #             if sink and seed in sink:
#         #                 return track
#         #
#         #         # If there is no neighborhood for the seed, return the track:
#         #         elif len(track) > 1:
#         #             return track
#         #         else:
#         #             return None
#
#         # If there is no neighborhood for the seed, return the track:
#         elif len(track) > 1:
#             return track
#         else:
#             return None
#
#     # Return the track:
#     if len(track) > 1:
#         return track
#     else:
#         return None


# ============================================================================
# Doctests
# ============================================================================
if __name__ == "__main__":
    import doctest

    doctest.testmod(verbose=True)  # py.test --doctest-modules
