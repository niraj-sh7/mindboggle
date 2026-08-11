#!/usr/bin/python
"""
Compute the Zernike moments of a collection of points.


Authors:
    - Arthur Mikhno, 2013, Columbia University (original MATLAB code)
    - Brian Rossa, 2013, Tank Think Labs, LLC (port to Python)
    - Arno Klein, 2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def zernike_moments(
    points,
    faces,
    order=10,
    scale_input=True,
    decimate_fraction=0,
    decimate_smooth=0,
    verbose=False,
):
    """
    Compute the Zernike moments of a surface patch of points and faces.

    Optionally decimate the input mesh.

    Note::
      Decimation sometimes leads to an error of "Segmentation fault: 11"
      (Twins-2-1 left label 14 gives such an error only when decimated.)

    Parameters
    ----------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex
    faces : list of lists of 3 integers
        each list contains indices to vertices that form a triangle on a mesh
    order : integer
        order of the moments being calculated
    scale_input : bool
        translate and scale each object so it is bounded by a unit sphere?
        (this is the expected input to zernike_moments())
    decimate_fraction : float
        fraction of mesh faces to remove for decimation (0 for no decimation)
    decimate_smooth : integer
        number of smoothing steps for decimation
    verbose : bool
        print statements?

    Returns
    -------
    descriptors : list of floats
        Zernike descriptors

    Examples
    --------
    >>> pass
    """
    import numpy as np

    from mindboggle.guts.mesh import decimate, reindex_faces_0to1
    from mindboggle.shapes.zernike.pipelines import DefaultPipeline as Pipeline

    # Convert 0-indices (Python) to 1-indices (Matlab) for all face indices:
    index1 = False  # already done elsewhere in the code
    if index1:
        faces = reindex_faces_0to1(faces)

    # Convert lists to numpy arrays:
    if isinstance(points, list):
        points = np.array(points)
    if isinstance(faces, list):
        faces = np.array(faces)

    # ------------------------------------------------------------------------
    # Translate all points so that they are centered at their mean,
    # and scale them so that they are bounded by a unit sphere:
    # ------------------------------------------------------------------------
    if scale_input:
        center = np.mean(points, axis=0)
        points = points - center
        maxd = np.max(np.sqrt(np.sum(points**2, axis=1)))
        points /= maxd

    # ------------------------------------------------------------------------
    # Decimate surface:
    # ------------------------------------------------------------------------
    if 0 < decimate_fraction < 1:
        points, faces, u1, u2 = decimate(
            points, faces, decimate_fraction, decimate_smooth, [], save_vtk=False
        )

        # Convert lists to numpy arrays:
        points = np.array(points)
        faces = np.array(faces)

    # ------------------------------------------------------------------------
    # Multiprocessor pipeline:
    # ------------------------------------------------------------------------
    pl = Pipeline()

    # ------------------------------------------------------------------------
    # Geometric moments:
    # ------------------------------------------------------------------------
    G = pl.geometric_moments_exact(points, faces, order)

    # ------------------------------------------------------------------------
    # ------------------------------------------------------------------------
    Z = pl.zernike(G, order)

    # ------------------------------------------------------------------------
    # Extract Zernike descriptors:
    # ------------------------------------------------------------------------
    descriptors = pl.feature_extraction(Z, order).tolist()

    if verbose:
        print(f"Zernike moments: {descriptors}")

    return descriptors


def zernike_moments_per_label(
    vtk_file,
    order=10,
    exclude_labels=[-1],
    scale_input=True,
    decimate_fraction=0,
    decimate_smooth=25,
    verbose=False,
):
    """
    Compute the Zernike moments per labeled region in a file.

    Optionally decimate the input mesh.

    Parameters
    ----------
    vtk_file : string
        name of VTK surface mesh file containing index scalars (labels)
    order : integer
        number of moments to compute
    exclude_labels : list of integers
        labels to be excluded
    scale_input : bool
        translate and scale each object so it is bounded by a unit sphere?
        (this is the expected input to zernike_moments())
    decimate_fraction : float
        fraction of mesh faces to remove for decimation (1 for no decimation)
    decimate_smooth : integer
        number of smoothing steps for decimation
    verbose : bool
        print statements?

    Returns
    -------
    descriptors_lists : list of lists of floats
        Zernike descriptors per label
    label_list : list of integers
        list of unique labels for which moments are computed

    Examples
    --------
    >>> pass
    """
    import numpy as np

    from mindboggle.guts.mesh import keep_faces
    from mindboggle.mio.vtks import read_vtk
    from mindboggle.shapes.zernike.zernike import zernike_moments

    min_points_faces = 4

    # ------------------------------------------------------------------------
    # Read VTK surface mesh file:
    # ------------------------------------------------------------------------
    points, indices, lines, faces, labels, scalar_names, npoints, input_vtk = read_vtk(
        vtk_file
    )

    # ------------------------------------------------------------------------
    # Loop through labeled regions:
    # ------------------------------------------------------------------------
    ulabels = [x for x in np.unique(labels) if x not in exclude_labels]
    label_list = []
    descriptors_lists = []
    for label in ulabels:
        # if label == 1022:  # 22:
        #    print("DEBUG: COMPUTE FOR ONLY ONE LABEL")

        # --------------------------------------------------------------------
        # Determine the indices per label:
        # --------------------------------------------------------------------
        Ilabel = [i for i, x in enumerate(labels) if x == label]
        if verbose:
            print(f"  {len(Ilabel)} vertices for label {label}")

        if len(Ilabel) > min_points_faces:
            # ----------------------------------------------------------------
            # Remove background faces:
            # ----------------------------------------------------------------
            pick_faces = keep_faces(faces, Ilabel)
            if len(pick_faces) > min_points_faces:
                # ------------------------------------------------------------
                # Compute Zernike moments for the label:
                # ------------------------------------------------------------
                descriptors = zernike_moments(
                    points,
                    pick_faces,
                    order,
                    scale_input,
                    decimate_fraction,
                    decimate_smooth,
                    verbose,
                )

                # ------------------------------------------------------------
                # Append to a list of lists of spectra:
                # ------------------------------------------------------------
                descriptors_lists.append(descriptors)
                label_list.append(int(label))

    return descriptors_lists, label_list


# ============================================================================
# Doctests
# ============================================================================
if __name__ == "__main__":
    import doctest

    doctest.testmod(verbose=True)  # py.test --doctest-modules
