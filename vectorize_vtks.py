from pathlib import Path

path = Path("src/mindboggle/mio/vtks.py")
text = path.read_text()

# Fix 1: read_lines
old1 = """    Data = Reader.GetOutput()
    Lns = Data.GetLines()

    lines = [
        [Lns.GetData().GetValue(j) for j in range(i * 3 + 1, i * 3 + 3)]
        for i in range(Data.GetNumberOfLines())
    ]"""
new1 = """    from vtk.util.numpy_support import vtk_to_numpy

    Data = Reader.GetOutput()
    Lns = Data.GetLines()

    if Data.GetNumberOfLines() > 0:
        lines_array = vtk_to_numpy(Lns.GetData())
        lines = lines_array.reshape(-1, 3)[:, 1:3].tolist()
    else:
        lines = []"""
assert text.count(old1) == 1, f"read_lines block match count: {text.count(old1)}"
text = text.replace(old1, new1)

# Fix 2: read_faces_points
old2 = """    Data = Reader.GetOutput()
    points = [
        list(Data.GetPoint(point_id)) for point_id in range(Data.GetNumberOfPoints())
    ]
    npoints = len(points)

    if Data.GetNumberOfPolys() > 0:
        faces = [
            [
                int(Data.GetPolys().GetData().GetValue(j))
                for j in range(i * 4 + 1, i * 4 + 4)
            ]
            for i in range(Data.GetPolys().GetNumberOfCells())
        ]
    else:
        faces = []

    return faces, points, npoints"""
new2 = """    from vtk.util.numpy_support import vtk_to_numpy

    Data = Reader.GetOutput()
    points = vtk_to_numpy(Data.GetPoints().GetData()).tolist()
    npoints = len(points)

    if Data.GetNumberOfPolys() > 0:
        polys_array = vtk_to_numpy(Data.GetPolys().GetData())
        faces = polys_array.reshape(-1, 4)[:, 1:4].tolist()
    else:
        faces = []

    return faces, points, npoints"""
assert text.count(old2) == 1, f"read_faces_points block match count: {text.count(old2)}"
text = text.replace(old2, new2)

# Fix 3: read_vtk main function (points + faces + lines)
old3 = """    Data = Reader.GetOutput()
    PointData = Data.GetPointData()
    points = [
        list(Data.GetPoint(point_id)) for point_id in range(0, Data.GetNumberOfPoints())
    ]
    npoints = len(points)

    if Data.GetNumberOfPolys() > 0:
        faces = [
            [
                int(Data.GetPolys().GetData().GetValue(j))
                for j in range(i * 4 + 1, i * 4 + 4)
            ]
            for i in range(Data.GetPolys().GetNumberOfCells())
        ]
    else:
        faces = []

    if Data.GetNumberOfLines() > 0:
        lines = [
            [Data.GetLines().GetData().GetValue(j) for j in range(i * 3 + 1, i * 3 + 3)]
            for i in range(Data.GetNumberOfLines())
        ]
    else:
        lines = []"""
new3 = """    from vtk.util.numpy_support import vtk_to_numpy

    Data = Reader.GetOutput()
    PointData = Data.GetPointData()
    points = vtk_to_numpy(Data.GetPoints().GetData()).tolist()
    npoints = len(points)

    if Data.GetNumberOfPolys() > 0:
        polys_array = vtk_to_numpy(Data.GetPolys().GetData())
        faces = polys_array.reshape(-1, 4)[:, 1:4].tolist()
    else:
        faces = []

    if Data.GetNumberOfLines() > 0:
        lines_array = vtk_to_numpy(Data.GetLines().GetData())
        lines = lines_array.reshape(-1, 3)[:, 1:3].tolist()
    else:
        lines = []"""
assert text.count(old3) == 1, f"read_vtk block match count: {text.count(old3)}"
text = text.replace(old3, new3)

path.write_text(text)
print("All 3 blocks replaced successfully.")
