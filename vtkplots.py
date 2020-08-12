import vtk

#Draw polygon grid
def plotPolyGrid(polys, savename):
    g = open("{}_polygrid_all.visit".format(savename), "w")
    g.write("!NBLOCKS {}\n".format(sum(list(map(lambda x : len(x), polys)))))
    print(len(polys))
    for j in range(len(polys)):
        poly = polys[j]
        for i in range(len(poly)):
            p1 = poly[i]
            p2 = poly[(i+1) % len(poly)]
            line = vtk.vtkLineSource()
            line.SetPoint1(p1[0], p1[1], 0)
            line.SetPoint2(p2[0], p2[1], 0)
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName("{}_polygrid_{}_{}.vtp".format(savename, j, i))
            writer.SetInputConnection(line.GetOutputPort())
            writer.Write()
            g.write("{}_polygrid_{}_{}.vtp\n".format(savename, j, i))
    g.close()

#Draw circular arcs
#centers = list of centers, intersects = list of pairs of intersects (x, y)
def plotArcFacets(centers, intersects, savename):
    assert len(centers) == len(intersects)
    g = open("{}_all.visit".format(savename), "w")
    g.write("!NBLOCKS {}\n".format(len(centers)))
    for i in range(len(centers)):
        arc = vtk.vtkArcSource()
        arc.SetPoint1(intersects[i][0][0], intersects[i][0][1], 0)
        arc.SetPoint2(intersects[i][1][0], intersects[i][1][1], 0)
        arc.SetCenter(centers[i][0], centers[i][1], 0)
        arc.SetResolution(8)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName("{}_arc_{}.vtp".format(savename, i))
        writer.SetInputConnection(arc.GetOutputPort())
        writer.Write()
        g.write("{}_arc_{}.vtp\n".format(savename, i))
    g.close()

#Draw linear facets
#intersects = list of pairs of intersects (x, y)
def plotLinearFacets(intersects, savename):
    g = open("{}_all.visit".format(savename), "w")
    g.write("!NBLOCKS {}\n".format(sum(list(map(lambda x : len(x)-1, intersects)))))
    for facetnum, facet in enumerate(intersects):
        for i in range(len(facet)-1):
            p1 = facet[i]
            p2 = facet[i+1]
            line = vtk.vtkLineSource()
            line.SetPoint1(p1[0], p1[1], 0)
            line.SetPoint2(p2[0], p2[1], 0)
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName("{}_line_{}_{}.vtp".format(savename, facetnum, i))
            writer.SetInputConnection(line.GetOutputPort())
            writer.Write()
            g.write("{}_line_{}_{}.vtp\n".format(savename, facetnum, i))
    g.close()

#Draw perturbed quad grid
#points = 2d-array, where points[x][y] = the point corresponding to the perturbation of (x, y) in Cartesian grid
def plotQuadGrid(points, savename):
    pointindices = [[-1] * len(points) for _ in range(len(points))]
    indexcounter = 0
    vtkpoints = vtk.vtkPoints()
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.Allocate(len(points)**2)

    for x in range(len(points)-1):
        for y in range(len(points)-1):
            neighbors = [[0, 0], [1, 0], [1, 1], [0, 1]]
            for neighbor in neighbors:
                if pointindices[x+neighbor[0]][y+neighbor[1]] == -1:
                    pointindices[x+neighbor[0]][y+neighbor[1]] = indexcounter
                    vtkpoints.InsertPoint(indexcounter, [points[x+neighbor[0]][y+neighbor[1]][0], points[x+neighbor[0]][y+neighbor[1]][1], 0])
                    indexcounter += 1
            ugrid.InsertNextCell(vtk.VTK_QUAD, 4, [pointindices[x][y], pointindices[x+1][y], pointindices[x+1][y+1], pointindices[x][y+1]])

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName("{}_grid.vtp".format(savename))
    writer.SetInputData(ugrid)
    writer.Write()