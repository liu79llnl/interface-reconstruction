import vtk
import os

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

#Draw perturbed quad grid
#points = 2d-array, where points[x][y] = the point corresponding to the perturbation of (x, y) in Cartesian grid
def plotQuadGrid(points, savename):
    try:
        os.mkdir("advection_vtk/{}".format(savename))
    except:
        pass

    g = open("advection_vtk/{}_grid.visit".format(savename), "w")
    facetnum = 0
    for x in range(len(points)-1):
        for y in range(len(points[0])-1):
            line = vtk.vtkLineSource()
            line.SetPoint1(points[x][y][0], points[x][y][1], 0)
            line.SetPoint2(points[x][y+1][0], points[x][y+1][1], 0)
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
            writer.SetInputConnection(line.GetOutputPort())
            writer.Write()
            facetnum += 1
            line = vtk.vtkLineSource()
            line.SetPoint1(points[x][y][0], points[x][y][1], 0)
            line.SetPoint2(points[x+1][y][0], points[x+1][y][1], 0)
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
            writer.SetInputConnection(line.GetOutputPort())
            writer.Write()
            facetnum += 1
    for x in range(len(points)-1):
        line = vtk.vtkLineSource()
        line.SetPoint1(points[x][len(points[0])-1][0], points[x][len(points[0])-1][1], 0)
        line.SetPoint2(points[x+1][len(points[0])-1][0], points[x+1][len(points[0])-1][1], 0)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
        writer.SetInputConnection(line.GetOutputPort())
        writer.Write()
        facetnum += 1
    for y in range(len(points[0])-1):
        line = vtk.vtkLineSource()
        line.SetPoint1(points[len(points)-1][y][0], points[len(points)-1][y][1], 0)
        line.SetPoint2(points[len(points)-1][y+1][0], points[len(points)-1][y+1][1], 0)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
        writer.SetInputConnection(line.GetOutputPort())
        writer.Write()
        facetnum += 1

    g.write("!NBLOCKS {}\n".format(facetnum))
    for i in range(facetnum):
        g.write("{}/facet_{}.vtp\n".format(savename, i))
    g.close()

#Draw facets, input in format of makeFacets
#Linear facet: ['linear', intersects]
#Corner facet: ['corner', intersects]
#Arc facet: ['arc', arccenter, arcradius, arcintersects]
def plotFacets(facets, savename):
    try:
        os.mkdir("advection_vtk/{}".format(savename))
    except:
        pass

    g = open("advection_vtk/{}_all.visit".format(savename), "w")

    facetnum = 0
    for facet in facets:
        if facet is not None:
            if facet[0] == 'linear' or facet[0] == 'corner':
                for i in range(len(facet[1])-1):
                    p1 = facet[1][i]
                    p2 = facet[1][i+1]
                    line = vtk.vtkLineSource()
                    line.SetPoint1(p1[0], p1[1], 0)
                    line.SetPoint2(p2[0], p2[1], 0)
                    writer = vtk.vtkXMLPolyDataWriter()
                    writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
                    writer.SetInputConnection(line.GetOutputPort())
                    writer.Write()
                    facetnum += 1
            elif facet[0] == 'arc':
                arc = vtk.vtkArcSource()
                arc.SetPoint1(facet[3][0][0], facet[3][0][1], 0)
                arc.SetPoint2(facet[3][-1][0], facet[3][-1][1], 0)
                arc.SetCenter(facet[1][0], facet[1][1], 0)
                arc.SetResolution(8)
                writer = vtk.vtkXMLPolyDataWriter()
                writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
                writer.SetInputConnection(arc.GetOutputPort())
                writer.Write()
                facetnum += 1
            elif facet[0] ==  'curvedcorner':
                if facet[1] is not None:
                    arc = vtk.vtkArcSource()
                    arc.SetPoint1(facet[-1][0][0], facet[-1][0][1], 0)
                    arc.SetPoint2(facet[-1][1][0], facet[-1][1][1], 0)
                    arc.SetCenter(facet[1][0], facet[1][1], 0)
                    arc.SetResolution(8)
                    writer = vtk.vtkXMLPolyDataWriter()
                    writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
                    writer.SetInputConnection(arc.GetOutputPort())
                else:
                    line = vtk.vtkLineSource()
                    line.SetPoint1(facet[-1][0][0], facet[-1][0][1], 0)
                    line.SetPoint2(facet[-1][1][0], facet[-1][1][1], 0)
                    writer = vtk.vtkXMLPolyDataWriter()
                    writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
                    writer.SetInputConnection(line.GetOutputPort())
                writer.Write()
                facetnum += 1
                if facet[2] is not None:
                    arc = vtk.vtkArcSource()
                    arc.SetPoint1(facet[-1][1][0], facet[-1][1][1], 0)
                    arc.SetPoint2(facet[-1][2][0], facet[-1][2][1], 0)
                    arc.SetCenter(facet[2][0], facet[2][1], 0)
                    arc.SetResolution(8)
                    writer = vtk.vtkXMLPolyDataWriter()
                    writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
                    writer.SetInputConnection(arc.GetOutputPort())
                else:
                    line = vtk.vtkLineSource()
                    line.SetPoint1(facet[-1][1][0], facet[-1][1][1], 0)
                    line.SetPoint2(facet[-1][2][0], facet[-1][2][1], 0)
                    writer = vtk.vtkXMLPolyDataWriter()
                    writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
                    writer.SetInputConnection(line.GetOutputPort())
                writer.Write()
                facetnum += 1

    g.write("!NBLOCKS {}\n".format(facetnum))
    for i in range(facetnum):
        g.write("{}/facet_{}.vtp\n".format(savename, i))
    g.close()