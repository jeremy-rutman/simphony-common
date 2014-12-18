""" foam_wrapper module

Wrapper module for OpenFOAM using pythonFlu -python interface
NOTE: at a moment pythonFlu supports only up to 2.0.1 -version of OpenFOAM
  
"""
import copy
from Foam import ref,man
from simphony.cuds.mesh import Mesh,Point,Face,Edge,Cell
import inspect
from Foam.OpenFOAM import IOobject,fileName,word
from collections import defaultdict
from Foam.finiteVolume import volScalarField


class Foam_wrapper(object):
    
    def __init__(self,model):
        self.model = model

    

    def has_common_vertex(self,facePoints0,facePoints1):
        for point in facePoints0:
                if point in facePoints1 :
                    return True
        return False
    

    def get_mesh(self, foamCaseName):
        """ Reads given OpenFoam -mesh from "foamCaseName" -directory and returns SimPhony mesh
        """ 
        argv=["test", "-case", foamCaseName]
        args = ref.setRootCase( len(argv), argv )
        runTime = man.createTime( args )
        foamMesh = man.createMesh( runTime )

        foamPoints = foamMesh.points()
        foamFaces = foamMesh.faces()
        foamCells = foamMesh.cells()

        

        simphonyMesh = Mesh()

        foamLabelToUuid=[]
 
        for i in range(foamPoints.size()):
            point = foamPoints[i]

            spoint = Point((point.x(),point.y(),point.z()))
            simphonyMesh.add_point(spoint)
            foamLabelToUuid.append(spoint.uuid)

                 
        print "Points added"

        for i in range(foamFaces.size()):
            face = foamFaces[i]
            points = []
           
            for pi in range(face.size()):
                pl=face[pi]
                
                points.append(foamLabelToUuid[pl])
            face = Face(points)
            simphonyMesh.add_face(face)
 
        print "Faces added"


        for i in range(foamCells.size()):
            faces = foamCells[i]
             
            nnodes = len(foamMesh.cellPoints(i))
           
            if faces.size()==6:
                if nnodes==8:
                    type="hex"
                else:
                    type="wedge"
            elif faces.size()==5:
                if nnodes==5:
                    type="pyr"
                else:
                    type="prism"
            elif faces.size()==4:
                if nnodes==4:
                    type="tet"
                else:
                    type="tetWedge"


            # make mapping from point to faces
            pointFaces =defaultdict(list)
            for fli in range(faces.size()):
                face = foamFaces[faces[fli]]
                for ii in range(face.size()):
                    point = face[ii]
                    if point in pointFaces:
                        if not fli in pointFaces[point]:
                            pointFaces[point].append(fli)
                    else:
                        pointFaces[point].append(fli)


            cellPoints=[]

            for fli in range(faces.size()):
                face = foamFaces[faces[fli]]

                points=[]
                for ii in range(face.size()):
                    points.append(face[ii])
                   
                if cellPoints:

                    if type=="hex" or type=="prism":

                        if not self.has_common_vertex(cellPoints,points):
                            # must be opposite face while has no common vertex points

                            # walk through first face points and find opposite point from opposite face. Opposite point has two common faces with the  point
                            oppositePoints=[]
                            for point in cellPoints:
                                
                                pointfaces = pointFaces[point]
                                for pf,ff in pointFaces.iteritems():
                                    
                                    if (not pf==point) & (pf not in cellPoints):
                                       
                                        commonFaces = set(ff) & set(pointfaces)

                                        if len(commonFaces)==2:
                                            oppositePoints.append(pf)
                                            break
                                    
                            break
                    elif type=="tet":
                       oppositePoints=[]
                       for point in points:
                           if point not in cellPoints:
                               oppositePoints.append(point)
                    else:
                        error_str = "Element type not supported yet: {}"
                        raise ValueError(error_str.format(type))
                else:
                    if type=="prism":
                        if face.size()==3:
                            for point in reversed(points):
                                cellPoints.append(point)
                    else:
                        for point in reversed(points):
                            cellPoints.append(point)

            
            cellPointLabels= []
            
            for point in cellPoints:
                cellPointLabels.append(foamLabelToUuid[point])
            for point in oppositePoints:
                cellPointLabels.append(foamLabelToUuid[point])

                
            simphonyMesh.add_cell(Cell(cellPointLabels))

            
        del foamLabelToUuid

        print "Cells added"

        return simphonyMesh


"""

# moves simphonyMesh to OpenFoam polyMesh structure
# to get this work pythonFlu must be extended


    def to_foam_mesh(self,simphonyMesh):

        foamPoints = []
        for point in simphonyMesh.iter_points():
            
            foamPoints.append(point.coordinates)

        print "Points added"

        foamFaces = []
        for face in simphonyMesh.iter_faces():
            foamFaces.append(face.points)

        print "Faces added"

        foamCells = []
        for cell in simphonyMesh.iter_cells():
            foamCells.append(cell.points)

        print "Cells added"


        print inspect.getmembers(Foam.src.OpenFOAM.meshes.polyMesh)
        foamMesh = Foam.src.OpenFOAM.meshes.polyMesh.polyMesh.new_polyMesh()
        
"""

"""
# Reads given OpenFoam -scalarvariable from "foamCaseName" -directory# and returns corresponding SimPhony variable 
# to get this to work needs extension of pythonFlu to get vertex based values. Tis can be done using OpenFoam interpolations and pointmesh
from Foam.src.OpenFOAM import *
    def get_scalarvariable(self, foamCaseName, variableName):
        
        argv=["test", "-case", foamCaseName]
        args = ref.setRootCase( len(argv), argv )
        runTime = man.createTime( args )
        foamMesh = man.createMesh( runTime )
        variable = volScalarField( IOobject( word( variableName ),
                        fileName( runTime.timeName() ),
                        foamMesh,
                        IOobject.MUST_READ,
                        IOobject.AUTO_WRITE ),
                        foamMesh)

        internalField = variable.internalField()
        
        for value in internalField:
            print value

        
        print "# variables ",len(internalField)
 
"""       