#import geompy

import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

sys.path.insert( 0, r'E:/Hull/salome meshes')

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS
import numpy
import os

geompy = geomBuilder.New(theStudy)




points = 8000
# num points per turn
#radius = tail_radius
#radius_list = [0.823, 5]
radius_list = radius_list_input


#pipeRadius = tail_pipeRadius
#pipeRadius_list = [0.27483, 0.2]
pipeRadius_list = pipeRadius_list_input

#nturns = num_turns
#nturns_list = [0.744, 0.2]
nturns_list = nturns_list_input
#fraction of complete circular arc to make curved rod

height_list = height_list_input
#list of cylinder heights for straight rod case - can use any placeholder numbers for curved rod entries


outputname_list = outputname_list_input

#maxsize_list = [0.2 0.2]
#minsize_list = [0.1 0.1]
#fineness_list = [2 2]

maxsize_list = maxsize_list_input
minsize_list = minsize_list_input
fineness_list = fineness_list_input

refine_crack_list = refine_crack_list_input
#limit_tol = 0.00025;




O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

def Helix(a, t):
	return [a*math.cos(t), 
	a *math.sin(t), 
	0]

def dHelix(a, t):
	return [-a*math.sin(t), 
	a*math.cos(t), 
	0]

for c in range(len(radius_list)):
	print('On iter ', c+1, '   out of ', len(radius_list) )
	radius = radius_list[c]
	pipeRadius = pipeRadius_list[c]
	nturns = nturns_list[c]
	outputname = outputname_list[c]
	maxsize = maxsize_list[c]
	minsize = minsize_list[c]
	fineness = fineness_list[c]
	height = height_list[c]
	refine_crack = refine_crack_list[c]
	
	if height == 0:  # there is no cylinder, this is the code for a spherical body
		Vertex_1 = geompy.MakeVertex(0, 0, 0)
		Sphere_1 = geompy.MakeSpherePntR(Vertex_1, pipeRadius)
		Body = Sphere_1
	elif math.isinf( radius ):  #infinite radius of curv, must have a straight rod
		Vertex_1 = geompy.MakeVertex(0, 0, 0)
		Vector_1 = geompy.MakeVectorDXDYDZ(1, 0, 0)
		Cylinder_1 = geompy.MakeCylinder(Vertex_1, Vector_1, pipeRadius, height)
		Cylinder_1 = geompy.MakeRotation( Cylinder_1, OX, 45*math.pi/180 )
		Translation_1 = geompy.MakeTranslation(Cylinder_1, -height/2, 0, 0)
		Translation_1_face_12 = geompy.GetSubShape(Translation_1, [12])
		Vertex_2 = geompy.MakeVertexOnSurface(Translation_1_face_12, 0.5, 0.5)
		Translation_1_face_10 = geompy.GetSubShape(Translation_1, [10])
		Vertex_3 = geompy.MakeVertexOnSurface(Translation_1_face_10, 0.5, 0.5)
		Sphere_1 = geompy.MakeSpherePntR(Vertex_3, pipeRadius)
		Sphere_2 = geompy.MakeSpherePntR(Vertex_2, pipeRadius)
		Sphere_1 = geompy.MakeRotation( Sphere_1, OX, 45*math.pi/180 )
		Sphere_2 = geompy.MakeRotation( Sphere_2, OX, 45*math.pi/180 )
		Fuse_1 = geompy.MakeFuseList([Translation_1, Sphere_1, Sphere_2], True)
		Body = Fuse_1
	else:
		# create vertices	
		pointList=[]
		samples = numpy.linspace(0, 2*math.pi * nturns, num = int(round(points * nturns)))
		for i in range(len(samples)):
			Helixx = Helix(radius, samples[i])
			pointList.append(geompy.MakeVertex(Helixx.pop(0), Helixx.pop(0), Helixx.pop(0)))
		
		HelixCurve = geompy.MakeInterpol(pointList)
		
		Helixx = tuple(Helix(radius, 0))
		dHelixx = dHelix(radius, 0)
		
		P1forVector = geompy.MakeVertex(Helixx[0], Helixx[1], Helixx[2])
		P2forVector = geompy.MakeVertex(Helixx[0]+dHelixx.pop(0), Helixx[1]+dHelixx.pop(0), Helixx[2]+dHelixx.pop(0))
		
		VectorForBaseCircle = geompy.MakeVector(P1forVector, P2forVector)
		
		BaseCircle = geompy.MakeCircle(P1forVector, VectorForBaseCircle, pipeRadius)
		
		CircleFace = geompy.MakeFace(BaseCircle, 1)
		
		Helix_Pipe = geompy.MakePipe(CircleFace, HelixCurve)
		
		Helix_Pipe_face_1 = geompy.GetSubShape(Helix_Pipe, [12])
		Vertex_1 = geompy.MakeVertexOnSurface(Helix_Pipe_face_1, 0.5, 0.5)
		Sphere_1 = geompy.MakeSpherePntR(Vertex_1, pipeRadius)
		Plane_1 = geompy.MakePlaneFace(Helix_Pipe_face_1, pipeRadius*4)
		Vector_Normal_1 = geompy.GetNormal(Plane_1, Vertex_1)
		Extrusion_1 = geompy.MakePrismVecH(Plane_1, Vector_Normal_1, -pipeRadius*4)
		hemisph_1 = geompy.MakeCutList(Sphere_1, [Extrusion_1], True)
		
		Helix_Pipe_face_2 = geompy.GetSubShape(Helix_Pipe, [3])
		Vertex_2 = geompy.MakeVertexOnSurface(Helix_Pipe_face_2, 0.5, 0.5)
		Sphere_2 = geompy.MakeSpherePntR(Vertex_2, pipeRadius)
		Plane_2 = geompy.MakePlaneFace(Helix_Pipe_face_2, pipeRadius*4)
		Vector_Normal_2 = geompy.GetNormal(Plane_2, Vertex_2)
		Extrusion_2 = geompy.MakePrismVecH(Plane_2, Vector_Normal_2, -pipeRadius*4)
		hemisph_2 = geompy.MakeCutList(Sphere_2, [Extrusion_2], True)
		
		
		
		geompy.addToStudy(Sphere_1,'Sphere_1')
		geompy.addToStudy(Sphere_2,'Sphere_2')
		geompy.addToStudy(hemisph_1,'hemisph_1')
		geompy.addToStudy(hemisph_2,'hemisph_2')
		geompy.addToStudy(Helix_Pipe,'Helix_Pipe')
		
		#Limit_tolerance_1 = geompy.LimitTolerance(Helix_Pipe, 0.0001)
		#Limit_tolerance_1 = geompy.LimitTolerance(Helix_Pipe, limit_tol)
		#Limit_tolerance_1 = Helix_Pipe;
		#Fuse_1 = geompy.MakeFuseList([Sphere_1, Sphere_2, Helix_Pipe], True)
		try:
			 
			 Fuse_1 = geompy.MakeFuseList([hemisph_1, hemisph_2, Helix_Pipe], True, True)
			 geompy.addToStudy(Fuse_1,'Fuse_1')
		except:
			try:
				hemisph_1_lt = geompy.LimitTolerance(hemisph_1, 1E-3)
				hemisph_2_lt = geompy.LimitTolerance(hemisph_2, 1E-3)
				Helix_Pipe_lt = geompy.LimitTolerance(Helix_Pipe, 1E-3)
				Fuse_1 = geompy.MakeFuseList([hemisph_1_lt, hemisph_2_lt, Helix_Pipe_lt], False, False)
				geompy.addToStudy(Fuse_1,'Fuse_1')
			except:
				try:
					hemisph_1_lt = geompy.LimitTolerance(hemisph_1, 1E-2)
					hemisph_2_lt = geompy.LimitTolerance(hemisph_2, 1E-2)
					Helix_Pipe_lt = geompy.LimitTolerance(Helix_Pipe, 1E-2)
					Fuse_1 = geompy.MakeFuseList([hemisph_1_lt, hemisph_2_lt, Helix_Pipe_lt], False, False)
					geompy.addToStudy(Fuse_1,'Fuse_1')
				except:
					print('Fuse failed, skipping')
					continue
		
		#Rotation_1 = geompy.MakeRotation(Fuse_1, OY, -90*math.pi/180.0)
		#Rotation_1 = geompy.MakeRotation(Rotation_1, OZ, -90*math.pi/180.0)
		Translation_1 = geompy.MakeTranslation(Fuse_1, -radius, 0, 0)
		Rotation_1 = geompy.MakeRotation(Translation_1, OZ, -90*math.pi/180.0)
		Body = Rotation_1
		
		if refine_crack:
			Fuse_1_edge_9 = geompy.GetSubShape(Fuse_1, [9])
			Vertex_1 = geompy.MakeVertexOnCurve(Fuse_1_edge_9, 0.5, True)
			Plane_1 = geompy.MakePlane(Vertex_1, OZ, 10)
			Translation_1 = geompy.MakeTranslation(Plane_1, 0, 0, -pipeRadius * 1/2)   # this is what determines extent of crack refinement (multiple by smaller fraction for wider region!) , used to be * 2/3
			Mirror_1 = geompy.MakeMirrorByAxis(Translation_1, OX)
			Extrusion_1 = geompy.MakePrismVecH(Translation_1, OZ, 10)
			Extrusion_2 = geompy.MakePrismVecH(Mirror_1, OZ, -10)
			
			Intersection_2 = geompy.MakeSection(Fuse_1, Extrusion_1, True)
			Intersection_3 = geompy.MakeSection(Fuse_1, Extrusion_2, True)
			Fuse_1_face_3 = geompy.GetSubShape(Fuse_1, [3])
			plane_pt_1 = geompy.MakeVertexOnSurface(Fuse_1_face_3, 1, 0.5)
			plane_pt_2 = geompy.MakeVertexOnSurface(Fuse_1_face_3, 0.5, 0.5)
			plane_pt_3 = geompy.MakeVertexOnSurface(Fuse_1_face_3, 0.75, 0.5)
			crack_plane = geompy.MakePlaneThreePnt(plane_pt_1, plane_pt_2, plane_pt_3, 5)
			body_crackplane_intersection = geompy.MakeSection(Fuse_1, crack_plane, True)
			Cut_1 = geompy.MakeCutList(body_crackplane_intersection, [Extrusion_1, Extrusion_2], True)
			[Edge_1,Edge_2] = geompy.SubShapeAllSortedCentres(Cut_1, geompy.ShapeType["EDGE"])
			
			centroid_1 = geompy.MakeCDG(Edge_1)
			centroid_2 = geompy.MakeCDG(Edge_2)
			coords_1 = geompy.PointCoordinates(centroid_1)
			coords_2 = geompy.PointCoordinates(centroid_2)
			
			if coords_1[1] < coords_2[1]: # lower y coord value should be correct one
				Crack = Edge_1
			else:
				Crack = Edge_2
			
			Crack = geompy.MakeTranslation(Crack, -radius, 0, 0)
			Crack = geompy.MakeRotation(Crack, OZ, -90*math.pi/180.0)
			geompy.addToStudy(Crack,'Crack')
			BodyCrack = geompy.MakePartition([Body], [Crack], [], [], geompy.ShapeType["SOLID"], 0, [], 1)
			geompy.addToStudy(BodyCrack,'BodyCrack')
			#[Edge_1,Edge_2,Edge_3,Edge_4,Crack_subshape,Edge_6,Edge_7,Edge_8,Edge_9,Edge_10] = geompy.SubShapeAllSortedCentres(BodyCrack, geompy.ShapeType["EDGE"])
			Crack_subshape = geompy.GetInPlace(BodyCrack, Crack)
			geompy.addToStudy(Crack_subshape, 'Crack_subshape')
	
	
	###
	### SMESH component
	###
	
	import  SMESH, SALOMEDS
	from salome.smesh import smeshBuilder
	smesh = smeshBuilder.New(theStudy)
	
	if not refine_crack:
		
		Mesh_1 = smesh.Mesh(Body)
		NETGEN_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D)
		NETGEN_2D_Parameters = NETGEN_2D.Parameters()
		NETGEN_2D_Parameters.SetMaxSize( maxsize )
		NETGEN_2D_Parameters.SetMinSize( minsize )
		NETGEN_2D_Parameters.SetSecondOrder( 1 )
		NETGEN_2D_Parameters.SetOptimize( 1 )
		NETGEN_2D_Parameters.SetFineness( fineness )
		#NETGEN_2D_Parameters.SetGrowthRate( 0.03 )
		#NETGEN_2D_Parameters.SetNbSegPerEdge( 15 )
		#NETGEN_2D_Parameters.SetNbSegPerRadius( 1 )
		
		NETGEN_2D_Parameters.SetUseSurfaceCurvature( 1 )
		NETGEN_2D_Parameters.SetFuseEdges( 1 )
		NETGEN_2D_Parameters.SetQuadAllowed( 0 )
		
		#NETGEN_2D_Parameters_sub = smesh.CreateHypothesis('NETGEN_Parameters_2D', 'NETGENEngine')
		#NETGEN_2D_Parameters_sub.SetSecondOrder( 1 )
		#NETGEN_2D_Parameters_sub.SetOptimize( 1 )
		#NETGEN_2D_Parameters_sub.SetUseSurfaceCurvature( 1 )
		#NETGEN_2D_Parameters_sub.SetFuseEdges( 1 )
		#NETGEN_2D_Parameters_sub.SetQuadAllowed( 0 )
		#NETGEN_2D_Parameters_sub.SetFineness( ends_fineness )
		#NETGEN_2D_Parameters_sub.SetMaxSize( ends_maxsize )
		#NETGEN_2D_Parameters_sub.SetMinSize( ends_minsize )
		
		#status = Mesh_1.AddHypothesis(NETGEN_2D,Auto_group_for_SubMesh_1)
		#status = Mesh_1.AddHypothesis(NETGEN_2D_Parameters_sub,Auto_group_for_SubMesh_1)
		
	else:
		Mesh_1 = smesh.Mesh(BodyCrack)
		NETGEN_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D)
		NETGEN_2D_Parameters = NETGEN_2D.Parameters()
		NETGEN_2D_Parameters.SetMaxSize( maxsize )
		NETGEN_2D_Parameters.SetMinSize( minsize )
		NETGEN_2D_Parameters.SetSecondOrder( 1 )
		NETGEN_2D_Parameters.SetOptimize( 1 )
		NETGEN_2D_Parameters.SetFineness( fineness )
		NETGEN_2D_Parameters.SetUseSurfaceCurvature( 1 )
		NETGEN_2D_Parameters.SetFuseEdges( 1 )
		NETGEN_2D_Parameters.SetQuadAllowed( 0 )
		Regular_1D = Mesh_1.Segment(geom=Crack_subshape)
		
		Max_Size_1 = Regular_1D.MaxSize( minsize )
		Max_Size_1.SetLength( minsize / 2 *3   )
		
		
	try:
		isDone = Mesh_1.Compute()
		Mesh_1.ConvertToQuadratic(0)
		if os.path.isfile(outputname):
			os.remove(outputname)
		Mesh_1.ExportDAT( outputname )
		
		## Set names of Mesh objects
		smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN_2D')
		smesh.SetName(NETGEN_2D_Parameters, 'NETGEN 2D Parameters')
		smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
		
	except:
		print('Meshing failed, continuing')
		
		
		


print('Done')
import killSalome
killSalome.killAllPorts()
sys.exit()
