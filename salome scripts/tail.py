#import geompy

import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

#sys.path.insert( 0, r'E:/Hull/salome meshes')

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS
import numpy
import os

geompy = geomBuilder.New(theStudy)
batch = batch_input  #if true, use "done" file and try to kill salome after
#points = points_per_turn # num points per turn
points = 4000

amp_list = amp_list_input  #radius of helix centerline i.e. helix amplitude
pipeRadius_list = pipeRadius_list_input  #radius of pipe around helix centerline
wavelength_list = wavelength_list_input   #wavelength of helix
nlambda_list = nlambda_list_input  #number of wavelengths
# kE = k in Shum et al
kE_list = kE_list_input  #parameter controlling exponential decay of helix amplitude to centered starting point

#limit_tol = 0.0025  #needed for boolean fuse to work between helix pipe and end spheres

#wavenumber = k = 2 * pi / wavelength
# t = k * xi
# xi = wavelength * t / 2 / pi
# 1/k = wavelength / 2 / pi
# k = 2 * pi / wavelength
outputname_list = outputname_list_input

#maxsize_list = [0.2 0.2]
#minsize_list = [0.1 0.1]
#fineness_list = [2 2]

maxsize_list = maxsize_list_input
minsize_list = minsize_list_input
fineness_list = fineness_list_input

maxsize_ends_list = maxsize_ends_list_input
minsize_ends_list = minsize_ends_list_input
fineness_ends_list = fineness_ends_list_input

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

if batch and os.path.isfile(done_filename):
	os.remove(done_filename)


for c in range(len(amp_list)):
	print('On iter ', c+1, '   out of ', len(amp_list) )
	amp = amp_list[c]
	pipeRadius = pipeRadius_list[c]
	nlambda = nlambda_list[c]
	wavelength = wavelength_list[c]
	kE = kE_list[c]
	
	outputname = outputname_list[c]
	maxsize = maxsize_list[c]
	minsize = minsize_list[c]
	fineness = fineness_list[c]
	maxsize_ends = maxsize_ends_list[c]
	minsize_ends = minsize_ends_list[c]
	fineness_ends = fineness_ends_list[c]
	
	try:
		def Helix(a, t):
			return [a*(1-math.exp(-(kE)**2 * (wavelength*t/2/math.pi)**2))*math.cos(t), 
			a*(1-math.exp(-(kE)**2 * (wavelength*t/2/math.pi)**2))*math.sin(t), 
			wavelength*t/(2*math.pi)]
		
		def dHelix(a, t):
			return [a*( 1/2*kE**2*wavelength**2*t/(math.pi**2)*math.exp(-kE**2*(wavelength*t/2/math.pi)**2)*math.cos(t) 
			- (1-math.exp(-(kE)**2 * (wavelength*t/2/math.pi)**2))*math.sin(t)), 
			a*( 1/2*kE**2*wavelength**2*t/(math.pi**2)*math.exp(-kE**2*(wavelength*t/2/math.pi)**2)*math.sin(t) 
			+ (1-math.exp(-(kE)**2 * (wavelength*t/2/math.pi)**2))*math.cos(t)), 
			wavelength/(2*math.pi)]
		
		
		# create vertices
		pointList=[]
		samples = numpy.linspace(0, 2*math.pi * nlambda, num = int(round(points * nlambda)))
		for i in range(len(samples)):
			Helixx = Helix(amp, samples[i])
			pointList.append(geompy.MakeVertex(Helixx.pop(0), Helixx.pop(0), Helixx.pop(0)))
		
		HelixCurve = geompy.MakeInterpol(pointList)
		
		Helixx = tuple(Helix(amp, 0))
		dHelixx = dHelix(amp, 0)
		
		P1forVector = geompy.MakeVertex(Helixx[0], Helixx[1], Helixx[2])
		P2forVector = geompy.MakeVertex(Helixx[0]+dHelixx.pop(0), Helixx[1]+dHelixx.pop(0), Helixx[2]+dHelixx.pop(0))
		
		VectorForBaseCircle = geompy.MakeVector(P1forVector, P2forVector)
		
		BaseCircle = geompy.MakeCircle(P1forVector, VectorForBaseCircle, pipeRadius)
		
		CircleFace = geompy.MakeFace(BaseCircle, 1)
		
		Helix_Pipe = geompy.MakePipe(CircleFace, HelixCurve)
		
		Helix_Pipe_face_12 = geompy.GetSubShape(Helix_Pipe, [12])
		Vertex_1 = geompy.MakeVertexOnSurface(Helix_Pipe_face_12, 0.5, 0.5)
		Helix_Pipe_face_3 = geompy.GetSubShape(Helix_Pipe, [3])
		Vertex_2 = geompy.MakeVertexOnSurface(Helix_Pipe_face_3, 0.5, 0.5)
		Sphere_1 = geompy.MakeSpherePntR(Vertex_2, pipeRadius)
		Sphere_2 = geompy.MakeSpherePntR(Vertex_1, pipeRadius)
		
		geompy.addToStudy(Sphere_1, "Sphere_1")
		geompy.addToStudy(Sphere_2, "Sphere_2")
		geompy.addToStudy(Helix_Pipe, "Helix Pipe")
		
		#Limit_tolerance_1 = geompy.LimitTolerance(Helix_Pipe, limit_tol)
		#Fuse_1 = geompy.MakeFuseList([Sphere_1, Sphere_2, Limit_tolerance_1], True, True)
		
		Fuse_1 = geompy.MakeFuseList([Sphere_1, Sphere_2, Helix_Pipe], True, False)
		Rotation_1 = geompy.MakeRotation(Fuse_1, OY, -90*math.pi/180.0)
		geompy.addToStudy(Rotation_1, "Rotation_1")
		[Face_1,Face_2,Face_3] = geompy.ExtractShapes(Rotation_1, geompy.ShapeType["FACE"], True)
		Auto_group_for_SubMesh_1 = geompy.CreateGroup(Rotation_1, geompy.ShapeType["FACE"])
		
		geompy.UnionList(Auto_group_for_SubMesh_1, [Face_1, Face_3])
		
		# add objects in the study
		#id_HelixCurve = geompy.addT oStudy(HelixCurve, "Helix Curve")
		#id_VectorForBase = geompy.addToStudy(VectorForBaseCircle, "Vector For Base Circle")
		#id_BaseCircle = geompy.addToStudy(BaseCircle, "Base Circle")
		#id_CircleFace = geompy.addToStudy(CircleFace, "Base Circle Face")
		# id_HelixPipe = geompy.addToStudy(Helix_Pipe, "Helix Pipe")
		
		# geompy.addToStudy( Rotation_1, 'Rotation_1' )
		
		# geompy.addToStudyInFather( Rotation_1, Face_1, 'Face_1' )
		# geompy.addToStudyInFather( Rotation_1, Face_2, 'Face_2' )
		# geompy.addToStudyInFather( Rotation_1, Face_3, 'Face_3' )
		# geompy.addToStudyInFather( Rotation_1, Auto_group_for_SubMesh_1, 'Auto_group_for_SubMesh_1' )
		
		
		###
		### SMESH component
		###
		
		import  SMESH, SALOMEDS
		from salome.smesh import smeshBuilder
		
		smesh = smeshBuilder.New(theStudy)
		Mesh_1 = smesh.Mesh(Rotation_1)
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
		
		NETGEN_2D_Parameters_sub = smesh.CreateHypothesis('NETGEN_Parameters_2D', 'NETGENEngine')
		NETGEN_2D_Parameters_sub.SetSecondOrder( 1 )
		NETGEN_2D_Parameters_sub.SetOptimize( 1 )
		NETGEN_2D_Parameters_sub.SetUseSurfaceCurvature( 1 )
		NETGEN_2D_Parameters_sub.SetFuseEdges( 1 )
		NETGEN_2D_Parameters_sub.SetQuadAllowed( 0 )
		NETGEN_2D_Parameters_sub.SetFineness( fineness_ends )
		NETGEN_2D_Parameters_sub.SetMaxSize( maxsize_ends )
		NETGEN_2D_Parameters_sub.SetMinSize( minsize_ends )
		
		status = Mesh_1.AddHypothesis(NETGEN_2D,Auto_group_for_SubMesh_1)
		status = Mesh_1.AddHypothesis(NETGEN_2D_Parameters_sub,Auto_group_for_SubMesh_1)
		
		
		isDone = Mesh_1.Compute()
		Mesh_1.ConvertToQuadratic(0)
		Mesh_1.ExportDAT( outputname )
		
		
		
		## Set names of Mesh objects
		smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN_2D')
		smesh.SetName(NETGEN_2D_Parameters, 'NETGEN 2D Parameters')
		smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
	
	except:
		if batch:
			f = open(done_filename, "w")
			f.close()
		
		continue
	
	
print('Done')


if ( not salome.sg.hasDesktop() ) and batch:
	f = open(done_filename, "w")
	f.close()
	
	import killSalome
	killSalome.killAllPorts()
	#sys.exit()
	
else:
	#we update the object browser
	salome.sg.updateObjBrowser(1)
