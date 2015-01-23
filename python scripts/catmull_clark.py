
class Vertex:
	def __init__(self, x_, y_, z_):
		self.x = x_
		self.y = y_
		self.z = z_

class Edge:
	def __init__(self, v_1, v_2, midpoint):
		self.v1 = v_1
		self.v2 = v_2
		self.middle_vertex = midpoint

class Triangle:
	def __init__(self, a, b, c):
		self.v1 = a
		self.v2 = b
		self.v3 = c 

class Quad:
	def __init__(self, a, b, c, d):
		self.v1 = a
		self.v2 = b
		self.v3 = c
		self.v4 = d

def parse_vertex(line_info, vertices):
	new_vert = Vertex(float(line_info[1]), float(line_info[2]), float(line_info[3]))
	vertices.append(new_vert)

def parse_triangle(line_info, tringles):
	new_triangle = Triangle(int(line_info[1]), int(line_info[2]), int(line_info[3]))
	triangles.append(new_triangle)

def parse_quad(line_info, quads):
	new_quad = Quad(int(line_info[1]), int(line_info[2]), int(line_info[3]), int(line_info[4]))
	quads.append(new_quad)

def print_file(vertices, quads):
	for v in vertices:
		print "v", v.x, v.y, v.z

	for q in quads:
		print "f", q.v1, q.v2, q.v3, q.v4

def find_edge_midpoint_vertex(edges, vertex1, vertex2):

	for edge in edges:
		if (vertex1 == edge.v1 and vertex2 == edge.v2) or (vertex1 == edge.v2 and vertex2 == edge.v1):
			return edge.middle_vertex
	return None

def add_midpoint(vertices, edges, vert1, vert2):
	
	midpoint_x = (vertices[vert1-1].x + vertices[vert2-1].x) / 2.0
	midpoint_y = (vertices[vert1-1].y + vertices[vert2-1].y) / 2.0
	midpoint_z = (vertices[vert1-1].z + vertices[vert2-1].z) / 2.0

	#print "Midpoint calculation"
	#print "  old x vals:", vertices[vert1].x, vertices[vert2].x

	#add the vertex to the vertex list
	new_vertex = Vertex(midpoint_x, midpoint_y, midpoint_z)
	vertices.append(new_vertex)

	#add the new edge midpoint to the edges list
	new_edge = Edge(vert1, vert2, vertices.index(new_vertex) + 1)
	edges.append(new_edge)

	return vertices.index(new_vertex) + 1

def add_centerpoint(vertices, vert1, vert2, vert3):
	
	middle_x = (vertices[vert1-1].x + vertices[vert2-1].x + vertices[vert3-1].x) / 3.0
	middle_y = (vertices[vert1-1].y + vertices[vert2-1].y + vertices[vert3-1].y) / 3.0
	middle_z = (vertices[vert1-1].z + vertices[vert2-1].z + vertices[vert3-1].z) / 3.0

	new_vertex = Vertex(middle_x, middle_y, middle_z)
	vertices.append(new_vertex)

	return vertices.index(new_vertex) + 1

def quadranglzie(quads, vertices, edges, triangle):
	
	#find all alredy existing midpoints
	midpoint1 = find_edge_midpoint_vertex(edges, triangle.v1, triangle.v2)
	midpoint2 = find_edge_midpoint_vertex(edges, triangle.v2, triangle.v3)
	midpoint3 = find_edge_midpoint_vertex(edges, triangle.v3, triangle.v1)
	
	#add midpoints that don't already exist
	if midpoint1 == None:
		midpoint1 = add_midpoint(vertices, edges, triangle.v1, triangle.v2)

	if midpoint2 == None:
		midpoint2 = add_midpoint(vertices, edges, triangle.v2, triangle.v3)

	if midpoint3 == None:
		midpoint3 = add_midpoint(vertices, edges, triangle.v3, triangle.v1)	
	

	#find center point for triangle
	center_point = add_centerpoint(vertices, triangle.v1, triangle.v2, triangle.v3)

	#add the new quads
	new_quad = Quad(triangle.v1, midpoint1, center_point, midpoint3)
	quads.append(new_quad)

	new_quad = Quad(triangle.v2, midpoint2, center_point, midpoint1)
	quads.append(new_quad)

	new_quad = Quad(triangle.v3, midpoint3, center_point, midpoint2)
	quads.append(new_quad)


#all vertices so far
vertices = []
#contains all edges that already have midpoints
edges = []
#all triangls -> need to be turned into quads
triangles = []
#all original and new quads
quads = []

file_name = 'diamond_3_converter.obj'
file_stream = open(file_name)
file_contents = file_stream.read()
lines = file_contents.split('\n')

#parse the file
for line in lines:
	line_info = line.split(' ')

	#different things to do based on the type of line
	if line_info[0] == 'v':
		parse_vertex(line_info, vertices)
	elif line_info[0] == 'f':
		if len(line_info) == 4:
			parse_triangle(line_info, triangles)
		elif len(line_info) == 5:
			parse_quad(line_info, quads)
		else:
			print "can't handle face with", len(line_info) -1, "vertices"

#turn each triangle into quads
for triangle in triangles:
	quadranglzie(quads, vertices, edges, triangle)

print_file(vertices, quads)