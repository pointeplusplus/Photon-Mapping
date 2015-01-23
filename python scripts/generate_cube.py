
class Vertex:
	def __init__(self, x_, y_, z_):
		self.x = x_
		self.y = y_
		self.z = z_

class Line:
	def __init__(self, v1, v2):
		self.start_point = v1
		self.end_point = v2

#Cube 
file_name = raw_input("Enter file name:  ")
center_x = float(raw_input("Enter center of cube, x coord:  "))
center_y = float(raw_input("Enter center of cube, y coord:  "))
center_z = float(raw_input("Enter center of cube, z coord:  "))
cube_width = float(raw_input("Enter cube width:  "))
starting_vertex_number = int(raw_input("Enter starting vertex number:  "))

draw_top_face = raw_input("Draw top face (y/n):  ")

cube_left = center_x - cube_width/2.0
cube_right = center_x + cube_width/2.0
cube_bottom = center_y - cube_width/2.0
cube_top = center_y + cube_width/2.0
cube_back = center_z - cube_width/2.0
cube_front = center_z + cube_width/2.0

#Label the vertices
right_back_bottom = starting_vertex_number
left_back_bottom = starting_vertex_number + 1
left_front_bottom = starting_vertex_number + 2
right_front_bottom = starting_vertex_number + 3
left_back_top = starting_vertex_number + 4
right_back_top = starting_vertex_number + 5
right_front_top = starting_vertex_number + 6
left_front_top = starting_vertex_number + 7

output_file = open(file_name, "w")

#Print all vertices needed for cube 
#Bottom plane
output_file.write("v " + str(cube_right) + " " + str(cube_bottom) + " " + str(cube_back) + '\n')
output_file.write("v " + str(cube_left) + " " + str(cube_bottom) + " " + str(cube_back) + '\n')
output_file.write("v " + str(cube_left) + " " + str(cube_bottom) + " " + str(cube_front) + '\n')
output_file.write("v " + str(cube_right) + " " + str(cube_bottom) + " " + str(cube_front) + '\n')
#Top plane 
output_file.write("v " + str(cube_left) + " " + str(cube_top) + " " + str(cube_back) + '\n')
output_file.write("v " + str(cube_right) + " " + str(cube_top) + " " + str(cube_back) + '\n')
output_file.write("v " + str(cube_right) + " " + str(cube_top) + " " + str(cube_front) + '\n')
output_file.write("v " + str(cube_left) + " " + str(cube_top) + " " + str(cube_front) + '\n')

#Print the faces
#Bottom face
output_file.write("f " + str(right_back_bottom) + " " + str(left_back_bottom) + " " + str(left_front_bottom) + " "  + str(right_front_bottom) + '\n')
#Top face
if draw_top_face == 'y':
	output_file.write("f " + str(left_back_top) +  " " + str(right_back_top) + " " + str(right_front_top) + " " + str(left_front_top) + '\n')
#Left face
output_file.write("f " + str(left_back_top) + " " + str(left_front_top) + " " + str(left_front_bottom) + " " + str(left_back_bottom) + '\n')
#Right face
output_file.write("f " + str(right_front_top) + " " + str(right_back_top) + " " + str(right_back_bottom) + " " + str(right_front_bottom) + '\n')
#Front face
output_file.write("f " + str(left_front_top) + " " + str(right_front_top) + " " + str(right_front_bottom) + " " + str(left_front_bottom) + '\n')
#Back face
output_file.write("f " + str(right_back_top) + " " + str(left_back_top) + " " + str(left_back_bottom) + " " + str(right_back_bottom) + '\n')

#Generate Light Source
#focus_x = float(raw_input("Enter light focus, x coord:  "))
#focus_y = float(raw_input("Enter light focus, y coord:  "))
#focus_z = float(raw_input("Enter light focus, z coord:  "))
light_size = float(raw_input("Enter light size:  ")) #length of side of square
#focus_point = Vertex(focus_x, focus_y, focus_z)

#Put light source 3/4 of the way between point of interest and the top left back corner (arbitrary choice)
corner_point = Vertex(cube_left, cube_top, cube_back)
'''plane_point = Vertex(corner_point.x * 0.75 + focus_point.x * 0.25,
					 corner_point.y * 0.75 + focus_point.y * 0.25,
					 corner_point.z * 0.75 + focus_point.z * 0.25)'''

output_file.write("\n\nv " + str(cube_left + light_size * 0.5) + " " + str(cube_top) + " " + str(cube_back + light_size * 0.5) + "\n")
output_file.write("v " + str(cube_left + light_size * 1.5) + " " + str(cube_top) + " " + str(cube_back + light_size * 0.5) + "\n")
output_file.write("v " + str(cube_left + light_size * 1.5) + " " + str(cube_top) + " " + str(cube_back + light_size * 1.5) + "\n")
output_file.write("v " + str(cube_left + light_size * 0.5) + " " + str(cube_top) + " " + str(cube_back + light_size * 1.5) + "\n")

output_file.write("f " + str(starting_vertex_number + 8) + " " + str(starting_vertex_number + 9) + " " + str(starting_vertex_number + 10) + " " + str(starting_vertex_number + 11))