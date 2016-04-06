# Aggregates 2 results files into 1

import sys

def main(argv):
	if len(argv) < 4:
		print("Usage: {} <input file 1> <input file 2> <output filename>".format(argv[0]));
		exit()
	
	inputfile1 = open(argv[1], 'r');
	inputfile2 = open(argv[2], 'r');
	outputfile = open(argv[3], 'w');

	# 1st line = Photons shot
	photons1 = int((inputfile1.readline().strip().split(' '))[-1]);
	photons2 = int((inputfile2.readline().strip().split(' '))[-1]);
	
	outputfile.write('Photons shot: {}\n'.format(photons1 + photons2));
	
	# 2nd line = input file
	tmp = inputfile1.readline()
	assert tmp == inputfile2.readline();
	outputfile.write(tmp);
	
	# 3rd line = num of bounces
	tmp = inputfile1.readline()
	assert tmp == inputfile2.readline();
	outputfile.write(tmp);
	
	# 4th line = num of shadow samples
	tmp = inputfile1.readline()
	assert tmp == inputfile2.readline();
	outputfile.write(tmp);
	
	# 5th line = column headers
	tmp = inputfile1.readline()
	assert tmp == inputfile2.readline();
	outputfile.write(tmp);
	
	# process faces
	line1 = inputfile1.readline().strip();
	line2 = inputfile2.readline().strip();
	while ":" in line1:
		assert ":" in line2
		
		part1 = line1.split(" ")
		part2 = line2.split(" ")
		
		# there are 8 parts - number, material, normal, area, and 4 ray counts
		# ray counts need to be added, the rest just output as normal
		# verify that number and material is equal - the others may have floating point error
		# so we don't bother checking if they are equal
		assert part1[0] == part2[0];
		assert part1[1] == part2[1];
		
		r1 = int(part1[4]) + int(part2[4])
		r2 = int(part1[5]) + int(part2[5])
		r3 = int(part1[6]) + int(part2[6])
		r4 = int(part1[7]) + int(part2[7])
		
		outputfile.write("{} {} {} {} {} {} {} {}\n".format(part1[0], part1[1], part1[2], part1[3], r1, r2, r3, r4));
		
		# read next lines
		line1 = inputfile1.readline().strip();
		line2 = inputfile2.readline().strip();
		
	print ("aggregation complete");
	
if __name__ == "__main__":
	main(sys.argv);

