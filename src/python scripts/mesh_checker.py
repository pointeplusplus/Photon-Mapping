import collections
for line in open("diamond.obj"):
	a = line.split(' ')
	if a[0] == 'f':
		b= [x for x, y in collections.Counter(a).items() if y > 1]
		if len(b) > 0:
			print a,b
	