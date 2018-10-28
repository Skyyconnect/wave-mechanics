from random import randint


chromaticScale = [440,466.16,493.88,523.25,554.37,587.33,622.25,659.26, 698.46,739.99,783.99,830.61,880]
 


#twelve root 2 = tr2
def octave(CS, multiplier):
	newCS = []
	for i in range(len(CS)):
		newCS.append(multiplier*CS[i])
	return newCS

def tr2(root):
	return (2**(root/12.0))

def tr2_Difference(CS):
	newCS = []
	for i in range(len(CS)):
		newCS.append((CS[i]*tr2(i)) - CS[i])
	return newCS


def tr2_All(CS):
	newCS = []
	for i in range(len(CS)):
		newCS.append(CS[i]*tr2(i))
	return newCS



for i in range(len(chromaticScale)):
	print "roots: " + str(tr2(i))

'''
print "all: " +  str(tr2_All(chromaticScale))
print "minus: "+ str(tr2_Difference(chromaticScale))
print "**********************************************************"
'''
for i in range(len(chromaticScale)):

	print i 
	print "octave : " + str(tr2_All(octave(chromaticScale,i)))
	print "octave-diff : " + str(tr2_Difference(octave(chromaticScale,i)))
