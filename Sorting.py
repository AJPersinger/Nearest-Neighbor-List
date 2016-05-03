#Sorting Algorithms
someList = [2082014598, 7045196965, 3900351543, 8152141763, 6837573588, 3957789636, 8783736452, 7320829674]
tempArray=[]
final = []

for x in xrange(0, len(someList)):
    try:
        tempArray[x] = x
        print tempArray[x]
    except IndexError:
        tempArray[x] = -1
        print tempArray[x]

print tempArray
