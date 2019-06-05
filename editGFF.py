import sys


def openGFF(filename):
    op = open(filename, "r")
    lines = op.readlines()
    return lines

def editGFF(lines):
    for line in lines:
        line = line.split("\t")
        line[8] = "ID={id}".format(id=line[8])
        line[2] = "super_enhancer"
        print("\t".join(line)[:-1])

        
editGFF(openGFF(sys.argv[1]))
