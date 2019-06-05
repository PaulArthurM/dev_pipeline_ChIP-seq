# This script was a experiment. No longer used for hockey stick plots generation.

import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import collections

print("Création d'un hockey stick plot à partir d'un gff et count htseq.")

def openFile(filename):
    print("Opening file: {}".format(filename))
    op = open(filename, "r")
    lines = op.readlines()
    return lines

def createDict(gff, *args):
    print('Creating dictionnary')
    gff = openFile(gff)
    def gffDict(gff):
        dict = collections.OrderedDict()
        for line in gff:
            line = line.split("\t")
            rank = int(line[5])
            ID = line[8][3:-2]
            dict[ID] = rank
        return dict

    def htseq_countDict(arg):
        print('Get counts')
        arg = openFile(arg)
        counts = {}
        for line in arg:
            line = line.split("\t")
            ID = line[0]
            count = line[1][:-1]
            counts[ID] = int(count)
        return counts

    def orderCounts(gff, counts):
        print('Put counts in order')
        ranking = [0] * len(gff)
        for count_id in counts:
            for gff_id in gff:
                if count_id == gff_id:
                    ranking[gff[gff_id]-1] = counts[count_id]
        return ranking

    def plot_counts(ranking):
        with PdfPages('multipage_pdf.pdf') as pdf:
            print('Creating Hockey stick plot')
            plt.plot(ranking)
            plt.ylabel('some numbers')
            plt.title('Hockey stick plot')
            pdf.savefig()
            plt.show()

    def write_ranking(ranking, args):
        print('Creating ranking file')
        args = "_".join([arg.split("/")[1].split(".")[0] for arg in args])
        filename = "htseq_count/" + args + "_ranking.txt"
        file = open(filename, 'w')
        for index, value in enumerate(ranking):
            value = str(index+1) + "\t" + str(value) + "\n"
            file.write(value)
        file.close()

    ranking = orderCounts(gffDict(gff), htseq_countDict(args[0]))[::-1]
    #plot_counts(ranking)
    write_ranking(ranking, args)


createDict(sys.argv[1], sys.argv[2])
