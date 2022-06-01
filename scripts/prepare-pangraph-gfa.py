# Converts a gfa to a visualization of linear blocks
import re
from random import randint
import sys
from random import seed


def readGFA(gfa_file, colour_seed=9999):
    """Reads in a gfa and returns useful format."""
    output_blocks=gfa_file+".blocks.csv"
    block_colour_file=gfa_file+".colours.csv"

    node_dict = {}
    genome_dict = {}
    block_count_dict = {}
    with open(gfa_file, 'r') as gfa:
        i = 1
        for line in gfa.readlines():
            entries = line.split()
            if entries[0] == "S":
                block_count_dict[entries[1]] = 0
                node_dict[entries[1]] = int(re.sub(".*:", "", entries[3]))
                i += 1
            if entries[0] == "P":
                blocks_in_path = entries[2].split(',')
                new_path = []
                i = 0
                for b in blocks_in_path:
                    block_count_dict[b[:-1]] += 1
                    new_path.append([b[:-1], b[-1], i, node_dict[b[:-1]]+i-1])
                    i += node_dict[b[:-1]]
                genome_dict[entries[1].strip('\n')] = new_path
    # Create random colours
    block_colour_dict = {}
    block_num = len([x for x in block_count_dict.values() if x>1])
    block_colors = []
    seed(colour_seed) # for random colours
    for i in range(block_num):
        block_colors.append("#%06X" % randint(0, 0xFFFFFF))
    i = 0
    for block, count in block_count_dict.items():
        if count==1:
            block_colour_dict[block] = "#BDBABA" # colour grey
        else:
            block_colour_dict[block] = block_colors[i]
            i += 1
    with open(output_blocks, 'w') as output_f:
        output_f.write("genome,block,strand,start,end,colour\n")
        for g, values in genome_dict.items():
            for v in values:
                block = v[0]
                output_f.write("%s\n" % (g+","+','.join([str(x) for x in v])+","+block_colour_dict[block]))
    with open(block_colour_file, 'w') as output_f_colours:
        output_f_colours.write("block,colour\n")
        for block, colour in block_colour_dict.items():
            output_f_colours.write("%s,%s\n" % (block, colour))

def rewriteGFA(gfa_file, colour_file, new_gfa_file):
    """Rewrites a GFA with colours."""
    block_colour_dict = {}
    with open(colour_file, 'r') as colours:
        for line in colours.readlines():
            block_colour_dict[line.split(',')[0]] = line.split(',')[1].strip('\n')
    with open(new_gfa_file, "w") as new_gfa:
        with open(gfa_file, 'r') as gfa:
            for line in gfa.readlines():
                entries = line.split()
                if entries[0] == "S":
                    entries.append("CL:z:"+block_colour_dict[entries[1]])
                new_gfa.write("%s\n" % "\t".join(entries))


input_gfa = sys.argv[1]
readGFA(input_gfa)
rewriteGFA(input_gfa,
        input_gfa+".colours.csv",
        input_gfa+".coloured.gfa")

