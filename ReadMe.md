# An Implementation of ESU (FanMod) Using OpenMP in C++

 This code decomposes a network (format edgelist -- see lfr_network.el as an example) onto 3-or 4-node graphlets, using the ESU algoritm (cf http://theinf1.informatik.uni-jena.de/~wernicke/motifs-wabi2005.pdf).

 * It builds all the Benson Matrices (ie a matrix W where W(i,j) is number of motif occurrences containing nodes and j) of each graphlet that appears in the network at least once (see https://arxiv.org/pdf/1802.06820.pdf Section 2.2.4 - Motif Adjacency Matrix).
 These matrices are saved in Results/BensonGraphs. motif_k_isomorphism_id.txt contains the edgelist of the Benson Matrix of the k-node graphlet identified by the integer id (see https://github.com/luleg/DiscriminantMotifs/blob/master/SupplementaryASONAM2020.pdf Section III for the identifiers)

 * The file motifs_occurrences.txt in the folder Results provides the number of occurrences of each graphlet.

 /!\ The nodes were relabelled during the run and the Benson Matrices are those of the relabelled network. The correspondance between initial node labels (column "label_node") and final node labels (column "id_node") is provided in the file Results/dict_graph.txt.

/!\ Names and locations of directories where the results are saved can be simply changed in Global.h.

Using the code : requires a c++ compiler (c++0x at least) and openMP to be installed/enabled.

-> in a Unix console, type "make"

-> then "./exec_esu -h" specifies how to use the code.
