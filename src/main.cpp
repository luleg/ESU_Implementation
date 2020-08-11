#include "Esu.h"


using namespace std;
using namespace std::chrono;

char *infile = NULL;

string global_graph_file="";
int global_size_motifs(3);
int global_nb_threads_to_decompose(4);
int global_nb_threads_to_write(2);


void usage(char *prog_name, const char *more) {
  cerr << more <<endl;
  cerr << "usage: " << prog_name << " -i input_file [-s motif_size] [-d nthread_decompose] [-w nthread_write]" << endl;
  cerr << "read an edgelist graph file and decompose it onto 3- or 4-node motifs using ESU algorithm." << endl <<endl;
  cerr << "Inputs :" << endl;
  cerr << "-s motif_size:\t\tsize of motifs (must be 3 or 4). Default is 3." << endl;
  cerr << "-d nthread_decompose:\tnumber of threads used by OpenMP to decompose the graph. Default is 4." << endl;
  cerr << "-w nthread_write:\tnumber of threads used by OpenMP to write the results. Default is 2." << endl;
  cerr << "-h\t\t\tshow this usage message" << endl<<endl;
  cerr << "Outputs :" << endl;
  cerr << "Results/motifs_occurrences.txt:\n\t\tnumber of occurrences of each k-node motif." << endl;
  cerr << "Results/BensonGraphs/motif_k_isomorphism_id.txt:\n\t\tthe k-id motif graph of the input graph on an edgelist format" << endl;
  cerr << "\t\tThis is a weighted graph W where w(i,j) is the number of co-occurrences\n\t\tof nodes i and j in the motif k-id." << endl;
  cerr << "Results/dict_graph:\n\t\tin the motif graphs, nodes have been relabelled.\n\t\tCorrespondance with the initial labels of nodes are in this file." <<endl;
  exit(0);
}

void parse_args(int argc, char **argv) {
  for (int i = 1; i < argc; i++) {
    if(argv[i][0] == '-') {
      switch(argv[i][1]) {
        case 'i':
	       if (i==argc-1)
	        usage(argv[0], "Infile missing\n");
	       infile = argv[i+1];
	       i++;
	       break;
        case 's':
	        if (i==argc-1)
	         usage(argv[0], "Size of motifs missing\n");
          global_size_motifs = atoi(argv[i+1]);
	        i++;
	        break;
        case 'd':
  	      if (i==argc-1)
  	        usage(argv[0], "Number of threads to decompose graph missing\n");
          global_nb_threads_to_decompose = atoi(argv[i+1]);
  	      i++;
  	      break;
        case 'w':
    	    if (i==argc-1)
    	      usage(argv[0], "Number of threads to write motif graph missing\n");
          global_nb_threads_to_write = atoi(argv[i+1]);
    	    i++;
    	    break;
          case 'h':
            usage(argv[0], "Help:");
            break;
        default:
	       usage(argv[0], "Unknown option\n");
      }
    } else {
      usage(argv[0], "More than one filename\n");
    }
  }
  if (infile==NULL ){
    usage(argv[0], "Infile missing\n");
  }else {
    string my_str(infile);
    global_graph_file = my_str;
  }
}


int main(int argc, char **argv)
{
  mkdir(global_outputs.c_str(),0777);
  mkdir((global_outputs+"/"+global_benson).c_str(),0777);
  // mkdir((global_outputs+"/"+global_louvain).c_str(),0777);

  parse_args(argc,argv);
  auto start = high_resolution_clock::now();
  rGraph G = rGraph(global_graph_file,true);
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop-start);
  printf("time spent to read the graph : %lu millisec.\n", duration.count());

  // G.display();



  // int size_motifs = 4;
  start = high_resolution_clock::now();
  try{
    Esu Sample = Esu(&G,global_size_motifs,global_nb_threads_to_decompose,global_nb_threads_to_write);
    auto temps1 = Sample.getTimeDecomp();
    auto temps2 = Sample.getTimeSaveFile();
    // printf("time spent to decompose the graph : %lu millisec.\n", temps1.count());
    // printf("time spent to write the results : %lu millisec.\n", temps2.count());


    // Sample.decompose();
    stop = high_resolution_clock::now();
    printf("time spent to decompose the graph : %lu millisec.\n",   duration_cast<milliseconds>(stop-start).count());
  }
  catch(exception const& e){
    cerr << "ERROR : "<<e.what()<<endl;
  }

  return 0;
}
