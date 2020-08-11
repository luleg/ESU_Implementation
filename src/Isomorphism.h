#ifndef _ISOMORPHISM_
#define _ISOMORPHISM_

#include "Esu_graph.h"


class Isomorphism
{
private:

  struct BensonNode{
    std::map<int,int> voisins;
    void clear();
  };

  struct BensonGraph{
    omp_lock_t lock_Benson;
    std::vector<BensonNode> lnodes;
    BensonGraph(){
      omp_init_lock(&lock_Benson);
    };
    void clear();
  };

  std::map<int,BensonGraph*> lBGraphs;
  std::vector<int> lid_iso;
  int nb_nodes;
  int global_size_motifs;
  int global_nb_threads_to_write;

public:


  Isomorphism(std::string path_to_file,int n,int s_mot, int thr_write);
  Isomorphism();
  ~Isomorphism();

  void update(int id_iso,std::set<int> *V_keep);
  void display(int id_iso);
  void save_in_files();
};



#endif
