#ifndef _ESU_
#define _ESU_

#include "Isomorphism.h"



class Esu{
  private :
  // Attributes :
    rGraph const *graph;
    Isomorphism *iso;
    int const global_size_motifs;

    std::chrono::milliseconds time_to_decompose;
    std::chrono::milliseconds time_to_write;

  // Private methods
    void decompose(int global_nb_threads_to_decompose); // ESU algorithm
    int rec_decompose(std::set<int> *V_keep, std::set<int> *V_test,int *id_init); // Recursive part of the ESU alg
    int compute_id(std::set<int> const *V_keep); // Finding the isomorphism of a motif occurrence


  public :
  // Public methods :
    Esu(rGraph const *one_graph,int s_mots,int thr_dec, int thr_write); // Constructor. It launchs the algorithm

    std::chrono::milliseconds getTimeDecomp();
    std::chrono::milliseconds getTimeSaveFile();


};

#endif
