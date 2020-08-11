#include "Esu.h"


using namespace std;
using namespace std::chrono;


// ESU constructor: it is provided a list of simple nodes and the targetted size of motifs. From here, it itsel builds its own isomorphism iso
Esu::Esu(rGraph const *one_graph,int s_mots,int thr_dec, int thr_write) : graph(one_graph), global_size_motifs(s_mots)//lnodes(&(graph->nodes))
{
  // Isomorphism reads a dictionary of isomorphisms in a file. Here, we only have these files for 3- and 4-node motifs, hence...
  string file_name;
  if (global_size_motifs == 3){
    file_name = dico_3;
  }
  else if(global_size_motifs == 4){
    file_name = dico_4;
  }
  else{
    throw invalid_argument("Esu::Constructor::The algorithm is implemented for 3 and 4 node motifs only (given : "+ to_string(global_size_motifs)+ ").");
  }
  iso = new Isomorphism(file_name,graph->nodes->size(),s_mots,thr_write);
  auto start = high_resolution_clock::now();
  this->decompose(thr_dec);
  auto stop = high_resolution_clock::now();
  time_to_decompose = duration_cast<milliseconds>(stop-start);

  start = high_resolution_clock::now();
  iso->save_in_files();
  stop = high_resolution_clock::now();
  time_to_write = duration_cast<milliseconds>(stop-start);
}

milliseconds Esu::getTimeDecomp(){
  return this->time_to_decompose;
}
milliseconds Esu::getTimeSaveFile(){
  return this->time_to_write;
}


void Esu::decompose(int global_nb_threads_to_decompose)
{

  omp_set_num_threads(global_nb_threads_to_decompose); // Number of threads working on parallel
  int tid(-1);
  #pragma omp parallel for private(tid) schedule(dynamic,1)
  for (int i=0; i < graph->nodes->size();++i) // for each node
  {
    tid = omp_get_thread_num();
    if(i%5000 == 0){
      tid = omp_get_thread_num();
      printf("Thread %d doing node %d\n",tid,i);
    }
    set<int> *V_keep = new set<int>(); // Set of nodes that generate a motif (or not)
    set<int>::iterator iter_V_keep = V_keep->begin();
    V_keep->insert(iter_V_keep,i);
    int id_init = i;

    set<int> *V_test = new set<int>();
    set<int>::iterator iter_V_test = V_test->begin();


    for (set<int>::iterator it_v = graph->nodes->at(i).voisins->begin();it_v!=graph->nodes->at(i).voisins->end();++it_v){
      if ((*it_v)>id_init){
        iter_V_test = V_test->insert(iter_V_test,*it_v); // Insertion in a set is log(N) (N=Card(V_test)) in worst case. With the iterator it is amortized constant
      }
    }
    Esu::rec_decompose(V_keep,V_test,&id_init);
  }
}

int Esu::rec_decompose(set<int> *V_keep, set<int> *V_test, int *id_init){

  int size_keep = V_keep->size();

  // terminal case 1 : V_keep is of size size_motifs
  if (size_keep == global_size_motifs)
  {
    int id_motif = Esu::compute_id(V_keep);
    iso->update(id_motif,V_keep);
    V_keep->clear();
    delete(V_keep);
    V_test->clear();
    delete(V_test);
    return 0;
  }
  // terminal case 2: this branch is dead, we cannot go further
  if (V_test->empty())
  {
    V_keep->clear();
    delete(V_keep);
    V_test->clear();
    delete(V_test);
    return 0;
  }
  // Non terminal case (recursive part) : each node v in V_test is added to V_keep (the motif is extended by 1 node), and removed from V_test.
  // Simultaneously, and independently, V_test is extended with the neighbours of v that
  // -have a label (here, an id) higher than the id of the root (the 1st node in V_keep)
  // -has no other neighbour than v in V_keep

  int size_Vt = V_test->size(); // constant
  while (size_Vt>0)
  {
    set<int>::iterator it_v = V_test->begin(); // Complexity constant
    int node_v = *it_v;
    V_test->erase(it_v); size_Vt = V_test->size(); // Complexity constant for both erase and size
    set<int> *V_keep_rec = new set<int>(*V_keep); // linear but very small set
    V_keep_rec->insert(node_v); // log(n), but very small set

    set<int> *V_test_rec = new set<int>();

    if (size_keep<global_size_motifs-1){ // If size_keep= size_motifs-1, it means V_keep_rec is a motif, and there is no need to update V_test as it will not be analysed in the next section
      V_test_rec = new set<int>(*V_test); // linear
      set<int>::iterator ite_Vtest_rec = V_test_rec->begin(); // constant
      int node_w;
      for (set<int>::iterator iter_w = graph->nodes->at(node_v).voisins->begin();iter_w!=graph->nodes->at(node_v).voisins->end();++iter_w){
        node_w = *iter_w;
        if (node_w<=*id_init){
          continue;
        }
        set<int>::iterator iter_Nw;
        bool no_neighb = true;
        for (set<int>::iterator iter_V_keep = V_keep->begin();iter_V_keep!=V_keep->end();++iter_V_keep){
          iter_Nw = graph->nodes->at(node_w).voisins->find(*iter_V_keep);
          if (iter_Nw != graph->nodes->at(node_w).voisins->end()) {
            no_neighb = false;
            break;
          }
        }
        if (no_neighb){
          ite_Vtest_rec = V_test_rec->insert(ite_Vtest_rec,node_w);
        }
      }
    }
    Esu::rec_decompose(V_keep_rec,V_test_rec,id_init);
  }
}

int Esu::compute_id(std::set<int> const *V_keep){
  int id_motif = 0;
  int cpt_u = 0;
  for (set<int>::iterator it_u = V_keep->begin();it_u != V_keep->end();it_u++)
  {
    rGraph::SimpleNode u = graph->nodes->at(*it_u);
    set<int>::iterator it_v(it_u);
    it_v++;
    int cpt_v = cpt_u+1;
    while (it_v!=V_keep->end()){
      map<int,char>::iterator ite_edge = u.ptr_neighb->find(*it_v);
      if (ite_edge==u.ptr_neighb->end()){
        it_v++;
        cpt_v++;
        continue;
      }
      switch(ite_edge->second){
        case '+':
          id_motif+=pow(2,global_size_motifs*(cpt_v)+cpt_u);
          break;
        case '-':
          id_motif+=pow(2,global_size_motifs*(cpt_u)+cpt_v);
          break;
        case '=':
          id_motif+=pow(2,global_size_motifs*(cpt_v)+cpt_u)+pow(2,global_size_motifs*(cpt_u)+cpt_v);
          break;
        default:
          throw domain_error("Esu::Compute_id::There is an invalid character to represent an edge.");
      }
      it_v++;
      cpt_v++;
    }
    cpt_u++;
  }
  return id_motif;
}
