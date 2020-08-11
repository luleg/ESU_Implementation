#include "Isomorphism.h"


using namespace std;


void Isomorphism::BensonNode::clear(){
  this->voisins.clear();
};

void Isomorphism::BensonGraph::clear(){
  omp_destroy_lock(&this->lock_Benson);
  for(vector<BensonNode>::iterator it = this->lnodes.begin();it!=this->lnodes.end();it++){
    it->clear();
  }
  this->lnodes.clear();
};

Isomorphism::Isomorphism(){
  nb_nodes = 0;
  global_nb_threads_to_write = 1;
}

Isomorphism::Isomorphism(string path_to_file,int n,int s_mot, int thr_write) : nb_nodes(n), global_size_motifs(s_mot), global_nb_threads_to_write(thr_write){
  // Initialisation of Benson matrices : One per isomorphism, all of the same size as the decomposed graph
  // this->nb_nodes = nb_nodes; // number of nodes in the graph that will be decomposed
  map<int,Isomorphism::BensonGraph*> visited_iso; // To know if an isomorphism has already been discovered, and have a trace of the corresponding BensonGraph

  ifstream myfile(path_to_file); // The isomorphisms are stored in dico_motif3 and dico_motif4 for 3-node and 4-node motifs, respectively. So we read the file and update the instancied object
  if (myfile){
    string id_src,id_tgt;
    map<int,Isomorphism::BensonGraph*>::iterator iter_bgr;
    while(!myfile.eof()){
      myfile >>id_src;
      if(myfile.eof()){
        break;
      }
      myfile >> id_tgt;
      iter_bgr = visited_iso.find(stoi(id_tgt));
      if (iter_bgr != visited_iso.end())// Si le Benson graphe existe deja, on ajoute juste le pointeur.
      {
        lBGraphs[stoi(id_src)] = iter_bgr->second;
      }else // Sinon, on cree le Benson Graphe: une liste de nb_nodes BensonNode.
      {
        lid_iso.push_back(stoi(id_tgt));
        Isomorphism::BensonGraph *crt_bensonGraph = new Isomorphism::BensonGraph();
        visited_iso[stoi(id_tgt)] = crt_bensonGraph;
        lBGraphs[stoi(id_src)] = crt_bensonGraph;
      }

    }
  }else{
    throw invalid_argument("Isomorphism::Constructor:: Unable to open the file \""+path_to_file+"\". Are you sure this file is in the code source directory?");
  }
}

Isomorphism::~Isomorphism(){
  for (vector<int>::iterator it = lid_iso.begin();it!=lid_iso.end();it++){
    lBGraphs[*it]->clear();
  }
  lBGraphs.clear();
}

void Isomorphism::update(int id_iso,set<int> *V_keep){
  map<int,Isomorphism::BensonGraph*>::iterator to_Benson = lBGraphs.find(id_iso);
  Isomorphism::BensonGraph* ptr_Benson = to_Benson->second;
  // Recupere le verrou !!
  omp_set_lock(&ptr_Benson->lock_Benson);
  if (ptr_Benson->lnodes.size() == 0){ // If this is the first time such a motif is encountered, the corresponding BensonGraph is built (The graph may be very large so we choose to built such graphs on the fly to avoid using too much memory)
    ptr_Benson->lnodes = vector<Isomorphism::BensonNode>(nb_nodes);
  }
  for (set<int>::iterator it_u = V_keep->begin();it_u != V_keep->end();++it_u)
  {
      set<int>::iterator it_v(it_u);
      it_v ++;
      while (it_v!=V_keep->end())
      {
        ptr_Benson->lnodes[*it_u].voisins[*it_v]+=1;
        it_v++;
      }
  }
  omp_unset_lock(&ptr_Benson->lock_Benson);
}


// Debugg function
void Isomorphism::display(int id_iso){
  map<int,Isomorphism::BensonGraph*>::iterator to_Benson = lBGraphs.find(id_iso);
  Isomorphism::BensonGraph* ptr_Benson = to_Benson->second;
  printf("-----------------------Benson Graph of the isomorphism %d\n",id_iso);
  for (int i = 0; i<nb_nodes; i++){
    printf("node %d : [",i);
    map<int,int> i_neighb = ptr_Benson->lnodes[i].voisins;
    for (map<int,int>::iterator it = i_neighb.begin();it!=i_neighb.end();it++){
      printf(" (%d : %d) ",it->first,it->second);
    }
    printf(" ]\n");
  }
}

void Isomorphism::save_in_files(){

  string to_save;
  string file_to;
  string name_folder(global_outputs+"/"+global_benson);
  string motif_decomp = "";

  omp_set_num_threads(global_nb_threads_to_write);
  //
  #pragma omp parallel for schedule(dynamic,1) private(to_save,file_to) reduction(merge:motif_decomp)
  for(int iter = 0;iter<lid_iso.size();iter++){
    to_save = "";
    file_to = "";
    Isomorphism::BensonGraph* crt_Benson = lBGraphs[lid_iso[iter]];
    if(crt_Benson->lnodes.size()==0){
      continue;
    }
    int nb_motifs = 0;
    for(int i = 0;i<nb_nodes;++i){
      for(map<int,int>::iterator it_v = crt_Benson->lnodes[i].voisins.begin();it_v != crt_Benson->lnodes[i].voisins.end();++it_v){
        to_save+=to_string(i)+" "+to_string(it_v->first)+" "+to_string(it_v->second)+"\n";
        nb_motifs+=it_v->second;
      }
    }
    file_to = name_folder+"/motif_"+to_string(global_size_motifs)+"_isomorphism_"+to_string(lid_iso[iter])+".txt";
    ofstream wfile(file_to.c_str());
    if(wfile){
      wfile << to_save;
    }
    wfile.close();
    crt_Benson->clear();

    nb_motifs = 2*nb_motifs/(global_size_motifs*(global_size_motifs-1));
    motif_decomp+=to_string(global_size_motifs)+"-"+to_string(lid_iso[iter])+" "+to_string(nb_motifs)+"\n";
  }

  motif_decomp = "id_motif occurrences_number\n"+motif_decomp;
  file_to = global_outputs+"/motifs_occurrences.txt";
  ofstream wfile(file_to.c_str());
  if(wfile){
    wfile << motif_decomp;
  }
  wfile.close();
}
