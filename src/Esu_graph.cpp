#include "Esu_graph.h"

using namespace std;


rGraph::Node::Node(string name)
{
  id = -1;
  label = name;
  degre = 0;
}

rGraph::Node::~Node()
{
  voisins.clear();
}



rGraph::rGraph() {
  vector<rGraph::SimpleNode> *nodes = new vector<rGraph::SimpleNode>();
}


rGraph::rGraph(string path_to_file,bool save_dict) {
  nodes = rGraph::readFromFile(path_to_file,save_dict);
}

rGraph::~rGraph(){
  for (vector<SimpleNode>::iterator it = nodes->begin();it!=nodes->end();++it){
    it->ptr_neighb->clear(); delete it->ptr_neighb;
    it->voisins->clear(); delete it->voisins;
  }
  nodes->clear();
  delete nodes;
}

vector<rGraph::SimpleNode>* rGraph::readFromFile(string path_to_file,bool save_dict){
  unordered_map<string,rGraph::Node*> visited_nodes;
  ifstream myfile(path_to_file);
  int cpt_edge(0);
  if (myfile)
  {
    while(!myfile.eof())
    {
      string src,tgt;
      myfile >>src;// here is the source node label
      if(myfile.eof()){
        break;
      }
      cpt_edge++;
      myfile >> tgt;// here is the target node label
      Node *node_src,*node_tgt;
      // If one node (or both) does not exist we need to create them
      // If they exist we have to update them
      unordered_map<string,rGraph::Node*>::iterator map_src = visited_nodes.find(src);
      unordered_map<string,rGraph::Node*>::iterator map_tgt = visited_nodes.find(tgt);
      bool is_src_vis = map_src!=visited_nodes.end(); if (is_src_vis){node_src = map_src->second;}
      bool is_tgt_vis = map_tgt!=visited_nodes.end(); if (is_tgt_vis){node_tgt = map_tgt->second;}
      //src and tgt exist:
      if(is_src_vis && is_tgt_vis){
        // is node_tgt already a neighbour of node_src (or the contrary)
        Node *node_from(node_src), *node_toward(node_tgt);
        bool src_is_from=true;
        if (node_src->degre > node_tgt->degre){ // to improve the computationality
          node_from =node_tgt;
          node_toward = node_src;
          src_is_from=false;
        }
        unordered_map<rGraph::Node*,char>::iterator voisin = node_from->voisins.find(node_toward);
        // Not neighbours
        if (voisin == node_from->voisins.end()){ // not neighbours : we add them
            node_src->voisins[node_tgt] = '-';
            node_src->degre+=1;
            node_tgt->voisins[node_src] = '+';
            node_tgt->degre+=1;
        }
        // Already neighbours
        else if((voisin->second == '+' && src_is_from)||(voisin->second == '-' && !src_is_from)){
          node_from->voisins[node_toward] = '=';
          node_toward->voisins[node_from] = '=';
        }
      }//src and/or tgt do not exist: we have to create them
      else {
        if(!is_src_vis){ // src does not exist
          node_src = new rGraph::Node(src); // creation of the pointer
          visited_nodes[src] = node_src; // we add it to visited map
        }
        if(!is_tgt_vis){ // tgt does not exist
          node_tgt = new rGraph::Node(tgt);
          visited_nodes[tgt] = node_tgt;
        }
        // Now the two nodes have been built, we add the corresponding edge
        node_src->voisins[node_tgt] = '-';
        node_src->degre+=1;
        node_tgt->voisins[node_src] = '+';
        node_tgt->degre+=1;
      }
    }
    printf("number of edges ? %d\n",cpt_edge);
  }
  else{
    throw invalid_argument("Graph::ReadFromFile:: Unable to open the file \""+path_to_file+"\". Are you sure this is a correct network name?");
  }
  // To sort nodes by ascending order of degrees (label does not matter)
  set<pair<int,Node*>> sort_nodes;

  for (unordered_map<string,Node*>::iterator it = visited_nodes.begin();it!=visited_nodes.end();++it)
  {
    sort_nodes.insert(pair<int,Node*>((it->second)->degre,it->second));
  }

  visited_nodes.clear();
  int cpt(0);
  for (set<pair<int,Node*>>::iterator it = sort_nodes.begin();it!=sort_nodes.end();++it)
  {
    it->second->id = cpt;
    cpt++;
  }
  vector<rGraph::SimpleNode> *snodes = new vector<rGraph::SimpleNode>();
  string to_save("id_node\tlabel_node");
  for (set<pair<int,Node*>>::iterator it = sort_nodes.begin();it!=sort_nodes.end();++it)
  {
    to_save+="\n"+to_string(it->second->id)+"\t"+it->second->label;
    rGraph::SimpleNode scrt_node;
    map<int,char> *scrt_neighb = new map<int,char>() ;
    set<int> *scrt_voisins = new set<int>();
    for (unordered_map<Node*,char>::iterator it_v = it->second->voisins.begin();it_v != it->second->voisins.end();++it_v){
      scrt_neighb->insert(pair<int,char>(it_v->first->id,it_v->second));
      scrt_voisins->insert(it_v->first->id);
    }
    scrt_node.ptr_neighb = scrt_neighb;
    scrt_node.voisins = scrt_voisins;
    snodes->push_back(scrt_node);
  }

  if(save_dict){
      ofstream myFile;
      myFile.open(global_outputs+"/dict_graph.txt");
      myFile << to_save;
    }
  return snodes;
}


void rGraph::display(){
  printf("--------Display of Graph : %lu nodes (only 20 displayed) ----------\n",nodes->size());
  int cpt(0);
  for(vector<rGraph::SimpleNode>::iterator it = nodes->begin();it!=nodes->end();++it){
    if (cpt<20){
      int cpt2(0);
      printf("node : %d -- %lu neighbours \t [",cpt,(it->voisins)->size());

      cpt++;
      for(map<int,char>::iterator it2 = it->ptr_neighb->begin();it2 != it->ptr_neighb->end();++it2){
        if (cpt2<5){
          cpt2++;
          string dir;
          if (it2->second == '-') dir ="-->";
          else if (it2->second == '+') dir = "<--";
          else dir = "<->";
          printf(" (%s , %d) ",dir.c_str(),it2->first);
        }
        else if(cpt2==5)
        {
          cpt2++;
          printf(" ... ");
        }
      }
      printf("]\n");
    }
    else if (cpt == 2){
      cpt++;
      printf(" ... \n");
    }
}
}
