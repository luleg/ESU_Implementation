#ifndef _GRAPH_M
#define _GRAPH_M

#include "../Global.h"


class rGraph
{
public:
  struct Node {
    int id;
    std::string label;
    std::unordered_map<Node*,char> voisins;
    int degre; // number of neighbours
    Node(std::string name);
    ~Node();
  };

  struct SimpleNode {
    std::map<int,char> *ptr_neighb;
    std::set<int> *voisins;
  };



  static std::vector<SimpleNode>* readFromFile(std::string path_to_file,bool save_dict);

  // Attributes :
  std::vector<SimpleNode> *nodes;
  // Constructors :
  rGraph();
  rGraph(std::string path_to_file,bool save_dict);
  // Destructor :
  ~rGraph();

  // Other
  void display();



};


#endif
