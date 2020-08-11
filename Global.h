#ifndef _GLOBAL_
#define _GLOBAL_

#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>

#include <omp.h>

#include <math.h>
#include <chrono>

#include <map>
#include <unordered_map>
#include <set>
#include <vector>
#include <list>

#pragma omp declare reduction(merge : std::string : omp_out+=omp_in)

const std::string global_outputs = "Results"; // where the outputs are saved
const std::string global_benson = "BensonGraphs"; // where the BensonGraphs are saved (in a folder : $global_outputs/global_benson)


const std::string dico_3 = "src/dico_motif3.txt";
const std::string dico_4 = "src/dico_motif4.txt";

#endif
