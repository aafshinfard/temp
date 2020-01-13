//  g++ -std=c++11 tests.cpp
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/biconnected_components.hpp>
//#include <boost/graph/graph_utility.hpp>
//#include <boost/graph/subgraph.hpp>

#if _OPENMP
#include <omp.h>
#endif

//static uint64_t;
using namespace std;


using adjacencyMatrix_t = vector<vector<int>>;

void func(adjacencyMatrix_t& temp){
    cout<<temp[1][1]<<endl;
}

int
main(int argc, char* argv[])
{
	adjacencyMatrix_t a(3, vector<int>(5,4));
	adjacencyMatrix_t b(a);
	//cout<<b.size()<<endl;
	//cout<<b[1][1];
    func(b);

	return 0;
}
