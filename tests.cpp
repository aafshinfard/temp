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

//static uint64_t
using namespace std;

int
main(int argc, char* argv[])
{
	vector< vector<int> > a(3, vector<int>(3,4));
	vector< vector<int> > b(a);
	cout<<b.size();

	return 0;
}
