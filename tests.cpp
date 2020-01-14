//  g++ -std=c++11 tests.cpp
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdint.h>
#include <chrono>
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/biconnected_components.hpp>
//#include <boost/graph/graph_utility.hpp>
//#include <boost/graph/subgraph.hpp>

#if _OPENMP
#include <omp.h>
#endif

//static uint64_t
using namespace std;
using namespace std::chrono;

using adjacencyMatrix_t = vector<vector<uint_fast32_t>>;
using adjacencyVector_t = std::vector<uint_fast32_t>;


void func(adjacencyMatrix_t& temp){
    cout<<temp[1][1]<<endl;
}


int test_fast_initialization1(){
    adjacencyVector_t tmp(50, 0);
    adjacencyMatrix_t mat(50, tmp);
    return mat[0][0]+1;
}
int test_fast_initialization2(){
    adjacencyMatrix_t mat(50, adjacencyVector_t(50, 0));
    return mat[0][0]+1;
}

adjacencyMatrix_t
square_matrix_ijk(
    adjacencyMatrix_t M,
    bool symmetric=true)
{
    int n = M.size();

    adjacencyVector_t tempVector(n, 0); // Fast initialization
    adjacencyMatrix_t M2(n, tempVector);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if ( j < i && symmetric )
                M2[i][j] = M2[j][i];
                continue;
            for (int k = 0; k < n; k++) {
                M2[i][j] += M[i][k] * M[k][j];
            }
        }
    }
    return M2;
}


adjacencyMatrix_t
square_matrix_ikj( // Might be faster than ijk, benchmark it
    adjacencyMatrix_t M,
    bool symmetric = true)
{
    // Fast initialization:
    int n = M.size();
    adjacencyVector_t tmp(n, 0);
    adjacencyMatrix_t M2(n, tmp);
    // Multiplication
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
//            if ( !M[i][k] )
//                continue;
            for (int j = 0; j < n; j++) {
                if ( j < i && symmetric ) {
                    M2[i][j] = M2[j][i];
                    continue;
                }
                M2[i][j] += M[i][k] * M[k][j];
            }
        }
    }
    return M2;
}

//boost::numeric::ublas::matrix<int>
//square_matrix_boost(
//    adjacencyMatrix_t M,
//)
//{
//    return boost::numeric::ublas::prod(M, M);
//}

int
main(int argc, char* argv[])
{
	//adjacencyMatrix_t a(3, adjacencyVector_t(5,4));
	//adjacencyMatrix_t b(a);
	//cout<<b.size()<<endl;
	//cout<<b[1][1];
    //func(b);
    adjacencyMatrix_t a(500, adjacencyVector_t(500,4));
    auto start = high_resolution_clock::now();
    for( int i=0 ; i < 10000 ; i++){
        adjacencyMatrix_t temp(square_matrix_ijk(a));
        if(temp[0][0]<0){
            cout<<"Failed";
            return 0;

        }
//        if(! test_fast_initialization1() ){
//            cout<<"Failed";
//            return 0;
//        }
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "\nElapsed time: "<<duration.count() << endl;

	return 0;
}
