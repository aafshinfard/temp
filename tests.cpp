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
#include <utility>
#include <tgmath.h>
#include <stdexcept>
#include <algorithm>
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
typedef std::chrono::high_resolution_clock::time_point TimeVar;


#define duration(a) std::chrono::duration_cast<std::chrono::microseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

TimeVar time_start;
long long total_duration=0;

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
square_matrix_bitwise(
    adjacencyMatrix_t M,
    bool symmetric=true)
{
    // transform the input matrix into a vector of integers
    // traverse over and apply bitwise and to achieve multiplication in a faster manner
    // transform the resulting vecot to adjacency matrix if needed ?! or continue using it with the same structure
}

adjacencyMatrix_t
square_matrix_ijk2(
    adjacencyMatrix_t M,
    bool symmetric=true)
{
    int n = M.size();

    adjacencyVector_t tempVector(n, 0); // Fast initialization
    adjacencyMatrix_t M2(n, tempVector);
    adjacencyMatrix_t::iterator M_iter = M.begin();
    adjacencyMatrix_t::iterator M2_iter = M2.begin();

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

vector<vector<double> >
square_matrix_ikj_d( // Might be faster than ijk, benchmark it
    vector<vector<double> > M,
    bool symmetric = true)
{
    // Fast initialization:
    int n = M.size();
    vector<double> tmp(n, 0);
    vector<vector<double> > M2(n, tmp);
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
inline double
cosine_similarity_vectors(
    adjacencyMatrix_t::iterator& row_i,
    adjacencyMatrix_t::iterator& row_j)
{
    // Input: 2 vectors (1D) as rows and columns of a Matrix
    // Output: Cosine similarity of the two vectors
    // (Cosine Similarity between 2 corresponding vertices)

    float mul = 0.0; // also test double
    float d_i = 0.0;
    float d_j = 0.0;

    if (row_i->size() != row_j->size())
    {
        throw logic_error("Vector A and Vector B are not the same size");
    }

    // Prevent Division by zero
    if (row_i->size() < 1)
    {
        throw logic_error("Input vectors for multiplication are empty");
    }

    adjacencyVector_t::iterator i_iter = row_i->begin();
    adjacencyVector_t::iterator j_iter = row_j->begin();
    for( ; i_iter != row_i->end(); i_iter++ , j_iter++ )
    {
        mul += *i_iter * *j_iter;
        d_i += *i_iter * *i_iter;
        d_j += *j_iter * *j_iter;
        // cout<<"\nDebug - mul:"<<mul<<" - d_i:"<<d_i<<" - d_j:"<<d_j<<endl;
    }
    if (mul == 0.0f)
    {
        return 0;
    }
    if (d_i == 0.0f || d_j == 0.0f)
    {
        return 0;
//        throw std::logic_error(
//                "cosine similarity is not defined whenever one or both "
//                "input vectors are zero-vectors.");
    }
    //return mul / (sqrt(d_a) * sqrt(d_b));
    return mul / sqrt(d_i * d_j);
}

inline
void
calculate_cosine_similarity_2d(
    adjacencyMatrix_t& adj_mat,
    vector<vector<double> >& cosimilarity)
{
    // NOT COMPLETE YET:
    // STRATEGY: NORMALIZE THEN SQUARE (instead of normalizing per vector while multip
    time_start = timeNow();
    int n = adj_mat.size();
    vector<double> temp(n, 0.0);
    vector<vector<double> > normalized(n, temp);
    double row_sum = 0;
    uint_fast32_t init = 0;

    adjacencyMatrix_t::iterator row_i;
    vector<vector<double> >::iterator normalized_row_i = normalized.begin();
    for (row_i = adj_mat.begin(); row_i != adj_mat.end(); ++row_i, ++normalized_row_i)
    {
        row_sum = 0;
        vector<uint_fast32_t>::iterator first = row_i->begin();
        vector<uint_fast32_t>::iterator last = row_i->end();
        while(first!=last){
            row_sum += *first * *first;
            ++first;
        }
        row_sum = sqrt(row_sum);
        //row_sum = accumulate(row_i->begin(), row_i->end(), init);

        first = row_i->begin();
        vector<double>::iterator first_normalized = normalized_row_i->begin();
        vector<double>::iterator last_normalized = normalized_row_i->end();
        while(first!=last){
            *first_normalized= (double)*first / row_sum;
            ++first;
            ++first_normalized;
        }
    }
//    time_duration = duration(timeNow() - time_start );
    total_duration += std::chrono::duration_cast<std::chrono::microseconds>( timeNow() - time_start ).count();
    time_start = timeNow();
    cosimilarity = square_matrix_ikj_d(normalized);
    total_duration += std::chrono::duration_cast<std::chrono::microseconds>( timeNow() - time_start ).count();
}

inline
void
calculate_cosine_similarity_2d_v2(
    adjacencyMatrix_t& adj_mat,
    vector<vector<double> >& cosimilarity)
{
    // Assumptions: the input matrix is symmetric and cubic
    // This function calculate the 2-dimensional cosine similarity of the input matrix
    // to itself, that is the similarity between vertices of the corresponding graph
    // for the input matrix (as adj matrix)
    adjacencyMatrix_t::iterator row_i;
    adjacencyMatrix_t::iterator row_j;
    vector<vector<double> >::iterator out_row = cosimilarity.begin();
    vector<double>::iterator out_col = out_row->begin();
    if (false/* inputs are not of same size*/)
        //Throw error
        ;

    int i = 0;
    int j = 0;
    //cosimilarity[0][0] = 4;
    for (row_i = adj_mat.begin(); row_i != adj_mat.end(); ++row_i)
//    for (int row_i = 0; row_i 4; row_i++)
    {
        j = 0;
        for (row_j = adj_mat.begin(); row_j != adj_mat.end(); ++row_j)
//        for (int row_j = 0; row_j 4; row_j++)
        {
//            cout<<" "<<i<<" "<<j<<" started!"<<endl;
            if (j < i)
            {
//                cosimilarity[i][j] = cosimilarity[j][i];
//                cout<<" "<<i<<" "<<j<<" copied!"<<endl;
            }
            else
            {
//                *out_col = cosine_similarity_vectors(row_i, row_j);
                cosimilarity[i][j] = cosine_similarity_vectors(row_i, row_j);
//                cout<<" "<<i<<" "<<j<<" done!"<<endl;
            }
            j += 1;
//            out_col++;
        }
        i += 1;
//        out_row++;
    }
}


int
main(int argc, char* argv[])
{
	//adjacencyMatrix_t a(3, adjacencyVector_t(5,4));
	//adjacencyMatrix_t b(a);
	//cout<<b.size()<<endl;
	//cout<<b[1][1];
    //func(b);

//    int aaa = 5;
//    int bbb = 3;
//    int ccc;
//    ccc = aaa & bbb;
//    cout<<"\n answer:"<<ccc<<endl;
//    return 0;

    //adjacencyMatrix_t a_small(3, adjacencyVector_t(3,4));
//    adjacencyMatrix_t a_small {
//				{ 0, 12, 4, 1 },
//				{ 12, 0, 1, 0 },
//				{ 4, 1, 0, 0 },
//				{ 1, 0, 0, 0 }
//			};
    adjacencyMatrix_t a_small {
				{ 0, 1, 1, 0 },
				{ 1, 0, 0, 1 },
				{ 1, 0, 0, 1 },
				{ 0, 1, 1, 0 }
			};

	int test_size = 50;
    adjacencyMatrix_t a(test_size, adjacencyVector_t(test_size,4));
    vector<double> tempVector(test_size, 0);
    vector<vector<double> > cosSimilarity2d(test_size, tempVector);
    vector<vector<double> > cos1(cosSimilarity2d);

    vector<double> tempVector_small(4, 0.0);
    vector<vector<double> > cosSimilarity2d_small(4, tempVector_small);


    calculate_cosine_similarity_2d_v2(a_small, cosSimilarity2d_small);
    calculate_cosine_similarity_2d(a_small, cos1);

    vector<vector<double> >::iterator rowrow1 = cosSimilarity2d_small.begin();

//    for ( ; rowrow1 != cosSimilarity2d_small.end(); ++rowrow1 )
//    {
//        for (vector<double>::iterator col = rowrow1->begin(); col != rowrow1->end(); ++col){
//            cout<<" "<<*col;
//        }
//        cout<<endl;
//    }
//    cout<<"NEXT:"<<endl;
//    vector<vector<double> >::iterator rowrow2 = cos1.begin();
//    for ( ; rowrow2!= cos1.end(); ++rowrow2 )
//    {
//        for (vector<double>::iterator col = rowrow2->begin(); col != rowrow2->end(); ++col){
//            cout<<" "<<*col;
//        }
//        cout<<endl;
//    }
//    return 0;
//
//    vector<vector<double> >::iterator row_d;
//    vector<double>::iterator col_d;
//    for (row_d = normalized.begin(); row_d != normalized.end(); ++row_d)
//    {
//        for (col_d = row_d->begin(); col_d != row_d->end(); ++col_d)
//        {
//            cout<<" "<<*col_d;
//        }
//        cout<<endl;
//    }
//    return 0;
//    adjacencyMatrix_t::iterator row_i;
//    uint_fast32_t ada = 0;
//    uint_fast32_t init = 0;
//
////    for (row_i = a_small.begin(); row_i != a_small.end(); ++row_i)
//    auto start = high_resolution_clock::now();
//    for (row_i = a.begin(); row_i != a.end(); ++row_i)
//    {
//        ada = accumulate(row_i->begin(), row_i->end(), init);
////        cout<<" sum is:"<<ada<<endl;
////        ada = 0;
////        for (vector<uint_fast32_t>::iterator col = row_i->begin(); col != row_i->end(); ++col){
////            ada += *col;
////        }
////        cout<<" sum is:"<<ada<<endl;
//    }
//    auto stop = high_resolution_clock::now();
//    auto duration = duration_cast<microseconds>(stop - start);
//    cout<<"duration:"<<duration.count()<<endl;
//    return 0;
//
//
//
//    calculate_cosine_similarity_2d_v2(a_small, cosSimilarity2d_small);
//    vector<vector<double> >::iterator row = cosSimilarity2d_small.begin();
//    for ( ; row != cosSimilarity2d_small.end(); ++row )
//    {
//        for (vector<double>::iterator col = row->begin(); col != row->end(); ++col){
//            cout<<" "<<*col;
//        }
//        cout<<endl;
//    }




//    start = timeNow();
    for( int i=0 ; i < 1000 ; i++){
        cout<<i<<endl;
        a[0][0] = a[0][0] + 1;
        adjacencyMatrix_t temp(square_matrix_ijk(a));
        if(temp[0][0]<0){
            cout<<"Failed";
            return 0;
        }
        calculate_cosine_similarity_2d(a,cosSimilarity2d);
//        calculate_cosine_similarity_2d(a,cosSimilarity2d);
//        if(! test_fast_initialization1() ){
//            cout<<"Failed";
//            return 0;
//        }
    }
//    auto stop = high_resolution_clock::now();
//    auto duration = duration_cast<microseconds>(stop - start);
//    duration = duration_cast<microseconds>(stop - start);
//    cout << "\nElapsed time: "<<duration.count() << endl;
//    cout << "\nElapsed time: "<<duration(timeNow() - start) << endl;
    cout << "\nElapsed time: "<<total_duration << endl;
	return 0;
}
