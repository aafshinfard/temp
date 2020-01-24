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
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/subgraph.hpp>

#if _OPENMP
#include <omp.h>
#endif


memory_usage()
{
	int mem = 0;
	std::ifstream proc("/proc/self/status");
	for (std::string s; std::getline(proc, s);) {
		if (s.substr(0, 6) == "VmSize") {
			std::stringstream convert(s.substr(s.find_last_of('\t'), s.find_last_of('k') - 1));
			if (!(convert >> mem)) {
				return 0;
			}
			return mem;
		}
	}
	return mem;
}

struct vertexProperties
{
	std::string name = "";
	int weight = 0;
	size_t indexOriginal = 0;
};

struct edgeProperties
{
	int weight = 0;
};

struct edgeComponent_t
{
	enum
	{
		num = INT_MAX
	};
	using kind = boost::edge_property_tag;
} edgeComponent;

//static uint64_t
using namespace std;
using namespace std::chrono;
using graph_t = boost::subgraph<boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::undirectedS,
    vertexProperties,
    boost::property<
        boost::edge_index_t,
        int,
        boost::property<edgeComponent_t, std::size_t, edgeProperties>>>>;
using vertex_t = graph_t::vertex_descriptor;
using edge_t = graph_t::edge_descriptor;
using barcodeToIndex_t = std::unordered_map<std::string, vertex_t>;
using indexToBarcode_t = std::unordered_map<vertex_t, std::string>;
using vertexSet_t = std::unordered_set<vertex_t>;
using componentToVertexSet_t = std::vector<vertexSet_t>;
using vertexToComponent_t = std::unordered_map<vertex_t, size_t>;
using vecVertexToComponent_t = std::vector<vertexToComponent_t>;
using vertexToIndex_t = std::unordered_map<vertex_t, size_t>; // wanna improve this? checkout boost::bimap
using indexToVertex_t = std::unordered_map<size_t, vertex_t>; // wanna improve this? checkout boost::bimap
using adjacencyMatrix_t = std::vector<vector<uint_fast32_t>>;
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


static void
printVersion()
{
	const char VERSION_MESSAGE[] =
	    PROGRAM " (" PACKAGE_NAME ") " GIT_REVISION "\n"
	            "Written by Johnathan Wong.\n"
	            "\n"
	            "Copyright 2019 Canada's Michael Smith Genome Science Centre\n";
	std::cerr << VERSION_MESSAGE << std::endl;
	exit(EXIT_SUCCESS);
}

static void
printErrorMsg(const std::string& progname, const std::string& msg)
{
	std::cerr << progname << ": " << msg
	          << "\nTry 'physlr-molecules --help' for more information.\n";
}

static void
printUsage(const std::string& progname)
{
	std::cout << "Usage:  " << progname
	          << "  [-s SEPARATION-STRATEGY] [-v] FILE...\n\n"
	             "  -v         enable verbose output\n"
	             "  -s --separation-strategy   \n"
	             "  SEPARATION-STRATEGY      `+` separated list of molecule separation strategies "
	             "[bc]\n"
	             "  --help     display this help and exit\n";
}

void
printGraph(const graph_t& g)
{
	std::cout << "U\tn" << std::endl;
	auto vertexItRange = boost::vertices(g);
	for (auto vertexIt = vertexItRange.first; vertexIt != vertexItRange.second; ++vertexIt) {
		auto& node1 = g[*vertexIt].name;
		auto& weight = g[*vertexIt].weight;
		std::cout << node1 << "\t" << weight << "\n";
	}
	std::cout << "\nU\tV\tn" << std::endl;
	auto edgeItRange = boost::edges(g);
	for (auto edgeIt = edgeItRange.first; edgeIt != edgeItRange.second; ++edgeIt) {
		auto& weight = g[*edgeIt].weight;
		auto& node1 = g[boost::source(*edgeIt, g)].name;
		auto& node2 = g[boost::target(*edgeIt, g)].name;
		std::cout << node1 << "\t" << node2 << "\t" << weight << "\n";
	}
}

void
readTSV(graph_t& g, const std::vector<std::string>& infiles, bool verbose)
{
	auto progname = "physlr-molecules";
	std::cerr << "Loading graph" << std::endl;
#if _OPENMP
	double sTime = omp_get_wtime();
#endif
	barcodeToIndex_t barcodeToIndex;
	indexToBarcode_t indexToBarcode;
	for (auto& infile : infiles) {
		infile == "-" ? "/dev/stdin" : infile;
		std::ifstream infileStream(infile);
		for (std::string line; std::getline(infileStream, line);) {
			if (line == "U\tn") {
				continue;
			}
			if (line.empty()) {
				break;
			}
			std::string node1;
			int weight;
			std::istringstream ss(line);
			if (ss >> node1 >> weight) {
				auto u = boost::add_vertex(g);
				g[u].name = node1;
				g[u].weight = weight;
				g[u].indexOriginal = u;
				barcodeToIndex[node1] = u;
				indexToBarcode[u] = node1;
			} else {
				printErrorMsg(progname, "unknown graph format");
				exit(EXIT_FAILURE);
			}
		}

		if (verbose) {
			std::cerr << "Loaded vertices to graph ";
#if _OPENMP
			std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
			sTime = omp_get_wtime();
#endif
		}
		for (std::string line; std::getline(infileStream, line);) {
			if (line == "U\tV\tn") {
				continue;
			}
			if (line.empty()) {
				printErrorMsg(progname, "unknown graph format");
				exit(EXIT_FAILURE);
			}
			std::string node1, node2;
			int weight;
			std::istringstream ss(line);
			if (ss >> node1 >> node2 >> weight) {
				auto E = boost::add_edge(barcodeToIndex[node1], barcodeToIndex[node2], g).first;
				g[E].weight = weight;
			} else {
				printErrorMsg(progname, "unknown graph format");
				exit(EXIT_FAILURE);
			}
		}

		if (verbose) {
			std::cerr << "Loaded edges to graph ";
		} else {
			std::cerr << "Loaded graph ";
		}
#if _OPENMP
		std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
		sTime = omp_get_wtime();
#endif
		std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576) << "GB"
		          << std::endl;
	}
}

/* Generate a molecule separated graph (molSepG) using component/community information from
molecule separation (vecVertexToComponent). The input graph (inG) is the barcode overlap graph
or a molecule separated graph from the previous round of molecule separation.*/
void
componentsToNewGraph(
    const graph_t& inG,
    graph_t& molSepG,
    vecVertexToComponent_t& vecVertexToComponent)
{
	barcodeToIndex_t molSepGBarcodeToIndex;
#if _OPENMP
	double sTime = omp_get_wtime();
#endif
	for (size_t i = 0; i < vecVertexToComponent.size(); ++i) {

		size_t maxVal = 0;
		if (!vecVertexToComponent[i].empty()) {
			maxVal =
			    std::max_element(
			        vecVertexToComponent[i].begin(),
			        vecVertexToComponent[i].end(),
			        [](const vertexToComponent_t::value_type& p1,
			           const vertexToComponent_t::value_type& p2) { return p1.second < p2.second; })
			        ->second;
		}

		for (size_t j = 0; j < maxVal + 1; ++j) {
			auto u = boost::add_vertex(molSepG);
			molSepG[u].name = inG[i].name + "_" + std::to_string(j);
			molSepG[u].weight = inG[i].weight;
			molSepG[u].indexOriginal = u;
			molSepGBarcodeToIndex[molSepG[u].name] = u;
		}
	}

	auto edgeItRange = boost::edges(inG);
	for (auto edgeIt = edgeItRange.first; edgeIt != edgeItRange.second; ++edgeIt) {
		auto& u = inG[boost::source(*edgeIt, inG)].indexOriginal;
		auto& v = inG[boost::target(*edgeIt, inG)].indexOriginal;

		if (vecVertexToComponent[u].find(v) == vecVertexToComponent[u].end() ||
		    vecVertexToComponent[v].find(u) == vecVertexToComponent[v].end()) {
			continue;
		}

		auto& uMolecule = vecVertexToComponent[u][v];
		auto& vMolecule = vecVertexToComponent[v][u];
		auto uName = inG[u].name + "_" + std::to_string(uMolecule);
		auto vName = inG[v].name + "_" + std::to_string(vMolecule);
		auto e =
		    boost::add_edge(molSepGBarcodeToIndex[uName], molSepGBarcodeToIndex[vName], molSepG)
		        .first;
		molSepG[e].weight = inG[*edgeIt].weight;
	}

	std::cerr << "Generated new graph ";
#if _OPENMP
	std::cerr << "in sec: " << omp_get_wtime() - sTime << std::endl;
	sTime = omp_get_wtime();
#endif

	std::cerr << "Memory usage: " << double(memory_usage()) / double(1048576) << "GB" << std::endl;
}

void
biconnectedComponents(graph_t& subgraph, vertexToComponent_t& vertexToComponent)
{
	// Find biconnected components
	boost::property_map<graph_t, edgeComponent_t>::type component =
	    boost::get(edgeComponent, subgraph);

	std::vector<vertex_t> artPointsVec;
	boost::biconnected_components(subgraph, component, std::back_inserter(artPointsVec));

	vertexSet_t artPoints(artPointsVec.begin(), artPointsVec.end());

	// Remove articulation points from biconnected components
	boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
	componentToVertexSet_t componentToVertexSet;

	for (boost::tie(ei, ei_end) = boost::edges(subgraph); ei != ei_end; ++ei) {
		size_t componentNum = component[*ei];
		if (componentNum + 1 > componentToVertexSet.size()) {
			componentToVertexSet.resize(componentNum + 1);
		}

		auto node1 = source(*ei, subgraph);
		auto node2 = target(*ei, subgraph);

		if (artPoints.find(node1) == artPoints.end()) {
			componentToVertexSet[componentNum].insert(subgraph[node1].indexOriginal);
		}
		if (artPoints.find(node2) == artPoints.end()) {
			componentToVertexSet[componentNum].insert(subgraph[node2].indexOriginal);
		}
	}

	size_t moleculeNum = 0;

	// Remove components with size less than 1
	for (auto&& vertexSet : componentToVertexSet) {
		if (vertexSet.size() <= 1) {
			continue;
		}
		for (auto&& vertex : vertexSet) {
			vertexToComponent[vertex] = moleculeNum;
		}
		++moleculeNum;
	}
}



template<typename K, typename V>
std::unordered_map<V,K> inverse_map(std::unordered_map<K,V> &map)
{
	std::unordered_map<V,K> inverse;
	for (const auto &p: map) {
		inverse.insert(std::make_pair(p.second, p.first));
	}
	return inverse;
}

vector<vector<double> >
square_matrix_ijk(
    vector<vector<double> > M,
    bool symmetric_output=true)
{
    // Compute M.M^T (not M.M)
    int n = M.size();

    vector<double> tempVector(n, 0.0); // Fast initialization
    vector<vector<double> > M2(n, tempVector);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if ( j < i && symmetric_output )
            {
                M2[i][j] = M2[j][i];
                continue;
             }
            for (int k = 0; k < n; k++)
            {
                // second argument is transposed implicitly
                M2[i][j] += M[i][k] * M[j][k];
            }
        }
    }
    return M2;
}

adjacencyMatrix_t
square_matrix_ijk(
    adjacencyMatrix_t M,
    bool symmetric_output=true)
{
    // Compute M.M^T, not M.M
    int n = M.size();

    adjacencyVector_t tempVector(n, 0); // Fast initialization
    adjacencyMatrix_t M2(n, tempVector);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if ( j < i && symmetric_output )
            {
                M2[i][j] = M2[j][i];
                continue;
            }
            for (int k = 0; k < n; k++)
            {
                // second argument is transposed implicitly
                M2[i][j] += M[i][k] * M[j][k];
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
square_matrix_ikj( // Might be faster than ijk, benchmark it
    vector<vector<double> > M,
    bool symmetric=true)
{
    // Square the input matrix iterating i, k, then j

    // Fast initialization:
    int n = M.size();
    vector<double> tempVector(n, 0.0); // Fast initialization
    vector<vector<double> > M2(n, tempVector);
    // Multiplication
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            if ( !M[i][k] )
                // one side of multip. is zero, so skip
                continue;
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
            row_sum += *first * *first;
            ++first;
        }

        first = row_i->begin();
        vector<double>::iterator first_normalized = normalized_row_i->begin();
        vector<double>::iterator last_normalized = normalized_row_i->end();
        while(first!=last){
            *first_normalized = *first / sqrt(1.0 * row_sum);
            ++first;
            ++first_normalized;
        }
    }
//    time_duration = duration(timeNow() - time_start );
    //total_duration += std::chrono::duration_cast<std::chrono::microseconds>( timeNow() - time_start ).count();
    //time_start = timeNow();
    cosimilarity = square_matrix_ijk(normalized);
//    cosimilarity = square_matrix_ikj(normalized);
//    cosimilarity = normalized;

//    total_duration += std::chrono::duration_cast<std::chrono::microseconds>( timeNow() - time_start ).count();
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
                cosimilarity[i][j] = cosimilarity[j][i];
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

void
community_detection_cosine_similarity(
    graph_t& subgraph, vertexToComponent_t& vertexToComponent,
    bool squaring = true, double threshold=0.7)
{
    // Detect communities using cosine similarity of vertices

    // 1- Calculate the cosine similarity:
    vertexToIndex_t vertexToIndex(num_vertices(subgraph));
    adjacencyMatrix_t adj_mat(convert_adj_list_adj_mat(subgraph, vertexToIndex));
    indexToVertex_t indexToVertex = inverse_map(vertexToIndex);

    size_t size_adj_mat = adj_mat.size();
    vector<double> tempVector(size_adj_mat, 0);
    vector<vector<double>> cosSimilarity2d(size_adj_mat, tempVector);
    calculate_cosine_similarity_2d(squaring ?
                                square_matrix_ikj(adj_mat, true) // may need some change
                                //square_matrix_ijk(adj_mat, true)
                                //square_matrix_boost(adj_mat)
                                :
                                adj_mat,
                        cosSimilarity2d);
    // 2- Determine the threshold:
    // not implemented yet; so use a predefined universal threshold.
    threshold = threshold;

    // 3- Filter out edges:
    for (int i = 0; i < adj_mat.size() ; i++)
    {
        for (int j = i+1; j < adj_mat.size() ; j++)
            {
                if (cosSimilarity2d[i][j] < threshold)
                {
                    adj_mat[i][j] = 0;
                    adj_mat[j][i] = 0;
                }
            }
    }
    // 4- Detect Communities (find connected components - DFS)
    //      Alternative implementation: convert to adjacency list and use boost to find cc

    // / use .reserve to set the capacity of the below 2d vector instead of initialization
    int max_communities = 30;
    vector<vector<uint_fast32_t>> communities(max_communities,vector<uint_fast32_t>(adj_mat.size(),-1));

    size_t community_id = 0;
    stack<int> toCheck;
    //unordered_map<int, int>;
    vector<int> zeros(adj_mat.size(),0);
    vector<int> isDetected(zeros);
    //vector<int> isVisited(zeros);
    bool isSingleton = true;

    for (int i = 0 ; i < adj_mat.size(); i++)
    {
        // DFS traversal
        isVisited = zeros;
        if (isDetected[i])
            continue; // this node is included in a community already.
        toCheck.push(i);
        isDetected[i] = 1;
        isSingleton = true;

        while(!toCheck.empty()){

            ii = toCheck.top();
            toCheck.pop();
            // /communities[community_id].push_back(ii);
            vertex_t vertex = indexToVertex_t.find(ii);
            vertexToComponent.insert (std::pair<vertex_t, size_t>(vertex, community_id));

            for (int j = 0 ; j < adj_mat.size(); j++)
            {
                if (isDetected[j])
                    continue; // this node is included in a community already.
                //if (isVisited[j])
                //    continue; // this node is included in this community already.
                if (adj_mat[ii][j] > 0){
                    toCheck.push(j);
                    isDetected[j]=1;
                    isSingleton = false;
                }
            }
        }
        if (isSingleton)
        {
            // /communities[community_id].pop_back();
            vertexToComponent.erase ( indexToVertex.find(i) );
        }
        else
            community_id++;
    }
}

void
Community_detection_k3_cliques(
    graph_t& subgraph, vertexToComponent_t& vertexToComponent,
    int k = 3)
{
    // k-cliques community detection in case of k=3
    // based on matrix multiplication
    if (k != 3)
    {
        cout<<" This implementation of k-cliques does not support any k other than 3.";
        exit (EXIT_FAILURE);
    }
    vertexToIndex_t vertexToIndex(num_vertices(subgraph));
    adjacencyMatrix_t adj_mat(convert_adj_list_adj_mat(subgraph, vertexToIndex));
    indexToVertex_t indexToVertex = inverse_map(vertexToIndex);

    size_t size_adj_mat = adj_mat.size();
    adjacencyMatrix_t squared_adj_mat(square_matrix_ijk(adj_mat));

    /// TEST WHICH IS FASTER:
    /// 1-MATRIX MULTIPLICATION TO FIND TRIANGLES?
//    int adj_mat_size = adj_mat.size();
//    for (int i = 0; i < adj_mat_size; i++)
//    {
//        for (int j = i+1; j < adj_mat_size; j++)
//        {
//            if ( adj_mat[i][j] > 0 && squared_adj_mat[i][j] > 0 )
//            {
//                // There is a Triangle of 3 vertices i,j, and ?
//
//            }
//        }
//    }
    /// 2-MATRIX TO VECTOR CONVERSION + BITWISE AND ON INTEGERS (compacted vectors)?

    /// 3-NORMAL K-CLIQUE DETECTION
}

int
main(int argc, char* argv[])
{
	//adjacencyMatrix_t a(3, adjacencyVector_t(5,4));
	//adjacencyMatrix_t b(a);
	//cout<<b.size()<<endl;
	//cout<<b[1][1];
    //func(b);


    //adjacencyMatrix_t a_small(3, adjacencyVector_t(3,4));
//    adjacencyMatrix_t a_small {
//				{ 0, 12, 4, 1 },
//				{ 12, 0, 1, 0 },
//				{ 4, 1, 0, 0 },
//				{ 1, 0, 0, 0 }
//			};
//    adjacencyMatrix_t a_small {
//				{ 0, 0, 8, 8 },
//				{ 0, 0, 8, 8 },
//				{ 8, 8, 0, 25 },
//				{ 8, 8, 25, 0 }
//			};
    adjacencyMatrix_t a_small {
				{ 0, 1, 1, 0 },
				{ 1, 0, 1, 1 },
				{ 1, 1, 0, 1 },
				{ 0, 1, 1, 0 }
			};
//    adjacencyMatrix_t a_small {
//				{ 0, 1, 1, 0 },
//				{ 1, 0, 0, 1 },
//				{ 1, 0, 0, 1 },
//				{ 0, 1, 1, 0 }
//			};

	int test_size = 50;
    adjacencyMatrix_t a(test_size, adjacencyVector_t(test_size,4));
    vector<double> tempVector(test_size, 0);
    vector<vector<double> > cosSimilarity2d(test_size, tempVector);


    vector<double> tempVector_small(4, 0.0);
    vector<vector<double> > cosSimilarity2d_small(4, tempVector_small);
    vector<vector<double> > cos1(cosSimilarity2d_small);

//    calculate_cosine_similarity_2d_v2(a_small, cosSimilarity2d_small);
    calculate_cosine_similarity_2d(a_small, cos1);

//    vector<vector<double> >::iterator rowrow1 = cosSimilarity2d_small.begin();
//
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




    auto start = timeNow();

//    uint_fast32_t num1;
//    uint_fast32_t num2;
//    uint_fast32_t res;

    for( int i=0 ; i < (10000*50) ; i++){
        //cout<<i<<endl;

        uint_fast32_t num1 = rand()%32;
        uint_fast32_t num2 = rand()%32;
//        uint_fast32_t res = num1 & num2;
//        if( i%100 != 0 )
//            cout<<res<<endl;

        vector <uint_fast32_t> vec1(50, num1);
        vector <uint_fast32_t> vec2(50, num2);
        int vecsize = vec1.size();
//        for( int i = 0 ; i < vecsize; i++){
//            vec1[i] = (vec1[i] !=0 && vec2[i] !=0 ? 1 : 0);
//        }
        vec1[0] = vec2[0] = 4;
        if( i%100 != 0 )
            cout<<vec1[0]<<endl;

//        a[0][0] = a[0][0] + 1;
//        adjacencyMatrix_t temp(square_matrix_ijk(a));
//        if(temp[0][0]<0){
//            cout<<"Failed";
//            return 0;
//        }
//        calculate_cosine_similarity_2d(a,cosSimilarity2d);
////        calculate_cosine_similarity_2d_v2(a,cosSimilarity2d);
//        if(! test_fast_initialization1() ){
//            cout<<"Failed";
//            return 0;
//        }
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    duration = duration_cast<microseconds>(stop - start);
    cout << "\nElapsed time: "<<duration.count() << endl;
//    cout << "\nElapsed time: "<<duration(timeNow() - start) << endl;
//    cout << "\nElapsed time: "<<total_duration << endl;
	return 0;
}
