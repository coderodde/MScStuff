#include "utils.h"
#include "rodde_current_time.h"
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <lemon/network_simplex.h>
#include <unordered_set>
#include <utility>
#include <fstream>

using std::runtime_error;
using namespace std;
using namespace lemon;
using rodde::current_time::milliseconds;

vector<string> get_vector_of_genome_file_names(string genome_list_file_name)
{
	ifstream genome_list_file(genome_list_file_name);
	vector<string> genome_file_name_vector;
	char line[1024];
	
	while (genome_list_file.good())
	{
		genome_list_file.getline(line, sizeof(line));
		string file_name(line);
		
		if (!file_name.empty())
		{
			genome_file_name_vector.push_back(string(line));			
		}
	}
	
	genome_list_file.close();
	return genome_file_name_vector;
}

string build_kmer(string genome_string, int start_index, int k)
{
	string ret;
	
	for (int i = 0; i < k; ++i)
	{
		ret.push_back(genome_string[(start_index + i) % genome_string.size()]);
	}
	
	return ret;
}

// Passes the test!
void test_build_kmer()
{
	string genome_string = "CGATATAG";
	vector<string> kmers;
	
	for (int start_index = 0; start_index < genome_string.size(); ++start_index)
	{
		kmers.push_back(build_kmer(genome_string, start_index, 3));
	}
	
	for (const auto& s : kmers)
	{
		cout << s << endl;
	}
	
	/*
	 // Lemon allows parallel arcs. :-(
	ListDigraph tmp;
	ListDigraph::Node n1 = tmp.addNode();
	ListDigraph::Node n2 = tmp.addNode();
	tmp.addArc(n1, n2);
	tmp.addArc(n1, n2);
	tmp.addArc(n1, n2);
	cout << "Arcs: " << countArcs(tmp) << endl;
	*/
}

string read_genome_file(const string& file_name)
{
	ifstream file(file_name);
	string genome_string;
	char line[1024];
	
	while (file.good())
	{
		file.getline(line, sizeof(line));
		
		if (line[0] != '>') // Comment lines begin with '>' Ignore them.
		{
			genome_string += line;
		}
	}
	
	file.close();
	return genome_string;
}

vector<string> read_genome_files(vector<string>& file_name_vector)
{
	vector<string> ret;
	
	for (const string& file_name : file_name_vector)
	{
		ret.push_back(read_genome_file(file_name));
	}
	
	return ret;
}

struct graph_result {
	StaticDigraph* p_graph;
	StaticDigraph::NodeMap<string>* p_nodeLabels;
};

graph_result construct_graph_from_genomes(vector<string>& genome_vector, int k)
{
	ListDigraph* list_digraph = new ListDigraph;
	ListDigraph::NodeMap<string>* list_digraph_node_labels = new ListDigraph::NodeMap<string>(*list_digraph);
	
	unordered_map<string, ListDigraph::Node> node_labels_to_list_digraph_nodes;
	
	cout << "SHIT: creating nodes" << endl;
	
	// Create all the nodes:
	for (const string& genome_string : genome_vector)
	{
		cout << "NEW GENOME OF LENGTH " << genome_string.length() << endl;
		
		for (int start_index = 0; start_index < genome_string.length(); ++start_index)
		{
			cout << start_index << endl;
			
			
			string node_label = build_kmer(genome_string, start_index, k);
			
			if (node_labels_to_list_digraph_nodes.find(node_label) ==
			    node_labels_to_list_digraph_nodes.end())
			{
				ListDigraph::Node new_node = list_digraph->addNode();
				node_labels_to_list_digraph_nodes[node_label] = new_node;
				(*list_digraph_node_labels)[new_node] = node_label;
			}
		}
	}
	
	unordered_map<int, unordered_map<int, bool>> arc_filter;
	
	cout << "SHIT: creating arcs" << endl;
	
	// Create all the arcs:
	for (const string& genome_string : genome_vector)
	{
		string previous_label = build_kmer(genome_string, genome_string.length() - 1, k);
		
		for (int start_index = 0; start_index < genome_string.length(); ++start_index)
		{
			string current_label = build_kmer(genome_string, start_index, k);
			
			ListDigraph::Node tail_node = node_labels_to_list_digraph_nodes[previous_label];
			ListDigraph::Node head_node = node_labels_to_list_digraph_nodes[current_label];
			
			int tail_node_id = list_digraph->id(tail_node);
			int head_node_id = list_digraph->id(head_node);
			
			if (!arc_filter[tail_node_id][head_node_id])
			{
				arc_filter[tail_node_id][head_node_id] = true;
				list_digraph->addArc(tail_node, head_node);
			}
			
			previous_label = current_label;
		}
	}
	
	StaticDigraph* output_digraph = new StaticDigraph;
	StaticDigraph::NodeMap<string>* output_digraph_node_labels =
		new StaticDigraph::NodeMap<string>(*output_digraph);
		
	DigraphCopy<ListDigraph, StaticDigraph> copy(*list_digraph, *output_digraph);
	copy.nodeMap(*list_digraph_node_labels, *output_digraph_node_labels);
	copy.run();
	
	graph_result result = { output_digraph, output_digraph_node_labels };
	return result;
}
	
void test_construct_graph_from_genomes()
{
	vector<string> genome_string_vector { "CGATATAG", "AGC" };
	graph_result result = construct_graph_from_genomes(genome_string_vector, 3);
	
	StaticDigraph* graph = result.p_graph;
	StaticDigraph::NodeMap<string>* labels = result.p_nodeLabels;
	
	for (StaticDigraph::NodeIt nodeit(*graph); nodeit != INVALID; ++nodeit)
	{
		cout << graph->id(nodeit) << ": " << (*labels)[nodeit];
		cout << "; incoming labels:";
		
		for (StaticDigraph::InArcIt arcit(*graph, nodeit); arcit != INVALID; ++arcit)
		{
			StaticDigraph::Node source = graph->source(arcit);
			cout << " " << (*labels)[source];
		}
		
		cout << "; outgoing labels:";
		
		for (StaticDigraph::OutArcIt arcit(*graph, nodeit); arcit != INVALID; ++arcit)
		{
			StaticDigraph::Node target = graph->target(arcit);
			cout << " " << (*labels)[target];
		}
		
		cout << endl;
	}
}

int N_THREADS;

//vector<vector<int>> get_node_covering_reconstruction(const StaticDigraph& graph, bool debug_print)
vector<pair<vector<StaticDigraph::Node>,
	    vector<StaticDigraph::Arc>>>
get_node_covering_reconstruction(const StaticDigraph& graph, bool debug_print)
{
	uint64_t start_time = milliseconds();
	
	int self_loop_count = 0;
	
	for (StaticDigraph::ArcIt arcit(graph); arcit != INVALID; ++arcit)
	{
		int tailNodeId = graph.id(graph.source(arcit));
		int headNodeId = graph.id(graph.target(arcit));
		
		if (tailNodeId == headNodeId)
		{
			self_loop_count++;
		}
	}
	
	cout << "[ALEXANDRU] Number of self-loops: " << self_loop_count << endl;
	
	if (debug_print)
	{
		cout << "[ALEXANDRU](get_node_covering_reconstruction) Computing the node covering reconstruction...\n";
	}
	
	unordered_map<int, unordered_map<int, StaticDigraph::Arc>> id_pair_to_static_graph_arc;
	
	for (StaticDigraph::ArcIt arcit(graph); arcit != INVALID; ++arcit)
	{
		StaticDigraph::Node tail = graph.source(arcit);
		StaticDigraph::Node head = graph.target(arcit);
		id_pair_to_static_graph_arc[graph.id(tail)][graph.id(head)] = arcit;
	}
	
	int nodes = graph.nodeNum();
	
	//// Subdivide the input graph:
	ListDigraph subdivided_graph;
	subdivided_graph.reserveNode(2 * nodes);
	subdivided_graph.reserveArc(nodes + countArcs(graph));
	
	unordered_map<int, int> list_node_id_to_static_node_id;
	
	//// Create the nodes of the subdivided graph:
	for (int id = 0; id < nodes; ++id)
	{
		ListDigraph::Node tail = subdivided_graph.addNode();
		ListDigraph::Node head = subdivided_graph.addNode();
	}
	
	if (debug_print)
	{
		cout << "[ALEXANDRU](get_node_covering_reconstruction) The size of the subdivided graph is "
		     << countNodes(subdivided_graph) << endl;	
	}
	
	//// Maps the indices of the tail and head nodes to the actual arc between them:
	unordered_map<int, unordered_map<int, ListDigraph::Arc>> arc_matrix;
	
	// The two maps are used for specifying the demand:
	ListDigraph::ArcMap<int64_t> lowerMap(subdivided_graph);
	ListDigraph::ArcMap<int64_t> upperMap(subdivided_graph);
	
	// This map stores the costs for all arcs in the subdivided graph:
	ListDigraph::ArcMap<int64_t> costMap (subdivided_graph);
	
	// This map maps each node in the subdivided graph to a node in the static graph that it represents:
	ListDigraph::NodeMap<StaticDigraph::Node> list_to_static_graph_node_map(subdivided_graph);
	
	// Create the subdivision arcs for the subdivided graph. We split each
	// node x into two nodes x_in and x_out, where x_in = x, and x_out is
	// a newly added node.
	for (int id = 0; id < nodes; ++id)
	{
		ListDigraph::Node tail = subdivided_graph.nodeFromId(id);
		ListDigraph::Node head = subdivided_graph.nodeFromId(id + nodes);
		ListDigraph::Arc arc   = subdivided_graph.addArc(tail, head);
		// Subdivided arcs: demand 1, cost 0:
		lowerMap[arc] = 1;
		upperMap[arc] = numeric_limits<int64_t>::max();
		costMap [arc] = 0;
		
		// Save the arc for further management. We will add the original
		// arcs into this map in the next for-loop.
		arc_matrix[id][id + nodes] = arc;
		
		// 'tail' and 'head' in the 'subdivided_graph' both correspond to a node with ID 'id' in the static input graph 'graph':
		list_to_static_graph_node_map[tail] = graph.nodeFromId(id);
		list_to_static_graph_node_map[head] = graph.nodeFromId(id);
	}
	
	for (auto p : list_node_id_to_static_node_id)
	{
		cout << p.first << ":" << p.second << endl;
	}
	
	if (debug_print)
	{
		cout << "[ALEXANDRU](get_node_covering_reconstruction) The number of divided arcs before copying the arcs is "
		     << countArcs(subdivided_graph)
		     << "\n";		
	}
	
	// Copy the original arcs to the subdivided graph.
	// For each arc (y, z) in G, add an arc (y_out, z_in) to G'.
	for (StaticDigraph::ArcIt arcit(graph); arcit != INVALID; ++arcit)
	{
		StaticDigraph::Node y = graph.source(arcit);
		StaticDigraph::Node z = graph.target(arcit);
		
		int y_index = graph.index(y) + nodes; // The index of y_out.
		int z_index = graph.index(z);         // The index of z_in.
		
		ListDigraph::Node new_arc_tail = subdivided_graph.nodeFromId(y_index);
		ListDigraph::Node new_arc_head = subdivided_graph.nodeFromId(z_index);
		ListDigraph::Arc new_arc = subdivided_graph.addArc(new_arc_tail, new_arc_head);
		
		// Original graph arcs: demand 0, cost 1:
		lowerMap[new_arc] = 0;
		upperMap[new_arc] = numeric_limits<int64_t>::max();
		costMap [new_arc] = 1;
		
		// Put in the a-matrix:
		arc_matrix[y_index][z_index] = new_arc;
	}
	
	if (debug_print)
	{
		cout << "[ALEXANDRU](get_node_covering_reconstruction) The number of divided arcs after copying the arcs is "
		     << countArcs(subdivided_graph)
		     << endl;		
	}
	
	//// Once here, we have the subdivided graph!
	ListDigraph::NodeMap<int64_t> supplyMap(subdivided_graph);
	
	// In 'supplyMap', map each node to zero:
	for (ListDigraph::NodeIt nodeit(subdivided_graph); nodeit != INVALID; ++nodeit)
	{
		supplyMap[nodeit] = 0;
	}
	
	NetworkSimplex<ListDigraph, int64_t> ns(subdivided_graph);
	ns.lowerMap(lowerMap).upperMap(upperMap).costMap(costMap).supplyMap(supplyMap);
	
	NetworkSimplex<ListDigraph, int64_t>::ProblemType min_flow_result =
	ns.run(NetworkSimplex<ListDigraph, int64_t>::CANDIDATE_LIST);
	
	if (min_flow_result == NetworkSimplex<ListDigraph, int64_t>::OPTIMAL)
	{
		cout << "[ALEXANDRU](get_node_covering_reconstruction) The min-cost flow is optimal!" << endl;
	}
	else
	{
		cerr << "[ALEXANDRU](get_node_covering_reconstruction) The min-cost flow is not found!" << endl;
		abort();
	}
	
	ListDigraph::ArcMap<int64_t> resultFlowMap(subdivided_graph);
	ns.flowMap(resultFlowMap);
	
	//cout << "GRAPH:" << endl;
	//cout << "Nodes: " << endl;
	
	/*
	for (ListDigraph::NodeIt nodeit(subdivided_graph); nodeit != INVALID; ++nodeit)
	{
		cout << subdivided_graph.id(nodeit) << endl;
	}
	cout << "Arcs: " << endl;
	for (ListDigraph::ArcIt arcit(subdivided_graph); arcit != INVALID; ++arcit)
	{
		cout << subdivided_graph.id(subdivided_graph.source(arcit)) << " -> "
		     << subdivided_graph.id(subdivided_graph.target(arcit))
		     << ": " << resultFlowMap[arcit] 
		     << endl;
	}*/
	
	//// Reconstruct the cycles:
	unordered_set<int> arc_id_set_with_nonzero_flows;
	
	for (ListDigraph::ArcIt arcit(subdivided_graph); arcit != INVALID; ++arcit)
	{
		if (resultFlowMap[arcit] > 0)
		{
			arc_id_set_with_nonzero_flows.insert(subdivided_graph.id(arcit));
		}
	}
	
	vector<vector<int>> cycles;
	int count = 0;
	while (true)
	{
		vector<int> cycle;
		/*ListDigraph::ArcIt target_arc;
		bool start_arc_found = false;*/
		
		if (arc_id_set_with_nonzero_flows.empty())
		{
			break;
		}
		
		auto iter = arc_id_set_with_nonzero_flows.begin();
		ListDigraph::Arc target_arc = subdivided_graph.arcFromId(*iter);
		
		// First find any arc with non-zero flow:
		/*for (ListDigraph::ArcIt arcit(subdivided_graph); arcit != INVALID; ++arcit)
		{
			if (resultFlowMap[arcit] > 0)
			{
				target_arc = arcit;
				start_arc_found = true;
				break;
			}
		}*/
		
		
		
		/*if (!start_arc_found)
		{
			// We are done with constructing cycles:
			break;
		}*/
		
		unordered_set<int> filter;
		
		int start_id = subdivided_graph.id(subdivided_graph.source(target_arc));
		ListDigraph::Node current_node = subdivided_graph.target(target_arc);
		int current_node_id = subdivided_graph.id(current_node);
		cycle.push_back(start_id);
		filter.insert(start_id);
		
		// Find the cycle:
		while (filter.find(current_node_id) == filter.end())
		{
			cycle.push_back(current_node_id);
			filter.insert(current_node_id);
			
			// Find next node to visit:
			for (ListDigraph::OutArcIt arcit(subdivided_graph, current_node); arcit != INVALID; ++arcit)
			{
				if (resultFlowMap[arcit] > 0)
				{
					ListDigraph::Node next_node = subdivided_graph.target(arcit);
					current_node = next_node;
					current_node_id = subdivided_graph.id(current_node);
					break;
				}
			}
		}
		
		// Prune the cycle:
		vector<int> pruned_cycle;
		size_t idx = 0;
		
		for (; cycle[idx] != current_node_id; ++idx) {}
		
		for (; idx < cycle.size(); ++idx)
		{
			pruned_cycle.push_back(cycle[idx]);
		}
		
		// Remove one unit of flow from each arc in the currently found cycle:
		
		for (size_t i = 0; i < pruned_cycle.size(); ++i)
		{
			int tail_node_id = pruned_cycle[i];
			int head_node_id = pruned_cycle[(i + 1) % pruned_cycle.size()];
			ListDigraph::Arc arc = arc_matrix[tail_node_id][head_node_id];
			resultFlowMap[arc]--;
			
			if (resultFlowMap[arc] == 0)
			{
				// Remove the arc by ID from the set 'arc_id_set_with_nonzero_flows':
				arc_id_set_with_nonzero_flows.erase(subdivided_graph.id(arc));
			}
		}
		
		cycles.push_back(pruned_cycle);
	}
	
	//// Next, convert the cycles in the subdivided graph into cycles in the input graph:
	vector<vector<int>> graph_cycles;
	
	vector<pair<vector<StaticDigraph::Node>,
	            vector<StaticDigraph::Arc>>> result;
	
	for (vector<int>& cycle : cycles)
	{
		vector<int> graph_cycle;
		
		for (int node_id : cycle)
		{
			if (node_id < nodes)
			{
				ListDigraph::Node tmp_node = subdivided_graph.nodeFromId(node_id);
				StaticDigraph::Node static_tmp_node = list_to_static_graph_node_map[tmp_node];
				graph_cycle.push_back(graph.id(static_tmp_node));
			}
		}
		
		// 1. Simply convert the node IDs to StaticDigraph::Node's.
		vector<StaticDigraph::Node> walk_as_node_vector;
		
		for (int id : graph_cycle)
		{
			walk_as_node_vector.push_back(graph.nodeFromId(id));
		}
		
		// 2. Convert to the list of StaticDigraph::Arc's.
		vector<StaticDigraph::Arc> walk_as_arc_vector;
		
		for (size_t i = 0; i < walk_as_node_vector.size(); ++i)
		{
			StaticDigraph::Node tail = walk_as_node_vector[i];
			StaticDigraph::Node head = walk_as_node_vector[(i + 1) % walk_as_node_vector.size()];
			
			int tail_id = graph.id(tail);
			int head_id = graph.id(head);
			
			walk_as_arc_vector.push_back(id_pair_to_static_graph_arc[tail_id][head_id]);
		}
		
		// Make a pair of the cycle (node version + arc version):
		result.push_back(pair<vector<StaticDigraph::Node>,
				      vector<StaticDigraph::Arc>>(walk_as_node_vector,
								  walk_as_arc_vector));
		//graph_cycles.push_back(graph_cycle);
	}
	
	//cout << "RESULT CYCLES: " << result.size() << endl;
	/*
	for (pair<vector<StaticDigraph::Node>, vector<StaticDigraph::Arc>> walk : result)
	{
		cout << "WALK BEGIN" << endl;
		cout << "Node walk: ";
		
		for (StaticDigraph::Node node : walk.first)
		{
			cout << graph.id(node) << " ";
		}
		
		cout << endl;
		
		cout << "Arc walk:  ";
		
		for (StaticDigraph::Arc arc : walk.second)
		{
			StaticDigraph::Node tail = graph.source(arc);
			StaticDigraph::Node head = graph.target(arc);
			cout << "(" << graph.id(tail) << " -> " << graph.id(head) << ") ";
		}
		
		cout << endl << "WALK END" << endl;
	}*/
	
	
	/*****************************************************************
	* Here we check that all input graph nodes belong to some cycle: *
	*****************************************************************/ 
	unordered_set<int> unvisited_graph_node_id_set;
	
	for (StaticDigraph::NodeIt nodeit(graph); nodeit != INVALID; ++nodeit)
	{
		unvisited_graph_node_id_set.insert(graph.id(nodeit));
	}
	
	for (pair<vector<StaticDigraph::Node>,
	          vector<StaticDigraph::Arc>> p : result)
	{
		for (StaticDigraph::Node node : p.first)
		{
			int id = graph.id(node);
			auto iterator = unvisited_graph_node_id_set.find(id);
			
			if (iterator != unvisited_graph_node_id_set.end())
			{
				unvisited_graph_node_id_set.erase(iterator);
			}
		}
	}
	
	cout << "[SANITY CHECK] Number of unvisited input graph nodes: " << unvisited_graph_node_id_set.size() << endl;
	
	/****************************************************
	* Here we check that all cycles are in fact cycles: *
	****************************************************/ 
	for (pair<vector<StaticDigraph::Node>,
		  vector<StaticDigraph::Arc>> p : result)
	{
		if (p.first.size() != p.second.size())
		{
			cerr << "[ERROR] Bad cycle." << endl;
			exit(1);
		}
		
		for (int i = 0; i < p.first.size(); ++i)
		{
			StaticDigraph::Node tailNode = p.first[i];
			StaticDigraph::Node headNode = p.first[(i + 1) % p.first.size()];
			
			StaticDigraph::Arc arc = p.second[i];
			
			int tailNodeId = graph.id(tailNode);
			int headNodeId = graph.id(headNode);
			
			if (tailNodeId != graph.id(graph.source(arc)))
			{
				cerr << "[ERROR] Cycle mismatch. (1)" << endl;
			}
			
			if (headNodeId != graph.id(graph.target(arc)))
			{
				cerr << "[ERROR] Cycle mismatch. (2)" << endl;
			}
		}
	}
	
	//// id_pair_to_static_graph_arc
	
	uint64_t end_time = milliseconds();
	cout << "[ALEXANDRU](get_node_covering_reconstruction) in "
	     << end_time - start_time << " milliseconds.\n";
	//return graph_cycles;
	return result;
}

// This function returns a circular node-covering walk in the input graph.
// A node-covering walk is a walk that visits each node at least once.
static pair<vector<StaticDigraph::Node>, vector<StaticDigraph::Arc>>
get_circular_walk(const StaticDigraph& graph, bool debug_print)
{
	uint64_t start_time = milliseconds();
	
	if (debug_print)
	{
		cout << "[ALEXANDRU](get_circular_walk) Computing the main circular walk...\n";	
	}
	
	vector<StaticDigraph::Node> main_walk;
	vector<StaticDigraph::Arc> main_walk_arcs;
	
	int nodes = graph.nodeNum();
	
	for (int node_id = 0; node_id < nodes; ++node_id)
	{
		int source_node_id = node_id;
		int target_node_id = (node_id + 1) % nodes;
		
		StaticDigraph::Node source_node = graph.nodeFromId(source_node_id);
		StaticDigraph::Node target_node = graph.nodeFromId(target_node_id);

		Bfs<StaticDigraph> bfs(graph); 
		
		if (!bfs.run(source_node, target_node))
		{
			// We get here only if 'target_node' is not reachable from 'source_node'.
			throw std::runtime_error("The input graph is not strongly connected!");
		}
		
		// Append the nodes to the node-covering walk under construction.
		Path<StaticDigraph> path = bfs.path(target_node);
		
		for (int i = 0; i < path.length(); ++i)
		{
			StaticDigraph::Arc arc = path.nth(i);
			main_walk.push_back(graph.source(arc));
			main_walk_arcs.push_back(arc);
		}
	}
	
	if (debug_print)
	{
		cout << "[ALEXANDRU](get_circular_walk) The length of the main circular walk is: " << main_walk.size() << "\n";	
	}
	
	uint64_t end_time = milliseconds();
	
	cout << "[ALEXANDRU] get_circular_walk() in " << (end_time - start_time) << " milliseconds.\n";
	return pair<vector<StaticDigraph::Node>, vector<StaticDigraph::Arc>>(main_walk, main_walk_arcs);
}

// This function computes a map. The map in question maps each node in the main input graph to
// an unordered (hash) set containing the IDs of the nodes that are the certificates of each mapped node.
static void find_certificate_sets(const StaticDigraph& graph,
				  StaticDigraph::NodeMap<unordered_set<int>>& map_node_to_certificate_set,
			  	  bool debug_print)
{
	uint64_t start_time = milliseconds();
	
	int nodes = graph.nodeNum();
	
	if (debug_print)
	{
		cout << "[ALEXANDRU](find_certificate_sets) Computing the certificate sets!\n";	
	}
	
	for (int id = 0; id < nodes; ++id)
	{
		StaticDigraph::Node node = graph.node(id);
		unordered_set<int> initial_certificate_set = { id };
		map_node_to_certificate_set[node] = initial_certificate_set;
	}
	
	//// Now subdivide the graph.
	ListDigraph subdivided_graph;
	subdivided_graph.reserveNode(2 * nodes);
	subdivided_graph.reserveArc(nodes + countArcs(graph));
	
	// Create the nodes of the subdivided graph.
	for (int id = 0; id < nodes; ++id)
	{
		subdivided_graph.addNode();
		subdivided_graph.addNode();
	}
	
	if (debug_print)
	{
		cout << "[ALEXANDRU](find_certificate_sets) The size of the subdivided graph is "
		     << countNodes(subdivided_graph) << endl;	
	}
	
	// Maps the indices of the tail and head nodes to the actual arc between them:
	unordered_map<int, unordered_map<int, ListDigraph::Arc>> arc_matrix;
	
	// Create the subdivision arcs for the subdivided graph. We split each
	// node x into two nodes x_in and x_out, where x_in = x, and x_out is
	// a newly added node.
	for (int id = 0; id < nodes; ++id)
	{
		ListDigraph::Node tail = subdivided_graph.nodeFromId(id);
		ListDigraph::Node head = subdivided_graph.nodeFromId(id + nodes);
		ListDigraph::Arc arc   = subdivided_graph.addArc(tail, head);
		// Save the arc for further management.
		arc_matrix[id][id + nodes] = arc;
	}
	
	if (debug_print)
	{
		cout << "[ALEXANDRU](find_certificate_sets) The number of divided arcs before copying the arcs is "
		     << countArcs(subdivided_graph)
		     << "\n";		
	}
	
	// Copy the original arcs to the subdivided graph.
	// For each arc (y, z) in G, add an arc (y_out, z_in) to G'.
	for (StaticDigraph::ArcIt arcit(graph); arcit != INVALID; ++arcit)
	{
		StaticDigraph::Node y = graph.source(arcit);
		StaticDigraph::Node z = graph.target(arcit);
		
		int y_index = graph.index(y) + nodes; // The index of y_out.
		int z_index = graph.index(z);         // The index of z_in.
		
		ListDigraph::Node new_arc_tail = subdivided_graph.nodeFromId(y_index);
		ListDigraph::Node new_arc_head = subdivided_graph.nodeFromId(z_index);
		subdivided_graph.addArc(new_arc_tail, new_arc_head);
	}
	
	if (debug_print)
	{
		cout << "[ALEXANDRU] The number of divided arcs after copying the arcs is " << countArcs(subdivided_graph) << endl;		
	}
	
	// Next, for each x in V(G) construct the graph G'_x.
	for (int x_node_id = 0; x_node_id < nodes; ++x_node_id)
	{
		int head_node_id = x_node_id + nodes;
		
		// Remove (x_in, x_out).
		ListDigraph::Node x_in  = subdivided_graph.nodeFromId(x_node_id);
		ListDigraph::Node x_out = subdivided_graph.nodeFromId(head_node_id);
		ListDigraph::Arc removed_arc = arc_matrix[x_node_id][head_node_id];
		subdivided_graph.erase(removed_arc);
		
		// Here, compute the strongly connected components.
		ListDigraph::NodeMap<int> scc(subdivided_graph);
		stronglyConnectedComponents(subdivided_graph, scc);
		
		for (int y_node_id = 0; y_node_id < nodes; ++y_node_id)
		{
			if (y_node_id == x_node_id)
			{
				// Here y is x, omit it.
				continue;
			}
			
			int y_in_node_id  = y_node_id;
			int y_out_node_id = y_node_id + nodes;
			
			ListDigraph::Node y_in  = subdivided_graph.nodeFromId(y_in_node_id);
			ListDigraph::Node y_out = subdivided_graph.nodeFromId(y_out_node_id);
			
			if (scc[y_in] != scc[y_out])
			{
				StaticDigraph::Node static_x_node = graph.node(x_node_id);	
				map_node_to_certificate_set[static_x_node].insert(y_node_id);
			}
		}
		
		// Return (x_in, x_out) to the graph and start the next iteration.
		subdivided_graph.addArc(x_in, x_out);
	}
	
	uint64_t end_time = milliseconds();
	
	cout << "[ALEXANDRU] find_certificate_sets() in " << (end_time - start_time) << " milliseconds.\n";
}


// Computes the m x m matrix of bools. For arcs (x -> y) and (z -> w), the matrix[x -> y][z -> w] tells
// whether there is a path x ---> w, with the first arc different from x -> y and the last arc different from z -> w.
static unordered_map<int, unordered_map<int, bool>> compute_a_matrix(const StaticDigraph& graph,
								     bool debug_print)
{
	uint64_t start_time = milliseconds();
	
	if (debug_print)
	{
		cout << "[ALEXANDRU](compute_a_matrix) Computing the a-matrix in O(m^2)...\n";
	}
	
	unordered_map<int, unordered_map<int, bool>> a_matrix;
	int nodes = graph.nodeNum();
	
	// Create a ListDigraph for manipulating the topology of the graph:
	ListDigraph work_graph;
	
	// Copy the input graph to the ListDigraph created above:
	DigraphCopy<StaticDigraph, ListDigraph> copy_graph(graph, work_graph);
	
	StaticDigraph::ArcMap<ListDigraph::Arc>   map_static_digraph_arcs_to_list_digraph_arcs(graph);
	StaticDigraph::NodeMap<ListDigraph::Node> map_static_digraph_nodes_to_list_digraph_nodes(graph);
	ListDigraph::ArcMap<StaticDigraph::Arc>   map_list_digraph_arcs_to_static_digraph_arcs(work_graph);
	
	copy_graph.arcRef(map_static_digraph_arcs_to_list_digraph_arcs);
	copy_graph.nodeRef(map_static_digraph_nodes_to_list_digraph_nodes);
	copy_graph.arcCrossRef(map_list_digraph_arcs_to_static_digraph_arcs);
	copy_graph.run();
	
	for (StaticDigraph::ArcIt arc(graph); arc != INVALID; ++arc)
	{
		// Get the corresponding arc from the list digraph (x_1, x_2)
		ListDigraph::Arc removed_arc = map_static_digraph_arcs_to_list_digraph_arcs[arc];
		ListDigraph::Node removed_arc_tail = work_graph.source(removed_arc); // x_1
		ListDigraph::Node removed_arc_head = work_graph.target(removed_arc); // x_2
		int removed_static_arc_id = graph.id(arc);
		int removed_arc_id = work_graph.id(removed_arc);
		work_graph.erase(removed_arc);
		
		// Run the DFS in order to find all the nodes reachable from the node x_1:
		Dfs<> dfs(work_graph);
		dfs.run(removed_arc_tail);
		
		// Consider each node z in V(G):
		for (StaticDigraph::NodeIt nodeit(graph); nodeit != INVALID; ++nodeit)
		{
			ListDigraph::Node z = map_static_digraph_nodes_to_list_digraph_nodes[nodeit];
			int d_z = 0;
			
			// Iterate over all in-neighbours of z:
			for (ListDigraph::InArcIt in_arc(work_graph, z); in_arc != INVALID; ++in_arc)
			{
				ListDigraph::Node incoming_node = work_graph.source(in_arc);
				
				if (dfs.reached(incoming_node))
				{
					d_z++;
				}
			}
			
			// Iterate over all in-neighbours of z once again:
			for (ListDigraph::InArcIt in_arc(work_graph, z); in_arc != INVALID; ++in_arc)
			{
				ListDigraph::Node w = work_graph.source(in_arc);
				int r_w = dfs.reached(w) ? 1 : 0;
				StaticDigraph::Arc corresponding_arc = map_list_digraph_arcs_to_static_digraph_arcs[in_arc];
				int corresponding_arc_id = graph.id(corresponding_arc);
				a_matrix[removed_static_arc_id][corresponding_arc_id] = d_z - r_w > 0;
			}
		}
		
		// Return the removed arc back to the work graph:
		work_graph.addArc(removed_arc_tail, removed_arc_head);
	}
	
	uint64_t end_time = milliseconds();
	cout << "[ALEXANDRU] compute_a_matrix() in " << (end_time - start_time) << " milliseconds.\n";
	return a_matrix;
}

// This function returns an unordered (hash) set of arc IDs that are strong
// bridges. Removing a strong bridge from the graph increases the number of
// strongly connected components.
static unordered_set<int> find_strong_bridges(const StaticDigraph& graph, bool debug_print)
{
	uint64_t start_time = milliseconds();
	
	if (debug_print)
	{
		cout << "[ALEXANDRU](find_strong_bridges) Computing strong bridges...\n";
	}
	
	unordered_set<int> ret;

	// Create a graph copy that we will manipulate:
	ListDigraph work_graph;
	DigraphCopy<StaticDigraph, ListDigraph> digraph_copy(graph, work_graph);
	
	// Create a map mapping the nodes from StaticDigraph to ListDigraph:
	StaticDigraph::NodeMap<ListDigraph::Node> input_to_work_graph_node_map(graph);
	digraph_copy.nodeRef(input_to_work_graph_node_map);
	
	// Create a map mapping the arcs of ListDigraph to the arcs of StaticDigraph:
	ListDigraph::ArcMap<StaticDigraph::Arc> work_graph_to_input_arc_map(work_graph);
	digraph_copy.arcCrossRef(work_graph_to_input_arc_map);
	
	digraph_copy.run();
	
	StaticDigraph::ArcMap<ListDigraph::Arc> map_input_graph_arcs_to_work_graph_arcs(graph);
	
	for (ListDigraph::ArcIt arcit(work_graph); arcit != INVALID; ++arcit)
	{
		map_input_graph_arcs_to_work_graph_arcs[work_graph_to_input_arc_map[arcit]] = arcit;
	}

	for (StaticDigraph::ArcIt arcit(graph); arcit != INVALID; ++arcit)
	{
		// Here I need to find an arc in the ListDigraph to which 'arcit' maps:
		ListDigraph::Arc current_arc = map_input_graph_arcs_to_work_graph_arcs[arcit];
		
		work_graph.erase(current_arc);	
		
		ListDigraph::NodeMap<int> scc(work_graph);
		int number_of_strongly_connected_components = stronglyConnectedComponents(work_graph, scc);
	
		if (number_of_strongly_connected_components != 1)
		{
			ret.insert(graph.id(arcit));
		}
		
		work_graph.addArc(work_graph.source(current_arc),
			          work_graph.target(current_arc));
	}
	
	if (debug_print)
	{
		cout << "[ALEXANDRU](find_strong_bridges) Number of strong bridges: " << ret.size() << "\n";
	}
	
	uint64_t end_time = milliseconds();
	
	cout << "[ALEXANDRU] find_strong_bridges() in " << (end_time - start_time) << " milliseconds.\n";
	
	return ret;
}

unordered_map<int, int> compute_funky_ell_indices(const StaticDigraph& graph,
						  vector<StaticDigraph::Node> main_walk,
						  StaticDigraph::NodeMap<unordered_set<int>>& map_node_to_certificate_set)
{
	uint64_t start_time = milliseconds();
	
	// The 'main_walk' consists of d nodes and d arcs, where n <= d <= n^2.
	// Iterate over all nodes of the walk:
	unordered_map<int, int> ell_map;
	
	// This map maps each StaticDigraph::Node ID to the number of times the
	// node appear in the certificate set intersection:
	unordered_map<int, int> counter_map;
	
	// Load the counts of the first certificate set:
	StaticDigraph::Node first_node = main_walk[0];
	unordered_set<int> first_node_certificate_set =
		map_node_to_certificate_set[first_node];
		
	for (const int& id: first_node_certificate_set)
	{
		counter_map[id] = 1;
	}
	
	int ell = 0;
	int index = 0;
	
	//const int stop_index = main_walk.size() + countNodes(graph);
	const int stop_index = main_walk.size() * 2;
	
	while (index < stop_index)
	{
		loop:
		while (index < stop_index)
		{
			ell_map[index++] = ell;
			StaticDigraph::Node node = main_walk[index % main_walk.size()];
			unordered_set<int> node_certificate_set = map_node_to_certificate_set[node];
			
			for (const int id : node_certificate_set)
			{
				counter_map[id]++;
			}
			
			const int minimum_required_count = index - ell + 1;
			
			for (const auto& entry : counter_map)
			{
				if (entry.second >= minimum_required_count)
				{
					goto loop;
				}
			}
			
			break;
		}
		
		while (index < stop_index && ell < index)
		{
			StaticDigraph::Node node = main_walk[ell % main_walk.size()];
			ell++;
			
			unordered_set<int> node_certificate_set =
				map_node_to_certificate_set[node];
				
			for (const int id : node_certificate_set)
			{
				counter_map[id]--;
			}
			
			const int minimum_required_count = index - ell + 1;
		
			for (const auto& entry : counter_map)
			{
				if (entry.second >= minimum_required_count)
				{
					goto loop;
				}
			}
		}
	}
	
	uint64_t end_time = milliseconds();
	
	//cout << "[ALEXANDRU] compute_funky_ell_indices() in " << (end_time - start_time) << " milliseconds.\n";
	
	return ell_map;
}

vector<contig> coderodde_project_algorithm(const StaticDigraph& graph,
					   const StaticDigraph::NodeMap<string>& nodeLabel,
					   const string& inputFileName,
					   const size_t kmersize,
					   const string& sequence,
					   const bool debug_print)
{
	uint64_t start_time = milliseconds();
	vector<contig> ret;
	size_t nodes = graph.nodeNum();
	
	if (debug_print)
	{
		cout << "[ALEXANDRU] Entered 'coderodde_project_algorithm'." << endl;	
		cout << "[ALEXANDRU] The size of the input graph is: " << nodes << endl;
		cout << "[ALEXANDRU] The number of arcs in the graph is: " << graph.arcNum() << endl;
	
		cout << "[ALEXANDRU] k-mer size: " << kmersize << endl;
		cout << "[ALEXANDRU] Sequence length: " << sequence.length() << endl;
		cout << "[ALEXANDRU] Input file name: " << inputFileName << endl;
	}
	
	//vector<pair<vector<StaticDigraph::Node>,
	//            vector<StaticDigraph::Arc>>> cycle_vector = get_node_covering_reconstruction(graph, debug_print);

	vector<pair<vector<StaticDigraph::Node>,
		    vector<StaticDigraph::Arc>>> cycle_vector;
		    
		    cycle_vector.push_back(get_circular_walk(graph, debug_print));

	size_t c = 0;
	cout << "[ALEXANDRU] Number of cycles: " << cycle_vector.size() << endl;
	
	for (pair<vector<StaticDigraph::Node>, vector<StaticDigraph::Arc>> p : cycle_vector)
	{
		c += p.first.size();
	}
	
	cout << "[ALEXANDRU] Average cycle length: " << 1.0 * c / cycle_vector.size() << endl;
	
	    ////////////////////////////////////////////////////
	  //// Computing a node-covering circular walk C. ////
	////////////////////////////////////////////////////
	/*pair<vector<StaticDigraph::Node>,
	     vector<StaticDigraph::Arc>> walk_data = get_circular_walk(graph, debug_print);*/
		  	
	    ////////////////////////////////////////////////////
	  //// Computing node certificate sets! Lemma 5.1 ////
	////////////////////////////////////////////////////
	StaticDigraph::NodeMap<unordered_set<int>> map_node_to_certificate_set(graph);
	find_certificate_sets(graph, map_node_to_certificate_set, debug_print);
	
	    /////////////////////////////////////////
	  //// Compute the a-matrix. Lemma 5.2 ////
	/////////////////////////////////////////
	unordered_map<int, unordered_map<int, bool>> a_matrix = compute_a_matrix(graph, debug_print);
	
	    ///////////////////////////////////////////////////////////////////////////////////
	  //// Computing the strong bridges. Lemma 5.3 with a relaxation O(m) -> O(m^2). ////
	///////////////////////////////////////////////////////////////////////////////////
	unordered_set<int> strong_bridge_id_set = find_strong_bridges(graph, debug_print);
	
	if (debug_print)
	{
		cout << "[ALEXANDRU] Number of strong bridges is " << strong_bridge_id_set.size() << endl;
	}
	
	uint64_t start_time_2 = milliseconds();
	const int n = graph.nodeNum();
		
	struct MyHashVector {
		size_t operator() (const vector<int> &vec) const {
			size_t hash = 0;
			size_t i = 0;
		
			for (int elem : vec)
			{
				hash += ++i * elem;
			}
		
			return hash;
		}
	};
	
	struct MyKeyEqualVector {
		bool operator()(vector<int> const& lhs, vector<int> const& rhs) const
		{
			size_t sz = lhs.size();
			
			if (sz != rhs.size())
			{
				return false;
			}
		
			for (size_t i = 0; i != sz; ++i)
			{
				if (lhs[i] != rhs[i])
				{
					return false;
				}
			}
		
			return true;
		}
	};
	
	unordered_set<vector<int>, MyHashVector, MyKeyEqualVector> filter;
	
	uint64_t main_alg_start_time = milliseconds();
	
	for (pair<vector<StaticDigraph::Node>, vector<StaticDigraph::Arc>>& pair : cycle_vector)
	{
		vector<StaticDigraph::Node> main_walk = pair.first;
		vector<StaticDigraph::Arc>  main_walk_arcs = pair.second;
		
		unordered_map<int, int> ell_map = compute_funky_ell_indices(graph, main_walk, map_node_to_certificate_set);
		
		const int n = main_walk.size();
		const int d = main_walk.size();
		
		vector<unordered_set<int>> S_k(n + 1);
		
		for (int k = 1; k <= n; ++k)
		{
			for (int i = 0; i < d; ++i)
			{
				if (k == 1)
				{
					StaticDigraph::Arc e_i = main_walk_arcs[i];
					unordered_set<int>::const_iterator iter =
						strong_bridge_id_set.find(graph.id(e_i));
					
					if (iter != strong_bridge_id_set.end())
					{
						int end_index = (i + 1) % d;
						int start_index = i;
						
						if (ell_map[end_index] <= start_index)
						{
							S_k[1].insert(i);
						}
					}
				}
				else /* (k > 1) */
				{
					// Checks that i \in S_{k - 1} and i + 1 \mod d \in S_{k - 1}:
					unordered_set<int>::const_iterator iter1 = S_k[k - 1].find(i);
					unordered_set<int>::const_iterator iter2 = S_k[k - 1].find((i + 1) % d);
					
					if (iter1 == S_k[k - 1].end() || iter2 == S_k[k - 1].end())
					{
						continue;
					}
					
					//// Check there is no v_{i + k - 1 mode d} - v_{i + 1 mod d} path with
					//// first edge different than e_{i + k - 1 mod d} and
					//// last edge different than e_{i}:
					
					// Get the edge e_{i + k - 1 mod d}:
					StaticDigraph::Arc first_arc  = main_walk_arcs[(i + k - 1) % d];
					
					// Get the edge e_i
					StaticDigraph::Arc second_arc = main_walk_arcs[i];
					
					if (a_matrix[graph.id(first_arc)][graph.id(second_arc)])
					{
						continue;
					}
					
					// Last check: Cert(v_i) \cap .. \cap Cert(v_{i + k mod d} not empty:
					const int end_index = (i + k) % d;
					const int start_index = i;
					
					if (ell_map[end_index] <= start_index)
					{
						S_k[k].insert(i);
					}
				}
			}
		}
		
		for (int k = 0; k < S_k.size(); ++k)
		{
			for (const auto i : S_k[k])
			{
				vector<int> pre_contig;
				
				for (int j = 0; j <= k; ++j)
				{
					pre_contig.push_back(graph.id(main_walk[(i + j) % main_walk.size()]));
				}
				
				filter.insert(pre_contig);
			}
		}
	}
	
	uint64_t main_alg_end_time = milliseconds();
	cout << "[CODERODDE] The main algorithm duration: " << (main_alg_end_time - main_alg_start_time) << " milliseconds." << endl;
	
	cout << "[ALEXANDRU] FILTER SIZE: " << filter.size() << endl;
	
	size_t count = 0;
	
	for (const vector<int>& vec : filter)
	{
		count += vec.size();
	}
	
	cout << "Average nodes per pre contig: " << 1.0 * count / filter.size() << endl;
	
	for (const vector<int>& pre_contig : filter)
	{
		contig current_contig;
				
		for (int id : pre_contig)
		{
			current_contig.nodes.push_back(id);
		}
				
		ret.push_back(current_contig);
	}
	
	uint64_t end_time = milliseconds();
	
	cout << "[ALEXNADRU] Inner algorithm duration: " << (end_time - start_time_2) << " milliseconds.\n";
	cout << "[ALEXANDRU] Total duration: " << (end_time - start_time) << " milliseconds.\n";
	cout << "[ALEXANDRU] Total number of contigs: " << ret.size() << endl;
	
	uint64_t sum = 0;
	
	for (contig& c : ret)
	{
		sum += c.nodes.size();
	}
	
	cout << "[ALEXANDRU] Average nodes per contig: " << (1.0 * sum) / ret.size() << endl;
	cout << "[ALEXANDRU] Populating strings..." << endl;
	
	populate_with_strings_from_node_labels(sequence, kmersize, graph, nodeLabel, ret);
	
	cout << "[ALEXANDRU] Populated! Printing to a file..." << endl;
	
	print_collection(ret, inputFileName + ".k" + std::to_string(kmersize), ".alexandru_omnitigs");
	
	cout << "[ALEXANDRU] Done!" << endl;
	return ret;
}

vector<contig> coderodde_project_algorithm(const StaticDigraph& graph,
					   const StaticDigraph::NodeMap<string>& nodeLabel,
					   const string& inputFileName,
					   const size_t kmersize,
					   const string& sequence)
{
	return coderodde_project_algorithm(graph,
					   nodeLabel,
					   inputFileName,
					   kmersize,
					   sequence,
					   false);
}

vector<contig> compute_unitigs(const StaticDigraph& graph, 
	const StaticDigraph::NodeMap<size_t>& length, 
	const StaticDigraph::NodeMap<size_t>& seqStart,
	const size_t kmersize,
	const string& sequence,
	const string inputFileName)
{
	vector<contig> ret;
	unordered_set<int> processedArcs;

	// first, merging the arcs incident to unary nodes
	for (StaticDigraph::NodeIt node(graph); node != INVALID; ++node)
	{
		if ((countInArcs(graph,node) == 1) and (countOutArcs(graph,node) == 1))
		{
			StaticDigraph::InArcIt in_arc(graph,node);
			StaticDigraph::OutArcIt out_arc(graph,node);

			contig entry;
			entry.nodes.push_back(graph.id(graph.source(in_arc)));
			entry.nodes.push_back(graph.id(node));
			entry.nodes.push_back(graph.id(graph.target(out_arc)));

			ret.push_back(entry);

			processedArcs.insert(graph.id(in_arc));
			processedArcs.insert(graph.id(out_arc));
		}
	}

	for (StaticDigraph::ArcIt a(graph); a != INVALID; ++a)
	{
		if (processedArcs.count(graph.id(a)) == 0)	
		{
			
			contig entry;
			entry.nodes.push_back(graph.id(graph.source(a)));
			entry.nodes.push_back(graph.id(graph.target(a)));
			ret.push_back(entry);	

		}
	}

	return ret;
}

vector<contig> compute_unitigs_by_traversal(const ListDigraph& graph)
{
	vector<contig> ret;

	cout << "Entered compute_unitigs_by_traversal" << endl;

	ListDigraph::Node current_node;

	// iterating over all arcs
	for (ListDigraph::ArcIt arc(graph); arc != INVALID; ++arc)	
	{
		ListDigraph::Node u;
		u = graph.source(arc);

		if ((countInArcs (graph,u) > 1) or (countOutArcs(graph,u) > 1))
		{
			// traversing the outtig (forward)
			contig entry;

			entry.nodes.push_back(graph.id(u));

			//cout << "creating a nice contig form arc " << graph.id(graph.source(arc)) << "," << graph.id(graph.target(arc)) << endl;
			current_node = graph.target(arc);
			while ((countOutArcs(graph, current_node) == 1) and (countInArcs(graph, current_node) == 1))
			{
				entry.nodes.push_back(graph.id(current_node));
				ListDigraph::OutArcIt out_arc(graph,current_node);
				if (current_node == graph.target(out_arc))
				{
					entry.nodes.push_back(graph.id(current_node));
					break;
				}
				current_node = graph.target(out_arc);
				//cout << "seen " << graph.id(current_node) << " forward" << endl;
			}
			entry.nodes.push_back(graph.id(current_node));

			ret.push_back(entry);
		}		
	}

	return ret;
}

vector<contig> compute_non_switching_contigs(const StaticDigraph& graph, 
	const StaticDigraph::NodeMap<size_t>& length, 
	const StaticDigraph::NodeMap<size_t>& seqStart,
	const size_t kmersize,
	const string& sequence,
	const string inputFileName)
{
	vector<contig> ret;

	cout << "Entered non_switching_contigs" << endl;

	StaticDigraph::Node current_node;
	string s;

	// iterating over all arcs
	for (StaticDigraph::ArcIt arc(graph); arc != INVALID; ++arc)	
	{
		StaticDigraph::Node u,v,z;
		u = graph.source(arc);
		v = graph.target(arc);
		z = INVALID;

		// OLD: if (((countOutArcs(graph,graph.source(arc)) > 1) and (countInArcs(graph,graph.target(arc)) > 1)) or 
		//  	((countOutArcs(graph,graph.source(arc)) == 1) and (countInArcs(graph,graph.source(arc)) == 1)))
		// NEW: as in the paper
		if (countInArcs(graph,u) == 1)
		{
			StaticDigraph::InArcIt in_arc(graph,u);
			z = graph.source(in_arc);
		}

		if ( ((countOutArcs(graph,u) >= 2) and (countInArcs(graph,v) >= 2)) or
			 ((z != INVALID) and (countOutArcs(graph,u) == 1) and 
			 	( ((countOutArcs(graph,z) >= 2) and (countInArcs (graph,v) >= 2)) or 
			 	  ((countInArcs (graph,z) >= 2) and (countOutArcs(graph,v) >= 2)) )
			 )
		   )
		{
			// traversing the outtig (forward)
			contig entry;

			//cout << "creating a nice contig form arc " << graph.id(graph.source(arc)) << "," << graph.id(graph.target(arc)) << endl;
			current_node = graph.target(arc);
			while (countOutArcs(graph, current_node) == 1)
			{
				entry.nodes.push_back(graph.id(current_node));
				StaticDigraph::OutArcIt out_arc(graph,current_node);
				if (current_node == graph.target(out_arc))
				{
					entry.nodes.push_back(graph.id(current_node));
					break;
				}
				current_node = graph.target(out_arc);
				//cout << "seen " << graph.id(current_node) << " forward" << endl;
			}
			entry.nodes.push_back(graph.id(current_node));

			// traversing the intig (backward)
			current_node = graph.source(arc);
			while (countInArcs(graph, current_node) == 1)
			{
				entry.nodes.push_front(graph.id(current_node));
				StaticDigraph::InArcIt in_arc(graph,current_node);
				if (current_node == graph.source(in_arc))
				{
					entry.nodes.push_front(graph.id(current_node));
					break;
				}
				current_node = graph.source(in_arc);
				//cout << "seen " << graph.id(current_node) << " backward" << endl;
			}
			entry.nodes.push_front(graph.id(current_node));

			ret.push_back(entry);
		}		
	}

	return ret;
}

vector<contig> compute_YtoV_contigs_and_reduce(StaticDigraph& graph_static, 
	StaticDigraph::NodeMap<size_t>& length_static, 
	StaticDigraph::NodeMap<size_t>& seqStart_static,
	StaticDigraph::NodeMap<string>& nodeLabel_static,
	const size_t kmersize,
	string& sequence,
	const string inputFileName,
	string genome_type,
	size_t seqLength
	)
{

	vector<contig> YtoV_contigs;
	cout << "Entered YtoV_contigs" << endl;

	// copying graph_static into a ListDigraph graph
	ListDigraph graph;
	ListDigraph::NodeMap<size_t> length(graph);
	ListDigraph::NodeMap<size_t> seqStart(graph);
	ListDigraph::NodeMap<string> nodeLabel(graph);

	DigraphCopy<StaticDigraph, ListDigraph> copy_graph(graph_static, graph);
	ListDigraph::NodeMap<StaticDigraph::Node> graph_to_graph_static_node(graph);
	// copy_graph.nodeCrossRef(graph_to_graph_static_node);
	copy_graph.nodeMap(length_static, length);
	copy_graph.nodeMap(seqStart_static, seqStart);
	copy_graph.nodeMap(nodeLabel_static, nodeLabel);
	copy_graph.run();
	// finished copying into graph

	vector<ListDigraph::Node> forwardYNodes, backwardYNodes;
	vector<ListDigraph::Arc> arcsToContract;
	size_t n_nodes = countNodes(graph);
	size_t n_arcs = countArcs(graph);
	forwardYNodes.reserve(n_nodes);
	backwardYNodes.reserve(n_nodes);
	arcsToContract.reserve(n_arcs);
	bool applied_reduction = true;
	size_t progress = 0;

	while (applied_reduction)
	{
		applied_reduction = false;

		///////////////////////////////////
		// now we try the forward Ynodes
		///////////////////////////////////
		for (ListDigraph::NodeIt node(graph); node != INVALID; ++node)
		{
			// forwardYNodes
			if ((countInArcs(graph,node) == 1) and (countOutArcs(graph,node) > 1))
			{
				if (node != graph.source( ListDigraph::InArcIt(graph,node) ) )
				{
					forwardYNodes.push_back(node);
					applied_reduction = true;	
					// break;				
				}
			}
		}
		// cout << "Found " << forwardYNodes.size() << " forwardYNodes" << endl;

		// we reduce the forward nodes
		for (int i = 0; i < forwardYNodes.size(); i++)
		{
			ListDigraph::Node node = forwardYNodes[i];

			if (not graph.valid(node))
			{
				cout << "node not valid" << endl;
			}

			ListDigraph::InArcIt in_arc(graph,node);
			ListDigraph::Node node_predecessor = graph.source(in_arc);

			ListDigraph::OutArcIt out_arc(graph,node);
			++out_arc;

			while (out_arc != INVALID)
			{				
				ListDigraph::OutArcIt next_out_arc = out_arc;
				++next_out_arc;

				ListDigraph::Node node_prime = graph.addNode();
				length[node_prime] = length[node];
				seqStart[node_prime] = seqStart[node];
				nodeLabel[node_prime] = nodeLabel[node];

				graph.addArc(node_predecessor, node_prime);
				graph.changeSource(out_arc, node_prime);

				out_arc = next_out_arc;
			}
			//graph.erase(node);
		}	

		///////////////////////////////////
		// now we try the backward Ynodes
		///////////////////////////////////
		for (ListDigraph::NodeIt node(graph); node != INVALID; ++node)
		{
			// backwardYNodes
			if ((countInArcs(graph,node) > 1) and (countOutArcs(graph,node) == 1))
			{
				if (node != graph.target( ListDigraph::OutArcIt(graph,node) ) )
				{
					backwardYNodes.push_back(node);
					applied_reduction = true;
					// break;	
				}
			}
		}
		// cout << "Found " << backwardYNodes.size() << " backwardYNodes" << endl;

		// we reduce the backward nodes
		progress = 1;
		for (int i = 0; i < backwardYNodes.size(); i++)
		{
			ListDigraph::Node node = backwardYNodes[i];

			if (not graph.valid(node))
			{
				cout << "node not valid" << endl;
			}

			ListDigraph::OutArcIt out_arc(graph,node);
			ListDigraph::Node node_successor = graph.target(out_arc);

			ListDigraph::InArcIt in_arc(graph,node);
			++in_arc;

			while (in_arc != INVALID)
			{
				ListDigraph::InArcIt next_in_arc = in_arc;
				++next_in_arc;

				ListDigraph::Node node_prime = graph.addNode();
				length[node_prime] = length[node];
				seqStart[node_prime] = seqStart[node];
				nodeLabel[node_prime] = nodeLabel[node];
				
				graph.changeTarget(in_arc, node_prime);
				graph.addArc(node_prime, node_successor);

				in_arc = next_in_arc;
			}
			//graph.erase(node);
		}	


		///////////////////////////////////
		// now we try the contract arcs
		///////////////////////////////////
		for (ListDigraph::ArcIt a(graph); a != INVALID; ++a)
		{
			if (progress % 100000 == 0)
			{
				cout << "Compacting phase, arc : #" << progress << "/" << n_arcs << " ";
			  	cout << "Time: " << currentDateTime();
			}
			progress++;

			ListDigraph::Node u,v;
			u = graph.source(a);
			v = graph.target(a);

			if ((countOutArcs(graph,u) == 1) and (countInArcs(graph,v) == 1))
			// if the arc (u,v) is the 'middle' arc of a unitig, then contract it
			// if ((countOutArcs(graph,u) == 1) and (countInArcs(graph,u) == 1) and
			// 	(countOutArcs(graph,v) == 1) and (countInArcs(graph,v) == 1))
			{
				arcsToContract.push_back(a);
				applied_reduction = true;
			}
		}

		for (auto arc : arcsToContract)
		{
			if (progress % 100000 == 0)
			{
				cout << "Compacting phase 2, arc : #" << progress << "/" << arcsToContract.size() << " ";
			  	cout << "Time: " << currentDateTime();
			}
			progress++;

			if (not graph.valid(arc))
			{
				cout << "arc not valid" << endl;
			}

			ListDigraph::Node u,v;
			u = graph.source(arc);
			v = graph.target(arc);
			nodeLabel[u] = plus_strings(nodeLabel[u],nodeLabel[v],kmersize);

			//graph.contract(u, v);
			for (ListDigraph::OutArcIt out_arc(graph,v); out_arc != INVALID;)
			{
				ListDigraph::OutArcIt next_out_arc = out_arc;
				++next_out_arc;
				graph.changeSource(out_arc,u);
				out_arc = next_out_arc;
			}
			graph.erase(v);

		}


		arcsToContract.clear();
		forwardYNodes.clear();
		backwardYNodes.clear();	
	}

	cout << "Resulting graph has " << countNodes(graph) << " nodes and " << countArcs(graph) << " arcs" << endl;
	cout << "Graph has " << countStronglyConnectedComponents(graph) << " strongly connected components" << endl;
	cout << "Applied all YtoV transformations. Reporting the unitigs of the resulting graph." << endl;

	// saving graph to graph_static
	ListDigraph::NodeMap<StaticDigraph::Node> nodeRef(graph);
	ListDigraph::ArcMap<StaticDigraph::Arc> arcRef(graph);
	graph_static.clear();
	graph_static.build(graph, nodeRef, arcRef);
	// copy maps
	for (ListDigraph::NodeIt node(graph); node != INVALID; ++node)
	{
		length_static[nodeRef[node]] = length[node];
		seqStart_static[nodeRef[node]] = seqStart[node];
		nodeLabel_static[nodeRef[node]] = nodeLabel[node];
	}

	YtoV_contigs = compute_unitigs(graph_static, 
		length_static, 
		seqStart_static, 
		kmersize,
		sequence,
		inputFileName);

	return YtoV_contigs;

}


vector<contig> compute_omnitigs_2(StaticDigraph& graph, 
	StaticDigraph::NodeMap<size_t>& length, 
	StaticDigraph::NodeMap<size_t>& seqStart,
	set_of_pairs& safe_pairs,
	const size_t kmersize,
	const string& sequence,
	const string inputFileName)
{
	vector<contig> ret;
	size_t n_nodes = countNodes(graph);
	ret.reserve(n_nodes * 100);

	cout << "Entered compute_omnitigs_2" << endl;


	StaticDigraph g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16;
	StaticDigraph::NodeMap<StaticDigraph::Node> nm1(graph),nm2(graph),nm3(graph),nm4(graph),nm5(graph),nm6(graph),nm7(graph),nm8(graph),nm9(graph),nm10(graph),nm11(graph),nm12(graph),nm13(graph),nm14(graph),nm15(graph),nm16(graph);


	switch (N_THREADS) {
		case 16: digraphCopy(graph,g16).nodeRef(nm16).run();
		case 15: digraphCopy(graph,g15).nodeRef(nm15).run();
		case 14: digraphCopy(graph,g14).nodeRef(nm14).run();
		case 13: digraphCopy(graph,g13).nodeRef(nm13).run();
		case 12: digraphCopy(graph,g12).nodeRef(nm12).run();
		case 11: digraphCopy(graph,g11).nodeRef(nm11).run();
		case 10: digraphCopy(graph,g10).nodeRef(nm10).run();
		case 9: digraphCopy(graph,g9).nodeRef(nm9).run();
		case 8: digraphCopy(graph,g8).nodeRef(nm8).run();
		case 7: digraphCopy(graph,g7).nodeRef(nm7).run();
		case 6: digraphCopy(graph,g6).nodeRef(nm6).run();
		case 5: digraphCopy(graph,g5).nodeRef(nm5).run();
		case 4: digraphCopy(graph,g4).nodeRef(nm4).run();
		case 3: digraphCopy(graph,g3).nodeRef(nm3).run();
		case 2: digraphCopy(graph,g2).nodeRef(nm2).run();
		case 1: digraphCopy(graph,g1).nodeRef(nm1).run();
	}

	cout << "Finished making copies" << endl;

	size_t progress = 0;

	// for (StaticDigraph::NodeIt node(graph); node != INVALID; ++node)
	// #pragma omp parallel for
	omp_set_dynamic(0);
	#pragma omp parallel for num_threads(N_THREADS) shared(ret,safe_pairs)
	for (uint64_t i = 0; i < n_nodes; i++)
	{
		StaticDigraph::Node node = graph.node(i);

		//************** progress ************** //
		#pragma omp critical
		{
			progress = progress + 1;
		}

		if (progress % 1000 == 0)
		{
			cout << "Node #" << progress << "/" << n_nodes << " ";
		  	cout << "Time: " << currentDateTime();
		  	cout.flush();
		}
		//************** progress ************** //

		// if one-in, then all right extensions are safe
		if (countInArcs(graph,node) == 1)
		{
			StaticDigraph::InArcIt in_arc(graph,node);
			for (StaticDigraph::OutArcIt out_arc(graph,node); out_arc != INVALID; ++out_arc)
			{
				if (out_arc != in_arc)
				{
					contig entry;
					entry.nodes.push_back(graph.id(graph.source(in_arc)));
					entry.nodes.push_back(graph.id(node));
					entry.nodes.push_back(graph.id(graph.target(out_arc)));

					#pragma omp critical
					{
						ret.push_back(entry);
						safe_pairs.insert(pair_of_ints(graph.id(in_arc),graph.id(out_arc)));	
					}	
				}
			}	
		} else 
		// if one-out, then all left extensions are safe
		if (countOutArcs(graph,node) == 1)
		{
			StaticDigraph::OutArcIt out_arc(graph,node);
			for (StaticDigraph::InArcIt in_arc(graph,node); in_arc != INVALID; ++in_arc)
			{
				if (in_arc != out_arc)
				{
					contig entry;
					entry.nodes.push_back(graph.id(graph.source(in_arc)));
					entry.nodes.push_back(graph.id(node));
					entry.nodes.push_back(graph.id(graph.target(out_arc)));

					#pragma omp critical
					{
						ret.push_back(entry);
						safe_pairs.insert(pair_of_ints(graph.id(in_arc),graph.id(out_arc)));	
					}
				}				
			}				
		} else
		{
			for (StaticDigraph::InArcIt in_arc(graph,node); in_arc != INVALID; ++in_arc)
			{
				for (StaticDigraph::OutArcIt out_arc(graph,node); out_arc != INVALID; ++out_arc)
				{
					if (out_arc == in_arc)
					{
						continue;
					}
					// checking whether we the walk in_arc,out_arc forms an edge centric contig
					int tid = omp_get_thread_num();
					tid++;
					//cout << tid << endl;
					StaticDigraph::NodeMap<bool> filterNodesMap(GET_GRAPH(tid), true);
					FilterNodes<StaticDigraph> subgraph(GET_GRAPH(tid), filterNodesMap);

					// StaticDigraph::NodeMap<bool> filterNodesMap(graph, true);
					// FilterNodes<StaticDigraph> subgraph(graph, filterNodesMap);
					
					// check whether some in_nbr is reachable in G - node from some out-neighbor of node different than target(out_arc)
					subgraph.disable((GET_NODE_MAP(tid))[node]);
					Bfs<FilterNodes<StaticDigraph>> visit(subgraph);
					visit.init();

					// adding sources of the visit
					for (StaticDigraph::OutArcIt out_arc2(graph,node); out_arc2 != INVALID; ++out_arc2)
					{
						if (graph.target(out_arc2) != graph.target(out_arc))
						{
							visit.addSource((GET_NODE_MAP(tid))[graph.target(out_arc2)]);
						}
					}
					assert(visit.queueSize() > 0);

					vector<StaticDigraph::Node> bad_guys;
					// adding the other in_nbrs as the other targets of visit
					FilterNodes<StaticDigraph>::NodeMap<bool> visitTargets(subgraph,false);
					// size_t n_bad = 0;
					for (StaticDigraph::InArcIt in_arc2(graph,node); in_arc2 != INVALID; ++in_arc2)
					{
						if (in_arc2 != in_arc)
						{
							visitTargets[(GET_NODE_MAP(tid))[graph.source(in_arc2)]] = true;
							bad_guys.push_back((GET_NODE_MAP(tid))[graph.source(in_arc2)]);
						}
					}

					// if no target vertex is reached, add these arcs as ec_contig_2 and safe pair
					if (visit.start(visitTargets) == INVALID)
					{
						// check that indeed all bad guys are not reached
						// LEMON seems to have a bug and the above check is not sufficient
						bool no_bad_inbr_reached = true;
						for (size_t i = 0; i < bad_guys.size(); i++)
						{
							if (visit.reached(bad_guys[i]))
							{
								no_bad_inbr_reached = false;
								break;
							}
						}
						if (no_bad_inbr_reached)
						{
							contig entry;
							entry.nodes.push_back(graph.id(graph.source(in_arc)));
							entry.nodes.push_back(graph.id(node));
							entry.nodes.push_back(graph.id(graph.target(out_arc)));

							#pragma omp critical
							{
								ret.push_back(entry);
								safe_pairs.insert(pair_of_ints(graph.id(in_arc),graph.id(out_arc)));
							}
							
						}
					}
				}
			}
			
		}
	}

	return ret;
}


void extend_pair(StaticDigraph& graph,
	StaticDigraph& thread_graph,
	StaticDigraph::NodeMap<StaticDigraph::Node>& originalNode2ThreadNode,
	const set_of_pairs& safe_pairs,
	StaticDigraph::Arc in_arc,
	list<int> current_walk,
	unordered_set<int> internal_nodes_set,
	vector<int> bad_in_nbrs,
	vector<contig>& collection_of_contigs,
	unordered_set<int>& reported_arcs
	)
{
	//cout << "Entered extend_pair " << endl;

	StaticDigraph::Node node = graph.target(in_arc);
	// iterate over all out-neighbors
	bool is_maximal = true;

	// checking first whether graph.target(in_arc) is repeated in the contig
	//if ( internal_nodes_set.count( graph.id(graph.target(in_arc)) ) == 0 )
	if ( true )
	{
		for (StaticDigraph::OutArcIt out_arc(graph,node); out_arc != INVALID; ++out_arc)
		{
			// if the pair (in_arc,out_arc) is safe
			if (safe_pairs.count(pair_of_ints(graph.id(in_arc),graph.id(out_arc))) > 0)
			{
				// check whether from 'node', which becomes an internal vertex, there is no path back to 
				// the in-nbrs of the internal vertices of current_walk
				StaticDigraph::NodeMap<bool> filterNodesMap(thread_graph,true);
				FilterNodes<StaticDigraph> subgraph(thread_graph, filterNodesMap);
				Bfs<FilterNodes<StaticDigraph>> visit(subgraph);
				visit.init();

				// adding the sources
				for (StaticDigraph::OutArcIt out_arc2(graph,node); out_arc2 != INVALID; ++out_arc2)
				{
					if (out_arc2 != out_arc)
					{
						visit.addSource(originalNode2ThreadNode[graph.target(out_arc2)]);
					}
				}

				// adding the targets of the search as the bad_in_nbrs
				FilterNodes<StaticDigraph>::NodeMap<bool> visitTargets(subgraph,false);
				// size_t n_added = 0;
				vector<StaticDigraph::Node> bad_guys;
				for (size_t i = 0; i < bad_in_nbrs.size(); i++)
				{
					visitTargets[originalNode2ThreadNode[graph.node(bad_in_nbrs[i])]] = true;
					bad_guys.push_back(originalNode2ThreadNode[graph.node(bad_in_nbrs[i])]);
				}

				subgraph.disable(originalNode2ThreadNode[node]);
				// if no bad_in_nbr is reached
				if ((visit.queueSize() == 0) or (visit.start(visitTargets) == INVALID))
				{
					// check that indeed all bad guys are not reached
					// LEMON seems to have a bug and the above check is not sufficient
					bool no_bad_inbr_reached = true;
					for (size_t i = 0; i < bad_guys.size(); i++)
					{
						if (visit.reached(bad_guys[i]))
						{
							no_bad_inbr_reached = false;
							break;
						}
					}

					// if yes, then add the target of out_arc to current_walk
					if (no_bad_inbr_reached)
					{
						list<int> new_walk = current_walk;
						unordered_set<int> new_internal_nodes_set = internal_nodes_set;

						new_walk.push_back(graph.id(graph.target(out_arc)));
						new_internal_nodes_set.insert(graph.id(graph.target(out_arc)));

						// add the in-nbrs of node not on current_walk to bad_in_nbrs
						vector<int> new_bad_in_nbrs = bad_in_nbrs;
						for (StaticDigraph::InArcIt in_arc2(graph,node); in_arc2 != INVALID; ++in_arc2)
						{
							if (in_arc2 != in_arc)
							{
								new_bad_in_nbrs.push_back(graph.id(graph.source(in_arc2)));
							}
						}
						// and recurse
						is_maximal = false;
						#pragma omp critical
						{
							reported_arcs.insert(graph.id(out_arc));
						}

						// if in_arc and out_arc are arcs between the same two vertices, but in opposite directions
						// then we must stop the walk here
						// because there is no other way to stop it
						if (graph.source(in_arc) == graph.target(out_arc))
						{
							contig entry;
							entry.nodes = new_walk;
							#pragma omp critical
							{
								collection_of_contigs.push_back(entry);
							}

						} 
						else
						{
							extend_pair(graph, thread_graph, originalNode2ThreadNode, safe_pairs, out_arc, new_walk, new_internal_nodes_set, new_bad_in_nbrs, collection_of_contigs, reported_arcs);		
						}
					}
				}
			}
		}

	}	
	
	// if we didn't manage to extend the contig, then it is right-maximal
	// and we add it to the collection of contigs
	// WARNING: it may not be left-maximal
	if (is_maximal)
	{
		contig entry;
		entry.nodes = current_walk;
		#pragma omp critical
		{
			collection_of_contigs.push_back(entry);
		}
	}

}


vector<contig> compute_omnitigs(StaticDigraph& graph, 
	StaticDigraph::NodeMap<size_t>& length, 
	StaticDigraph::NodeMap<size_t>& seqStart,
	const set_of_pairs& safe_pairs
	)
{
	cout << "Entered omnitigs" << endl;
	vector<contig> ret;
	unordered_set<int> reported_arcs;
	size_t progress = 0;

	vector<pair_of_ints> safe_pairs_v(safe_pairs.begin(),safe_pairs.end());
	int n_safe_pairs = safe_pairs_v.size();

	ret.reserve(n_safe_pairs * 100);

	StaticDigraph g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16;
	StaticDigraph::NodeMap<StaticDigraph::Node> nm1(graph),nm2(graph),nm3(graph),nm4(graph),nm5(graph),nm6(graph),nm7(graph),nm8(graph),nm9(graph),nm10(graph),nm11(graph),nm12(graph),nm13(graph),nm14(graph),nm15(graph),nm16(graph);

	switch (N_THREADS) {
		case 16: digraphCopy(graph,g16).nodeRef(nm16).run();
		case 15: digraphCopy(graph,g15).nodeRef(nm15).run();
		case 14: digraphCopy(graph,g14).nodeRef(nm14).run();
		case 13: digraphCopy(graph,g13).nodeRef(nm13).run();
		case 12: digraphCopy(graph,g12).nodeRef(nm12).run();
		case 11: digraphCopy(graph,g11).nodeRef(nm11).run();
		case 10: digraphCopy(graph,g10).nodeRef(nm10).run();
		case 9: digraphCopy(graph,g9).nodeRef(nm9).run();
		case 8: digraphCopy(graph,g8).nodeRef(nm8).run();
		case 7: digraphCopy(graph,g7).nodeRef(nm7).run();
		case 6: digraphCopy(graph,g6).nodeRef(nm6).run();
		case 5: digraphCopy(graph,g5).nodeRef(nm5).run();
		case 4: digraphCopy(graph,g4).nodeRef(nm4).run();
		case 3: digraphCopy(graph,g3).nodeRef(nm3).run();
		case 2: digraphCopy(graph,g2).nodeRef(nm2).run();
		case 1: digraphCopy(graph,g1).nodeRef(nm1).run();
	}	

	cout << "Finished making copies" << endl;

	//for (set_of_pairs::const_iterator itr = safe_pairs.begin(); itr != safe_pairs.end(); ++itr) 
	omp_set_dynamic(0);
	#pragma omp parallel for num_threads(N_THREADS) shared(ret)
	for (int i = 0; i < n_safe_pairs; i++)
	{

		//************** progress ************** //
		#pragma omp critical
		{
			progress = progress + 1;
		}

		if (progress % 1000 == 0)
		{
			cout << "Safe pair #" << progress << " / " << n_safe_pairs << " ";
		  	cout << "Time: " << currentDateTime();
		  	cout.flush();
		}
		//************** progress ************** //

		StaticDigraph::Arc first_arc = graph.arc(safe_pairs_v[i].first);
		StaticDigraph::Arc second_arc = graph.arc(safe_pairs_v[i].second);
		
		list<int> current_walk;
		current_walk.push_back(graph.id(graph.source(first_arc)));
		current_walk.push_back(graph.id(graph.source(second_arc)));
		current_walk.push_back(graph.id(graph.target(second_arc)));

		unordered_set<int> internal_nodes_set;
		internal_nodes_set.insert(graph.id(graph.source(second_arc)));

		// cout << graph.id(graph.source(first_arc)) << ", " << graph.id(graph.source(second_arc)) << ", " << graph.id(graph.target(second_arc)) << endl;

		#pragma omp critical
		{
			reported_arcs.insert(safe_pairs_v[i].first);
			reported_arcs.insert(safe_pairs_v[i].second);	
		}
		
		// problematic case of a 2-cycle a,b,a
		// if so, insert a,b,a as contig and continue
		if (graph.source(first_arc) == graph.target(second_arc))
		{
			contig entry;
			entry.nodes = current_walk;
			#pragma omp critical
			{
				ret.push_back(entry);
			}
			continue;
		}


		vector<int> bad_in_nbrs;
		for (StaticDigraph::InArcIt in_arc2(graph,graph.target(first_arc)); in_arc2 != INVALID; ++in_arc2)
		{
			if (in_arc2 != first_arc)
			{
				bad_in_nbrs.push_back(graph.id(graph.source(in_arc2)));
			}
 		}

		int tid = omp_get_thread_num();
		tid++;
    	extend_pair(graph,GET_GRAPH(tid),GET_NODE_MAP(tid), safe_pairs, second_arc, current_walk, internal_nodes_set, bad_in_nbrs, ret, reported_arcs);
	}

	// including also the arcs not reported in some ec-contig
	for (StaticDigraph::ArcIt a(graph); a != INVALID; ++a)
	{
		if (reported_arcs.count(graph.id(a)) == 0)	
		{
			contig entry;
			entry.nodes.push_back(graph.id(graph.source(a)));
			entry.nodes.push_back(graph.id(graph.target(a)));
			ret.push_back(entry);				
		}
	}

	return ret;
}

static void test_cycle_reconstruction()
{
	ListDigraph my_graph;
	
	ListDigraph::Node a = my_graph.addNode();
	ListDigraph::Node b = my_graph.addNode();
	ListDigraph::Node c = my_graph.addNode();
	ListDigraph::Node d = my_graph.addNode();
	
	ListDigraph::Arc ab = my_graph.addArc(a, b);
	ListDigraph::Arc bc = my_graph.addArc(b, c);
	ListDigraph::Arc ca = my_graph.addArc(c, a);
	ListDigraph::Arc cd = my_graph.addArc(c, d);
	ListDigraph::Arc dc = my_graph.addArc(d, c);
	
	ListDigraph::NodeMap<StaticDigraph::Node> graph_nodes_to_static_graph_nodes(my_graph);
	ListDigraph::ArcMap<StaticDigraph::Arc>   graph_arcs_to_static_graph_arcs(my_graph);
	
	StaticDigraph static_graph;
	static_graph.build(my_graph, graph_nodes_to_static_graph_nodes, graph_arcs_to_static_graph_arcs);
	
	cout << "a: " << static_graph.id(graph_nodes_to_static_graph_nodes[a]) << endl;
	cout << "b: " << static_graph.id(graph_nodes_to_static_graph_nodes[b]) << endl;
	cout << "c: " << static_graph.id(graph_nodes_to_static_graph_nodes[c]) << endl;
	cout << "d: " << static_graph.id(graph_nodes_to_static_graph_nodes[d]) << endl;
	
	/*
	vector<vector<int>> result = get_node_covering_reconstruction(static_graph, true);
	
	for (const vector<int>& cycle : result)
	{
		for (int id : cycle)
		{
			cout << id << " ";
		}
		cout << endl;
	}*/
	
	vector<pair<vector<StaticDigraph::Node>,
		    vector<StaticDigraph::Arc>>> result = get_node_covering_reconstruction(static_graph, true);
		    
	for (pair<vector<StaticDigraph::Node>, vector<StaticDigraph::Arc>>& pair : result)
	{
		for (StaticDigraph::Node& node : pair.first)
		{
			cout << static_graph.id(node) << " ";
		}
		
		cout << endl;
	}
}

static void test_cycle_reconstruction_2()
{
	ListDigraph my_graph;
	
	ListDigraph::Node a = my_graph.addNode();
	ListDigraph::Node b = my_graph.addNode();
	ListDigraph::Node c = my_graph.addNode();
	ListDigraph::Node d = my_graph.addNode();
	ListDigraph::Node e = my_graph.addNode();
	
	ListDigraph::Arc ab = my_graph.addArc(a, b);
	ListDigraph::Arc bc = my_graph.addArc(b, c);
	ListDigraph::Arc ca = my_graph.addArc(c, a);
	ListDigraph::Arc cd = my_graph.addArc(c, d);
	ListDigraph::Arc de = my_graph.addArc(d, e);
	ListDigraph::Arc ec = my_graph.addArc(e, c);
	
	ListDigraph::NodeMap<StaticDigraph::Node> graph_nodes_to_static_graph_nodes(my_graph);
	ListDigraph::ArcMap <StaticDigraph::Arc>  graph_arcs_to_static_graph_arcs(my_graph);
	
	StaticDigraph static_graph;
	static_graph.build(my_graph, graph_nodes_to_static_graph_nodes, graph_arcs_to_static_graph_arcs);
	
	cout << "a: " << static_graph.id(graph_nodes_to_static_graph_nodes[a]) << endl;
	cout << "b: " << static_graph.id(graph_nodes_to_static_graph_nodes[b]) << endl;
	cout << "c: " << static_graph.id(graph_nodes_to_static_graph_nodes[c]) << endl;
	cout << "d: " << static_graph.id(graph_nodes_to_static_graph_nodes[d]) << endl;
	cout << "e: " << static_graph.id(graph_nodes_to_static_graph_nodes[e]) << endl;
	
	vector<pair<vector<StaticDigraph::Node>,
		    vector<StaticDigraph::Arc>>> result = get_node_covering_reconstruction(static_graph, true);
		    
	for (pair<vector<StaticDigraph::Node>, vector<StaticDigraph::Arc>>& pair : result)
	{
		for (StaticDigraph::Node& node : pair.first)
		{
			cout << static_graph.id(node) << " ";
		}
		
		cout << endl;
	}
}

static void test_strong_bridges()
{
	ListDigraph my_graph;
	
	ListDigraph::Node a = my_graph.addNode();
	ListDigraph::Node b = my_graph.addNode();
	ListDigraph::Node c = my_graph.addNode();
	ListDigraph::Node d = my_graph.addNode();
	ListDigraph::Node e = my_graph.addNode();
	
	// Strong bridges:
	ListDigraph::Arc ab = my_graph.addArc(a, b);
	ListDigraph::Arc bd = my_graph.addArc(b, d);
	ListDigraph::Arc dc = my_graph.addArc(d, c);
	ListDigraph::Arc ca = my_graph.addArc(c, a);
	ListDigraph::Arc be = my_graph.addArc(b, e);
	ListDigraph::Arc ed = my_graph.addArc(e, d);
	
	// Other arcs:
	ListDigraph::Arc db = my_graph.addArc(d, b);
	ListDigraph::Arc ac = my_graph.addArc(a, c);
	
	ListDigraph::NodeMap<StaticDigraph::Node> graph_nodes_to_static_graph_nodes(my_graph);
	ListDigraph::ArcMap<StaticDigraph::Arc>   graph_arcs_to_static_graph_arcs(my_graph);
	
	StaticDigraph static_graph;
	static_graph.build(my_graph, graph_nodes_to_static_graph_nodes, graph_arcs_to_static_graph_arcs);
	
	StaticDigraph::Arc bridge_ab = graph_arcs_to_static_graph_arcs[ab];
	StaticDigraph::Arc bridge_bd = graph_arcs_to_static_graph_arcs[bd];
	StaticDigraph::Arc bridge_dc = graph_arcs_to_static_graph_arcs[dc];
	StaticDigraph::Arc bridge_ca = graph_arcs_to_static_graph_arcs[ca];
	StaticDigraph::Arc bridge_be = graph_arcs_to_static_graph_arcs[be];
	StaticDigraph::Arc bridge_ed = graph_arcs_to_static_graph_arcs[ed];
	
	cout << "Static graph nodes: " << countNodes(static_graph) << "\n";
	cout << "Static graph arcs:  " << countArcs(static_graph)  << "\n";
	
	cout << "a -> b: " << static_graph.id(bridge_ab) << "\n";
	cout << "d -> c: " << static_graph.id(bridge_dc) << "\n";
	cout << "c -> a: " << static_graph.id(bridge_ca) << "\n";
	cout << "b -> e: " << static_graph.id(bridge_be) << "\n";
	cout << "e -> d: " << static_graph.id(bridge_ed) << "\n";
	
	unordered_set<int> strong_bridge_id_set = find_strong_bridges(const_cast<const StaticDigraph&>(static_graph), false);
	
	cout << "Number of strong bridges: " << strong_bridge_id_set.size() << "\n";
	cout << "Contents (IDs): ";
	
	for (const auto id : strong_bridge_id_set)
	{
		cout << id << " ";
	}
	
	cout << "\n";
}

static void test_certificate_preprocessing()
{
	// rodde: this is my test I verified with Alexandru:
	ListDigraph graph;
	
	ListDigraph::Node l1 = graph.addNode();
	ListDigraph::Node l2 = graph.addNode();
	ListDigraph::Node m1 = graph.addNode();
	ListDigraph::Node m2 = graph.addNode();
	ListDigraph::Node m3 = graph.addNode();
	ListDigraph::Node r1 = graph.addNode();
	ListDigraph::Node r2 = graph.addNode();
	
	ListDigraph::Arc l1l2 = graph.addArc(l1, l2);
	ListDigraph::Arc l2m3 = graph.addArc(l2, m3);
	ListDigraph::Arc m3m2 = graph.addArc(m3, m2);
	ListDigraph::Arc m2m1 = graph.addArc(m2, m1);
	ListDigraph::Arc m1l1 = graph.addArc(m1, l1);
	ListDigraph::Arc m1r1 = graph.addArc(m1, r1);
	ListDigraph::Arc r1r2 = graph.addArc(r1, r2);
	ListDigraph::Arc r2m3 = graph.addArc(r2, m3);
	
	StaticDigraph static_graph;
	
	ListDigraph::NodeMap<StaticDigraph::Node> graph_nodes_to_static_graph_nodes(graph);
	ListDigraph::ArcMap<StaticDigraph::Arc>   graph_arcs_to_static_graph_arcs(graph);
	
	static_graph.build(graph, graph_nodes_to_static_graph_nodes, graph_arcs_to_static_graph_arcs);
	StaticDigraph::NodeMap<unordered_set<int>> map_node_to_certificate_set(static_graph);
	
	find_certificate_sets(static_graph, map_node_to_certificate_set, false);
	
	cout << "List of StaticDigraph node IDs:\n";
	
	cout << "l1: " << static_graph.id(graph_nodes_to_static_graph_nodes[l1]) << "\n";
	cout << "l2: " << static_graph.id(graph_nodes_to_static_graph_nodes[l2]) << "\n";
	cout << "m1: " << static_graph.id(graph_nodes_to_static_graph_nodes[m1]) << "\n";
	cout << "m2: " << static_graph.id(graph_nodes_to_static_graph_nodes[m2]) << "\n";
	cout << "m3: " << static_graph.id(graph_nodes_to_static_graph_nodes[m3]) << "\n";
	cout << "r1: " << static_graph.id(graph_nodes_to_static_graph_nodes[r1]) << "\n";
	cout << "r2: " << static_graph.id(graph_nodes_to_static_graph_nodes[r2]) << "\n";
	
	// Certificate sets of m1, m2, and m3:
	unordered_set<int> cert_set_of_m1 = map_node_to_certificate_set[graph_nodes_to_static_graph_nodes[m1]];
	unordered_set<int> cert_set_of_m2 = map_node_to_certificate_set[graph_nodes_to_static_graph_nodes[m2]];
	unordered_set<int> cert_set_of_m3 = map_node_to_certificate_set[graph_nodes_to_static_graph_nodes[m3]];
	
	cout << "First test passed: " << std::boolalpha
				      << (cert_set_of_m1 == cert_set_of_m2 &&
					  cert_set_of_m2 == cert_set_of_m3)
				      << "\n";
	
	cout << "The node ID's are:\n";
	
	for (auto i : cert_set_of_m1)
	{
		cout << i << " ";
	}
	
	cout << "\n";
	
	// The following test setting is provided by Alex. Thanky you, Alexandru! :-)
	// Test (a):
	graph.clear();
	
	ListDigraph::Node a = graph.addNode();
	ListDigraph::Node b = graph.addNode();
	ListDigraph::Node c = graph.addNode();
	ListDigraph::Node d = graph.addNode();
	ListDigraph::Node e = graph.addNode();
	
	ListDigraph::Arc ab = graph.addArc(a, b);
	ListDigraph::Arc bc = graph.addArc(b, c);
	ListDigraph::Arc cd = graph.addArc(c, d);
	ListDigraph::Arc de = graph.addArc(d, e);
	ListDigraph::Arc ea = graph.addArc(e, a);
	ListDigraph::Arc ca = graph.addArc(c, a);
	ListDigraph::Arc db = graph.addArc(d, b);
	
	StaticDigraph static_graph_2;
	
	ListDigraph::NodeMap<StaticDigraph::Node> graph_nodes_to_static_graph_nodes_2(graph);
	ListDigraph::ArcMap<StaticDigraph::Arc>   graph_arcs_to_static_graph_arcs_2(graph);
	
	static_graph_2.build(graph, graph_nodes_to_static_graph_nodes_2, graph_arcs_to_static_graph_arcs_2);
	StaticDigraph::NodeMap<unordered_set<int>> map_node_to_certificate_set_2(static_graph_2);
	
	find_certificate_sets(static_graph_2, map_node_to_certificate_set_2, false);
	
	cout << "List of StaticDigraph node IDs in Test (a):\n";
	
	cout << "a: " << static_graph_2.id(graph_nodes_to_static_graph_nodes_2[a]) << "\n";
	cout << "b: " << static_graph_2.id(graph_nodes_to_static_graph_nodes_2[b]) << "\n";
	cout << "c: " << static_graph_2.id(graph_nodes_to_static_graph_nodes_2[c]) << "\n";
	cout << "d: " << static_graph_2.id(graph_nodes_to_static_graph_nodes_2[d]) << "\n";
	cout << "e: " << static_graph_2.id(graph_nodes_to_static_graph_nodes_2[e]) << "\n";
	
	cout << "Cert(a) = {";
	
	for (auto i : map_node_to_certificate_set_2[graph_nodes_to_static_graph_nodes_2[a]])
	{
		cout << " " << i;
	}
	
	cout << " }\n";
	cout << "Cert(b) = {";
	
	for (auto i : map_node_to_certificate_set_2[graph_nodes_to_static_graph_nodes_2[b]])
	{
		cout << " " << i;
	}
	
	cout << " }\n";
	cout << "Cert(c) = {";
	
	for (auto i : map_node_to_certificate_set_2[graph_nodes_to_static_graph_nodes_2[c]])
	{
		cout << " " << i;
	}
	
	cout << " }\n";
	cout << "Cert(d) = {";
	
	for (auto i : map_node_to_certificate_set_2[graph_nodes_to_static_graph_nodes_2[d]])
	{
		cout << " " << i;
	}
	
	cout << " }\n";
	
	graph.erase(e);
	
	StaticDigraph static_graph_3;
	
	ListDigraph::NodeMap<StaticDigraph::Node> graph_nodes_to_static_graph_nodes_3(graph);
	ListDigraph::ArcMap<StaticDigraph::Arc>   graph_arcs_to_static_graph_arcs_3(graph);
	
	static_graph_3.build(graph, graph_nodes_to_static_graph_nodes_3, graph_arcs_to_static_graph_arcs_3);
	StaticDigraph::NodeMap<unordered_set<int>> map_node_to_certificate_set_3(static_graph_3);
	
	find_certificate_sets(static_graph_3, map_node_to_certificate_set_3, false);
	
	cout << "List of StaticDigraph node IDs in Test (b):\n";
	
	cout << "a: " << static_graph_3.id(graph_nodes_to_static_graph_nodes_3[a]) << "\n";
	cout << "b: " << static_graph_3.id(graph_nodes_to_static_graph_nodes_3[b]) << "\n";
	cout << "c: " << static_graph_3.id(graph_nodes_to_static_graph_nodes_3[c]) << "\n";
	cout << "d: " << static_graph_3.id(graph_nodes_to_static_graph_nodes_3[d]) << "\n";
	
	cout << "Cert(a) = {";
	
	for (auto i : map_node_to_certificate_set_3[graph_nodes_to_static_graph_nodes_3[a]])
	{
		cout << " " << i;
	}
	
	cout << " }\n";
	cout << "Cert(b) = {";
	
	for (auto i : map_node_to_certificate_set_3[graph_nodes_to_static_graph_nodes_3[b]])
	{
		cout << " " << i;
	}
	
	cout << " }\n";
	cout << "Cert(c) = {";
	
	for (auto i : map_node_to_certificate_set_3[graph_nodes_to_static_graph_nodes_3[c]])
	{
		cout << " " << i;
	}
	
	cout << " }\n";
	cout << "Cert(d) = {";
	
	for (auto i : map_node_to_certificate_set_3[graph_nodes_to_static_graph_nodes_3[d]])
	{
		cout << " " << i;
	}
	
	cout << " }\n";
	
	{		
		graph.clear();
	
		ListDigraph::Node a = graph.addNode();
		ListDigraph::Node b = graph.addNode();
		ListDigraph::Node c = graph.addNode();
		ListDigraph::Node d = graph.addNode();
		ListDigraph::Node e = graph.addNode();
	
		ListDigraph::Arc ab = graph.addArc(a, b);
		ListDigraph::Arc bc = graph.addArc(b, c);
		ListDigraph::Arc cd = graph.addArc(c, d);
		ListDigraph::Arc de = graph.addArc(d, e);
		ListDigraph::Arc ea = graph.addArc(e, a);
		ListDigraph::Arc cb = graph.addArc(c, b);
		ListDigraph::Arc db = graph.addArc(d, b);
		ListDigraph::Arc ca = graph.addArc(c, a);
		
		StaticDigraph static_graph;
	
		ListDigraph::NodeMap<StaticDigraph::Node> graph_nodes_to_static_graph_nodes(graph);
		ListDigraph::ArcMap<StaticDigraph::Arc>   graph_arcs_to_static_graph_arcs(graph);
	
		static_graph.build(graph, graph_nodes_to_static_graph_nodes, graph_arcs_to_static_graph_arcs);
		StaticDigraph::NodeMap<unordered_set<int>> map_node_to_certificate_set(static_graph);
		
		find_certificate_sets(static_graph, map_node_to_certificate_set, false);
		
		cout << "List of StaticDigraph node IDs in Test (c):\n";
	
		cout << "a: " << static_graph.id(graph_nodes_to_static_graph_nodes[a]) << "\n";
		cout << "b: " << static_graph.id(graph_nodes_to_static_graph_nodes[b]) << "\n";
		cout << "c: " << static_graph.id(graph_nodes_to_static_graph_nodes[c]) << "\n";
		cout << "d: " << static_graph.id(graph_nodes_to_static_graph_nodes[d]) << "\n";
		cout << "e: " << static_graph.id(graph_nodes_to_static_graph_nodes[e]) << "\n";
	
		cout << "Cert(a) = {";
	
		for (auto i : map_node_to_certificate_set[graph_nodes_to_static_graph_nodes[a]])
		{
			cout << " " << i;
		}
	
		cout << " }\n";
		cout << "Cert(b) = {";
		
		for (auto i : map_node_to_certificate_set[graph_nodes_to_static_graph_nodes[b]])
		{
			cout << " " << i;
		}
	
		cout << " }\n";
		cout << "Cert(c) = {";
		
		for (auto i : map_node_to_certificate_set[graph_nodes_to_static_graph_nodes[c]])
		{
			cout << " " << i;
		}
	
		cout << " }\n";
		cout << "Cert(d) = {";
		
		for (auto i : map_node_to_certificate_set[graph_nodes_to_static_graph_nodes[d]])
		{
			cout << " " << i;
		}
	
		cout << " }\n";
		cout << "Cert(e) = {";
		
		for (auto i : map_node_to_certificate_set[graph_nodes_to_static_graph_nodes[d]])
		{
			cout << " " << i;
		}
	
		cout << " }\n";	
	}
}

static void test_list_digraph_node_ids()
{
	ListDigraph graph;
	graph.addNode();
	graph.addNode();
	graph.addNode();
	graph.addNode();
	
	for (ListDigraph::NodeIt node(graph); node != INVALID; ++node)
	{
		cout << graph.id(node) << " ";
	}
	
	cout << "\n";
}

static void crash()
{
	cout << "CRASH: A test failed!\n";
	std::exit(1);
}

static void test_a_matrix_algo()
{
	ListDigraph list_graph;
	
	ListDigraph::Node x = list_graph.addNode();
	ListDigraph::Node y = list_graph.addNode();
	ListDigraph::Node z = list_graph.addNode();
	ListDigraph::Node w = list_graph.addNode();
	
	ListDigraph::Arc a = list_graph.addArc(x, y);
	ListDigraph::Arc b = list_graph.addArc(y, z);
	ListDigraph::Arc c = list_graph.addArc(x, w);
	ListDigraph::Arc d = list_graph.addArc(w, z);
	
	ListDigraph::NodeMap<StaticDigraph::Node> graph_nodes_to_static_graph_nodes(list_graph);
	ListDigraph::ArcMap<StaticDigraph::Arc>   graph_arcs_to_static_graph_arcs(list_graph);
	
	StaticDigraph static_graph;
	static_graph.build(list_graph, graph_nodes_to_static_graph_nodes, graph_arcs_to_static_graph_arcs);
	unordered_map<int, unordered_map<int, bool>> a_matrix = compute_a_matrix(static_graph, false);
	
	StaticDigraph::Arc sa = graph_arcs_to_static_graph_arcs[a];
	StaticDigraph::Arc sb = graph_arcs_to_static_graph_arcs[b];
	StaticDigraph::Arc sc = graph_arcs_to_static_graph_arcs[c];
	StaticDigraph::Arc sd = graph_arcs_to_static_graph_arcs[d];
	
	int isa = static_graph.id(sa);
	int isb = static_graph.id(sb);
	int isc = static_graph.id(sc);
	int isd = static_graph.id(sd);
	
	
	
	cout << "     (x,a) (x,b) (x,c) (x,d)\n";
	cout << "(a,x) " << a_matrix[isa][isa] << "     " << a_matrix[isa][isb] << "     " << a_matrix[isa][isc] << "     " << a_matrix[isa][isd] << "\n";
	cout << "(b,x) " << a_matrix[isb][isa] << "     " << a_matrix[isb][isb] << "     " << a_matrix[isb][isc] << "     " << a_matrix[isb][isd] << "\n";
	cout << "(c,x) " << a_matrix[isc][isa] << "     " << a_matrix[isc][isb] << "     " << a_matrix[isc][isc] << "     " << a_matrix[isc][isd] << "\n";
	cout << "(d,x) " << a_matrix[isd][isa] << "     " << a_matrix[isd][isb] << "     " << a_matrix[isd][isc] << "     " << a_matrix[isd][isd] << "\n";
	
	cout << "HELLO: " << a_matrix[isa].size() << " and " << a_matrix.size() << "\n";

	{
		ListDigraph list_graph;
		
		ListDigraph::Node x1 = list_graph.addNode();
		ListDigraph::Node y1 = list_graph.addNode();
		ListDigraph::Node x2 = list_graph.addNode();
		ListDigraph::Node y2 = list_graph.addNode();
		ListDigraph::Node z1 = list_graph.addNode();
		ListDigraph::Node z2 = list_graph.addNode();
		
		ListDigraph::Arc x1y1 = list_graph.addArc(x1, y1);
		ListDigraph::Arc y1x2 = list_graph.addArc(y1, x2);
		ListDigraph::Arc x2y2 = list_graph.addArc(x2, y2);
		ListDigraph::Arc x1z1 = list_graph.addArc(x1, z1);
		ListDigraph::Arc z1z2 = list_graph.addArc(z1, z2);
		ListDigraph::Arc z2x2 = list_graph.addArc(z2, x2);
		ListDigraph::Arc z2y2 = list_graph.addArc(z2, y2);
		
		ListDigraph::NodeMap<StaticDigraph::Node> graph_nodes_to_static_graph_nodes(list_graph);
		ListDigraph::ArcMap<StaticDigraph::Arc>   graph_arcs_to_static_graph_arcs(list_graph);
	
		StaticDigraph static_graph;
		static_graph.build(list_graph, graph_nodes_to_static_graph_nodes, graph_arcs_to_static_graph_arcs);
		unordered_map<int, unordered_map<int, bool>> a_matrix = compute_a_matrix(static_graph, false);
	
		StaticDigraph::Arc sx1y1 = graph_arcs_to_static_graph_arcs[x1y1];
		StaticDigraph::Arc sy1x2 = graph_arcs_to_static_graph_arcs[y1x2];
		StaticDigraph::Arc sx2y2 = graph_arcs_to_static_graph_arcs[x2y2];
		StaticDigraph::Arc sx1z1 = graph_arcs_to_static_graph_arcs[x1z1];
		StaticDigraph::Arc sz1z2 = graph_arcs_to_static_graph_arcs[z1z2];
		StaticDigraph::Arc sz2x2 = graph_arcs_to_static_graph_arcs[z2x2];
		StaticDigraph::Arc sz2y2 = graph_arcs_to_static_graph_arcs[z2y2];
		
		int isx1y1 = static_graph.id(sx1y1); 
		int isy1x2 = static_graph.id(sy1x2); 
		int isx2y2 = static_graph.id(sx2y2);
		int isx1z1 = static_graph.id(sx1z1); 
		int isz1z2 = static_graph.id(sz1z2); 
		int isz2x2 = static_graph.id(sz2x2);
		int isz2y2 = static_graph.id(sz2y2); 
		
		cout << "(x1 -> y1): " << isx1y1 << "\n";
		cout << "(y1 -> x2): " << isy1x2 << "\n";
		cout << "(x2 -> y2): " << isx2y2 << "\n";
		cout << "(x1 -> z1): " << isx1z1 << "\n";
		cout << "(z1 -> z2): " << isz1z2 << "\n";
		cout << "(z2 -> x2): " << isz2x2 << "\n";
		cout << "(z2 -> y2): " << isz2y2 << "\n";
		
		cout << "HWA\n";
		for (const auto& n : a_matrix)
		{
			cout << "size: " << n.second.size() << ": ";
			for (const auto& nn : n.second)
			{
				cout << "(" << n.first << ", " << nn.first << ")=" << nn.second << " ";
			}
			
			cout << "\n";
		}
		
		if (a_matrix[isx1y1][isx2y2] == false) crash();
		if (a_matrix[isx1y1][isy1x2] == false) crash();
		if (a_matrix[isy1x2][isx2y2] == true) crash();
		if (a_matrix[isx2y2][isx1y1] == true) crash();
		if (a_matrix[isx1y1][isy1x2] == false) crash();
		if (a_matrix[isy1x2][isx2y2] == true) crash();
		if (a_matrix[isx1z1][isz2y2] == false) crash();
		if (a_matrix[isx1z1][isz1z2] == true) crash();
	}
}

void coderodde_processing(string& genome_list_file_name, size_t kmersize)
{
	vector<string> genome_file_name_vector = get_vector_of_genome_file_names(genome_list_file_name);
	cout << "[INFO] Number of genome files in the list: " << genome_file_name_vector.size() << endl;

	/*
	for (string& s : genome_file_name_vector)
	{
		cout << s << endl;
	}
	
	cout << kmersize << endl;*/
	
	cout << "[INFO] Reading the genome files..." << endl;
	vector<string> genome_string_vector = read_genome_files(genome_file_name_vector);
	cout << "[INFO] Constructing the graph..." << endl;
	graph_result g_r = construct_graph_from_genomes(genome_string_vector, kmersize);
	
	StaticDigraph& graph = *g_r.p_graph;
	StaticDigraph::NodeMap<string>& nodeLabels = *g_r.p_nodeLabels;
	
	cout << "[INFO] The result graph has " << countNodes(graph) << " and " << countArcs(graph)  << endl;
}


int main(int argc, char **argv)
{
	//test_certificate_preprocessing();
	//test_list_digraph_node_ids();
	//test_strong_bridges();
	//test_a_matrix_algo();
	//test_cycle_reconstruction();
	//test_cycle_reconstruction_2();
	//test_construct_graph_from_genomes();
	//test_build_kmer();
	
	/*
	vector<string> genome_file_name_vector =
	get_vector_of_genome_file_names("TARGET_LIST");
	cout << "Total files: " << genome_file_name_vector.size() << endl;
	
	for (const auto& s : genome_file_name_vector)
	{
		cout << s << endl;
	}
	
	exit(0);*/
	
	//////////////////////////////////////////
	//////////////////////////////////////////
	//////////////////////////////////////////
	//////////////////////////////////////////
	StaticDigraph sd;
	vector<string> vs = {string("hello"), string("and"), string("world")};
	load_data(vs, sd, string("fds"));
	exit(0);

	StaticDigraph graph;
	StaticDigraph::NodeMap<size_t> length(graph);
	StaticDigraph::NodeMap<size_t> seqStart(graph);
	StaticDigraph::NodeMap<string> nodeLabel(graph);
	set_of_pairs safe_pairs(1 << 17,hash_pair);

	string sequence;
	size_t seqLength;
	size_t kmersize;
	string inputFileName;
    string genome_type = "linear";
    bool circular_genome = false;
	bool build_only = false;
	bool input_from_reads = false;
	bool do_not_contract_arcs;
	bool do_not_compute_omnitigs;
	int abundance = 1;
	N_THREADS = 1;
	struct timespec start_clock, finish_clock;


	// command line argument parser
	string usage = "\n  %prog OPTIONS"
		"\n\nBrief example:"
		"\n  %prog -i <input file> -k 31 [-a 2] [-t 4]";
	const string version = "%prog 0.2\nCopyright (C) 2014-2015 Alexandru Tomescu & Paul Medvedev\n"
		"License GPLv3+: GNU GPL version 3 or later "
		"<http://gnu.org/licenses/gpl.html>.\n"
		"This is free software: you are free to change and redistribute it.\n"
		"There is NO WARRANTY, to the extent permitted by law.";
	const string desc = "";
	const string epilog = "";
	
	optparse::OptionParser parser = optparse::OptionParser()
    	.usage(usage)
    	.version(version)
    	.description(desc)
    	.epilog(epilog);

	parser.add_option("-i", "--input") .type("string") .dest("i") .set_default("") .help("input file (.fasta/.fa/.fastq)");
	parser.add_option("-k", "--kmersize") .type("int") .dest("k") .action("store") .set_default(31) .help("kmer size (default: %default)");
	parser.add_option("-a", "--abundance") .type("int") .dest("a") .action("store") .set_default(1) .help("minimum abundance (default: %default)");
	parser.add_option("-t", "--threads") .type("int") .dest("t") .action("store") .set_default(1) .help("number of threads, in [1..16] (default: %default)");
	parser.add_option("-g", "--genome-type") .action("store") .dest("g") .type("string") .set_default("circular") .help("genome type: linear|circular (default: %default)");
	parser.add_option("-b", "--build-only") .action("store_true") .dest("build_only") .help("build the de bruijn graph and then exit");
	parser.add_option("-c", "--nocontract") .action("store_true") .set_default(false) .dest("nocontract") .help("do not contract arcs");
	parser.add_option("-x", "--noomnitigs") .action("store_true") .set_default(false) .dest("noomnitigs") .help("do not compute omnitigs");

	parser.add_option("-l", "--list-file").type("string").dest("l").set_default("").help("genome list file");

	optparse::Values& options = parser.parse_args(argc, argv);

	inputFileName = (string) options.get("i");
	input_from_reads = (inputFileName.substr(inputFileName.find_last_of(".") + 1) == "fastq") ? true : false;
	input_from_reads = (inputFileName.substr(inputFileName.find_last_of(".") + 1) == "FASTQ") ? true : false;

	kmersize = (size_t) options.get("k");
	abundance = (int) options.get("a");
	N_THREADS = (int) options.get("t");
	build_only = (options.get("build_only") ? true : false);
	do_not_contract_arcs = (options.get("nocontract") ? true : false);
	do_not_compute_omnitigs = (options.get("noomnitigs") ? true : false);
	genome_type = (string) options.get("g");
	
	// If a genome file list is given, call coderodde_processing and then exit:
	string genome_list_file_name = (string) options.get("l");
	
	if (!genome_list_file_name.empty())
	{
		cout << "[INFO]: Genome list file name = " << genome_list_file_name << endl;
		coderodde_processing(genome_list_file_name, kmersize);
		cout << "[INFO] Done dealing with this file list." << endl;
		return 0;
	}

	if (inputFileName == "")
	{
		cerr << "Parameter -i | --input is required" << endl;
		cerr << "Use --help to see OPTIONS" << endl;
		return EXIT_FAILURE;
	}

	
	if (genome_type.compare("circular") == 0) 
	{
		circular_genome = true;
	}
	else if (genome_type.compare("linear") == 0) 
	{
		circular_genome = false;
	}
	else 
	{
		cerr << "Genome type (option -g) must be either \'linear\' or \'circular\'" << endl;
		return EXIT_FAILURE;
	}

	if ((N_THREADS < 0) or (N_THREADS > 16))
	{
		cerr << "The number of threads must be between 1 and 16" << endl;
		return EXIT_FAILURE;
	}

	if (input_from_reads)
	{
		if (EXIT_SUCCESS != load_data_from_reads(sequence, seqLength, inputFileName, kmersize, circular_genome, do_not_contract_arcs, abundance, graph, length, seqStart, safe_pairs))
		{
			return EXIT_FAILURE;
		}
	}
	else
	{
		if (EXIT_SUCCESS != load_data(sequence, seqLength, inputFileName, kmersize, circular_genome, do_not_contract_arcs, graph, length, seqStart, nodeLabel, safe_pairs))
		{
			return EXIT_FAILURE;
		}		
	}


	/*stats*/ ofstream fileStats;
	/*stats*/ fileStats.open(inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".stats.csv");
	/*stats*/ fileStats << "algorithm,number of nodes,number of edges,runtime (sec)" << endl;

	std::cout.setf(std::ios::boolalpha);
	cout << "Graph has " << countNodes(graph) << " nodes" << endl;
	cout << "Graph has " << countArcs(graph) << " arcs" << endl;
	cout << "Graph has " << countStronglyConnectedComponents(graph) << " strongly connected components" << endl;
	cout << "Graph has unary nodes: " << count_unary_nodes(graph) << endl;
	if (parallelFree(graph) == false)
	{
		// cerr << "The graph still has parallel arcs. Contact a developer." << endl;
		// return EXIT_FAILURE;
	}
	if (count_unary_arcs(graph) != 0)
	{
		// cerr << "Graph still has unary arcs. Contact a developer." << endl;	
		// return EXIT_FAILURE;
	}

	if (build_only)
	{
		return EXIT_SUCCESS;
	}

	/*stats*/ clock_gettime(CLOCK_MONOTONIC, &start_clock);
	cout << "Time: " << currentDateTime();
	vector<contig> unitigs;
	unitigs = compute_unitigs(graph, length, seqStart, kmersize, sequence, inputFileName);
	populate_with_strings(sequence, kmersize, graph, seqStart, length, unitigs);
	cout << "----------------------" << endl;
	cout << "Statistics on unitigs:" << endl;
	cout << "----------------------" << endl;
	compute_statistics(unitigs, seqLength);
	print_collection(unitigs, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type, ".unitigs");
	cout << "Time: " << currentDateTime();
	/*stats*/ clock_gettime(CLOCK_MONOTONIC, &finish_clock);
	/*stats*/ fileStats << "unitigs," << countNodes(graph) << "," << countArcs(graph) << "," << (finish_clock.tv_sec - start_clock.tv_sec) + (finish_clock.tv_nsec - start_clock.tv_nsec) / 1000000000.0 << endl;

	// /*stats*/ clock_gettime(CLOCK_MONOTONIC, &start_clock);
	// vector<contig> non_switching_contigs;
	// non_switching_contigs = compute_non_switching_contigs(graph, length, seqStart, kmersize, sequence, inputFileName);
	// populate_with_strings(sequence, kmersize, graph, seqStart, length, non_switching_contigs);
	// cout << "----------------------" << endl;
	// cout << "Statistics on non-switching contigs:" << endl;
	// cout << "----------------------" << endl;
	// compute_statistics(non_switching_contigs, seqLength);
	// print_collection(non_switching_contigs, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type, ".non-switching-contigs");
	// cout << "Time: " << currentDateTime();
	// /*stats*/ clock_gettime(CLOCK_MONOTONIC, &finish_clock);
	// /*stats*/ fileStats << "non-switching contigs," << countNodes(graph) << "," << countArcs(graph) << "," << (finish_clock.tv_sec - start_clock.tv_sec) + (finish_clock.tv_nsec - start_clock.tv_nsec) / 1000000000.0 << endl;

	/*stats*/ clock_gettime(CLOCK_MONOTONIC, &start_clock);
	cout << "Time: " << currentDateTime();
	vector<contig> YtoV_contigs;
	YtoV_contigs = compute_YtoV_contigs_and_reduce(graph, length, seqStart, nodeLabel, kmersize, sequence, inputFileName, genome_type, seqLength);
	populate_with_strings_from_node_labels(sequence, kmersize, graph, nodeLabel, YtoV_contigs);
	cout << "----------------------" << endl;
	cout << "Statistics on YtoV contigs:" << endl;
	cout << "----------------------" << endl;
	compute_statistics(YtoV_contigs, seqLength);
	print_collection(YtoV_contigs, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type, ".YtoV-contigs");
	cout << "Time: " << currentDateTime();
	/*stats*/ clock_gettime(CLOCK_MONOTONIC, &finish_clock);
	/*stats*/ fileStats << "YtoV contigs," << countNodes(graph) << "," << countArcs(graph) << "," << (finish_clock.tv_sec - start_clock.tv_sec) + (finish_clock.tv_nsec - start_clock.tv_nsec) / 1000000000.0 << endl;

	// vector<contig> omnitigs_2;
	// omnitigs_2 = compute_omnitigs_2(graph, length, seqStart, safe_pairs, kmersize, sequence, inputFileName);
	// populate_with_strings(sequence, kmersize, graph, seqStart, length, omnitigs_2);
	// cout << "----------------------" << endl;
	// cout << "Statistics on omnitigs_2:" << endl;
	// cout << "----------------------" << endl;
	// compute_statistics(omnitigs_2, seqLength);
	// print_collection(omnitigs_2, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type, ".omnitigs_2");
	
	if (! do_not_compute_omnitigs)
	{
		/*stats*/ clock_gettime(CLOCK_MONOTONIC, &start_clock);
		cout << "safe_pairs.size(): " << safe_pairs.size() << endl;
		if (true or (safe_pairs.size() == 0))
		{
			safe_pairs.clear();
			cout << "Will compute safe pairs" << endl;
			vector<contig> omnitigs_2;
			omnitigs_2 = compute_omnitigs_2(graph, length, seqStart, safe_pairs, kmersize, sequence, inputFileName);
			save_safe_pairs(safe_pairs, inputFileName, kmersize, genome_type);
			//populate_with_strings(sequence, kmersize, graph, seqStart, length, omnitigs_2);
			//print_collection(omnitigs_2, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type, ".ec_contig_2");
		}

		vector<contig> omnitigs;
		try_loading_omnitigs(omnitigs,inputFileName,kmersize,genome_type);
		if (omnitigs.size() == 0)
		{
			omnitigs = compute_omnitigs(graph, length, seqStart, safe_pairs);	
		}
		
		// remove_non_maximal_contigs(omnitigs); 
		//populate_with_strings(sequence, kmersize, graph, seqStart, length, omnitigs);
		populate_with_strings_from_node_labels(sequence, kmersize, graph, nodeLabel, omnitigs);
		cout << "----------------------" << endl;
		cout << "Statistics on omnitigs:" << endl;
		cout << "----------------------" << endl;
		compute_statistics(omnitigs, seqLength);
		print_collection(omnitigs, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type, ".omnitigs");
		cout << "Time: " << currentDateTime();	
		/*stats*/ clock_gettime(CLOCK_MONOTONIC, &finish_clock);
		/*stats*/ fileStats << "omnitigs," << countNodes(graph) << "," << countArcs(graph) << "," << (finish_clock.tv_sec - start_clock.tv_sec) + (finish_clock.tv_nsec - start_clock.tv_nsec) / 1000000000.0 << endl;
	}
	
	//////////////////////////////////
	//// MY ENTRY POINT IS BELOW! ////
	//////////////////////////////////
	cout << "[CODERODDE] Steps into the room..." << endl;
	//coderodde_project_algorithm(graph, length, seqStart, kmersize, sequence, inputFileName, true);
	coderodde_project_algorithm(graph, nodeLabel, inputFileName, kmersize, sequence, true);
	cout << "[CODERODDE] Exited the funky algorithm." << endl;
	
	fileStats.close();
	
	get_node_covering_reconstruction(graph, true);

	return EXIT_SUCCESS;
}