#include "utils.h"
#include <stdexcept>

#define OUT

struct tm * timeinfo;

size_t hash_pair( const pair_of_ints& p )
{
    return p.first ^ p.second;
}

typedef unordered_set<pair_of_ints,function<decltype(hash_pair)>> set_of_pairs;


void make_upper_case(string& s)
{
	for (size_t i = 0; i < s.length(); i++)
	{
		if (s[i] == 'a') s[i] = 'A';
		else if (s[i] == 'c') s[i] = 'C';
		else if (s[i] == 'g') s[i] = 'G';
		else if (s[i] == 't') s[i] = 'T';
		else if (s[i] == 'n') s[i] = 'N';
		else if (s[i] == 'A') s[i] = 'A';
		else if (s[i] == 'C') s[i] = 'C';
		else if (s[i] == 'G') s[i] = 'G';
		else if (s[i] == 'T') s[i] = 'T';
		else s[i] = 'N';
	}
}

// Check if a file is readable
bool is_readable( const std::string & file ) 
{ 
    std::ifstream f( file.c_str() ); 
    return !f.fail(); 
} 

size_t count_unary_nodes(StaticDigraph& graph)
{
	size_t n_unary_nodes = 0;
	for (StaticDigraph::NodeIt node(graph); node != INVALID; ++node)
	{
		if ((countInArcs(graph,node) == 1) and (countOutArcs(graph,node) == 1))
		{
			n_unary_nodes++;
		}
	}
	
	return n_unary_nodes;
}

size_t count_unary_arcs(StaticDigraph& graph)
{
	size_t n_unary_arcs = 0;
	for (StaticDigraph::ArcIt arc(graph); arc != INVALID; ++arc)
	{
		if ((countOutArcs(graph,graph.source(arc)) == 1) and (countInArcs(graph,graph.target(arc)) == 1))
		{
			n_unary_arcs++;
		}
	}
	
	return n_unary_arcs;
}

size_t count_unary_arcs_ld(ListDigraph& graph)
{
	size_t n_unary_arcs = 0;
	for (ListDigraph::ArcIt arc(graph); arc != INVALID; ++arc)
	{
		if ((countOutArcs(graph,graph.source(arc)) == 1) and (countInArcs(graph,graph.target(arc)) == 1))
		{
			n_unary_arcs++;
		}
	}
	
	return n_unary_arcs;
}

void populate_with_strings(const string& sequence, 
	size_t kmersize, 
	const StaticDigraph& graph,
	const StaticDigraph::NodeMap<size_t>& seqStart,
	const StaticDigraph::NodeMap<size_t>& length,
	vector<contig>& collection)
{
	for (size_t i = 0; i < collection.size(); i++)
	{
		string contig = "";
		for (list<int>::const_iterator itr = collection[i].nodes.begin(); itr != collection[i].nodes.end(); ++itr)
		{	
			StaticDigraph::Node node = graph.node(*itr);
			string s = sequence.substr(seqStart[node],length[node]);
			contig = plus_strings(contig,s,kmersize);
		}
		collection[i].str = contig;
	}
}

void populate_with_strings_from_node_labels(const string& sequence, 
	size_t kmersize, 
	const StaticDigraph& graph,
	const StaticDigraph::NodeMap<string>& nodeLabel,
	vector<contig>& collection)
{
	for (size_t i = 0; i < collection.size(); i++)
	{
		string contig = "";
		for (list<int>::const_iterator itr = collection[i].nodes.begin(); itr != collection[i].nodes.end(); ++itr)
		{	
			StaticDigraph::Node node = graph.node(*itr);
			string s = nodeLabel[node]; //sequence.substr(seqStart[node],length[node]);
			contig = plus_strings(contig,s,kmersize);
		}
		collection[i].str = contig;
	}
}

void populate_with_strings_list_digraph(const string& sequence, 
	size_t kmersize, 
	const ListDigraph& graph,
	const ListDigraph::NodeMap<string>& nodeLabel,
	vector<contig>& collection)
{
	for (size_t i = 0; i < collection.size(); i++)
	{
		string contig = "";
		for (list<int>::const_iterator itr = collection[i].nodes.begin(); itr != collection[i].nodes.end(); ++itr)
		{	
			ListDigraph::Node node = graph.nodeFromId(*itr);
			string s = nodeLabel[node]; //sequence.substr(seqStart[node],length[node]);
			contig = plus_strings(contig,s,kmersize);
		}
		collection[i].str = contig;
	}
}

void print_collection(vector<contig>& collection, string inputFileName, string suffix) 
{

	ofstream unitigsFile;
	unitigsFile.open(inputFileName + suffix);

	for (size_t i = 0; i < collection.size(); i++)
	{
		unitigsFile << "> ";
		std::list<int>::const_iterator iterator = collection[i].nodes.begin();
		unitigsFile << *iterator;

		for (++iterator; iterator != collection[i].nodes.end(); ++iterator) 
		{
    		unitigsFile << "," << *iterator;
		}
		unitigsFile << endl;
		unitigsFile << collection[i].str << endl;
	}

	unitigsFile.close();

}

void compute_statistics(vector<contig>& collection, size_t seqLength)  
{
	size_t total_size = 0;
	size_t n_strings = collection.size();
	size_t sum_of_squares = 0;
	size_t max_length = 0;
	size_t n_nodes = 0;
	size_t max_n_nodes = 0;

	for (size_t i = 0; i < collection.size(); i++)
	{
		total_size += collection[i].str.length();
		sum_of_squares += collection[i].str.length() * collection[i].str.length();
		if (collection[i].str.length() > max_length)
		{
			max_length = collection[i].str.length();
		}
		n_nodes += collection[i].nodes.size();
		if (collection[i].nodes.size() > max_n_nodes)
		{
			max_n_nodes = collection[i].nodes.size();
		}
	}

	cout << "Number of strings: " << n_strings << endl;
	cout << "Total length of the strings: " << total_size << endl;
	cout << "Average length: " << total_size / (double)n_strings << endl;
	cout << "Maximum length: " << max_length << endl;
	cout << "Average n. nodes: " << n_nodes / (double)n_strings << endl;
	cout << "Maximum n. nodes: " << max_n_nodes << endl;
	cout << "E-size normalized to assembly length: " << sum_of_squares / (double)total_size << endl;
	//cout << "E-size normalized to genome length: " << sum_of_squares / (double)seqLength << endl;

}

void save_safe_pairs(set_of_pairs& safe_pairs, 
	const string inputFileName, 
	const size_t kmersize,
	const string genome_type)
{
	ofstream spFile;
	spFile.open(inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".safepairs");
	for (set_of_pairs::iterator itr = safe_pairs.begin(); itr != safe_pairs.end(); ++itr) 
	{
		spFile << itr->first << " " << itr->second << endl;
	}
	spFile.close();
}

void try_loading_omnitigs(vector<contig>& omnitigs,
	const string inputFileName, 
	const size_t kmersize,
	const string genome_type)
{
	string ecFileName = inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".omnitigs.nonmaximal";
	ifstream ecFile;
}


void print_graph_in_dot(StaticDigraph& graph, 
	const string inputFileName, 
	const size_t kmersize,
	const string genome_type)
{
	ofstream dotFile;
	dotFile.open(inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".dot");
	dotFile << "digraph DB {" << endl;

	for (StaticDigraph::ArcIt arc(graph); arc != INVALID; ++arc)
	{
		dotFile << graph.id(graph.source(arc)) << " -> " << graph.id(graph.target(arc)) << ";" << endl;
	}
	
	dotFile << "}" << endl;
	dotFile.close();
}

void contract_arcs(ListDigraph& graph, 
	ListDigraph::NodeMap<size_t>& length, 
	ListDigraph::NodeMap<size_t>& seqStart,
	const size_t kmersize)
{
	vector<ListDigraph::Arc> arcsToContract;


	size_t n_arcs = countArcs(graph);
	arcsToContract.reserve(n_arcs);

	for (ListDigraph::ArcIt a(graph); a != INVALID; ++a)
	{
		ListDigraph::Node u,v;
		u = graph.source(a);
		v = graph.target(a);

		// OLD: if ((countOutArcs(graph,u) == 1) and (countInArcs(graph,v) == 1))
		// NEW: if the arc (u,v) is the 'middle' arc of a unitig, then contract it
		if ((countOutArcs(graph,u) == 1) and (countInArcs(graph,u) == 1) and
			(countOutArcs(graph,v) == 1) and (countInArcs(graph,v) == 1))
		{
			arcsToContract.push_back(a);
		}
	}

	for (auto arc : arcsToContract)
	{

		if (not graph.valid(arc))
		{
		    
		}

		ListDigraph::Node u,v;
		u = graph.source(arc);
		v = graph.target(arc);
		length[u] = length[u] + length[v] - (kmersize - 1);

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
}

void contract_arcs_from_reads(ListDigraph& graph, 
	ListDigraph::NodeMap<size_t>& length, 
	ListDigraph::NodeMap<size_t>& seqStart,
	const size_t kmersize,
	string& sequence)
{
	vector<ListDigraph::Arc> arcsToContract;
	ListDigraph::NodeMap<string> nodeLabels(graph);

	for (ListDigraph::NodeIt node(graph); node != INVALID; ++node)
	{
		nodeLabels[node] = sequence.substr(seqStart[node],length[node]);
	}
	

	size_t progress = 0;
	size_t n_arcs = countArcs(graph);
	arcsToContract.reserve(n_arcs);

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

		// OLD: if ((countOutArcs(graph,u) == 1) and (countInArcs(graph,v) == 1))
		// NEW: if the arc (u,v) is the 'middle' arc of a unitig, then contract it
		if ((countOutArcs(graph,u) == 1) and (countInArcs(graph,u) == 1) and
			(countOutArcs(graph,v) == 1) and (countInArcs(graph,v) == 1))
		{
			arcsToContract.push_back(a);
		}
	}

	cout << "arcsToContract has size: " << arcsToContract.size() << endl;
	cout << "count_unary_arcs gives: " << count_unary_arcs_ld(graph) << endl;

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
		//length[u] = length[u] + length[v] - (kmersize - 1);
		nodeLabels[u] = plus_strings(nodeLabels[u], nodeLabels[v], kmersize);

		// the following code is doing "graph.contract(u, v);" (this lemon function doesn't behave as expected)
		for (ListDigraph::OutArcIt out_arc(graph,v); out_arc != INVALID;)
		{
			ListDigraph::OutArcIt next_out_arc = out_arc;
			++next_out_arc;
			graph.changeSource(out_arc,u);
			out_arc = next_out_arc;
		}
		graph.erase(v);

	}

	cout << "after compacting, count_unary_arcs gives: " << count_unary_arcs_ld(graph) << endl;

	// creating the new sequence and changing
	// length[] and seqStart[] maps to be relative to this new sequence
	sequence = "";
	for (ListDigraph::NodeIt node(graph); node != INVALID; ++node)
	{
		seqStart[node] = sequence.length();
		length[node] = nodeLabels[node].length();
		sequence += nodeLabels[node];
	}

	cout << "Compacted the graph" << endl;

}

void construct_graph_from_multiple_sequences(ListDigraph& graph,
					     ListDigraph::NodeMap<size_t>& length,
					     ListDigraph::NodeMap<size_t>& seqStart,
					     const size_t kmersize,
					     vector<string>& sequence_vector,
					     OUT string& output_total_sequence,
					     bool debugPrint)
{
    if (debugPrint)
    {
	cout << "[DEBUG] construct_graph_from_multiple_sequences is here!" << endl;
    }
	
    output_total_sequence.clear();
    ListDigraph::Node current_node, previous_node = INVALID;
    ListDigraph::NodeMap<string> nodes_to_kmers_map(graph);
    
    // maps the k-mer to the node ID.
    unordered_map<string, int> node_map;
    
    // maps the k-mer to the set of IDs of the parent nodes.
    unordered_map<string, unordered_set<int>> arc_map;
    
    int char_index = 0;
    
    for (string sequence : sequence_vector)
    {
	size_t kmers_limit = sequence.length();
	//sequence = sequence + sequence.substr(0, kmersize - 1);
	sequence = sequence + sequence;
	output_total_sequence += sequence;
	
	string initial_kmer = sequence.substr(0, kmersize);
	
	for (size_t i = 0; i != kmers_limit; ++i, ++char_index)
	{
	    string current_kmer = sequence.substr(i, kmersize);
	    
	    if (current_kmer.find("#") != std::string::npos)
	    {
		previous_node = INVALID;
		throw std::runtime_error{"Found '#' in a k-mer!"};
	    }
	    else
	    {
		//cout << "Current k-mer: " << current_kmer << endl;
		
		if (node_map.find(current_kmer) != node_map.end())
		{
		    //cout << "... already exists." << endl;
		    current_node = graph.nodeFromId(node_map[current_kmer]);
		}
		else
		{
		    //cout << "... does not exist." << endl;
		    current_node = graph.addNode();
		    nodes_to_kmers_map[current_node] = current_kmer;
		    node_map[current_kmer] = graph.id(current_node);
		    length[current_node] = kmersize;
		    seqStart[current_node] = char_index;
		}
		
		if (previous_node != INVALID)
		{
		    if (arc_map[current_kmer].find(graph.id(previous_node)) == arc_map[current_kmer].end())
		    {
			//cout << "ARC!" << endl;
			graph.addArc(previous_node, current_node);
			arc_map[current_kmer].insert(graph.id(previous_node));
		    }
		}
		
		previous_node = current_node;
	    }    
	}
	
	// Deal with the closing arc of the circular genome:
	string first_kmer = sequence.substr(0, kmersize);
	string last_kmer = sequence.substr(kmers_limit - 1, kmersize);
	//string last_kmer = sequence.substr(sequence.length() - kmersize, kmersize);
	
	unordered_set<int>& parent_node_id_set_of_first_node = arc_map[first_kmer];
	int last_kmer_node_id = node_map[last_kmer];
	
	if (parent_node_id_set_of_first_node.find(last_kmer_node_id) == parent_node_id_set_of_first_node.end())
	{
	    //cout << "Adding the missing cycle closing arc." << endl;
	    // The closing arc is not yet in the graph, add it:
	    ListDigraph::Node first_kmer_node = graph.nodeFromId(node_map[first_kmer]);
	    ListDigraph::Node last_kmer_node  = graph.nodeFromId(node_map[last_kmer]);
	    
	    graph.addArc(last_kmer_node, first_kmer_node);
	    arc_map[first_kmer].insert(graph.id(last_kmer_node));
	}
	else
	{
	    //cout << "No need for closing the cycle." << endl;
	}
	
	char_index += kmers_limit;
	previous_node = INVALID;
    }
    
    cout << "[DEBUG]: Before contracting, the graph has " << countNodes(graph) << " nodes "
         << "and " << countArcs(graph) << " arcs." << endl;
	 
    /*cout << "Uncontracted graph:" << endl;
    
    for (ListDigraph::NodeIt nodeit(graph); nodeit != INVALID; ++nodeit)
    {
	cout << "Current node: \"" << nodes_to_kmers_map[nodeit] << "\"" << endl;
	cout << "Children:" << endl;
	
	for (ListDigraph::OutArcIt outArcIt(graph, nodeit); outArcIt != INVALID; ++outArcIt)
	{
	    ListDigraph::Node child = graph.target(outArcIt);
	    cout << "    " << nodes_to_kmers_map[child] << endl;
	}
	
	cout << endl;
    }*/
    
    contract_arcs(graph,
		  length,
		  seqStart,
		  kmersize);
    
    cout << "[DEBUG]: After contracting, the graph has " << countNodes(graph) << " nodes "
         << "and " << countArcs(graph) << " arcs." << endl;
	 
    /*cout << "Contracted graph:" << endl;
    
    for (ListDigraph::NodeIt nodeit(graph); nodeit != INVALID; ++nodeit)
    {
	// Compute the node label:
	string current_kmer = output_total_sequence.substr(seqStart[nodeit], length[nodeit]);
	
	cout << "Current node: \"" << current_kmer << "\"" << endl;
	cout << "Children:" << endl;
	
	for (ListDigraph::OutArcIt outArcIt(graph, nodeit); outArcIt != INVALID; ++outArcIt)
	{
	    ListDigraph::Node child = graph.target(outArcIt);
	    string child_kmer = output_total_sequence.substr(seqStart[child], length[child]);
	    cout << "    " << child_kmer << endl;
	}
	
	cout << endl;
    }*/
    
    //cout << "[DEBUG] The result graph has " << countNodes(graph) << " nodes." << endl;
    //cout << "[DEBUG] The result graph has " << countArcs(graph) << " arcs." << endl;
    cout << "[DEBUG] construct_graph_from_multiple_sequences is done!" << endl;
}

void construct_graph(ListDigraph& graph, 
	ListDigraph::NodeMap<size_t>& length, 
	ListDigraph::NodeMap<size_t>& seqStart,
	const size_t kmersize, 
	const bool circular_genome,
	const string& sequence,
	const size_t seqLength,
	const int abundance)
{
	cout << "[DEBUG] construct_graph is here!" << endl;
	unordered_map<std::string,int> kmersToNodes;
	unordered_map<int,int> nodeAbundance;
	unordered_map<std::string,unordered_set<int>> in_nbrs;

	ListDigraph::Node currentNode, previousNode = INVALID;
	std::string currentKmer;
	// std::string currentKmer = sequence.substr(0,kmersize);
	// // cout << "current kmer is: " << currentKmer << endl;
	// previousNode = graph.addNode();
	// kmersToNodes[currentKmer] = graph.id(previousNode);
	// length[previousNode] = kmersize;
	// seqStart[previousNode] = 0;

	size_t kmerLimit;
	if (circular_genome)
	{
		kmerLimit = seqLength;
	} else
	{
		kmerLimit = seqLength - kmersize + 1;
	}

	size_t progress = 0;

	// traversing the sequence
	for (size_t i = 0; i <= kmerLimit; i++)
	{

		if (progress % 1000000 == 0)
		{
			cout << "Sequence position: #" << progress << "/" << kmerLimit << " ";
		  	cout << "Time: " << currentDateTime();
		}
		progress++;

		// if new kmer, then add a node to the graph
		currentKmer = sequence.substr(i,kmersize);

		// if current kmer contains # then skip it
		if (currentKmer.find("#") != std::string::npos)
		{
			previousNode = INVALID;
		}
		else 
		{
			// current kmer does not contain #
			// cout << "current kmer is: " << currentKmer << endl;
			if (kmersToNodes.count(currentKmer) > 0) 
			{
				currentNode = graph.nodeFromId(kmersToNodes[currentKmer]);
				nodeAbundance[kmersToNodes[currentKmer]]++;
			}
			else
			{
				currentNode = graph.addNode();
				//kmersToNodes.insert(std::make_pair<std::string,StaticDigraph::Node>(currentKmer, currentNode));
				kmersToNodes[currentKmer] = graph.id(currentNode);
				nodeAbundance[kmersToNodes[currentKmer]] = 1;
				length[currentNode] = kmersize;
				seqStart[currentNode] = i;
				// cout << "inserted seq start" << seqStart[kmersToNodes[currentKmer]] << endl;
			}

			if (previousNode != INVALID)
			{
				if (in_nbrs[currentKmer].count(graph.id(previousNode)) == 0)
				{
					graph.addArc(previousNode, currentNode);
					in_nbrs[currentKmer].insert(graph.id(previousNode));
				}	
			}
			previousNode = currentNode;
		}
	}	

	// traversing the graph and removing the nodes with abundance strictly lower than abundance
	vector<ListDigraph::Node> nodesToRemove;

	// finding the nodes to remove
	for (ListDigraph::NodeIt node(graph); node != INVALID; ++node)
	{
		if (nodeAbundance[graph.id(node)] < abundance)
		{
			nodesToRemove.push_back(node);
		}
	}

	// removing the nodes from nodesToRemove
	for (auto node : nodesToRemove)
	{
		graph.erase(node);
	}


	cout << "Constructed the graph" << endl;

}

int load_data(string& sequence, 
	size_t& seqLength, 
	const string& inputFileName,
	const size_t kmersize,
	const bool circular_genome,
	const bool do_not_contract_arcs,
	StaticDigraph& graph,
	StaticDigraph::NodeMap<size_t>& length,
	StaticDigraph::NodeMap<size_t>& seqStart,
	StaticDigraph::NodeMap<string>& nodeLabel,
	set_of_pairs& safe_pairs
)
{
	cout << "[DEBUG] sequence.length(): " << sequence.length() << endl;
	cout << "[DEBUG] load_data() is here!" << endl;
	
	string genome_type;
	if (circular_genome)
	{
		genome_type = "circular";
	}
	else
	{
		genome_type = "linear";
	}

	try 
	{
		ifstream sequenceFile;
		string line;
        sequenceFile.open(inputFileName);
        getline(sequenceFile, line); // reading the header
        while (getline(sequenceFile, line))
        {
        	sequence += line;
        }

        seqLength = sequence.length();
        make_upper_case(sequence);

        if (circular_genome)
        {
        	sequence = sequence + sequence; // + sequence.substr(0,kmersize);
        }
        sequenceFile.close();  
	} catch (Exception& error) 
	{ // check if there was any error
		std::cerr << "Error: " << error.what() << std::endl;
		return EXIT_FAILURE;
	}

	// NOTE: Program control continues to else branch.
	if (! is_readable(inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".lgf"))
	{
		cout << "[DEBUG] !is_readable" << endl;
		ListDigraph temporary_graph;
		ListDigraph::NodeMap<size_t> temporary_length(temporary_graph);
		ListDigraph::NodeMap<size_t> temporary_seqStart(temporary_graph);

		construct_graph(temporary_graph, temporary_length, temporary_seqStart, kmersize, circular_genome, sequence, seqLength, 1);
		if (not do_not_contract_arcs)
		{
			contract_arcs(temporary_graph, temporary_length, temporary_seqStart, kmersize);
		}
		//contract_arcs(temporary_graph, temporary_length, temporary_seqStart, kmersize);
		// saving graph to file
		try 
		{
			// copy temporary_graph into graph
			ListDigraph::NodeMap<StaticDigraph::Node> nodeRef(temporary_graph);
			ListDigraph::ArcMap<StaticDigraph::Arc> arcRef(temporary_graph);
			graph.build(temporary_graph, nodeRef, arcRef);
			// copy lists
			for (ListDigraph::NodeIt node(temporary_graph); node != INVALID; ++node)
			{
				length[nodeRef[node]] = temporary_length[node];
				seqStart[nodeRef[node]] = temporary_seqStart[node];
			}
			digraphWriter(graph, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".lgf").
				nodeMap("length", length).
				nodeMap("seqStart", seqStart).
				run();
			// print_graph_in_dot(graph, inputFileName, kmersize, genome_type);

		} catch (Exception& error) 
		{ // check if there was any error
			std::cerr << "Error: " << error.what() << std::endl;
			return EXIT_FAILURE;
		}
		
	}
	// else
	// THERE IS SOME BUG WHEN NOT LOADING THE GRAPH FROM FILE
	// SO WE ALWAYS LOAD IT FROM FILE
	{
		// NOTE: Always entered branch.
		cout << "[DEBUG] is_readable" << endl;
		try 
		{
			cout << "[DEBUG] sequence.length(): " << sequence.length() << endl;
			ListDigraph temporary_graph;
			ListDigraph::NodeMap<size_t> temporary_length(temporary_graph);
			ListDigraph::NodeMap<size_t> temporary_seqStart(temporary_graph);
    
			cout << "[DEBUG] Temporary graph size before: " << countNodes(temporary_graph) << endl;
			digraphReader(temporary_graph, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".lgf").
				nodeMap("length", temporary_length).
				nodeMap("seqStart", temporary_seqStart).
				run();
			cout << "[DEBUG] Temporary graph size after:  " << countNodes(temporary_graph) << endl;
				
			// copy temporary_graph into graph
			ListDigraph::NodeMap<StaticDigraph::Node> nodeRef(temporary_graph);
			ListDigraph::ArcMap<StaticDigraph::Arc> arcRef(temporary_graph);
			
			cout << "[DEBUG] StaticGraph size: " << countNodes(graph) << endl;
			graph.build(temporary_graph, nodeRef, arcRef);
			cout << "[DEBUG] StaticGraph size: " << countNodes(graph) << endl;
			
			// copy lists
			for (ListDigraph::NodeIt node(temporary_graph); node != INVALID; ++node)
			{
				length[nodeRef[node]] = temporary_length[node];
				seqStart[nodeRef[node]] = temporary_seqStart[node];
			}
			// attach the node labels
			for (StaticDigraph::NodeIt node(graph); node != INVALID; ++node)
			{
				nodeLabel[node] = sequence.substr(seqStart[node],length[node]);
			}


		}
		catch (Exception& error) 
		{ // check if there was any error
			std::cerr << "Error: " << error.what() << std::endl;
			return EXIT_FAILURE;
		}
	}
	
	cout << "[DEBUG] sequence.length(): " << sequence.length() << endl;

	if (is_readable(inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".safepairs"))
	{
		int first, second;
		ifstream spFile;
		spFile.open(inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".safepairs");
		while (spFile >> first)
		{
			spFile >> second;
			safe_pairs.insert(pair_of_ints(first,second));
		}
		spFile.close();
		cout << safe_pairs.size() << endl;
	}

	return EXIT_SUCCESS;
}

int load_data_from_reads(string& sequence, 
	size_t& seqLength, 
	const string& inputFileName,
	const size_t kmersize,
	const bool circular_genome,
	const bool do_not_contract_arcs,
	const int abundance,
	StaticDigraph& graph,
	StaticDigraph::NodeMap<size_t>& length,
	StaticDigraph::NodeMap<size_t>& seqStart,
	set_of_pairs& safe_pairs
)
{
	cout << "[DEBUG] load_data_from_reads is here!" << endl;
    
	string genome_type;
	if (circular_genome)
	{
		genome_type = "circular";
	} else
	{
		genome_type = "linear";
	}

	if (! is_readable(inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".lgf"))
	{
		try 
		{
			ifstream sequenceFile;
			string line;
	        sequenceFile.open(inputFileName);
	        
	        while (getline(sequenceFile, line)) // reading the header
	        {
	        	getline(sequenceFile, line); // the read itself
	        	make_upper_case(line);
	        	sequence += line + "#";
	        	sequence += reverse_complement(line) + "#";
	        	getline(sequenceFile, line);
	        	getline(sequenceFile, line);
	        }

	        seqLength = sequence.length();
	        make_upper_case(sequence);

	        sequenceFile.close();  
		} catch (Exception& error) 
		{ // check if there was any error
			std::cerr << "Error: " << error.what() << std::endl;
			return EXIT_FAILURE;
		}

		ListDigraph temporary_graph;
		ListDigraph::NodeMap<size_t> temporary_length(temporary_graph);
		ListDigraph::NodeMap<size_t> temporary_seqStart(temporary_graph);

		construct_graph(temporary_graph, temporary_length, temporary_seqStart, kmersize, circular_genome, sequence, seqLength, abundance);
		if (not do_not_contract_arcs)
		{
			contract_arcs_from_reads(temporary_graph, temporary_length, temporary_seqStart, kmersize, sequence);	
		}
		//contract_arcs(temporary_graph, temporary_length, temporary_seqStart, kmersize);
		// saving graph to file
		try 
		{
			// copy temporary_graph into graph
			ListDigraph::NodeMap<StaticDigraph::Node> nodeRef(temporary_graph);
			ListDigraph::ArcMap<StaticDigraph::Arc> arcRef(temporary_graph);
			graph.build(temporary_graph, nodeRef, arcRef);
			// copy lists
			for (ListDigraph::NodeIt node(temporary_graph); node != INVALID; ++node)
			{
				length[nodeRef[node]] = temporary_length[node];
				seqStart[nodeRef[node]] = temporary_seqStart[node];
			}
			digraphWriter(graph, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".lgf").
				nodeMap("length", length).
				nodeMap("seqStart", seqStart).
				attribute("sequence", sequence).
				run();
			// print_graph_in_dot(graph, inputFileName, kmersize, genome_type);

		} catch (Exception& error) 
		{ // check if there was any error
			std::cerr << "Error: " << error.what() << std::endl;
			return EXIT_FAILURE;
		}
		
	}
	else
	{
		try 
		{
			ListDigraph temporary_graph;
			ListDigraph::NodeMap<size_t> temporary_length(temporary_graph);
			ListDigraph::NodeMap<size_t> temporary_seqStart(temporary_graph);

			digraphReader(temporary_graph, inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".lgf").
				nodeMap("length", temporary_length).
				nodeMap("seqStart", temporary_seqStart).
				attribute("sequence", sequence).
				run();

			// copy temporary_graph into graph
			ListDigraph::NodeMap<StaticDigraph::Node> nodeRef(temporary_graph);
			ListDigraph::ArcMap<StaticDigraph::Arc> arcRef(temporary_graph);
			graph.build(temporary_graph, nodeRef, arcRef);
			// copy lists
			for (ListDigraph::NodeIt node(temporary_graph); node != INVALID; ++node)
			{
				length[nodeRef[node]] = temporary_length[node];
				seqStart[nodeRef[node]] = temporary_seqStart[node];
			}

			print_graph_in_dot(graph, inputFileName, kmersize, genome_type);

		}
		catch (Exception& error) 
		{ // check if there was any error
			std::cerr << "Error: " << error.what() << std::endl;
			return EXIT_FAILURE;
		}
	}

	if (is_readable(inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".safepairs"))
	{
		int first, second;
		ifstream spFile;
		spFile.open(inputFileName + ".k" + std::to_string(kmersize) + "." + genome_type + ".safepairs");
		while (spFile >> first)
		{
			spFile >> second;
			safe_pairs.insert(pair_of_ints(first,second));
		}
		spFile.close();
		cout << safe_pairs.size() << endl;
	}

	return EXIT_SUCCESS;
}

