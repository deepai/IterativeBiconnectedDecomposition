#include <cstdio>
#include <iostream>
#include <omp.h>
#include <climits>
#include <set>
#include <list>
#include <unordered_map>
#include <utility>

#include "FileReader.h"
#include "Files.h"
#include "utils.h"
#include "Host_Timer.h"
#include "CsrGraph.h"

debugger dbg;

char *input_file_name;

struct NodeLists
{
	std::list<int> vertices;
};

struct NodeTrack
{
	std::list<int>::iterator track;
};

int main(int argc,char *argv[])
{

	if(argc < 3)
	{
		printf("Ist Argument should indicate the Input Graph\n");
		printf("2nd argument should indicate the number of threads.(Optional) \n");
		exit(1);
	}

	int num_threads = 1;

	if(argc == 3)
		num_threads = atoi(argv[2]);

    	omp_set_num_threads(num_threads);

        input_file_name = argv[1];

        //Read the Inputfile.
	FileReader Reader(input_file_name);

	int v1,v2,Initial_Vertices;;

	int nodes,edges;

	Reader.get_nodes_edges(nodes,edges);

	csr_graph graph;
	graph.Nodes = nodes;

	for(int i=0;i<edges;i++)
	{
		Reader.read_edge(v1,v2);
		graph.insert(v1,v2,false);
	}

	Reader.fileClose();

	//sort and calculate degree offset.
	graph.calculateDegreeandRowOffset();

	debug("Graph Construction Complete!");

	//Node Buckets
	NodeLists *list_edge_by_size = new NodeLists[nodes];
	NodeTrack *elements_track = new NodeTrack[nodes];

	//Ordering in Degeneracy number.
	int *degeneracy_ordering = new int[nodes];
	int *degree_vertex = new int[nodes];

	//Fill the Buckets
	for(int i=0; i<nodes; i++)
	{
		int start_offset = graph.rowOffsets->at(i);
		int end_offset   = graph.rowOffsets->at(i+1);

		int degree = end_offset - start_offset;

		list_edge_by_size[degree].vertices.push_front(i);

		elements_track[i].track = list_edge_by_size[degree].vertices.begin();

		degree_vertex[i] = degree;

	}

	debug("Buckets build complete");


	int degeneracy_ordering_count = 0;


	debug("Initialization Complete!");

	//This loop is run until all the nodes are ordered.
	while(degeneracy_ordering_count != nodes)
	{
		//obtain the first non_empty bucket.
		for(int i=0; i<nodes; i++ )
		{
			if(list_edge_by_size[i].vertices.size() == 0)
				continue;

			int first_vtx = list_edge_by_size[i].vertices.front();

			list_edge_by_size[i].vertices.erase(elements_track[first_vtx].track);

			for(int j = graph.rowOffsets->at(first_vtx); j<graph.rowOffsets->at(first_vtx + 1); j++)
			{
				int dest_vtx = graph.columns->at(j);
				if(degree_vertex[dest_vtx] < 0)
					continue;

				list_edge_by_size[degree_vertex[dest_vtx]].vertices.erase(elements_track[dest_vtx].track);

				degree_vertex[dest_vtx]--;

				if(degree_vertex[dest_vtx] >= 0)
				{
					list_edge_by_size[degree_vertex[dest_vtx]].vertices.push_front(dest_vtx);
					elements_track[dest_vtx].track = 
						list_edge_by_size[degree_vertex[dest_vtx]].vertices.begin();
				}

			}

			degree_vertex[first_vtx] = -1;

			degeneracy_ordering[first_vtx] = degeneracy_ordering_count++;

			break;

		}
	}

	debug("degeneracy ordering complete");

	int degeneracy = 0;

	// To obtain degeneracy value of a node,count the number of neighbours of a node 
	// lying after it in the degeneracy ordering.
	#pragma omp parallel for reduction(max:degeneracy)
	for(int i=0; i<nodes; i++)
	{
		int start = graph.rowOffsets->at(i);
		int end   = graph.rowOffsets->at(i + 1);

		int _local_degeneracy = 0;

		int curr_position = degeneracy_ordering[i];

		for(int j=start; j<end; j++)
		{
			if(degeneracy_ordering[graph.columns->at(j)] > curr_position)
				_local_degeneracy += 1;
		}

		degeneracy = std::max(degeneracy,_local_degeneracy);
	}

	printf("Degeneracy of the graph is %d\n",degeneracy);

	return 0;

}