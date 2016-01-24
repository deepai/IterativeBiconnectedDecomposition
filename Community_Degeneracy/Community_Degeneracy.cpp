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

struct Intersections
{
	int count;

	int edge_index;

	int ordering;

	std::vector<int> Head;

	std::list<int>::iterator tracker;
};

struct NodeLists
{
	std::list<int> edge_values;
};

inline unsigned long long merge(unsigned long long upper,unsigned long long lower)
{
	unsigned long long result = 0;
	result = ((upper << 32) | lower);
	return result;
}

std::unordered_map<unsigned long long,int> edge_map;


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

	Intersections *edge_lists = new Intersections[edges];

	for(int i=0,cnt = 0;i < 2*edges ;i++)
	{
		unsigned src_vtx = graph.rows->at(i);
		unsigned dest_vtx = graph.columns->at(i);

		if(src_vtx < dest_vtx)
		{
			edge_lists[cnt].edge_index = i;
			edge_lists[cnt].ordering = -1;

			edge_map.insert(std::make_pair(merge(src_vtx,dest_vtx),cnt));

			cnt++;
		}
	}

	//Build the Intersection Edge List
	#pragma omp parallel for
	for(int i = 0; i < edges; ++i)
	{
		unsigned src_vtx = graph.rows->at(edge_lists[i].edge_index);
		unsigned dest_vtx = graph.columns->at(edge_lists[i].edge_index);

		int start_src = graph.rowOffsets->at(src_vtx);
		int start_dest = graph.rowOffsets->at(dest_vtx);

		int end_src = graph.rowOffsets->at(src_vtx + 1);
		int end_dest = graph.rowOffsets->at(dest_vtx + 1);

		int j = start_src;
		int k = start_dest;

		int count = 0;

		while( (j<end_src) && (k<end_dest))
		{
			if( graph.columns->at(j) < graph.columns->at(k) )
				j++;
			else if( graph.columns->at(j) > graph.columns->at(k) )
				k++;
			else
			{
				unsigned common_vtx = graph.columns->at(j);

				int edge_src_common = edge_map[merge(std::min(src_vtx,common_vtx),
					std::max(src_vtx,common_vtx))];

				int edge_dest_common = edge_map[merge(std::min(dest_vtx,common_vtx),
					std::max(dest_vtx,common_vtx))];

				//debug( src_vtx + 1,dest_vtx + 1,graph.columns->at(j) + 1);
				edge_lists[i].Head.push_back(edge_src_common);
				edge_lists[i].Head.push_back(edge_dest_common);

				j++;
				k++;

				count++;
			}
		}
		edge_lists[i].count = count;

		//debug("size is ",edge_lists[i].Head.size(),count);
	}

	debug("Intersection List build complete.");

	//Node Buckets
	NodeLists *list_edge_by_size = new NodeLists[nodes];

	debug("Allocated Buckets",edges);

	//Fill the Buckets
	for(int i=0; i<edges ;i++)
	{
		int count_common_neighbours = edge_lists[i].count;
		list_edge_by_size[count_common_neighbours].edge_values.push_front(i); 

		edge_lists[i].tracker = list_edge_by_size[count_common_neighbours].edge_values.begin();
		//debug("insert ",i,count_common_neighbours,nodes);
	}

	debug("Buckets build complete");

	int edge_ordering_count = 0;

	//Continue till all the edges are ordered.
	while(edge_ordering_count != edges)
	{
		//Find the first non_empty bucket.
		for(int i=0; i< (nodes - 1) ;i++)
		{
			//continue if bucket[i] is empty.
			if(list_edge_by_size[i].edge_values.size() == 0)
				continue;

			int first_edge_in_list = list_edge_by_size[i].edge_values.front();

			list_edge_by_size[i].edge_values.erase(edge_lists[first_edge_in_list].tracker);

			int csr_edge = edge_lists[first_edge_in_list].edge_index;

			///debug("Edges are ",graph.rows->at(csr_edge) + 1, graph.columns->at(csr_edge) + 1);

			//Iterate through the common nodes of the first_edge_in_the_list end points.
			for(int j = 0; j < edge_lists[first_edge_in_list].Head.size()/2 ; j++)
			{
				//each common vertex
				int edge_index_first  = edge_lists[first_edge_in_list].Head[2*j];
				int edge_index_second = edge_lists[first_edge_in_list].Head[2*j + 1];

				if((edge_lists[edge_index_first].count >= 0) && 
					(edge_lists[edge_index_second].count >=0) )
				{

					int current_count_first  = edge_lists[edge_index_first].count--;
					int current_count_second = edge_lists[edge_index_second].count--;

					list_edge_by_size[current_count_first].edge_values.
						erase(edge_lists[edge_index_first].tracker);
					list_edge_by_size[current_count_second].edge_values.
						erase(edge_lists[edge_index_second].tracker);
					
					if(current_count_first >= 1)
					{
						list_edge_by_size[current_count_first - 1].edge_values.
							push_front(edge_index_first);

						edge_lists[edge_index_first].tracker = 
							list_edge_by_size[current_count_first -1].edge_values.begin();
						
					}

					if(current_count_second >= 1)
					{
						list_edge_by_size[current_count_second - 1].edge_values.
							push_front(edge_index_second);

						edge_lists[edge_index_second].tracker = 
							list_edge_by_size[current_count_second -1].edge_values.begin();
					}

				}

			}

			//remove this edge.
			edge_lists[first_edge_in_list].count = -1;
			edge_lists[first_edge_in_list].ordering = edge_ordering_count++;

			if(edge_ordering_count%1000 == 0)
			{
				//printf("\redges_completed : %d",edge_ordering_count);
			}

			break;
		}
	}

	debug("\nEdge_Ordering Complete");

	int community_degeneracy = 0;

	#pragma omp parallel for reduction(max:community_degeneracy)
	for(int i = 0; i < edges; ++i)
	{

		int current_ordering = edge_lists[i].ordering;

		unsigned _src_vtx_current = graph.rows->at(edge_lists[i].edge_index);
		unsigned _dest_vtx_current = graph.columns->at(edge_lists[i].edge_index);

		//debug("EdgeList:",_src_vtx_current + 1,_dest_vtx_current + 1,"order:",current_ordering);

		std::set<unsigned> unique_nodes;

		for(int j=0 ;j<edge_lists[i].Head.size()/2 ; j++)
		{
			//debug("Size of list is ",edge_lists[i].Head.size());
			
			int edge_index_first  = edge_lists[i].Head[2*j];
			int edge_index_second = edge_lists[i].Head[2*j + 1];

			//debug("Indexes:",i,edge_index_first,edge_index_second);

			if( (edge_lists[edge_index_first].ordering > current_ordering) &&
				(edge_lists[edge_index_second].ordering > current_ordering) )
			{

				//Vertices of the edges.

				unsigned _src_vtx_first = graph.rows->at(edge_lists[edge_index_first].edge_index);
				unsigned _dest_vtx_first = graph.columns->at(edge_lists[edge_index_first].edge_index);

				//debug("Edge :",_src_vtx_first + 1,_dest_vtx_first + 1);

				unsigned _src_vtx_second = graph.rows->at(edge_lists[edge_index_second].edge_index);
				unsigned _dest_vtx_second = graph.columns->at(edge_lists[edge_index_second].edge_index);

				//debug("Edge :",_src_vtx_second + 1,_dest_vtx_second + 1);

				//Unique Nodes. 
				unsigned common_node = ((_src_vtx_first + _src_vtx_second + _dest_vtx_first + _dest_vtx_second)
					- (_src_vtx_current + _dest_vtx_current))/2;

				unique_nodes.insert(common_node);

			}
		}
		//debug(unique_nodes.size());

		community_degeneracy = std::max(community_degeneracy,(int)(unique_nodes.size()));

		unique_nodes.clear();
	}

	printf("Community Degeneracy is %d\n",community_degeneracy);

	return 0;

}