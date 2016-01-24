#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <utility>
#include <map>
#include <set>
#include <string>

#include "FileReader.h"
#include "Files.h"
#include "CsrGraph.h"

#include "KernelPruning.h"
#include "utils.h"

int main(int argc,char *argv[])
{
	if(argc < 5)
	{
		printf("There are 3 argument\n");
		printf("The first argument is Input directory path.\n");
		printf("The second argument is the degree_threshold value.\n");
		printf("The third argument should give the output directory name.\n");
		printf("The 4th argument should give the number of nodes.\n");
		exit(1);
	}

	//Create a global csr_graph
	csr_graph *global_graph = new csr_graph();

	//directory Name
	std::string InputDirectoryName = argv[1];
	std::string OutputDirectoryName = argv[3];

	int degree_threshold = atoi(argv[2]);
	int global_nodes_count = atoi(argv[4]);

	//fileList
	std::vector<std::string> fileList = openDirectory(InputDirectoryName);

	/** Create the following vector to hold the global mapping.
	*	This global mapping contains as Keys the index and values as a pair<string,int>
	* The first element of the pair indicates the filename and the second indicates the value of the original vertex.
	**/
	std::vector<std::pair<int,int> > *globalMapping=new std::vector<std::pair<int,int> >();

	//This offset is used to distinguish between the vertices of the graph.
	int GraphOffset = 0; 
	//iterate through the fileList
	for(int i=0;i<fileList.size();i++)
	{	
		//std::cout << fileList[i] << std::endl;

		std::string filePath = InputDirectoryName + "/" + fileList[i];

		//Make an instance of a file reader
		FileReader Reader(filePath.c_str());

		//Construct a local graph from the current file.
		csr_graph *localGraph = new csr_graph();

		int Nodes,Edges;

		Reader.get_nodes_edges(Nodes,Edges);

		#ifdef INFO
			printf("Nodes : %d,Edges: %d\n",Nodes,Edges);
		#endif

		int u,v;

		//Insert an edge between the nodes
		for(int t=0;t<Edges;t++)
		{
			Reader.read_edge(u,v);

			localGraph->insert(u,v,false);
		}

		Reader.fileClose();

		//Calculate the degree and rowoffset of the localgraph
		localGraph->calculateDegreeandRowOffset();

		//This map is used to hold the Vertex --> Relabelled Vertex Information which is used while processing the column vertices
		std::map<int,int > forwardPairs;

		//We iterate through the rowOffsets of the individual graph and add all the unique Nodes in the forwardPairs first
		for(int j=0;j<localGraph->rowOffsets->size();j++)
		{
			//relabelledIndex is j + graphOffset
			int relabelledIndex = j + GraphOffset;

			//Add an entry in the globalMapping vector
			globalMapping->push_back(std::make_pair(std::stoi(fileList[i]),localGraph->rowOffsets->at(j)));

			//Add an entry in the forwardpair
			forwardPairs[localGraph->rowOffsets->at(j)] = relabelledIndex;

			//Push the relabelled Vertex and the degree in the global CSR_GRAPH
			global_graph->rowOffsets->push_back(relabelledIndex);
			global_graph->degree->push_back(localGraph->degree->at(j));

		}

		//We iterate through the rows and columns in the local_graph and add them to the CSR_Graph. We take the relabelled Vertex from the
		//forwardPairs Mapping.
		for(int j=0;j<localGraph->rows->size();j++)
		{
			global_graph->rows->push_back(forwardPairs[localGraph->rows->at(j)]);
			global_graph->columns->push_back(forwardPairs[localGraph->columns->at(j)]);
		}

		forwardPairs.clear();

		GraphOffset += localGraph->rowOffsets->size();

		localGraph = NULL;

	}

	//Call this method to reset the device and free all resources
	DeviceReset();

	int intialNodeCount = global_graph->rowOffsets->size();
	int finalNodeCount  = 0;
	
	//This method is used for pruning the graph. i.e. Removing Vertices having degree less than degree_threshold.
	float totalTime = KernelPruningWrapper(*global_graph,degree_threshold);

	//If atleast one such graph exist
	if(global_graph->rows->size() > 0)
	{
		int graph_count = 1;  //Count of graph.

		int prev_bcc_no = globalMapping->at(global_graph->rows->at(0)).first; //get the first bcc no

		csr_graph *localgraph = new csr_graph();

		for(int index = 0; index < global_graph->rows->size(); index++)
		{
			int curr_bcc_no = globalMapping->at(global_graph->rows->at(index)).first;

			assert(curr_bcc_no == globalMapping->at(global_graph->columns->at(index)).first);

			if(curr_bcc_no == prev_bcc_no)
			{
				//Insert the edge into the local_graph
				localgraph->insert(globalMapping->at(global_graph->rows->at(index)).second,globalMapping->at(global_graph->columns->at(index)).second,true);
			}
			else
			{
				//clear the previous local_graph
				localgraph->calculateDegreeandRowOffset();

				//Output fileName
				std::string fileName = OutputDirectoryName + "/" + std::to_string(graph_count)+".mtx";

				//print the graph to the above file.
				localgraph->PrintToFile(fileName,global_nodes_count);

				//Count of Nodes
				finalNodeCount += localgraph->rowOffsets->size();

				//Clear the localgraph
				localgraph = NULL;

				//Make a new Csr_graph
				localgraph = new csr_graph();

				localgraph->insert(globalMapping->at(global_graph->rows->at(index)).second,globalMapping->at(global_graph->columns->at(index)).second,true);

				//increment the count of the graph.
				graph_count++;

				prev_bcc_no = curr_bcc_no;
			}

		}

		//Print the file at the end.
		localgraph->calculateDegreeandRowOffset();

				//Output fileName

		std::string fileName = OutputDirectoryName + "/" + std::to_string(graph_count)+".mtx";

				//print the graph to the above file.
		localgraph->PrintToFile(fileName,global_nodes_count);

				//Clear the localgraph
		localgraph = NULL;

	}

	//print the total time taken.
	printf("%d\n",intialNodeCount - finalNodeCount);
	printf("%f\n",(totalTime/1000));



	#ifdef INFO
		printf("global rowoffset size = %d,global columns size = %d\n",global_graph->rowOffsets->size(),global_graph->rows->size());
	#endif

	return 0;
}