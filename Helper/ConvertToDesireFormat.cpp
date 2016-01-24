#include <cstdio>
#include <cstdlib>

#include "FileReader.h"
#include "FileWriter.h"

#include <map>
#include <set>
#include <fstream>
#include <assert.h>

//Method used to insert into the adjacency list
void insert(std::map<int,std::set<int>*> &adjList,int u,int v,bool direction)
{
    if(adjList.find(u)==adjList.end())
    {
        std::set<int> *temp = new std::set<int>();
        temp->insert(v);
        adjList[u] = temp;
    }
    else
    {
        adjList[u]->insert(v);
    }

    if(!direction)
        insert(adjList,v,u,true);
}

int main(int argc,char *argv[])
{
    if(argc < 5)
    {
        printf("There are 4 arguments. 1) Input file and 2) Outputfile, 3) Nodes. 4) mtxformat_output 1 for true and 0 for false\n");
        exit(1);
    }

    int global_node_count = atoi(argv[3]);
    bool mtxformat_output = atoi(argv[4]);

    FileReader Reader(argv[1]);

    std::map<int,std::set<int>*> adjList;

    int Nodes,Edges;

    Reader.get_nodes_edges(Nodes,Edges);

    int v1,v2;

    for(int i=0;i<Edges;i++)
    {
      Reader.read_edge(v1,v2);
      if(v1!=v2)
        insert(adjList,v1,v2,false);
    }

    Reader.fileClose();

    if(!mtxformat_output)
    {

        std::ofstream fout(argv[2],std::ios::out);

        //count edges
        Edges = 0;
        for(std::map<int,std::set<int> *>::iterator it = adjList.begin();it!=adjList.end() ;it++)
        {  
            Edges += it->second->size();
        }
        
        fout << global_node_count << std::endl; //write the nodes
        fout << Edges << std::endl; //write the edges

        for(std::map<int,std::set<int> *>::iterator it = adjList.begin();it!=adjList.end() ;it++)
        {
            int src = it->first;
            for(std::set<int>::iterator ij=it->second->begin();ij!=it->second->end();ij++)
            {
              int dest = *ij;

              fout << src << " " << dest << std::endl;
            }
        }

        fout.close();

    }

    else
    {
        Edges = 0;
        for(std::map<int,std::set<int> *>::iterator it = adjList.begin();it!=adjList.end() ;it++)
        {  
            Edges += it->second->size();
        }
        
        FileWriter Writer(argv[2],global_node_count,Edges/2);
        for(std::map<int,std::set<int> *>::iterator it = adjList.begin();it!=adjList.end() ;it++)
        {
            int src = it->first;
            for(std::set<int>::iterator ij=it->second->begin();ij!=it->second->end();ij++)
            {
                  int dest = *ij;

                  if(src > dest)
                  {
                          Writer.write_edge(src,dest);
                  }
            }
        }

        Writer.fileClose();

    }

    
}