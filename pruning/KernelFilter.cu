#include "KernelPruning.h"

#include "cub.cuh"
#include "utils.h"

#include <thrust/device_vector.h>
#include <vector>

#include <cstdio>
#include "gputimer.h"

#define CAST(dv) (thrust::raw_pointer_cast(dv.data()))

debugger dbg;


/**
 * @brief [This is a device method which returns a boolean value depending on the whether the condition is satisfied]
 * @details [In this method, the condition is whether the degree is greater than the threshold. If yes, it returns true else false.]
 * 
 * @param degree [degree of the current node]
 * @param degreeThreshold [threshold condition]
 * 
 * @return [nothing is returned]
 */
__device__
bool condition(unsigned degree,unsigned degreeThreshold)
{
	return (degree > degreeThreshold);
}

/**
 * @brief [This Kernel is used to filter the vertices having degree > degreeThreshold]
 * @details [We first check the degree of the row vertices , if the result is 0, then set the result in the current position,
 * 			else we check the column vertices and then set the corresponding result.]
 * 
 * @param d_rows [This array contains the "source" vertices.]
 * @param d_cols [This array contains the "destination" vertices] 
 * @param d_degree [This array contains the degree of each of the vertex]
 * @param d_filterArray [This array is used to output the result i.e 0 if condition fails else 1.]
 * @param degreeThreshold [Condition threshold]
 * @param count [Count of number of columns.]
 */
__global__ 
void KernelFilterColOffset(unsigned *d_rows,unsigned *d_cols,unsigned *d_degree,unsigned *d_filterArray,unsigned degreeThreshold,int count)
{
	//thread Id of the current thread;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	//bounds check
	if(tid >= count)
		return;
	
	//obtain the row first
	unsigned currentRowValue = d_rows[tid];

	//degree corresponding to the row.
	unsigned degree_RowValue = d_degree[currentRowValue];


	int result = condition(degree_RowValue , degreeThreshold);

	if(result == 0)
	{
		d_filterArray[tid] = result;
	}
	else
	{

		//obtain the column value
		unsigned currentColValue = d_cols[tid];

		//obtain the degree corresponding to the column.
		unsigned degree_ColValue = d_degree[currentColValue];

		//Evaluate the condition and obtain the results
		int result = condition(degree_ColValue , degreeThreshold);

		//set the filter value.
		d_filterArray[tid] = result;
	}

}

/**
 * @brief [This Kernel involved each thread updating the degree of the remaining unfiltered Vertices]
 * @details [We update the d_degree[d_uniqueOut[tid]] = d_countsOut[tid]]
 * 
 * @param d_degree [original degree array]
 * @param d_uniqueOut [unique rows offset Array]
 * @param d_countsOut [unique Counts array corresponding to each RowOffset]
 * @param count_elements [Number of element, used for boundary condition]
 * @return [description]
 */
__global__
void KernelUpdateDegree(unsigned *d_degree,unsigned *d_uniqueOut,unsigned *d_countsOut,int count_elements)
{
		int tid = threadIdx.x + blockDim.x*blockIdx.x;
		if(tid >= count_elements)
				return;
		d_degree[d_uniqueOut[tid]] = d_countsOut[tid];
}


void DeviceReset()
{
	CudaError(cudaDeviceReset());
}
/**
 * @brief [This method is used to invoke the cub::select method]
 * @details [We first determine the buffer size in the first invocation. We then pass either NULL or the buffer array as argument]
 * 
 * @param d_temp_storage1 [contains NULL]
 * @param d_temp_storage2 [buffer array]
 * @param d_in [Input array]
 * @param d_flags [Input flags array]
 * @param d_out [Output array]
 * @param d_num_selected_out [number of values selected]
 * @param num_items [number of items]
 */
static void SelectMethodHelper(void *d_temp_storage1,unsigned *d_temp_storage2,unsigned *d_in,unsigned *d_flags,unsigned *d_out,int *d_num_selected_out,int num_items)
{
		size_t buffer_size = 0;

		unsigned temp;
		//Ist invocation to determine the buffer size required.
		CudaError(cub::DeviceSelect::Flagged(d_temp_storage1,buffer_size,d_in,d_flags,d_out,d_num_selected_out,num_items));

		void *buffer =  d_temp_storage2;

		//Actual invocation of the flag array.
		CudaError(cub::DeviceSelect::Flagged(buffer,buffer_size,d_in,d_flags,d_out,d_num_selected_out,num_items));

		CudaError(cudaDeviceSynchronize());

}

/**
 * @brief [This method is used to invoke the cub::RunLengthEncode method]
 * @details [We first determine the buffer size in the first invocation. We then pass either NULL or the buffer array as argument]
 * 
 * @param d_temp_storage1 [contains NULL]
 * @param d_temp_storage2 [contains Buffer Array]
 * @param d_in [Input Array]
 * @param d_uniqueOut [Output_array]
 * @param d_countsOut [Counts for each element]
 * @param d_num_selected_out [number of unique element]
 * @param num_items [number of input items]
 */
static void RunLengthEncodingHelper(void *d_temp_storage1,unsigned *d_temp_storage2,unsigned *d_in,unsigned *d_uniqueOut,unsigned *d_countsOut,int *d_num_runs_out,int num_items)
{
		size_t buffer_size = 0;

		unsigned temp;

		//Ist invocation of the method to determine the space required
		CudaError(cub::DeviceRunLengthEncode::Encode(d_temp_storage1,buffer_size,d_in,d_uniqueOut,d_countsOut,d_num_runs_out,num_items));

		void *buffer = d_temp_storage2;  //determing the buffer 

		//debug("buffer_size ", buffer_size);

		CudaError(cub::DeviceRunLengthEncode::Encode(buffer,buffer_size,d_in,d_uniqueOut,d_countsOut,d_num_runs_out,num_items));

}

/**
 * @brief [brief description]
 * @details [long description]
 * 
 * @param graph [description]
 * @param degreeThreshold [description]
 */
float KernelPruningWrapper(csr_graph &graph,int degreeThreshold)
{
		float totalTime = 0;

		GpuTimer gpuTimer;

		#define BLOCKSIZE 128    //BlockSize 
		#define GRIDSIZE(n) (ceil((double)n/BLOCKSIZE)) //GridSize 
		#define max(a,b) ((a < b)? b : a)
		//Declare the variables to hold the corresponding sizes of the auxilliary array.

		size_t aux_filter_rows;
		size_t aux_filter_run_length;

		void *d_temp_storage = NULL;

		bool toBePruned = true;

		bool flag_type; //flagtype is either 0 or 1. 0 indicates no vertices left, 1 indicates all remaining vertices are left

		#ifdef INFO
			int num_vertices_pruned = 0;
		#endif

		//Get the dimensions of the elements
		int count_rowoffset = graph.degree->size();
		int count_columns   = graph.rows->size();


		//Declare the corresponding containers for the thrust vectors.
		thrust::device_vector<unsigned> d_rows(graph.rows->begin(),graph.rows->end());
		thrust::device_vector<unsigned> d_cols(graph.columns->begin(),graph.columns->end());
		thrust::device_vector<unsigned> d_degree(graph.degree->begin(),graph.degree->end());

		thrust::device_vector<unsigned> d_filter(max(count_rowoffset,count_columns));
		thrust::device_vector<unsigned> d_uniqueOut(max(count_rowoffset,count_columns));
		thrust::device_vector<unsigned> d_countsOut(max(count_rowoffset,count_columns));

		thrust::device_vector<int> d_num_selected_out(1);

		//Get tempBufferSize for d_rowOffset
		CudaError(cub::DeviceSelect::Flagged(d_temp_storage,aux_filter_rows,CAST(d_cols),CAST(d_filter),CAST(d_cols),CAST(d_num_selected_out),count_columns));
		CudaError(cub::DeviceRunLengthEncode::Encode(d_temp_storage,aux_filter_run_length,CAST(d_rows),CAST(d_uniqueOut),CAST(d_countsOut),CAST(d_num_selected_out),count_columns));

		//Allocate buffer size
		size_t allocate_temporary_buffer_size = max(aux_filter_rows,aux_filter_run_length);

		//This device_vector is used as a buffer in the cub runtime calls.
		thrust::device_vector<unsigned> d_temp_buffer(allocate_temporary_buffer_size/4);

		//debug("max_size_allocated =",allocate_temporary_buffer_size , "bytes");
		//debug("Unique Vertexes count:", count_rowoffset);
		//debug("Columns count :" , count_columns);

		if(count_rowoffset == 0 && count_columns ==0)
			toBePruned = false;

		gpuTimer.Start();

		//Run until this flag is false.
		while(toBePruned)
		{

			//Step 1: Apply  the filtering Kernel to store valyues 1 and 0 in the d_filter array corresponding to the columns
			KernelFilterColOffset<<<GRIDSIZE(count_columns),BLOCKSIZE>>>(CAST(d_rows),CAST(d_cols),CAST(d_degree),CAST(d_filter),degreeThreshold,count_columns);

			CudaError(cudaDeviceSynchronize());

			//Step 2: Use the filter array to Select the Columns 
			SelectMethodHelper(d_temp_storage,CAST(d_temp_buffer),CAST(d_cols),CAST(d_filter),CAST(d_cols),CAST(d_num_selected_out), count_columns);

			//Step 3: Use the previous filter array to Select the rows too.
			SelectMethodHelper(d_temp_storage,CAST(d_temp_buffer),CAST(d_rows),CAST(d_filter),CAST(d_rows),CAST(d_num_selected_out), count_columns);


			//get the new count of the rows
			int new_columns_count = d_num_selected_out[0];

			if((new_columns_count == 0) || (new_columns_count == count_columns))
			{
				//debug("Finished : Final Columns Count",count_columns , " Final Unique Rows count = ", count_rowoffset);

				if(new_columns_count == 0)
				{
					flag_type = 0;
					//debug("All vertices are removed!!!");
				}
				else
				{
					flag_type = 1;
					//debug("Some vertices remain.")
				}

				break;
			}

			//Step 4: Do a runLengthEncoding on the Rows Array.
			RunLengthEncodingHelper(d_temp_storage,CAST(d_temp_buffer),CAST(d_rows),CAST(d_uniqueOut),CAST(d_countsOut),CAST(d_num_selected_out),new_columns_count);

			count_rowoffset = d_num_selected_out[0]; //update the count of rowoffset


			//Step 5: Update the degree of the degree Array based on the d_uniqueOut and d_countsOut
			KernelUpdateDegree<<<GRIDSIZE(new_columns_count),BLOCKSIZE>>>(CAST(d_degree),CAST(d_uniqueOut),CAST(d_countsOut),count_rowoffset);

			CudaError(cudaDeviceSynchronize());

			count_columns = new_columns_count;

		}

		totalTime += gpuTimer.StopGetTime();


		//Do something with the values here

		if(flag_type == 1)
		{
			//Make host_vectors from device vectors
			thrust::host_vector<unsigned> prunedRows(d_rows.begin(),d_rows.begin() + count_columns);
			thrust::host_vector<unsigned> prunedCols(d_cols.begin(),d_cols.begin() + count_columns);
			thrust::host_vector<unsigned> prunedDegree(d_degree.begin(),d_degree.end());

			//Clear the current host graph configuration
			graph.rows->clear();
			graph.columns->clear();
			graph.degree->clear();

			//Fill the pruned graph into the graph.(This will update the graph with the pruned version from the GPU)
			graph.rows->insert(graph.rows->end(),prunedRows.begin(),prunedRows.end());
			graph.columns->insert(graph.columns->end(),prunedCols.begin(),prunedCols.end());
			graph.degree->insert(graph.degree->end(),prunedDegree.begin(),prunedDegree.end());

			
			//Free the arrays

			prunedRows.clear();
			prunedCols.clear();
			prunedDegree.clear();

		}
		else
		{
			graph.rows->clear();
			graph.columns->clear();
			graph.degree->clear();
		}



		d_rows.clear();
		d_cols.clear();
		d_degree.clear();
		d_filter.clear();
		d_uniqueOut.clear();
		d_countsOut.clear();
		d_num_selected_out.clear();
		d_temp_buffer.clear();

		return totalTime;
		

}
