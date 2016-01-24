#ifndef __KERNEL_WRAPPER
#define __KERNEL_WRAPPER

#include "CsrGraph.h"

float KernelPruningWrapper(csr_graph &graph,int degreeThreshold);

void DeviceReset();

#endif