This is the introduction document of proposed method for calculating survival signature under INID

It contains three MATLAB codes.

"example-Tsinghua.m" and "example-series-parallel-system.m" are codes for Tsinghua water supply system and series-parallel system in our paper.

"model.m" is the codes that generates the hydraulic model of Tsinghua water supply system. It contains 1366 components, while there are 261 pipes in mainstream pipeline system.

In "example-series-parallel-system.m":
Users can change the number of components (variable "comnum" in nodes), the failure probability of components (matrix "reliability" in nodes) and number of samples in weighted random sampling (variable "MCnum1" in nodes) at will. There are 4 steps in our divide-and-conquer-based method.1. calculation of division num and division unit; 2. probabilisty structure calculation of each minimum unit; 3. recursive for the origin system probability structure; 4.weighted random sampling from top to bottom.

In "example-Tsinghua.m":
The general process is consistent with the series parallel system code. Replace the variables in the above code with the parameters in the real problem."comnum=261, MCnum1=10000"