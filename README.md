# SCC-LAST
In this project, we study the problem of calculating the strongly connected components
of a graph in the presence of failures. The solution we implement is the one that
Surender, Keerti and Liam proposed and it is based on the relation between
strongly connected components (SCCs) and reachability. The solution under study
uses different techniques, such as the heavy path decomposition and the creation
of subgraphs that maintain the reachability between vertices under the presence of
failures. More specifically, we consider a restricted variant of the problem in which we
only compute strongly connected components that intersect a certain path. Restricting
our attention to a path allows us to implicitly compute reachability between the path
vertices and the rest of the graph in time that depends logarithmically rather than
linearly in the size of the path. Therefore, we need to find an efficient way to represent
the strongly connected components using paths.
