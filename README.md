# Strongly Conected Components in the Fault Tolerant Model

In this project, we study the problem of calculating the strongly connected components of a graph in the presence of failures. 

## Intro

There are several classical algorithms for computing the SCCs in O(m + n) time that are taught in any standard undergraduate algorithms course. One of these algorithms is Τarjan's approach of the problem. In this project we study the following natural variant of the problem in dynamic graphs.

The solution we implement is the one that **Surender**, **Keerti** and **Liam** proposed and it is based on the relation between strongly connected components (SCCs) and reachability. 

Different techniques are used, such as:
- the heavy path decomposition 
- and the creation of subgraphs that maintain the reachability between vertices under the presence of failures. 

More specifically, we consider a restricted variant of the problem in which we only compute strongly connected components that intersect a certain path. 

Restricting our attention to a path allows us to implicitly compute reachability between the path vertices and the rest of the graph in time that depends logarithmically rather than linearly in the size of the path. 

Therefore, we need to find an efficient way to represent the strongly connected components using paths.


## Tarjan's Algorithn for SCCs

If we change a bit Tarjan's algorithm to support k edge failures then the computational time becomes O(k*m + n). In a dense graph the maximum number of edges is n * (n - 1). That is because if you have n nodes, there are n - 1 directed edges than can lead from it (going to every other node). 

Therefore, the time complexity becomes <img src="https://latex.codecogs.com/svg.image?\bg_white&space;O(k*\frac{(n^2-n)}{2}&space;&plus;&space;n)" title="\bg_white O(k*\frac{(n^2-n)}{2} + n)" />, the question is can we achieve something better than this?

 


