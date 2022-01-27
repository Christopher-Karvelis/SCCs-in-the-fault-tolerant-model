# Strongly Conected Components in the Fault Tolerant Model

In this project, we study the problem of calculating the strongly connected components of a graph in the presence of failures. 

## Intro

There are several classical algorithms for computing the SCCs in <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\bg_white&space;\inline&space;O(k*m&space;&plus;&space;n)" title="\bg_white \inline O(m + n)" /> time that are taught in any standard undergraduate algorithms course. One of these algorithms is Τarjan's approach of the problem. In this project we study the following natural variant of the problem in dynamic graphs.

## Tarjan's Algorithn for SCCs
Tarjanm's algorithm is based on simple DFS traversal of the graph hence the <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\bg_white&space;\inline&space;O(k*m&space;&plus;&space;n)" title="\bg_white \inline O(m + n)" /> time.

If we change a bit Tarjan's algorithm to support k edge failures then the computational time becomes <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\bg_white&space;\inline&space;O(k*m&space;&plus;&space;n)" title="\bg_white \inline O(k*m + n)" />. In a dense graph the maximum number of edges is n * (n - 1). That is because if you have n nodes, there are n - 1 directed edges than can lead from it (going to every other node). 

Therefore, the time complexity becomes <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\bg_white&space;\inline&space;O(k*\frac{(n^2-n)}{2}&space;&plus;&space;n)" title="\bg_white \inline O(k*\frac{(n^2-n)}{2} + n)" />, the question is can we achieve something better than this?

## Main Idea
The solution we implement is the one that **Surender**, **Keerti** and **Liam** proposed and it is based on the relation between strongly connected components (SCCs) and reachability. 

More specifically, we consider a restricted variant of the problem in which we only compute strongly connected components that intersect a certain path. 

Restricting our attention to a path allows us to implicitly compute reachability between the path vertices and the rest of the graph in time that depends logarithmically rather than linearly in the size of the path.

Different techniques are used, such as:
- the heavy path decomposition 
- and the creation of subgraphs that maintain the reachability between vertices under the presence of failures using known techniques such as flow algorithms. 

 
Therefore, we need to find an efficient way to represent the strongly connected components using paths.

 


