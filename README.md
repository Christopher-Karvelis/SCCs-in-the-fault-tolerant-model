# Strongly Conected Components in the Fault Tolerant Model

## Intro
Networks in most real life applications are prone to failures. Such a network can be modeled as a graph
where edges may change their status from active to failed, and vice versa. These failures,
though unpredictable, are small in numbers and are transient due to some simultaneous repair process that is
undertaken in these applications. This aspect can be captured by associating a parameter k with the network
such that there are at most k edges that are failed at any stage, where k is much smaller than the
number of vertices in the underlying graph.

In this project, we study the problem of calculating the strongly connected components of a graph in the presence of failures. 

There are several classical algorithms for computing the SCCs in <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\bg_white&space;\inline&space;O(m&space;&plus;&space;n)" title="\bg_white \inline O(m + n)" /> time that are taught in any standard undergraduate algorithms course. One of these algorithms is Τarjan's approach of the problem. In this project we study the following natural variant of the problem in dynamic graphs.

## Tarjan's Algorithn for SCCs
Tarjanm's algorithm is based on simple DFS traversal of the graph hence the <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\bg_white&space;\inline&space;O(km&space;&plus;&space;n)" title="\bg_white \inline O(m + n)" /> time.

If we change a bit Tarjan's algorithm to support k edge failures then the computational time becomes <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\bg_white&space;\inline&space;O(km&space;&plus;&space;n)" title="\bg_white \inline O(km + n)" />. In a dense graph the maximum number of edges is n * (n - 1). That is because if you have n nodes, there are n - 1 directed edges than can lead from it (going to every other node). 

Therefore, the time complexity becomes <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\bg_white&space;\inline&space;O(k\frac{(n^2-n)}{2}&space;&plus;&space;n)" title="\bg_white \inline O(k\frac{(n^2-n)}{2} + n)" />, the question is can we achieve something better than this?

## Main Idea
The solution we implement is the one that **Surender**, **Keerti** and **Liam** proposed and it is based on the relation between strongly connected components (SCCs) and reachability. 

More specifically they showed that there is an algorithm that computes the SCCs of G\F, for any set F of k edges or vertices, in  <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\bg_white&space;\inline&space;O(2^knlog^{2}n)" title="\bg_white \inline O(2^knlog^{2}n)" /> time. The algorithm uses a data structure of size O(2^kn^2) computed in O(2^kn^2m) time for G during a preprocessing phase.

The solution is based on a restricted variant of the problem in which we only compute strongly connected components that intersect a certain path. 

By doing this we implicitly compute reachability between the path vertices and the rest of the graph in time that depends logarithmically rather than linearly in the size of the path.

Different techniques are used, such as:
- the heavy path decomposition for creating the paths on which we will compute the intersecting SCCs.
- and the creation of subgraphs that maintain the reachability between vertices under the presence of failures using known techniques such as flow algorithms. 

 
Therefore, we need to find an efficient way to represent the strongly connected components using paths.

## Heavy Path Decomposition

Designed by Sleator and Tarjan Given any rooted tree T , this decomposition splits T into a set P of vertex disjoint paths.

Any path from the root to a leaf node in T can be expressed as a concatenation of at most 1 + log n subpaths of paths in P. 

This decomposition is carried out as follows:
- Starting from the root, we follow the path downward choosing the next node whose subtree is of maximum size.
- We terminate upon reaching a leaf node 
- Let P be the path obtained in this manner. If we remove P from T, we are left with a collection of subtrees each of size at most n/2. 
- Each of these trees hangs from P through an edge in T. 
- We carry out the decomposition of these trees recursively.
 
![hld](https://user-images.githubusercontent.com/25777650/151389548-f858166e-8379-42cc-9152-79525b0531bc.png)


