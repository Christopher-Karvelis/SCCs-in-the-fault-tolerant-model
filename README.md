# Strongly Conected Components in the Fault Tolerant Model

## Table of Contents
1.  [Purpose](https://github.com/Christopher-Karvelis/SCC-LAST#purpose)

2.  [Preliminaries](https://github.com/Christopher-Karvelis/SCC-LAST#preliminaries)

3.  [Tarjan's Algorithn for SCCs](https://github.com/Christopher-Karvelis/SCC-LAST#tarjans-algorithn-for-sccs)

4.  [Main Idea](https://github.com/Christopher-Karvelis/SCC-LAST#main-idea)

5.  [Structure of the Project](https://github.com/Christopher-Karvelis/SCC-LAST#structure-of-the-project)

6.  [Heavy Path Decomposition](https://github.com/Christopher-Karvelis/SCC-LAST#heavy-path-decomposition)

7.  [k-Fault Tolerant Reachability Subgraph (k-FTRS)](https://github.com/Christopher-Karvelis/SCC-LAST#k-fault-tolerant-reachability-subgraph-k-ftrs)

8.  [Computation of SCCs Intersecting a Given Path](https://github.com/Christopher-Karvelis/SCC-LAST#computation-of-SCCs-intersecting-a-given-path)

9.  [Main Algorithm](https://github.com/Christopher-Karvelis/SCC-LAST#main-algorithm)

10. [Examples](https://github.com/Christopher-Karvelis/SCC-LAST#examples)

## Purpose
In this project, we study the problem of calculating the strongly connected components of a graph in the presence of failures. It is important to study problems like this as we come across more and more networks in real life like the internet, powerline and telephone networks, etc. This rise in demand for networks in our daily lives has spurred the  need  for tools thatcan handdle any issues occurred when working with them.

Fault Tolerant Subgraph for Single Source Reachability: Generic and Optimal(Baswana S., Choudhary K., and Roditty L., 2018):
> Networks in most real life applications are prone to failures. Such a network can be modeled as a graph
where edges may change their status from active to failed, and vice versa. These failures,
though unpredictable, are small in numbers and are transient due to some simultaneous repair process that is
undertaken in these applications. This aspect can be captured by associating a parameter k with the network
such that there are at most k edges that are failed at any stage, where k is much smaller than the
number of vertices in the underlying graph. 

![fault tolerance](https://user-images.githubusercontent.com/25777650/151781733-a97d70fc-8bde-4ce5-86c8-429ec5511f82.png)

There are several classical algorithms for computing the SCCs in <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\bg_white&space;\inline&space;O(m&space;&plus;&space;n)" title="\bg_white \inline O(m + n)" /> time.

One of these algorithms is Τarjan's approach of the problem. In this project we study the following natural variant of the problem in dynamic graphs.



## Preliminaries
- Disclaimers:
  - This is a project that I worked on as my diploma thesis.
  - In this README I provide some basic information to get a better idea of the project.
  - **For more information about the theory behind this algorithm please read the pappers:**
       1. Pappers/1. An Efficient Strongly Connected Components.pdf
       2. Pappers/2. Fault Tolerant Subgraph for Single Source Reachability.pdf
- **Assumption:** The out-degree of all vertices in the input graph G is at most two. If that's not the case then the input graph can be transformed to a graph where the out-degree of all vertices is at most two as described in page 15 of "2. Fault Tolerant Subgraph for Single Source Reachability.pdf".

- The algorithm for computing SCCs in a fault tolerant environment crucially uses the concept of a k-fault tolerant reachability subgraph (k-FTRS) which is a sparse
subgraph that preserves reachability from a given source vertex even after the failure of at most k edges in G.

- In this project we construct a k-FTRS with respect to edge failures only. Vertex failures can be handled by simply splitting a vertex v into an edge (vin, vout), where the incoming and outgoing edges of v are respectively directed into vin and directed out of vout.
 

 
## Tarjan's Algorithn for SCCs
Tarjanm's algorithm is based on simple DFS traversal of the graph hence the <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\bg_white&space;\inline&space;O(km&space;&plus;&space;n)" title="\bg_white \inline O(m + n)" /> time.

If we change a bit Tarjan's algorithm to support k edge failures then the computational time becomes <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\bg_white&space;\inline&space;O(km&space;&plus;&space;n)" title="\bg_white \inline O(km + n)" />. In a dense graph the maximum number of edges is n * (n - 1). That is because if you have n nodes, there are n - 1 directed edges than can lead from it (going to every other node). 

Therefore, the time complexity becomes <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\bg_white&space;\inline&space;O(k\frac{(n^2-n)}{2}&space;&plus;&space;n)" title="\bg_white \inline O(k\frac{(n^2-n)}{2} + n)" />, the question is can we achieve something better than this?



## Main Idea
The solution we implement is the one that **Surender**, **Keerti** and **Liam** proposed in "An Efficient Strongly Connected Components
Algorithm in the Fault Tolerant Mode" and it is based on the relation between strongly connected components (SCCs) and reachability. 

More specifically they showed that there is an algorithm that computes the SCCs of G\F, for any set F of k edges or vertices, in  <img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\bg_white&space;\inline&space;O(2^knlog^{2}n)" title="\bg_white \inline O(2^knlog^{2}n)" /> time. The algorithm uses a data structure of size O(2^kn^2) computed in O(2^kn^2m) time for G during a preprocessing phase.

The solution is based on a restricted variant of the problem in which we only compute strongly connected components that intersect a certain path. 

By doing this we implicitly compute reachability between the path vertices and the rest of the graph in time that depends logarithmically rather than linearly in the size of the path.

Therefore, we need an efficient way to represent the strongly connected components using paths.

Different techniques are used, such as:
- the heavy path decomposition for creating the paths on which we will compute the intersecting SCCs.
- and the creation of subgraphs that maintain the reachability between vertices under the presence of failures using known techniques such as flow algorithms. 





## Structure of the Project
The algorithm for computing the Strongly Connected Components of a given strongly connected graph G consists of two phases:
1. a preprocessing phase where:
      - We contsruct a heavy path ecomposition P of G by inspecting the DFS tree of G.
      - For every path in P we compute a k-Fault Tolerant Reachability Subgraph(k-FTRS) of every node on the path as source.
2. and the computation of SCCs phase where:
      -  Given k failed egdes we compute the final paths of P.
      -  We computate the SCCs intersecting a given path of P using the k-FTRSs computed in the preprocessing phase.
      -  Combining all the above we return the SCCs of the G under k failed edges. 

![structure](https://user-images.githubusercontent.com/25777650/151577246-cbf2474f-e722-438e-815c-292cd146698a.png)



## Heavy Path Decomposition

Designed by Sleator and Tarjan Given any rooted tree T, this decomposition splits T into a set P of vertex disjoint paths.

Any path from the root to a leaf node in T can be expressed as a concatenation of at most 1 + log n subpaths of paths in P. 

This decomposition is carried out as follows:
- Starting from the root, we follow the path downward choosing the next node whose subtree is of maximum size.
- We terminate upon reaching a leaf node 
- Let P be the path obtained in this manner. If we remove P from T, we are left with a collection of subtrees each of size at most n/2. 
- Each of these trees hangs from P through an edge in T. 
- We carry out the decomposition of these trees recursively.
 
![hld](https://user-images.githubusercontent.com/25777650/151389548-f858166e-8379-42cc-9152-79525b0531bc.png)



## k-Fault Tolerant Reachability Subgraph (k-FTRS)

(k-FTRS): A sparse subgraph that preserves the reachability from a given fixed source s even after k failures. 
So a vertex v is reachable from s in the original graph in the preif and only if it is reachable from the s in the subgraph.

The algorithm for computing a k-FTRS involves the concepts of max-flow, min-cut, and edge disjoint paths. So we will visualize the same graph G as a network with unit edge capacities. 


### Main Idea of k-FTRS
The following well known result from Graph theory shows the connection between max-flow and edge-disjoint paths.

For any positive integer α, there is a flow from a source set S to a destination vertex t of value α if and only if there are α edge disjoint paths originating from set S and terminating at t.

So if we have α disjoint paths originating from a node s and terminating at t, then we have α different ways to reach t starting from s.



### Farthest Min-Cut (FMC)

Ford and Fullkerson gave an algorithm for constructing the farthest (S, t)-min-cut and also established its uniqueness.

A set of edges is said to be an (S, t)-cut if each path from any s ∈ S to t must pass through at least one edge from C. The size of a cut C is the number of edges present in C. An (S, t)-cut of smallest size is called (S, t)-min-cut.

Let S be a source set and t be a destination vertex. Any (S, t)-min-cut C partitions the vertex set into two sets: A(C) containing S, and B(C) containing t. An (S, t)-min-cut C
is said to be the farthest min-cut if A(C ) ) A(C) for any (S, t)-min-cut C other than C. We denote C∗ with FMC(G, S, t)

![fig1](https://user-images.githubusercontent.com/25777650/151565804-c8f699b3-8997-41a3-a879-f0210768bc7e.png)



### Computing k-FTRS(t)

**Definition:** Given a vertex v ∈ V and an integer k ≥ 1, a subgraph G0 = (V, E0), E0 ⊆ E is said to be k-FTRS(v) if for any set F of k edge failures, the following condition holds: v is reachable from s in G\F if and only if v is reachable from s in G0\F

The algorithm implemented runs in O(2kmn) time and for any given integer k ≥ 1, and any given directed graph G on n vertices, m edges and a designated source vertex s, computes a k-FTRS for G with at most 2 kn edges. Moreover, the in-degree of each vertex in this k-FTRS is bounded by 2 k.

![kftrs_alg](https://user-images.githubusercontent.com/25777650/151420102-046bbf84-4fa3-4092-b360-0bf5be608a1c.png) 

![σδφ](https://user-images.githubusercontent.com/25777650/151563009-b7884b01-8eb1-4677-bea7-499b3a8030ec.png)



## Computation of SCCs Intersecting a Given Path

For each v ∈ V, let Xin(v) be the vertex of X of minimum index (if exists) that is reachable from v in G\F. Similarly, let Xout(v) be the vertex of X of maximum index (if exists) that has a path to v in G\F.

![xinsxouts](https://user-images.githubusercontent.com/25777650/151846120-06974b24-ac6b-4a68-8d05-83ebcca8be1d.png)

For any vertex w ∈ V , the SCC that contains w in G\F intersects X if
and only if the following two conditions are satisfied.
(i) Both Xin(w) and Xout(w) are defined, and
(ii) Either Xin(w) = Xout(w), or Xin(w) appears before Xout(w) on X.

So when two vertices belong to the sanme SCC?

Let a, b be any two vertices in V whose SCCs intersect X. Then a and b
lie in the same SCC if and only if Xin(a) = Xin(b) and Xout(a) = Xout(b).



### Compute Xin and Xout
In order to compute the SCCs in G\F that intersect with X, it suffices to compute Xin(·) and Xout(·) for all vertices in V.

It suffices to focus on computation of Xout(·) for all the vertices of V, since Xin(·) can be computed in an analogous manner by just looking at graph G reversed.

One trivial approach to achieve this goal is to compute the set Vi consisting of all vertices reachable from each xi by performing a BFS or DFS traversal of graph G(xi)\F, for 1 ≤ i ≤ t = |X|.

Using this straightforward approach it takes O(2knt) time to complete the task of computing Xout(v) for every v ∈ V, while our target is to do so in O(2kn log n) time.

Observe the nested structure underlying Vi’s, that is, V1 ⊇ V2 ⊇···⊇ Vt . Consider any vertex xl, 1 <l< t. The nested structure implies that for every v ∈ Vt, Xout(v)
must be on the portion (xl,..., xt) of X.

Similarly, it implies that for every v ∈ V1\Vl, Xout(v)must be on the portion (x1,..., xl−1) of X. 

This suggests a divide and conquer approach to efficiently compute Xout(·).

We first compute the sets V1 and Vt in O(2kn) time each.

For each v ∈ V\V1, we assign NULL to Xout(v) as it is not reachable from any vertex on X; and for each v ∈ Vt we set Xout(v) to xt . For vertices in set V1 \ Vt , Xout(·) is computed by calling the function Binary-Search(1, t − 1, V1\Vt).

![bfs](https://user-images.githubusercontent.com/25777650/151580197-9e5e803c-2e18-43ef-9d0f-a88630fa4fea.png)



### Implementation of Function Reach
It suffices to do a traversal from xmid in the graph GA, the induced subgraph of A in G(x)\F, that has O(2k |A|) edges.

![reach](https://user-images.githubusercontent.com/25777650/151580955-e0142471-8f91-4be1-88c3-08e6782b8ba4.png)



## Main Algorithm

-Let C denote the collection of SCCs in G\F initialized to ∅.

-We process the paths from P in non-decreasing order of their depths.

-Let Path(a, b) be any path in P and let A be the set of vertices belonging to T (a).

-We use the data structure Da,b to compute SCCs of G(A)\F intersecting Path(a, b). Let these be S1,..., St.

-Note that some of these SCCs might be a part of some bigger SCC computed earlier.

-We can detect it by keeping a set W of all vertices for which we have computed their SCCs.

-So if Si ⊆ W, then we can discard Si , else we add Si to collection C.


![scc_algo](https://user-images.githubusercontent.com/25777650/151698368-256f4f9f-f25e-454a-b54b-45cf63e3323a.png)




## Examples
coming soon!
