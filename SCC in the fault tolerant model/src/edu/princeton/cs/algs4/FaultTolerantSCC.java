package edu.princeton.cs.algs4;
/******************************************************************************
 *  Compilation:  javac FaultTolerantSCC.java
 *  Execution:    Java FaultTolerantSCC V E
 *  Dependencies: edu.princeton.cs.algs4.Digraph.java StdOut.java
 *  Data files:   https://algs4.cs.princeton.edu/42digraph/tinyDG.txt
 *                https://algs4.cs.princeton.edu/42digraph/mediumDG.txt
 *                https://algs4.cs.princeton.edu/42digraph/largeDG.txt
 *
 *  Compute the strongly-connected components of a digraph
 *
 *  % java FaultTolerantSCC tinyDG.txt
 *  5 components
 *  1
 *  0 2 3 4 5
 *  9 10 11 12
 *  6 8
 *  7
 *
 ******************************************************************************/

import java.util.*;

public class FaultTolerantSCC {
    private int Vsize;
    private int pathSize;
    private int[] path;
    private int[] Xins;
    private int[] Xouts;
    private int[] subtreeOfTSize;
    private int[] subtreeOfTrSize;
    private boolean[] marked;                                         // marked[v] = has v been visited?
    private Digraph dfsTreeG;                                         // DFS Tree of G rooted at vertex 0
    private Digraph dfsTreeGr;                                        // DFS Tree of Gr rooted at vertex 0
    private ArrayList<int[]> P;                                       //heavy path decomposition of T, where paths are sorted in non-decreasing order of their depths
    private ArrayList<int[]> failedEdges;
    private HashSet<?>[][] InEdgesG;                                       //set of all incoming edges to v in graph G(for now, in the future G-> g(mid)\F)
    private HashSet<?>[][] InEdgesGr;                                      //set of all incoming edges to v in graph Gr(for now, in the future Gr-> g(mid)\F)
    private Map<Integer, HashSet<Integer>> subtreeVerticesT;          //vertices that belong in the subtree of v rooted at vertex v in T
    private Map<Integer, HashSet<Integer>> subtreeVerticesTr;         //vertices that belong in the subtree of v rooted at vertex v in Tr
    private FordFulkerson maxflow;
    private Map<Integer, Map<Integer, DynamicDigraph>> kFTRSs;
    private Map<Integer, Map<Integer, DynamicDigraph>> kFTRSsR;
    private Digraph[] kFTRSes;                                        //Contains the kFTRS(vi) for each vi as source
    private Digraph[] kFTRSesR;                                        //Contains the kFTRS(vi) for each vi as source
    private FlowNetwork Gf;
    private Digraph Gglob;
    /**
     * Computes the strong components of the digraph {@code G}.
     * @param G the digraph
     */
    public FaultTolerantSCC(Digraph G, Digraph Gr, ArrayList<int[]> failedEdges, int k){
        Vsize = G.V();
        maxflow = null;
        Gglob = Gr.reverse();
        this.failedEdges = failedEdges;
        P = new ArrayList<int[]>();
        Gf = new FlowNetwork(G.V());
        for (int v = 0; v < G.V(); v++) {
            for (int w : G.adj(v)) {
                FlowEdge e = new FlowEdge(v, w, 1);
                Gf.addEdge(e);
            }
        }
    }

    public void preprocessing(Digraph G, Digraph Gr, int k){
        System.out.println("--------------------------------------------------- Preprocessing  --------------------------------------------------");
        marked = new boolean[G.V()];
        kFTRSes = new Digraph[G.V()];
        kFTRSesR = new Digraph[G.V()];
        dfsTreeG = new Digraph(G.V());
        dfsTreeGr = new Digraph(G.V());
        subtreeOfTSize = new int[G.V()];
        subtreeOfTrSize = new int[G.V()];
        InEdgesG = new HashSet<?>[G.V()][G.V()];
        InEdgesGr = new HashSet<?>[G.V()][G.V()];
        kFTRSs = new HashMap<Integer, Map<Integer, DynamicDigraph>>();
        kFTRSsR = new HashMap<Integer, Map<Integer, DynamicDigraph>>();
        subtreeVerticesT = new HashMap<Integer, HashSet<Integer>>();
        subtreeVerticesTr = new HashMap<Integer, HashSet<Integer>>();

        //Calculate DFS tree T of G
        System.out.println("-Calculating a DFS tree of G...");
        for (int v = 0; v < G.V(); v++) {
            if (!marked[v]) dfs(G, v, dfsTreeG, subtreeOfTSize);
        }

        //Calculate DFS tree Tr of Gr
        marked = new boolean[Gr.V()];
        for (int v = 0; v < G.V(); v++) {
            if (!marked[v]) dfs(G, v, dfsTreeGr, subtreeOfTrSize);
        }

        /* print DFS tree of G
        System.out.println(dfsTreeG.toString());
        for(int i = 0; i < G.V(); i++){
            System.out.println("Node: " + i + ", |T("+ i + ")| = " + subtreeOfTSize[i]);
        }
        */

        //Calculate heavy path decomposition of G
        /*
        marked = new boolean[G.V()];
        pathSize = 0;
        int startingNode  = ThreadLocalRandom.current().nextInt(0, G.V());
        heavyPathDecomposition(dfsTreeG, startingNode, pathSize);
        path[pathSize] = startingNode;
        P.add(path);
        System.out.print("Array Path = [");
        for (int x: path) {
            System.out.print(x + ", ");
        }
        System.out.println("]");
         */
        marked = new boolean[G.V()];
        heavyPathDecomposition(G, dfsTreeG);
        /*
        marked = new boolean[G.V()];
        for (int v = 0; v < G.V(); v++) {
            if (!marked[v]){
                pathSize = 0;
                heavyPathDecomposition(dfsTreeG, v, pathSize);
                path[pathSize] = v;

                System.out.print("Array Path = [");
                for (int x: path) {
                    System.out.print(x + ", ");
                }
                System.out.println("]");


                P.add(path);
            }
        }
        */

        //Calculate Subtrees Of DFS tree T

        for (int v = 0; v < G.V(); v++) {
            marked = new boolean[G.V()];
            subtreeVerticesT.put(v, new HashSet<Integer>());
            int subPointer = 0;
            calculateSubtreesOfT(dfsTreeG, v, subtreeVerticesT.get(v));
        }
        //Calculate Subtrees Of DFS tree Tr
        marked = new boolean[G.V()];
        for (int v = 0; v < Gr.V(); v++) {
            if(!marked[v]){
                subtreeVerticesTr.put(v, new HashSet<Integer>());
                int subPointer = 0;
                calculateSubtreesOfT(dfsTreeGr, v, subtreeVerticesTr.get(v));
            }

        }
        //Calculate Data Structures of Theorem 2
        //Calcalculate the kFTRSes for each v in G
        calculateKFTRSvi(G, Gr, InEdgesG, kFTRSes , k);
        //Calcalculate the kFTRSes for each v in Gr
        //calculateKFTRSvi(Gr, InEdgesGr, kFTRSesR , k);
        System.out.println("------------------------------------------------- Preprocessing Done ------------------------------------------------\n\n");
    }

    private void dfs(Digraph G, int v, Digraph dfsTree, int[] subtreeSize) {
        subtreeSize[v]++;
        marked[v] = true;
        for (int w : G.adj(v)) {
            if (!marked[w]){
                dfsTree.addEdge(v, w);                      //Add (v,w) to DFS tree
                dfs(G, w, dfsTree, subtreeSize);
                subtreeSize[v] += subtreeSize[w];
            }
        }
    }

    private void calculateSubtreesOfT(Digraph T, int v, HashSet<Integer>  subtreeVertices){
        subtreeVertices.add(v);
        marked[v] = true;
        for (int w : T.adj(v)) {
            if (!marked[w]){
                calculateSubtreesOfT(T, w, subtreeVertices);
            }
        }
    }

    private ArrayList<int[]> calculateFinalHeavyPathDecomposition() {
        ArrayList<ArrayList<Integer>> Pnew = new ArrayList<ArrayList<Integer>>();
        for(int[] Pab: P){
            int a = Pab[0];
            int previous = -1;
            int current = -1;
            ArrayList<Integer> newPath = new ArrayList<Integer>();
            for(int v: Pab){
                current = v;
                if(previous != -1){
                    if(isFailedEdge(previous, current)){
                        Pnew.add(newPath);
                        newPath = new ArrayList<>();
                    }
                }
                newPath.add(v);
                previous = v;
            }
            Pnew.add(newPath);
        }
        ArrayList<int[]> Paths = new ArrayList<>();

        for (ArrayList<Integer> path: Pnew) {
            int[] pathNew = new int[path.size()];
            int i = 0;
            for(int v: path){
                pathNew[i] = v;
                i++;
            }
            for(int[] Pab: P){
                int a = Pab[0];
                for(int v: Pab){
                    if(pathNew[0] == v){
                        if(!kFTRSs.containsKey(v)){
                            kFTRSs.put(v, kFTRSs.get(a));
                            kFTRSsR.put(v, kFTRSsR.get(a));
                        }
                    }
                }
            }
            //System.out.println("NEW START ======= " + pathNew[0]);
            Paths.add(pathNew);
        }
        return Paths;
    }

    private  void heavyPathDecomposition(Digraph G, Digraph T){
        System.out.println("-Calculating Heavy Path Decomposition of G...");
        int[] heavyD = new int[G.V()];
        Arrays.fill(heavyD, -1);
        for (int v = 0; v < G.V(); v++) {
            int heaviestChild = -1;
            int heaviestSize = -1;
            for(int w : T.adj(v)){
                if(subtreeOfTSize[w] >= heaviestSize){
                    heaviestChild = w;
                    heaviestSize = subtreeOfTSize[w];
                }
            }
            heavyD[v] = heaviestChild;
        }

        for (int v = 0; v < G.V(); v++) {
            int pathSize = 0;

            if(heavyD[v] != -1 && !marked[v]) {
                marked[v] = true;
                recursivePathBuilding(heavyD, v, pathSize);
                path[pathSize] = v;
                /*
                System.out.print("Array Path = [");
                for (int x: path) {
                    System.out.print(x + ", ");
                }
                System.out.println("]");
                */
                P.add(path);
            }else if(!marked[v] && heavyD[v] == -1){
                //System.out.println("Array Path = ["+v+"]");
                int [] p = {v};
                P.add(p);
            }
        }
    }

    private void recursivePathBuilding(int[] heavyD, int v, int pathSize){
        marked[v] = true;
        if(heavyD[v] != -1 && !marked[heavyD[v]]){
            pathSize++;
            int next = heavyD[v];
            heavyD[v] = -1;
            recursivePathBuilding(heavyD, next, pathSize);
            path[pathSize] = next;
        }else{
            pathSize++;
            path = new int[pathSize];
        }
    }

    private void calculateKFTRSvi(Digraph G, Digraph Gr, HashSet[][] InEdges, Digraph[] kFTRSes,  int k){
        System.out.print("-Preprocessing paths...");
        for(int[] Pab: P){
            int a = Pab[0];
            System.out.println("\n..................................................... Preprocessing path : " + Arrays.toString(Pab) + "\n");
            kFTRSs.put(a, new HashMap<Integer, DynamicDigraph>());
            kFTRSsR.put(a, new HashMap<Integer, DynamicDigraph>());
            if(subtreeVerticesT.get(a).size() == 0){
                subtreeVerticesT.get(a).add(a);
            }
            System.out.println("-Creating subgraph of G induced by vertices in T(" + a +")...");
            DynamicDigraph Ga = new DynamicDigraph(subtreeVerticesT.get(a));
            for(int v: Ga.vertices()) {
                for (int w : G.adj(v)) {
                    Ga.addEdge(v, w);
                }
            }

            DynamicDigraph Gar = new DynamicDigraph(subtreeVerticesT.get(a));
            for(int v: Ga.vertices()) {
                for (int w : Gr.adj(v)) {
                    Gar.addEdge(v, w);
                }
            }
            for(int v: Pab) {
                System.out.println("-Creating " + k + "-FTRS with s = " + v);
                DynamicDigraph kFTRS = Ga;
                for (int t: Ga.vertices()) {
                    if(v != t){
                        kFTRS = kFTRS(kFTRS, InEdges, k, v, t);
                    }
                }
                System.out.println("---> " + k + "-FTRS with s = " + v + " has " + kFTRS.E() + "(< 2^k*n = " + "2^"+ k + "*"+ kFTRS.vertices().size() + " = " + Math.pow(2, k)*(kFTRS.vertices().size()) + ") edges <---");
                //System.out.println(kFTRS.toString());
                kFTRSs.get(a).put(v, kFTRS);

                System.out.println("-Creating " + k + "-FTRS^R with s = " + v);
                DynamicDigraph kFTRSr = Gar;
                for (int t: Gar.vertices()) {
                    if(v != t){
                        kFTRSr = kFTRS(kFTRSr, InEdgesGr, k, v, t);
                    }
                }
                System.out.println("---> " + k + "-FTRS^R with s = " + v + " has " + kFTRSr.E() + "(< 2^k*n = " + "2^"+ k + "*"+ kFTRSr.vertices().size() + " = " + Math.pow(2, k)*(kFTRSr.vertices().size()) + ") edges <---");
                //System.out.println(kFTRSr.toString());
                kFTRSsR.get(a).put(v, kFTRSr);
            }
        }
        /*
        for(int v = 0; v < G.V(); v++){                 //for each node v as source calculate the coresponding kFTRS
            System.out.println("=====================PROCESSING " + k + "-FTRS with s = " + v + " =====================");
            Digraph Gi = G;
            for (int t = 0; t < G.V(); t++) {           //Theorem
                if(v != t){
                    //System.out.println("===================== PROCESSING " + k + "-FTRS with s = " + v + " and t = " + t + "=====================" );
                    Gi = kFTRS(Gi, InEdges, k, v, t);
                    //System.out.println(Gi.toString());
                }
            }
            kFTRSes[v] = Gi;
            System.out.println(k + "-FTRS with s = " + v + " has " + kFTRSes[v].E() + "(< 2^k*n = " + "2^"+ k + "*"+ Gi.V() + " = " + Math.pow(2, k)*Gi.V() + ") edges ");
            //System.out.println(Gi.toString());
        }
         */
        //System.out.println(kFTRSes[1].toString());
    }

    //Algorithm 1
    private DynamicDigraph kFTRS(DynamicDigraph G, HashSet[][] InEdges, int k, int s, int t){
        HashSet<Integer> Si = new HashSet<Integer>();
        Si.add(s);

        /*PRINTPRINTPRINTPRINTPRINTPRINTPRINT
        System.out.print("Si = [");
        for (int v: Si) {
            System.out.print(v + ", ");
        }
        System.out.println("]");

         */
        //PRINTPRINTPRINTPRINTPRINTPRINTPRINT

        for (int i = 1; i <= k; i++) {
            //System.out.println("Si = " + Si);
            HashSet[] FMC = calculateFMC(G, Si, t);
            HashSet<int []> Ci = new HashSet<int[]>();
            Ci = FMC[0];

            //PRINTPRINTPRINTPRINTPRINTPRINTPRINT
            /*
            System.out.print("Ci =[ ");
            for (int[] x: Ci) {
                System.out.print("[ ");
               for (int w: x){
                   System.out.print(w+ ", ");
               }
                System.out.print("]");
            }
            System.out.println("]");

             */
            //PRINTPRINTPRINTPRINTPRINTPRINTPRINT

            HashSet<Integer> Ai = new HashSet<Integer>();
            Ai = FMC[1];

            //PRINTPRINTPRINTPRINTPRINTPRINTPRINT
            /*
            System.out.print("Ai =[ ");
            for (int x: Ai) {
                System.out.print(x+ ", ");
            }
            System.out.println("]");
            //PRINTPRINTPRINTPRINTPRINTPRINTPRINT

             */

            HashSet<Integer> Bi = new HashSet<Integer>();
            Bi = FMC[2];

            //PRINTPRINTPRINTPRINTPRINTPRINTPRINT
            /*
            System.out.print("Bi =[ ");
            for (int x: Bi) {
                System.out.print(x+ ", ");
            }
            System.out.println("]");

             */
            //PRINTPRINTPRINTPRINTPRINTPRINTPRINT

            //S(i+1) <- (Ai U OUT(Ai))\{t}
            for (int v: Ai) {
                if(G.vertices().contains(s) && s != t) {
                    Si.add(v);
                    for (int w : G.adj(G.vertices().indexOf(v))) {
                        if (w != t) {
                            Si.add(w);
                        }
                    }
                }
            }
            /*
            System.out.print("Si = [");
            for (int v: Si) {
                System.out.print(v + ", ");
            }
            System.out.println("]");

             */
        }
        /*
        Gf = new FlowNetwork(G.V());
        for (int v = 0; v < G.V(); v++) {
            for (int w : G.adj(v)) {
                FlowEdge e = new FlowEdge(v, w, 1);
                Gf.addEdge(e);
            }
        }
        */
        int fakeNode = Vsize;
        for (int node: Si) {
            if(G.vertices().contains(s) && s != t) {
                FlowEdge e = new FlowEdge(fakeNode, node, 1);
                Gf.addEdge(e);
            }
        }
        maxflow = new FordFulkerson(Gf, fakeNode, t);
        //System.out.println(maxflow.value());
        //System.out.println(Gf.toString());
        HashSet<Integer> InEdgesT = new HashSet<Integer>();
        HashSet<Integer> DeletedEdgesT = new HashSet<Integer>();
        for (int v: G.vertices()){
            for(FlowEdge edge: Gf.adj(v)){
                //System.out.println(edge.toString());
                if(v == edge.from() && edge.to() == t && edge.flow() > 0){
                    InEdgesT.add(v);
                    //System.out.println("Added to InEdges --> " + v );
                }else if(v == edge.from() && edge.to() == t){
                    DeletedEdgesT.add(v);
                }
            }
        }
        InEdges[s][t] = InEdgesT;
        DynamicDigraph kFTRS = new DynamicDigraph(G.vertices());
        for(int v: G.vertices()){
            for(int w: G.adj(G.vertices().indexOf(v))){
                if(w != t && w != s){
                    kFTRS.addEdge(v, w);
                }
            }
        }
        for (int v: InEdgesT) {
            //System.out.println("Edge " + v + " -> " + t + " added!");
            kFTRS.addEdge(v, t);
        }
        /*
        System.out.print("Deleted Edges =[ ");
        for (Integer e: DeletedEdgesT) {
            System.out.print("[" + e + ", " + t + "], ");
        }
        System.out.println("]");

         */
        return kFTRS;
    }

    private HashSet[] calculateFMC(DynamicDigraph G, HashSet<Integer> Si, int t){
        //----------------------------Calculate max-flow, residual graph Gf and FMC-------------------------------------
        int fakeNode = Vsize;
        Gf = new FlowNetwork(Vsize+1);
        for (int v : G.vertices()) {
            for (int w : G.adj(G.vertices().indexOf(v))) {
                FlowEdge e = new FlowEdge(v, w, 1);
                Gf.addEdge(e);
            }
        }
        for (int s: Si) {
            if(G.vertices().contains(s) && s != t){
                FlowEdge e = new FlowEdge(fakeNode, s, 1);
                Gf.addEdge(e);
            }
        }
        maxflow = new FordFulkerson(Gf, fakeNode, t);

        //StdOut.println("Max flow from " + s + " to " + t);
        /*
        for (int v = 0; v < Gf.V(); v++) {
            for (FlowEdge e : Gf.adj(v)) {
                if ((v == e.from()) && e.flow() > 0){
                    StdOut.println("   " + e);
                }
            }
        }
        */

        HashSet<Integer> B = new HashSet<Integer>();               // B = vertices from which there is a path to t in Gf
        for (int v = 0; v < Gf.V(); v++) {
            if(maxflow.hasAugmentingPath(Gf, v, t)){
                B.add(v);
            }
        }
        //System.out.println("B = " + B);
        HashSet<Integer> V = new HashSet<Integer>();               // V is the set of vertices of Ga
        for (int i: G.vertices()) {
            V.add(i);
        }
        //System.out.println("V = " + V);
        HashSet<Integer> A = new HashSet<Integer>();
        A = calculateSetSubtraction(V, B);                         // A = V \ B

        for (int v: Si){
            if(!A.contains(v) && G.vertices().contains(v) && v != t){
                A.add(v);
            }
        }

        //System.out.println("A = V-B = " + A);
        HashSet<int []> C = new HashSet<int []>();                 // C will be the FMC
        for (int v: A) {
            for(int w: G.adj(G.vertices().indexOf(v))) {
                if (B.contains(w)) {
                    C.add(new int[]{v, w});
                }
            }
        }
        for (int[] edge: C) {
            //System.out.println(edge[0] + " -> " + edge[1] );
        }
        return new HashSet[]{C, A, B};
    }

    public ArrayList<ArrayList<Integer>> ComputeSCC(Digraph G, Digraph Gr, int[] F){
        System.out.println("-----------------------------------------------  Computing SCCs of G\\F ----------------------------------------------");
        ArrayList<ArrayList<Integer>> C = new ArrayList<>();    //Collection of SCCs
        boolean[] W = new boolean[G.V()];                       // subset of V whose SCC have been computed

        HashSet<Integer> A;
        ArrayList<int[]> Paths;
        Paths = calculateFinalHeavyPathDecomposition();
        /*
        System.out.println("Paths = ");
        for(int[] Pab: Paths){
            System.out.println(Arrays.toString(Pab));
        }
        */
        for(int[] Pab: Paths){
            System.out.println("Processing path : " + Arrays.toString(Pab));
            int a = Pab[0];                 //Get the the starting vertex of the path
            int b = Pab[Pab.length-1];      //Get the the ending vertex of the path

            Xouts = new int[G.V()];
            Arrays.fill(Xouts, -1);  //Fill Xouts with: -1 == Null
            Xins = new int[G.V()];
            Arrays.fill(Xins, -1);   //Fill Xins with: -1 == Null

            /*
            System.out.println("=====================================================");
            System.out.println("path: " +  Arrays.toString(Pab) + " size = " + Pab.length);
            System.out.println("=====================================================");
             */

            //############################### BINARY SEARCH FOR Xouts GRAPH ##############################################
            System.out.println("Calculating Xouts...");
            //-----------------------------------------Calculate V1-----------------------------------------------------
            marked = new boolean[G.V()];
            HashSet<Integer> V1 = new HashSet<Integer>();
            getReachableVertices(a, kFTRSs.get(a).get(a), V1, false);

            //System.out.println(kFTRSs.get(a).get(a).toString());

            //-----------------------------------------Calculate Vt-----------------------------------------------------
            marked = new boolean[G.V()];
            HashSet<Integer> Vt = new HashSet<Integer>();
            getReachableVertices(b, kFTRSs.get(a).get(b), Vt, false);

            //System.out.println(V1+ ", " + Vt);

            //We need A = V1\Vt
            A = calculateSetSubtraction(V1, Vt);

            if(!A.isEmpty()){
                BinarySearch(0, Pab.length-2, A, InEdgesG, Xouts, Pab, false);
            }
            //Assign all vertices in Vt to Xout = t
            for (int v: Vt) {
                Xouts[v] = Pab.length-1;
            }
            //System.out.println("Xouts = " +Arrays.toString(Xouts));
            //############################### BINARY SEARCH FOR Xins IN THE REVERSED GRAPH ###############################
            System.out.println("Calculating Xins...");
            //------------------------------------Calculate V1r in Gr(reversed)------------------------------------------
            marked = new boolean[G.V()];
            HashSet<Integer> V1r = new HashSet<Integer>();
            getReachableVertices(a, kFTRSsR.get(a).get(a), V1r, true);

            //------------------------------------Calculate Vtr in Gr(reversed)------------------------------------------
            marked = new boolean[G.V()];
            HashSet<Integer> Vtr = new HashSet<Integer>();
            getReachableVertices(b, kFTRSsR.get(a).get(b), Vtr, true);

            //System.out.println(V1r+ ", " + Vtr);

            //We need A = Vtr\V1r
            A = calculateSetSubtraction(Vtr, V1r);

            if(!A.isEmpty() && Pab.length > 1 ){
                int[] PabR = new int[Pab.length];
                for (int i = 0; i < Pab.length; i++) {
                    PabR[i] = Pab[(Pab.length-1)-i];
                }
                BinarySearch(0, Pab.length-2, A, InEdgesGr, Xins, PabR, true);
            }
            //Assign all vertices in V1 to Xin = 0
            for (int v: V1r) {
                Xins[v] = 0;
            }
            //System.out.println("Xins = " + Arrays.toString(Xins));

            //################################### FIND THE SCCS INTESECTING P(a,b) #######################################
            Map<String, ArrayList<Integer>> SCCs = calculateSCCsFromXinfo(Xins, Xouts, G, Pab);

            //------------ Check if we have already calculated these SCCs or a bigger one that contains them -----------
            for (Map.Entry<String, ArrayList<Integer>> scc: SCCs.entrySet()) {
                ArrayList<Integer> vertices = scc.getValue();
                if(W[vertices.get(0)] == false){
                    for (int v: vertices){
                        W[v] = true;
                    }
                    C.add(vertices);
                }
            }
        }
        System.out.println("---------------------------------------------  Computing SCCs of G\\F Done -------------------------------------------");
        return C;
    }

    private void BinarySearch(int i, int j, HashSet<Integer> A, HashSet<?>[][] InEdges, int[] Xinfo, int[] Pab, boolean inReverse){
        HashSet<Integer> AdeepCopy = new HashSet<Integer>();
        AdeepCopy.addAll(A);
        if (i == j){
            for (int v: A) {
                if(inReverse){
                    Xinfo[v] = (Pab.length-1) - i;  //Xins
                    //System.out.println("Xinfo["+ v +"] = " + i);
                }else{
                    Xinfo[v] = i;                   //Xouts
                    //System.out.println("Xinfo["+ v +"] = " + i);
                }
            }
        }else {
            //System.out.println(i + ", " + j);
            int mid = (int) Math.ceil((double)(i+j)/2);
            //System.out.println("Pab[mid] = " + Pab[mid]);
            if(!AdeepCopy.isEmpty()){
                HashSet<Integer> B = Reach(Pab[mid], AdeepCopy, InEdges, inReverse);
                //System.out.println("A == " + A);
                //System.out.println("B == " + B);
                BinarySearch(i, mid-1, calculateSetSubtraction(AdeepCopy, B), InEdges, Xinfo, Pab, inReverse);
                BinarySearch(mid, j, B, InEdges, Xinfo, Pab, inReverse);
            }
        }
    }

    private HashSet<Integer> Reach(int mid, HashSet<Integer> A, HashSet<?>[][] InEdges, boolean inReverse){
        DynamicDigraph Ga = new DynamicDigraph(A);
        for (int v: A) {
            if(InEdges[mid][v] != null){
                for (Object y: InEdges[mid][v]){
                    if (A.contains(y)){
                        //System.out.println("Adding " + y + " -> " + v);
                        Ga.addEdge((Integer) y, v);
                    }
                }
            }
        }
        // B <- Vertices reachable from Xmid obtained by a BFS or DFS traversal of graph Ga
        HashSet<Integer> B = new HashSet<>();
        marked = new boolean[Vsize];
        getReachableVertices(mid, Ga, B, inReverse);
        return B;
    }

    private HashSet<Integer> calculateSetSubtraction(HashSet<Integer> A, HashSet<Integer> B){
        HashSet<Integer> result = new HashSet<Integer>();
        for (int v: A) {
            if(!B.contains(v)){
                result.add(v);
            }
        }
        return result;
    }

    private void getReachableVertices(int v, DynamicDigraph Ga, HashSet<Integer> Vi, boolean inReverse) {
        //System.out.println(v);
        Vi.add(v);
        marked[v] = true;
        if(Ga.vertices().contains(v)){
            for (int w : Ga.adj(Ga.vertices().indexOf(v))) {
                //System.out.println("inspecting edge = (" + v + "->" + w + ")");
                if(inReverse){
                    if (!marked[w]){
                        if((isEdge(w, v) && !isFailedEdge(w,v))){
                            getReachableVertices(w, Ga, Vi, true);
                        }
                    }
                }else{
                    if (!marked[w] && !isFailedEdge(v, w)){
                        getReachableVertices(w, Ga, Vi, false);
                    }
                }
            }
        }

    }

    private Map<String, ArrayList<Integer>> calculateSCCsFromXinfo(int[] Xins, int[] Xouts, Digraph G, int[] Pab) {
        ArrayList<Integer> componentsID;
        Map<String, ArrayList<Integer>> connectedComponents = new HashMap<String, ArrayList<Integer>>();
        for (int i = 0; i < G.V(); i++) {
            //System.out.println("Xouts["+ i + "] = " + Xouts[i] + ",  Xins["+ i + "] = " + Xins[i]);
            if((Xouts[i] != -1 && Xins[i] != -1) && (Xins[i] <= Xouts[i])){ //&& Pab[Xins[i]] <= Pab[Xouts[i]] ){
                String ID = Integer.toString(Pab[(Xins[i])]) + "-" + Integer.toString(Pab[(Xouts[i])]);
                if(!connectedComponents.containsKey(ID)){
                    connectedComponents.put(ID, new ArrayList<Integer>());
                }
                connectedComponents.get(ID).add(i);
            }
        }
        return connectedComponents;
    }

    public boolean isFailedEdge(int v, int w){
        boolean hasFailed = false;
        for (int[] edge: failedEdges) {
            if(edge[0] == v && edge[1] == w){
                hasFailed = true;
            }
        }
        return hasFailed;
    }

    public boolean isEdge(int v, int w) {
        for (int x : Gglob.adj(v)) {
            if (x == w){
                return true;
            }
        }
        return false;
    }

    public void printSCCS(ArrayList<ArrayList<Integer>> C){
        System.out.println();
        System.out.println(C.size() + " components");
        //PRINT OUT ALL THE SCCS
        for(ArrayList<Integer> SCC: C){
            for (int v: SCC) {
                System.out.print(v + " ");
            }
            System.out.println("");
        }
    }

    public static void main(String[] args) {
        In in = new In(args[0]);
        Digraph G = new Digraph(in);
        Digraph Gr = G.reverse();

        ArrayList<int[]> failedEdges = null;
        FaultTolerantSCC scc = new FaultTolerantSCC(G, Gr, failedEdges, 2);
    }

}