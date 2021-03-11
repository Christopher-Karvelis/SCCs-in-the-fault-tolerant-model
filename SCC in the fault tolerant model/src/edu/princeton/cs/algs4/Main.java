package edu.princeton.cs.algs4;

import java.util.ArrayList;

public class Main {
    public static void main(String[] args) {
        In in = new In(args[0]);
        Digraph G = new Digraph(in);
        Digraph Gr = G.reverse();
        int k = 2;
        ArrayList<int[]> failedEdges = new ArrayList<>();

        //exampleSimple
        //failedEdges.add(new int[]{2, 1});
        //failedEdges.add(new int[]{1, 0});

        //failedEdges.add(new int[]{4, 2});
        //failedEdges.add(new int[]{1, 0});

        //exampleNice
        failedEdges.add(new int[]{5, 1});
        failedEdges.add(new int[]{6, 7});

        //exampleTree
        //failedEdges.add(new int[]{0, 1});
        //failedEdges.add(new int[]{2, 6});

        //exampleDense
        //failedEdges.add(new int[]{5, 0});
        //failedEdges.add(new int[]{6, 2});

        //exampleDenseToThin
        //failedEdges.add(new int[]{6, 3});
        //failedEdges.add(new int[]{10, 3});

        //example2
        //failedEdges.add(new int[]{4, 5});
        //failedEdges.add(new int[]{3, 2});

        //exampleCycle
        //failedEdges.add(new int[]{1, 0});
        //failedEdges.add(new int[]{2, 3});

        FaultTolerantSCC kSCC = new FaultTolerantSCC(G, Gr, failedEdges, k);
        kSCC.preprocessing(G, Gr, k);
        ArrayList<ArrayList<Integer>> C = kSCC.ComputeSCC(G, Gr,null); //F will be the set failedEdges in the future

        System.out.print("\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ k-Fault Tolerant Algorithm Result $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$");
        kSCC.printSCCS(C);
        System.out.println("\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Tarjan Result $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$");
        TarjanSCC scc = new TarjanSCC(G, failedEdges);

        // number of connected components
        int m = scc.count();
        StdOut.println(m + " components");

        // compute list of vertices in each strong component
        Queue<Integer>[] components = (Queue<Integer>[]) new Queue[m];
        for (int i = 0; i < m; i++) {
            components[i] = new Queue<Integer>();
        }
        for (int v = 0; v < G.V(); v++) {
            components[scc.id(v)].enqueue(v);
        }

        // print results
        for (int i = 0; i < m; i++) {
            for (int v : components[i]) {
                StdOut.print(v + " ");
            }
            StdOut.println();
        }
    }
}