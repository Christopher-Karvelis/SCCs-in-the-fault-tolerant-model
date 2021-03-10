package edu.princeton.cs.algs4;

import java.util.ArrayList;
import java.util.HashSet;

public class DynamicDigraph {
    private static final String NEWLINE = System.getProperty("line.separator");
    private final int V;           // number of vertices in this digraph
    private int E;                 // number of edges in this digraph
    private ArrayList<Integer> vertices;
    private ArrayList<ArrayList<Integer>> adj;    // adj[v] = adjacency list for vertex v
    private int[] indegree;        // indegree[v] = indegree of vertex v
    private int nextFreePosition;   //Indicates the next free position on the array idsToLabesl;
    private int[] labesOfVertexIDs;      //Contains vertex names for valid ids in the [0,V)

    /**
     * Initializes an empty digraph with <em>V</em> vertices.
     *
     * @param  V the number of vertices
     * @throws IllegalArgumentException if {@code V < 0}
     */
    public DynamicDigraph (HashSet<Integer> vertices) {
        if (vertices.size() == 0) throw new IllegalArgumentException("Number of vertices in a edu.princeton.cs.algs4.Digraph must be nonnegative");
        this.V = vertices.size();
        this.E = 0;
        this.vertices = new ArrayList<Integer>();
        adj = new ArrayList<ArrayList<Integer>>();
        for (int v: vertices) {
            this.vertices.add(v);
            adj.add(new ArrayList<Integer>());
        }
        nextFreePosition = 0;
        indegree = new int[V];
        labesOfVertexIDs = new int[V];
    }

    public DynamicDigraph (ArrayList<Integer> vertices) {
        if (vertices.size() == 0) throw new IllegalArgumentException("Number of vertices in a edu.princeton.cs.algs4.Digraph must be nonnegative");
        this.V = vertices.size();
        this.E = 0;
        this.vertices = new ArrayList<Integer>();
        adj = new ArrayList<ArrayList<Integer>>();
        for (int v: vertices) {
            this.vertices.add(v);
            adj.add(new ArrayList<Integer>());
        }
        nextFreePosition = 0;
        indegree = new int[V];
        labesOfVertexIDs = new int[V];
    }

    /**
     * Initializes a digraph from the specified input stream.
     * The format is the number of vertices <em>V</em>,
     * followed by the number of edges <em>E</em>,
     * followed by <em>E</em> pairs of vertices, with each entry separated by whitespace.
     *
     * @param  in the input stream
     * @throws IllegalArgumentException if {@code in} is {@code null}
     * @throws IllegalArgumentException if the endpoints of any edge are not in prescribed range
     * @throws IllegalArgumentException if the number of vertices or edges is negative
     * @throws IllegalArgumentException if the input stream is in the wrong format
     */
    /*
    public Digraph(In in) {
        if (in == null) throw new IllegalArgumentException("argument is null");
        try {
            this.V = in.readInt();
            if (V < 0) throw new IllegalArgumentException("number of vertices in a edu.princeton.cs.algs4.Digraph must be nonnegative");
            indegree = new int[V];
            adj = (Bag<Integer>[]) new Bag[V];
            for (int v = 0; v < V; v++) {
                adj[v] = new Bag<Integer>();
            }
            int E = in.readInt();
            if (E < 0) throw new IllegalArgumentException("number of edges in a edu.princeton.cs.algs4.Digraph must be nonnegative");
            for (int i = 0; i < E; i++) {
                int v = in.readInt();
                int w = in.readInt();
                addEdge(v, w);
            }
        }
        catch (NoSuchElementException e) {
            throw new IllegalArgumentException("invalid input format in edu.princeton.cs.algs4.Digraph constructor", e);
        }
    }
    */
    /**
     * Initializes a new digraph that is a deep copy of the specified digraph.
     *
     * @param  G the digraph to copy
     * @throws IllegalArgumentException if {@code G} is {@code null}
     */
    /*
    public Digraph(Digraph G) {
        if (G == null) throw new IllegalArgumentException("argument is null");

        this.V = G.V();
        this.E = G.E();
        if (V < 0) throw new IllegalArgumentException("Number of vertices in a edu.princeton.cs.algs4.Digraph must be nonnegative");

        // update indegrees
        indegree = new int[V];
        for (int v = 0; v < V; v++)
            this.indegree[v] = G.indegree(v);

        // update adjacency lists
        adj = (Bag<Integer>[]) new Bag[V];
        for (int v = 0; v < V; v++) {
            adj[v] = new Bag<Integer>();
        }

        for (int v = 0; v < G.V(); v++) {
            // reverse so that adjacency list is in same order as original
            java.util.Stack<Integer> reverse = new Stack<Integer>();
            for (int w : G.adj[v]) {
                reverse.push(w);
            }
            for (int w : reverse) {
                adj[v].add(w);
            }
        }
    }
    */
    /**
     * Returns the number of vertices in this digraph.
     *
     * @return the number of vertices in this digraph
     */
    public int V() {
        return V;
    }

    /**
     * Returns the number of edges in this digraph.
     *
     * @return the number of edges in this digraph
     */
    public int E() {
        return E;
    }


    // throw an IllegalArgumentException unless {@code 0 <= v < V}
    private void validateVertex(int v) {
        if (v < 0 || v >= V)
            throw new IllegalArgumentException("vertex " + v + " is not between 0 and " + (V-1));
    }

    /**
     * Adds the directed edge v→w to this digraph.
     *
     * @param  v the tail vertex
     * @param  w the head vertex
     * @throws IllegalArgumentException unless both {@code 0 <= v < V} and {@code 0 <= w < V}
     */
    public void addEdge(int v, int w) {
        //System.out.println("for v = " + v);
        int index = vertices.indexOf(v);
        //System.out.println(index);
        adj.get(index).add(w);
        E++;
    }

    /**
     * Returns the vertices adjacent from vertex {@code v} in this digraph.
     *
     * @param  v the vertex
     * @return the vertices adjacent from vertex {@code v} in this digraph, as an iterable
     * @throws IllegalArgumentException unless {@code 0 <= v < V}
     */
    public ArrayList<Integer> adj(int v) {
        return adj.get(v);
    }

    public ArrayList<Integer> vertices() {
        return vertices;
    }

    /**
     * Returns the number of directed edges incident from vertex {@code v}.
     * This is known as the <em>outdegree</em> of vertex {@code v}.
     *
     * @param  v the vertex
     * @return the outdegree of vertex {@code v}
     * @throws IllegalArgumentException unless {@code 0 <= v < V}
     */
    public int outdegree(int v) {
        validateVertex(v);
        return adj.get(v).size();
    }

    /**
     * Returns the number of directed edges incident to vertex {@code v}.
     * This is known as the <em>indegree</em> of vertex {@code v}.
     *
     * @param  v the vertex
     * @return the indegree of vertex {@code v}
     * @throws IllegalArgumentException unless {@code 0 <= v < V}
     */
    public int indegree(int v) {
        validateVertex(v);
        return indegree[v];
    }

    /**
     * Returns the reverse of the digraph.
     *
     * @return the reverse of the digraph
     */
    public Digraph reverse() {
        Digraph reverse = new Digraph(V);
        for (int v: vertices) {
            for (int w : adj.get(v)) {
                reverse.addEdge(w, v);
            }
        }
        return reverse;
    }

    /**
     * Returns a string representation of the graph.
     *
     * @return the number of vertices <em>V</em>, followed by the number of edges <em>E</em>,
     *         followed by the <em>V</em> adjacency lists
     */
    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append(V + " vertices, " + E + " edges " + NEWLINE);
        for (int v: vertices) {
            s.append(String.format("%d: ", v));
            for (int w : adj.get(vertices.indexOf(v))) {
                s.append(String.format("%d ", w));
            }
            s.append(NEWLINE);
        }
        return s.toString();
    }

    /**
     * Unit tests the {@code edu.princeton.cs.algs4.Digraph} data type.
     *
     * @param args the command-line arguments
     */
    public static void main(String[] args) {
        In in = new In(args[0]);
        Digraph G = new Digraph(in);
        StdOut.println(G);
    }
}
