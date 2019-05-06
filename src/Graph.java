import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.Set;
import java.util.HashSet;
import java.util.Iterator;

public class Graph {

    // Keep a fast index to nodes in the map
    private Map<Integer, Vertex> vertexNames;

    /**
     * Construct an empty Graph with a map. The map's key is the name of a vertex
     * and the map's value is the vertex object.
     */
    public Graph() {
        vertexNames = new HashMap<>();
    }

    /**
     * Adds a vertex to the graph. Throws IllegalArgumentException if two vertices
     * with the same name are added.
     *
     * @param v
     *          (Vertex) vertex to be added to the graph
     */
    public void addVertex(Vertex v) {
        if (vertexNames.containsKey(v.name))
            throw new IllegalArgumentException("Cannot create new vertex with existing name.");
        vertexNames.put(v.name, v);
    }

    /**
     * Gets a collection of all the vertices in the graph
     *
     * @return (Collection<Vertex>) collection of all the vertices in the graph
     */
    public Collection<Vertex> getVertices() {
        return vertexNames.values();
    }

    /**
     * Gets the vertex object with the given name
     *
     * @param name
     *          (String) name of the vertex object requested
     * @return (Vertex) vertex object associated with the name
     */
    public Vertex getVertex(String name) {
        return vertexNames.get(name);
    }

    /**
     * Adds a directed edge from vertex u to vertex v
     *
     * @param nameU
     *          (String) name of vertex u
     * @param nameV
     *          (String) name of vertex v
     * @param cost
     *          (double) cost of the edge between vertex u and v
     */
    public void addEdge(int nameU, int nameV, Double cost) {
        if (!vertexNames.containsKey(nameU))
            throw new IllegalArgumentException(nameU + " does not exist. Cannot create edge.");
        if (!vertexNames.containsKey(nameV))
            throw new IllegalArgumentException(nameV + " does not exist. Cannot create edge.");
        Vertex sourceVertex = vertexNames.get(nameU);
        Vertex targetVertex = vertexNames.get(nameV);
        Edge newEdge = new Edge(sourceVertex, targetVertex, cost);
        sourceVertex.addEdge(newEdge);
    }

    /**
     * Adds an undirected edge between vertex u and vertex v by adding a directed
     * edge from u to v, then a directed edge from v to u
     *
     * @param name
     *          (String) name of vertex u
     * @param name2
     *          (String) name of vertex v
     * @param cost
     *          (double) cost of the edge between vertex u and v
     */
    public void addUndirectedEdge(int name, int name2, double cost) {
        addEdge(name, name2, cost);
        addEdge(name2, name, cost);
    }


    /**
     * Computes the euclidean distance between two points as described by their
     * coordinates
     *
     * @param ux
     *          (double) x coordinate of point u
     * @param uy
     *          (double) y coordinate of point u
     * @param vx
     *          (double) x coordinate of point v
     * @param vy
     *          (double) y coordinate of point v
     * @return (double) distance between the two points
     */
    public double computeEuclideanDistance(double ux, double uy, double vx, double vy) {
        return Math.sqrt(Math.pow(ux - vx, 2) + Math.pow(uy - vy, 2));
    }

    /**
     * Computes euclidean distance between two vertices as described by their
     * coordinates
     *
     * @param u
     *          (Vertex) vertex u
     * @param v
     *          (Vertex) vertex v
     * @return (double) distance between two vertices
     */
    public double computeEuclideanDistance(Vertex u, Vertex v) {
        return computeEuclideanDistance(u.x, u.y, v.x, v.y);
    }

    /**
     * Calculates the euclidean distance for all edges in the map using the
     * computeEuclideanCost method.
     */
    public void computeAllEuclideanDistances() {
        for (Vertex u : getVertices())
            for (Edge uv : u.adjacentEdges) {
                Vertex v = uv.target;
                uv.distance = computeEuclideanDistance(u.x, u.y, v.x, v.y);
            }
    }



    // STUDENT CODE STARTS HERE

    public void generateRandomVertices(int n) {

            vertexNames = new HashMap<>(); // reset the vertex hashmap

            Random r = new Random ();

            for (int i=0; i<n; i++) {

                Vertex v = new Vertex (i, r.nextInt (100), r.nextInt (100));

                addVertex (v);

                // add this new Vertex as an edge to each one of the previous edges
                for (int j=i-1; j>=0; j--) {

                    //addEdge (i, j, 0.0);
                    addUndirectedEdge(i, j, 0.0);

                }
            }

            computeAllEuclideanDistances(); // compute distances
    }

    public List<Edge> nearestNeighborTsp() {
        Random r = new Random ();
        Vertex startVertex = vertexNames.get (r.nextInt (vertexNames.size()));

        ArrayList<Edge> edges = new ArrayList<> ();

        HashSet<Vertex> visited = new HashSet <> ();

        Double distance = 0.0;
        System.out.println ("Start: " + startVertex + ", distance : " + distance);
        Vertex nextVertex = startVertex;
        visited.add (nextVertex);

        while (visited.size() < vertexNames.size()) {

            Iterator i = nextVertex.adjacentEdges.iterator();
            Edge closest = null;

            while (i.hasNext()) {

                Edge nextEdge = (Edge) i.next();

                if (visited.contains (nextEdge.target)) {
                    continue;
                }
                if (closest == null) {
                    closest = nextEdge;
                    continue;
                }

                if (closest.distance > nextEdge.distance) {
                    closest = nextEdge;
                }
            }


            // add closest to edges and move nextVertex to
            // target of the closest edge and mark it as visited
            distance += closest.distance;
            edges.add (closest);
            nextVertex = closest.target;
            visited.add (nextVertex);
            System.out.println ("Next : " + nextVertex + ", distance : " + distance);
        }

        // add the return from the last edge to starting point
        Edge nextEdge = getEdge (nextVertex, startVertex);
        edges.add (nextEdge);
        distance += nextEdge.distance;
        System.out.println ("Back : " + nextEdge.target + ", distance : " + distance);

        return edges;
    }

    // Return the edge that connects Source to Target
    private Edge getEdge (Vertex source, Vertex target) {
        Iterator i = source.adjacentEdges.iterator();

        while (i.hasNext()) {
            Edge nextEdge = (Edge) i.next();
            if (nextEdge.target == target) {
                return nextEdge;
            }
        }
        System.err.println("Could not find an edge on " + source + " for " + target);
        return null;

    }


    //
    // Generating permutation using Heap Algorithm for a list if names which are int's
    // https://www.geeksforgeeks.org/heaps-algorithm-for-generating-permutations/
    //
    static void heapPermutation(int a[], int size, ArrayList<int []> output)
    {
        // if size becomes 1 then prints the obtained
        // permutation
        if (size == 1) {
            int b[] = a.clone();
            // add the first one to the last slot for the return
            b[b.length-1] = b[0];
            output.add(b);
        }

        for (int i=0; i<size; i++)
        {
            heapPermutation(a, size-1, output);

            // if size is odd, swap first and last
            // element
            if (size % 2 == 1)
            {
                int temp = a[0];
                a[0] = a[size-1];
                a[size-1] = temp;
            }

            // If size is even, swap ith and last
            // element
            else
            {
                int temp = a[i];
                a[i] = a[size-1];
                a[size-1] = temp;
            }
        }
    }
    static private void printPerm (int perm []) {
        for (int i=0; i<perm.length; i++)
            System.out.print(perm[i] + " ");
        System.out.println();
    }

    // used for testing the heap Algorithm
    public static void main (String args[]) {
        int [] a = {1, 2, 3, 0};
        ArrayList<int []> output = new ArrayList<> ();
        heapPermutation (a, a.length-1, output);
        for (int [] one: output) {
            printPerm (one);
        }
    }

    public List<Edge> bruteForceTsp() {
        // make a list of names (int's) to be used
        // in the haapAlgorithm to build all permuations
        // of vertexes
        Set<Integer> keys = vertexNames.keySet();
        // +1 since need to return to the starting point.
        // need one extra to add the starting vertex name
        int names [] = new int [vertexNames.size()+1];
        Iterator it = keys.iterator();
        int index = 0;
        while (it.hasNext()) {
            names[index] = ((Integer)it.next()).intValue();
            index++;
        }

        ArrayList<int []> perms = new ArrayList<> ();
        heapPermutation(names, names.length-1, perms);

        ArrayList<Edge> shortest = new ArrayList<> ();
        double shortDist = Double.MAX_VALUE;
        // Iterate over each permutation to find total distance
        // and save shortest.
        for (int one[] : perms) {
            double currDist = 0;
            ArrayList<Edge> currEdges = new ArrayList<>();
            // calculate the distance for this permutation of Vertexes.
            for (int i=1; i < one.length; i++) {
                Vertex source = vertexNames.get (one[i-1]);
                Vertex target = vertexNames.get (one[i]);
                Edge edge = getEdge (source, target);
                currDist += edge.distance;
                currEdges.add (edge);
            }
            printPerm (one);
            System.out.println (shortDist + ", CurrDist : " + currDist);
            if (currDist < shortDist) {
                shortDist = currDist;
                shortest = currEdges;
            }
        }

        return shortest;
    }

    // STUDENT CODE ENDS HERE



    /**
     * Prints out the adjacency list of the graph for debugging
     */
    public void printAdjacencyList() {
        for (int u : vertexNames.keySet()) {
            StringBuilder sb = new StringBuilder();
            sb.append(u);
            sb.append(" -> [ ");
            for (Edge e : vertexNames.get(u).adjacentEdges) {
                sb.append(e.target.name);
                sb.append("(");
                sb.append(e.distance);
                sb.append(") ");
            }
            sb.append("]");
            System.out.println(sb.toString());
        }
    }
}
