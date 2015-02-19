"""
SImple Network Python Library for Education (SINPLE)

SINPLE is a Python 3 library of functions to create and analyse networks,
i.e. graphs, of social relations, transport links, communication channels, etc.

SINPLE is a pedagogical (i.e. simple, readable, and documented) implementation
to help introduce Python, data structures and algorithms, and network/graph
theory. SINPLE includes almost 50 exercises and a small real-life network.
SINPLE presents network concepts through executable definitions (functions)
and concrete examples (unit tests).

SINPLE handles undirected graphs with nodes represented by any hashable value,
e.g. an integer, a string of the node's name, or a tuple with additional
information about the node.

SINPLE's functions assume their arguments are correct, i.e. are of the right
type and satisfy any additional stated conditions. All functions
that take as arguments a graph and some edges (or nodes),
except `network()` and `add_edges()`,
assume that the graph contains those edges or nodes.

SINPLE is available under a permissive license from http://tiny.cc/sinple.
"""

# Coding conventions
# ------------------
# Functions named with adjectives or `is_`*noun* return booleans.
# Functions named with nouns return network entities or metrics.
# The code passes the [PEP8 style checker](https://github.com/jcrocholl/pep8).

# Exercises
# ---------
# SINPLE contains many exercises, of varying type and difficulty.
# Some exercises are more relevant for introductory Python courses,
# some for data structures and algorithms courses, others for network courses.

# Exercises aren't ordered by type or difficulty: they appear as soon as they
# can be solved, using the preceding functions.
# Exercises at the end of SINPLE require changes throughout the library.
# Exercises are grouped in bullet lists; there is no heading introducing them.
#
# A general, implicit, exercise is to study and understand the code.
# Some exercises focus on explaining particular functions.
#
# Some exercises ask to change existing functions, others to add new code.
# When writing a new function, think about its name and parameters,
# any preconditions the user needs to know about,
# and add the corresponding unit tests.
#
# Please don't share exercise solutions publicly.

# Network representation
# ----------------------
# The mathematical basis of network science is graph theory.
# Graph terminology varies among authors.
# We follow the Wikipedia [glossary of graph
# theory](http://en.wikipedia.org/wiki/Glossary_of_graph_theory).
# except in using 'node' instead of 'vertex'.
# If a term has multiple meanings in the literature, we chose one of them.
#
# In SINPLE, a **network** (also called a **graph**)
# is a set of **nodes** connected by a set of undirected **edges**.
# Each pair of nodes is connected by one edge or it is not connected.
# A node may be connected to itself.
#
# SINPLE represents a network by a tuple `(nodes, edges)`,
# where `nodes` is a set of hashable values (numbers, strings, tuples, etc.)
# and `edges` is a set of edges.
# Each edge is a set of nodes `{node1, node2}`,
# because the edge has no direction.

# Core functions
# --------------
# These functions hide the network representation
# from the rest of the library.

# Start with the **empty graph**, which has no nodes and hence no edges.
EMPTY = (set(), set())


def edge_set(graph):
    """Return the set of the graph's edges."""
    (nodes, edges) = graph
    return edges

assert edge_set(EMPTY) == set()


def node_set(graph):
    """Return the set of the graph's nodes."""
    (nodes, edges) = graph
    return nodes

assert node_set(EMPTY) == set()


def edge(node1, node2):
    """Create an edge connecting the two nodes, which may be the same."""
    return frozenset({node1, node2})

# The order of nodes doesn't matter.
assert edge(1, 2) == edge(2, 1)

# - Consult the Python documentation about sets.
#   Why is an edge represented by a frozen, i.e. immutable, set?
#   Hint: read again the definition of graph at the start of this section.


# The function to deconstruct a previously constructed edge returns a pair
# instead of a set to allow obtaining each node with the assignment
# `(n1, n2) = endpoints(e)`.
def endpoints(edge):
    """
    Return a pair with the edge's **endpoints** (the nodes it connects).
    """
    if len(edge) == 2:
        return tuple(edge)
    else:
        return tuple(edge) + tuple(edge)

assert endpoints(edge(2, 1)) == (2, 1) or endpoints(edge(2, 1)) == (1, 2)
assert endpoints(edge(1, 1)) == (1, 1)


def add_edges(edges, graph):
    """
    Create a new graph by adding the set of edges and their endpoints to graph.
    """
    # Start with the empty set. Note that `{}` is the empty dictionary.
    nodes = set()
    for edge in edges:
        # Add both endpoints to the set.
        nodes.update(endpoints(edge))
    # The new graph's nodes/edges is the union of the old and new sets.
    return (node_set(graph) | nodes, edge_set(graph) | edges)

assert node_set(add_edges({edge(2, 1)}, EMPTY)) == {1, 2}


# A convenience function for the common case of adding a single edge.
def add_edge(node1, node2, graph):
    """
    Create a new graph by adding to graph the nodes and an edge between them.

    The two nodes can be the same.
    """
    return add_edges({edge(node1, node2)}, graph)

assert node_set(add_edge(1, 1, EMPTY)) == {1}


def network(edges, nodes=set()):
    """
    Create a graph with the given sets of nodes, edges, and their endpoints.
    """
    return add_edges(edges, (nodes, set()))

assert network(set()) == EMPTY
assert network({edge(1, 2)}) == network({edge(2, 1)}, {1, 2})


def delete_edges(edges, graph):
    """
    Create a new graph by removing the set of edges from graph.

    Edges that don't exist in the graph are ignored.
    """
    return (node_set(graph), edge_set(graph) - edges)

assert delete_edges({edge(1, 2)}, EMPTY) == EMPTY


# A convenience function for the common case of deleting a single edge.
def delete_edge(node1, node2, graph):
    """
    Create a new graph by removing from graph the edge between the nodes.
    """
    return delete_edges({edge(node1, node2)}, graph)

assert delete_edge(1, 2, network({edge(2, 1)})) == network(set(), {1, 2})

# - Write functions `delete_nodes()` and `delete_node()`.
#   Don't forget to remove the attached edges.


# ### Example networks
#
# These networks are used in further unit tests.

# A **null graph** has no edges.
N1 = network(set(), {1})
N2 = network(set(), {1, 2})

# The smallest non-simple graph has a single loop.
LOOP = network({edge(1, 1)})

# A **complete graph** K*n* connects each node to all other *n-1* nodes.
K2 = network({edge(1, 2)})
K3 = network({edge(1, 2), edge(1, 3), edge(2, 3)})

# - Define K1.
# - Write a function to create a complete graph, given the number of nodes.
#   Write tests comparing the output to the 'manually' created K1, K2, K3.

# A **cycle graph** C*n* has *n* nodes connected in a 'round-robin' fashion.
# Cycle graphs can be drawn as regular shapes: triangle, square, pentagon, etc.
C3 = network({edge(1, 2), edge(2, 3), edge(3, 1)})
C4 = add_edges({edge(3, 4), edge(4, 1)}, delete_edge(3, 1, C3))
C5 = add_edges({edge(4, 5), edge(5, 1)}, delete_edge(4, 1, C4))

# - Define C2.
# - Explain how C4 is constructed from C3.
# - Write a function to create a cycle graph, given the number of nodes.

# The **Petersen graph** is a (counter-)example of various graph concepts.
PETERSEN = add_edges({
    # It is usually drawn as an inner pentagram (C5)
    edge(6, 8), edge(8, 10), edge(10, 7), edge(7, 9), edge(9, 6),
    # connected to
    edge(1, 6), edge(2, 7), edge(3, 8), edge(4, 9), edge(5, 10)},
    # an outer pentagon (C5).
    C5)

# - Draw the Petersen graph on paper and label the nodes as in `PETERSEN`.

# The [utility graph](http://en.wikipedia.org/wiki/Water,_gas,_and_electricity)
# connects 3 houses to 3 utilities.
UTILITY = network({
    edge(house, utility)
    for house in [1, 2, 3] for utility in ["gas", "water", "power"]
    })

# In December 1970, the Arpanet (the precursor of the Internet) had 13 nodes.
# Source: http://som.csudh.edu/cis/lpress/history/arpamaps.
ARPANET = network({
    edge("SRI", "Stanford"), edge("SRI", "UCSB"),
    edge("SRI", "UCLA"), edge("SRI", "Utah"),
    edge("Stanford", "UCLA"), edge("UCLA", "UCSB"), edge("UCLA", "RAND"),
    edge("RAND", "SDC"), edge("SDC", "Utah"), edge("Utah", "MIT"),
    edge("RAND", "BBN"), edge("MIT", "BBN"), edge("BBN", "Harvard"),
    edge("Harvard", "Carnegie"), edge("Carnegie", "CASE"),
    edge("CASE", "Lincoln"), edge("Lincoln", "MIT")
    })


# The 'kite network' by David Krackhardt illustrates different centralities.
# Source: http://www.orgnet.com/sna.html
KITE = network({
    edge("Andre", "Beverly"), edge("Andre", "Carol"),
    edge("Andre", "Diane"), edge("Andre", "Fernando"),
    edge("Beverly", "Diane"), edge("Beverly", "Ed"),
    edge("Carol", "Diane"), edge("Carol", "Fernando"),
    edge("Diane", "Ed"), edge("Diane", "Garth"),
    edge("Fernando", "Garth"), edge("Fernando", "Heather"),
    edge("Garth", "Heather"), edge("Heather", "Ike"), edge("Ike", "Jane")})


# Wayne Zachary obtained in the 1970s the ties between a karate club's members.
# Source: http://www-personal.umich.edu/~mejn/netdata/
KARATE = network({
    edge(2, 1), edge(3, 1), edge(3, 2), edge(4, 1), edge(4, 2), edge(4, 3),
    edge(5, 1), edge(6, 1), edge(7, 1), edge(7, 5), edge(7, 6),
    edge(8, 1), edge(8, 2), edge(8, 3), edge(8, 4), edge(9, 1), edge(9, 3),
    edge(10, 3), edge(11, 1), edge(11, 5), edge(11, 6), edge(12, 1),
    edge(13, 1), edge(13, 4),
    edge(14, 1), edge(14, 2), edge(14, 3), edge(14, 4),
    edge(17, 6), edge(17, 7), edge(18, 1),  edge(18, 2),
    edge(20, 1), edge(20, 2), edge(22, 1), edge(22, 2),
    edge(26, 24), edge(26, 25), edge(28, 3), edge(28, 24), edge(28, 25),
    edge(29, 3), edge(30, 24), edge(30, 27), edge(31, 2), edge(31, 9),
    edge(32, 1), edge(32, 25), edge(32, 26), edge(32, 29),
    edge(33, 3), edge(33, 9), edge(33, 15), edge(33, 16), edge(33, 19),
    edge(33, 21), edge(33, 23), edge(33, 24), edge(33, 30), edge(33, 31),
    edge(33, 32),
    edge(34, 9), edge(34, 10), edge(34, 14), edge(34, 15), edge(34, 16),
    edge(34, 19), edge(34, 20), edge(34, 21), edge(34, 23), edge(34, 24),
    edge(34, 27), edge(34, 28), edge(34, 29), edge(34, 30), edge(34, 31),
    edge(34, 32), edge(34, 33)})


# The 'Southern Women' dataset about attendance of social events by women
# in Natchez MS, USA, was collected in the 1930s by Davis et al.
# Source: Figure 1 of http://moreno.ss.uci.edu/86.pdf
# Note that other online datasets,
# e.g. [igraph](http://nexus.igraph.org/api/dataset_info?id=23&format=html),
# use the probably erroneous Figure 2.

WOMEN = network({
    edge("Evelyn", 1), edge("Evelyn", 2), edge("Evelyn", 3), edge("Evelyn", 4),
    edge("Evelyn", 5), edge("Evelyn", 6), edge("Evelyn", 8), edge("Evelyn", 9),
    edge("Laura", 1), edge("Laura", 2), edge("Laura", 3), edge("Laura", 5),
    edge("Laura", 6), edge("Laura", 7), edge("Laura", 8),
    edge("Theresa", 2), edge("Theresa", 3), edge("Theresa", 4),
    edge("Theresa", 5), edge("Theresa", 6), edge("Theresa", 7),
    edge("Theresa", 8), edge("Theresa", 9),
    edge("Brenda", 1), edge("Brenda", 3), edge("Brenda", 4), edge("Brenda", 5),
    edge("Brenda", 6), edge("Brenda", 7), edge("Brenda", 8),
    edge("Charlotte", 3), edge("Charlotte", 4),
    edge("Charlotte", 5), edge("Charlotte", 7),
    edge("Frances", 3), edge("Frances", 5),
    edge("Frances", 6), edge("Frances", 8),
    edge("Eleanor", 5), edge("Eleanor", 6),
    edge("Eleanor", 7), edge("Eleanor", 8),
    edge("Pearl", 6), edge("Pearl", 8), edge("Pearl", 9),
    edge("Ruth", 5), edge("Ruth", 7), edge("Ruth", 8), edge("Ruth", 9),
    edge("Verne", 7), edge("Verne", 8), edge("Verne", 9), edge("Verne", 12),
    edge("Myra", 8), edge("Myra", 9), edge("Myra", 10), edge("Myra", 12),
    edge("Katherine", 8), edge("Katherine", 9), edge("Katherine", 10),
    edge("Katherine", 12), edge("Katherine", 13), edge("Katherine", 14),
    edge("Sylvia", 7), edge("Sylvia", 8), edge("Sylvia", 9),
    edge("Sylvia", 10), edge("Sylvia", 12),
    edge("Sylvia", 13), edge("Sylvia", 14),
    edge("Nora", 6), edge("Nora", 7), edge("Nora", 9), edge("Nora", 10),
    edge("Nora", 11), edge("Nora", 12), edge("Nora", 13), edge("Nora", 14),
    edge("Helen", 7), edge("Helen", 8), edge("Helen", 10),
    edge("Helen", 11), edge("Helen", 12),
    edge("Dorothy", 8), edge("Dorothy", 9),
    edge("Olivia", 9), edge("Olivia", 11),
    edge("Flora", 9), edge("Flora", 11)
})


# Basics
# ------


def size(graph):
    """Return the graph's **size** (number of edges)."""
    return len(edge_set(graph))

assert size(EMPTY) == 0
assert size(N1) == 0
assert size(C3) == 3
assert size(PETERSEN) == 15
assert size(WOMEN) == 89
assert size(KARATE) == 78


def order(graph):
    """Return the graph's **order** (number of nodes)."""
    return len(node_set(graph))

assert order(EMPTY) == 0
assert order(N1) == 1
assert order(C3) == order(K3)
assert order(PETERSEN) == 10
assert order(KARATE) == 34
# The Southern women dataset has 14 events and 18 women.
assert order(WOMEN) == 14 + 18


def is_loop(edge):
    """Check if edge is a **loop** (connects a node to itself)."""
    (node1, node2) = endpoints(edge)
    return node1 == node2

assert is_loop(edge("A", "A"))
assert not is_loop(edge("A", "B"))


def simple(graph):
    """Check if the graph is **simple** (has no loops)."""
    for edge in edge_set(graph):
        if is_loop(edge):
            return False
    return True

assert simple(EMPTY)
assert simple(N1)
# Any complete graph is simple.
assert simple(K3)
assert simple(PETERSEN)
assert not simple(LOOP)


def density(graph):
    """
    Return the graph's **density** (ratio of actual to potential edges).

    It's the size divided by the size of a complete graph of the same order.
    The density of the graphs with one or zero nodes is set to zero.
    """
    o = order(graph)
    return size(graph) / (o * (o - 1) / 2) if o > 1 else 0

assert density(EMPTY) == 0
assert density(LOOP) == 0
# Every null graph has density 0.
assert density(N2) == 0
# Every complete graph has density 1.
assert density(K3) == 1
# A graph with loops can have density > 1.
assert density(add_edge(1, 1, K3)) > 1

# Further graph construction
# --------------------------


def null(n):
    """Create a **null graph** with n nodes numbered 1 to n and no edges."""
    return network(set(), {x for x in range(1, n+1)})

assert null(0) == EMPTY
assert null(2) == N2


def renumber(graph, n=1):
    """Create a new graph by renaming the graph's nodes as n, n+1, n+2, etc."""
    number = n
    map = {}
    for node in node_set(graph):
        map[node] = number
        number = number + 1
    nodeset = set(range(n, number))
    edgeset = set()
    for e in edge_set(graph):
        (node1, node2) = endpoints(e)
        edgeset.add(edge(map[node1], map[node2]))
    return network(edgeset, nodeset)

assert renumber(EMPTY, 1) == EMPTY
assert renumber(N2, -1) == network(set(), {-1, 0})
assert renumber(network({edge("A", "B")})) == network({edge(1, 2)})


def subgraph(nodes, graph):
    """
    Return the graph's **subgraph induced** by the set of nodes.

    The induced subgraph has the graph's edges that connect the given nodes.
    """
    edges = {edge(n1, n2) for n1 in nodes for n2 in nodes}
    return network(edges & edge_set(graph), nodes)

assert subgraph({1}, N2) == N1
assert subgraph({1, 2, 3, 4, 5}, PETERSEN) == C5

# - Explain how the function works. Hint: `&` denotes set intersection.
# - In which cases is the above implementation efficient?
#   Implement the function in a way that is more efficient for other cases.
#   Hint: what happens for sparse/dense graphs when many/few nodes are kept?
# - Write a function to check if *g1* is a **subgraph** of *g2*,
#   i.e. if the nodes and edges of *g1* are also in *g2*.
# - Rewrite the `delete_nodes` function using `subgraph`.


# Import a random number generator with uniform distribution
from random import uniform


def random(n, p):
    """
    Return a simple **random graph** with n nodes and edge probability p.

    The two numbers can't be negative, and p can't be larger than 1.
    """
    edges = set()
    for node1 in range(1, n+1):
        for node2 in range(node1+1, n+1):
            if uniform(0, 1) < p:
                edges.add(edge(node1, node2))
    return add_edges(edges, null(n))

assert random(0, 1) == EMPTY
assert random(1, 1) == N1
assert random(3, 1) == K3
assert order(random(4, 0)) == 4
assert size(random(4, 0)) == 0
assert size(random(4, 1)) == 6
assert simple(random(4, 0.5))

# - Write a function to create a simple random graph of given order and size.


# Adjacency
# ---------
# These functions relate to immediate adjacency.


def adjacent(node1, node2, graph):
    """Check if the nodes are **adjacent** (an edge connects them)."""
    return edge(node1, node2) in edge_set(graph)

assert adjacent(1, 2, K2)
assert adjacent(1, 4, C4)
assert adjacent(1, 1, LOOP)
assert not adjacent(1, 1, K2)
assert not adjacent(1, 3, C4)


def incident(edge, node):
    """
    Check if edge and node are **incident** to each other (edge connects node).
    """
    return node in endpoints(edge)

assert incident(edge(1, 2), 1)
assert incident(edge(1, 2), 2)
assert not incident(edge(1, 2), 3)


def degree(node, graph):
    """Return the node's **degree** (number of incident edge tips).

    Loops count twice, so that each edge contributes 2 to the degree sum.
    """
    degree = 0
    for edge in edge_set(graph):
        if incident(edge, node):
            degree = degree + 1
            if is_loop(edge):
                degree = degree + 1
    return degree

assert degree(1, N1) == 0
assert degree(1, LOOP) == 2
assert degree(1, C4) == 2
assert degree(6, PETERSEN) == 3

# - Explain why all nodes in a cycle graph have degree 2.
# - What is the degree of each node in K*n*? Write some unit tests to confirm.
# - Which woman or women attended the most events?
#   Which event or events were attended by most women?


def k_regular(graph):
    """
    Return k if graph is **k-regular** (all nodes have degree k), otherwise -1.
    """
    k = -1
    for node in node_set(graph):
        if k == -1:
            k = degree(node, graph)
        elif k != degree(node, graph):
            return -1
    return k

assert k_regular(EMPTY) == -1
assert k_regular(N2) == 0
assert k_regular(LOOP) == 2
assert k_regular(C3) == 2

# - Write a function to check if a graph is a null graph.
# - Write a function to check if a graph is **cubic**, i.e. 3-regular.
#   Include a unit test showing that the Petersen graph is cubic.


def degree_sequence(graph):
    """
    Return the graph's **degree sequence** (descending list of node degrees).
    """
    degrees = [degree(node, graph) for node in node_set(graph)]
    degrees.sort(reverse=True)
    return degrees

assert degree_sequence(EMPTY) == []
assert degree_sequence(N2) == [0, 0]
assert degree_sequence(K3) == [2, 2, 2]
assert degree_sequence(ARPANET) == [4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2]

# - Write a function to check if a graph is **regular**,
#   i.e. if all its nodes have the same degree. Implement 2 versions,
#   using `k_regular()` and `degree_sequence()`.


def degree_distribution(graph):
    """
    Return a list: the n-th value is the probability that a node has degree n.
    """
    degrees = degree_sequence(graph)
    if degrees:
        order = len(degrees)
        # n-th value is number of nodes with degree n by total number of nodes
        return [degrees.count(n) / order for n in range(max(degrees)+1)]
    else:
        return []

assert degree_distribution(EMPTY) == []
assert degree_distribution(N2) == [1]
assert degree_distribution(LOOP) == [0, 0, 1]

# - Rewrite `degree_distribution()` to make it more efficient.


def mean_degree(graph):
    """Return the graph's **average degree** (mean of all node degrees)."""
    degrees = degree_sequence(graph)
    return sum(degrees) / len(degrees) if degrees else 0

assert mean_degree(EMPTY) == 0
assert mean_degree(N2) == 0
assert mean_degree(C3) == 2

# - Rewrite `mean_degree()` to make it much more efficient.
#   Hint: each edge contributes 2 to the total degree count.


def neighbours(node, graph):
    """Return the set of the node's **neighbours** (adjacent nodes)."""
    return {node2 for node2 in node_set(graph) if adjacent(node, node2, graph)}

assert neighbours(1, LOOP) == {1}
assert neighbours(1, PETERSEN) == {2, 5, 6}


def neighbourhood(node, graph):
    """
    Return the node's **open neighbourhood** (subgraph induced by neighbours).
    """
    return subgraph(neighbours(node, graph), graph)

assert neighbourhood(1, N2) == EMPTY
assert neighbourhood(3, K3) == K2

# - Write a function that returns a node's **closed neighbourhood**,
#   which includes the node itself.


def local_clustering_coefficient(node, graph):
    """
    Return node's **local clustering coefficient** (density of neighbourhood).
    """
    return density(neighbourhood(node, graph))

assert local_clustering_coefficient(1, K3) == 1
assert local_clustering_coefficient(3, PETERSEN) == 0


def global_clustering_coefficient(graph):
    """
    Return the graph's **clustering coefficient**, a measure of transitivity.

    It measures the mean probability that A and C are connected directly
    if they are connected through a common node.
    """
    # Count the number of paths of length 2 (ABC) and those 'closed' (ABCA).
    all_paths_2 = 0
    closed_paths_2 = 0
    for a in node_set(graph):
        for b in node_set(graph) - {a}:
            for c in node_set(graph) - {a, b}:
                if adjacent(a, b, graph) and adjacent(b, c, graph):
                    all_paths_2 += 1
                    if adjacent(a, c, graph):
                        closed_paths_2 += 1
    return closed_paths_2 / all_paths_2 if all_paths_2 else 0

# Any cycle graph has coefficient 0.
assert global_clustering_coefficient(C5) == 0
# Any complete graph has coefficient 1.
assert global_clustering_coefficient(K3) == 1
# The following graph has 6 closed paths of length 2 in the K3 triangle,
# and 4 open paths of length 2 due to the additional edge.
assert global_clustering_coefficient(add_edge(3, 4, K3)) == 0.6


# Independence
# ------------
# These functions relate to non-adjacency.

# - Write a function to check if a node is **isolated**
#   (not connected to any node). Implement 4 versions,
#   using `incident()`, `adjacent()`, `degree()`, and `neighbours()`.
#   Which version is the most efficient?
# - Write a function to check if a graph is **bipartite**, i.e.
#   if its nodes can be separated into two disjoint subsets so that
#   each edge connects nodes in different subsets.
#   Write tests to confirm that the utility and women graphs are bipartite.
# - Explain why bipartite graphs are simple and have clustering coefficient 0.


# Paths and distances
# -------------------
# These functions relate to paths within the graph and their distances.


def geodesics(node1, node2, graph):
    """
    Return a list of all geodesics `[node1,...,node2]` between the two nodes.

    A **geodesic** is a shortest path.
    A **path** is a node sequence n_1, n_2, ... with n_i adjacent to n_i+1.
    Return each path as a list of nodes. Return `[]` if there is no path.
    Return `[[node1]]` if both nodes are the same.
    """
    paths = []
    candidates = [[node1]]
    while candidates:
        path = candidates.pop(0)
        if paths and len(path) > len(paths[0]):
            break
        last = path[-1]
        if last == node2:
            paths.append(path)
        else:
            for neighbour in neighbours(last, graph):
                if neighbour not in path:
                    candidates.append(path + [neighbour])
    return paths

assert geodesics(1, 2, N2) == []
assert geodesics(1, 1, LOOP) == [[1]]
assert geodesics(3, 1, K3) == [[3, 1]]
assert geodesics(1, 4, C5) == [[1, 5, 4]]
# Paths can be listed in any order: sort them to compare with the solution.
assert sorted(geodesics(1, 3, C4)) == [[1, 2, 3], [1, 4, 3]]
# Example of 2 geodesics going through same node.
# Taken from Newman's "Networks: an introduction", page 187.
assert sorted(geodesics(3, 5, add_edge(1, 5, C4))) \
    == [[3, 2, 1, 5], [3, 4, 1, 5]]

# - Explain how `geodesics` works.
# - Write a simplified `geodesic` (singular!) function that returns
#   the first shortest path found, or `[]` if there is none.


def distance(node1, node2, graph):
    """
    Return the **distance** (length of shortest path) between the nodes.

    The length of a path is the number of its edges.
    If there's no path, the distance is infinite.
    """
    paths = geodesics(node1, node2, graph)
    return len(paths[0]) - 1 if paths else float("infinity")

assert distance(2, 4, C4) == 2
# In a complete graph, the distance between different nodes is always 1.
assert distance(3, 2, K3) == 1
# In a null graph, the distance between different nodes is always infinite.
assert distance(1, 2, N2) == float("infinity")
# The distance between the same node is always 0.
assert distance(1, 1, N2) == 0

# - Make `distance` more efficient by calling `geodesic`, if you've written it.


def diameter(graph):
    """Return the graph's **diameter** (longest of all pairwise distances)."""
    if order(graph) == 0:
        return float("infinity")
    nodes = node_set(graph)
    return max([distance(n1, n2, graph) for n1 in nodes for n2 in nodes])

assert diameter(EMPTY) == float("infinity")
assert diameter(N1) == 0
assert diameter(N2) == float("infinity")
assert diameter(K3) == 1
assert diameter(C5) == 2

# - Explain how `diameter()` works.
# - Make `diameter()` more efficient.
#   Hint: the distance from A to B is the distance from B to A.
# - Write a function to return a node's **eccentricity**, i.e.
#   its maximum distance to any other node.
# - Rewrite `diameter()` using `eccentricity()`.
# - Write a function to return a graph's **center**, the set of
#   nodes with minimum eccentricity.


# Connectivity
# ------------
# These functions relate to connectivity, which extends adjacency.


def connected(graph):
    """
    Check if graph is **connected** (there's a path between any pair of nodes).
    """
    for node1 in node_set(graph):
        for node2 in node_set(graph):
            if not geodesics(node1, node2, graph):
                return False
    return True

assert connected(EMPTY)
assert connected(N1)
assert not connected(N2)
assert connected(C5)

# - Make `connected()` more efficient.
# - Add a function to check if a network is resilient to single edge failure,
#   i.e. if after removing any one edge, the remaining graph is connected.
#   Write a unit test showing that ARPANET is resilient.
# - Add a function to check if a network is resilient to single node failure,
#   i.e. if after removing any one node, the remaining graph is connected.
#   Show that ARPANET is also resilient in this sense.
#   Hint: first complete the earlier node removal exercise.
# - Write a function to check if a graph is a complete graph.
#   Include a unit test with a 2-regular graph of order 3 that is not K3.
# - Write a function to check if a graph is a cycle graph.
#   Include a unit test with a 2-regular graph of order 2 that is not C2.


def components(graph):
    """
    Return a list of the graph's **components** (maximal connected subgraphs).
    """
    components = []
    visited = set()
    for node in node_set(graph):
        component = set()
        to_visit = {node}
        while to_visit:
            node = to_visit.pop()
            if node not in visited:
                visited.add(node)
                component.add(node)
                to_visit = to_visit | neighbours(node, graph)
        if component:
            components.append(subgraph(component, graph))
    return components

assert components(EMPTY) == []
assert components(K3) == [K3]
assert [order(c) for c in components(N2)] == [1, 1]
assert [size(c) for c in components(N2)] == [0, 0]

# - Explain how `components()` works. Explain its unit tests.
# - Rewrite `connected()` using `components()`. Which is more efficient?


def giant_component(graph):
    """
    Return the **giant component** (component with majority of graph's nodes).

    Return the empty graph if there is no giant component.
    """
    half = order(graph) / 2
    for component in components(graph):
        if order(component) > half:
            return component
    return null(0)

assert giant_component(EMPTY) == EMPTY
assert giant_component(C3) == C3
assert giant_component(N2) == EMPTY

# - Rewrite `giant_component()` to return the largest component (most nodes).


def k_core(k, graph):
    """
    Return the *k*-**core** (maximal subgraph of nodes with degree k or more).
    """
    # Keep removing nodes with degree < k until no longer possible
    while True:
        nodes = {node for node in node_set(graph) if degree(node, graph) >= k}
        if nodes == node_set(graph):
            return graph
        graph = subgraph(nodes, graph)

# Any k-core of the empty graph is empty.
assert k_core(-1, EMPTY) == EMPTY
# The 0-core of any graph is itself.
assert k_core(0, C5) == C5
# If k is larger than the maximum degree, the k-core is empty.
assert k_core(3, C5) == EMPTY
# If the graph is k-regular, its k-core is the whole graph.
assert k_core(k_regular(K3), K3) == K3
# In the early Arpanet, no node had 3 neighbours with 3 neighbours themselves.
assert k_core(3, ARPANET) == EMPTY

# - Rewrite `k_core` using the `delete_node` function.


# Centrality
# ----------
# These functions relate to how 'important' a node is.
# The degree is one such metric.
#
# - Which node or nodes in the kite and Arpanet graphs
#   have the highest degree centrality?


def betweenness(node, graph):
    """Return the **betweenness centrality** of the node.

    It captures the extent to which shortest paths go through the node.
    """
    centrality = 0
    for node1 in node_set(graph):
        for node2 in node_set(graph):
            paths = geodesics(node1, node2, graph)
            if paths:
                on_paths = sum(1 for path in paths if node in path)
                centrality = centrality + on_paths / len(paths)
    return centrality

# In K3, a node A is on all 9 paths except on B-B, C-C, B-C and C-B.
assert betweenness(1, K3) == 5
# In a null graph, any node is only on the path with itself.
assert betweenness(1, N2) == 1
# In C5, a node is on the 4+4 paths that start and end at the node,
# on the 1 path to itself, and on the 2 paths between its neighbours.
assert betweenness(1, C5) == 11
# Utah is more central to connecting East and West coasts than CASE.
assert betweenness("Utah", ARPANET) > betweenness("CASE", ARPANET)

# - Which node or nodes in the kite and Arpanet graphs
#   have the highest betweenness centrality?
# - Some alternative definitions, e.g. on Wikipedia,
#   don't consider geodesics that start or end at the node.
#   Change the function and unit tests accordingly.


def ego_graph(ego, d, graph):
    """
    Return the subgraph induced by all nodes up to distance d of the ego node.
    """
    nodes = {n for n in node_set(graph) if distance(ego, n, graph) <= d}
    return subgraph(nodes, graph)

# If the distance is set to the diameter, the ego graph is the whole graph.
assert ego_graph(5, diameter(PETERSEN), PETERSEN) == PETERSEN
# For any simple graph, an ego graph with distance zero is just the ego node.
assert ego_graph(5, 0, PETERSEN) == renumber(N1, 5)
# In a complete graph, any ego graph with distance > 0 is the whole graph.
assert ego_graph(1, 1, K3) == K3


# Input / Output
# --------------
# These functions help read and write graphs
# in GML (Graph Modelling Language) format. For more details, see the GML
# [specification](https://www.fim.uni-passau.de/fileadmin/files/lehrstuhl/brandenburg/projekte/gml/gml-technical-report.pdf).


def to_gml(graph):
    """Return a string representing the graph in GML format.

    The nodes of the graph must be represented by integers.
    """
    # Concatenating many strings is slow; join a list of strings instead.
    lines = ["graph [", "  directed 0"]  # in GML, 0 means false
    for node in node_set(graph):
        lines.append("  node [id {}]".format(node))
    for edge in edge_set(graph):
        (node1, node2) = endpoints(edge)
        lines.append("  edge [source {} target {}]".format(node1, node2))
    lines.append("]")
    # Join the lines, separated by newlines.
    return "\n".join(lines)

assert to_gml(EMPTY) == """\
graph [
  directed 0
]"""
assert to_gml(LOOP) == """\
graph [
  directed 0
  node [id 1]
  edge [source 1 target 1]
]"""
assert to_gml(N2) == """\
graph [
  directed 0
  node [id 1]
  node [id 2]
]"""


def from_gml(text):
    """
    Return the list of graphs read from `text`, a string in GML format.

    This function assumes that the input follows the GML specification,
    provides unique integer ids even for isolated nodes, and
    defines one or more graphs.
    This function ignores anything other than node ids and edge endpoints.
    This means directed graphs are read as undirected graphs,
    node labels and edge weights are discarded, etc.
    If an edge endpoint (integer) is an unknown node id, the node is created.
    """
    # Define the grammar with [pyparsing](http://pyparsing.wikispaces.com).
    # Don't use `from pyparsing import *` as it adds many constants
    # to the generated documentation.

    from pyparsing import (
        srange, oneOf, Forward, Optional, Suppress, Word, ZeroOrMore,
        dblQuotedString, pythonStyleComment
    )

    digit = srange("[0-9]")
    sign = Optional(oneOf("+ -"))
    mantissa = Optional("E" + sign + digit)

    # `Word(x)` is a sequence of one or more characters from the set x.
    digits = Word(digit)
    integer = sign + digits
    real = sign + Optional(digits) + "." + Optional(digits) + mantissa

    # For simplicity, use pyparsing's string with double-quotes,
    # hoping that it is a generalisation of GML's definition of a string.
    string = dblQuotedString

    # A GML file is a list of key-value pairs, where a value may be a list.
    # To handle this recursive definition, we delay what a pair is.
    pair = Forward()
    list = ZeroOrMore(pair)
    # A file may have comments, which are as in Python. Ignore them.
    list.ignore(pythonStyleComment)

    # `Word(x, y)` is 1 character from x followed by 0 or more from y.
    key = Word(srange("[a-zA-Z]"), srange("[a-zA-Z0-9]"))
    # `Suppress(x)` matches x but doesn't put it in the list of parsed tokens.
    listValue = Suppress("[") + list + Suppress("]")
    value = real | integer | string | listValue

    # The mandatory key-value pairs for graphs are as follows.
    graph = Suppress("graph") + listValue
    node = Suppress("node") + listValue
    anEdge = "edge" + listValue     # to avoid conflict with edge() function
    id = Suppress("id") + integer
    source = Suppress("source") + integer
    target = Suppress("target") + integer
    # First try to parse graph-specific key-value pairs; otherwise ignore pair.
    pair <<= graph | node | anEdge | id | source | target | Suppress(key+value)

    # The above suppressions lead to the GML string
    # `'graph [ node [id 1 label "ego"] edge [source 1 target 1 weight 0.5] ]'`
    # being parsed into the list of tokens
    # `["1", "edge", "1", "1"]`,
    # which is converted by the following functions into a graph.

    def to_int(text, position, tokens):
        # Convert parsed integer tokens to integers, e.g. `["1"]` to `1`.
        return int(tokens[0])

    def to_edge(text, position, tokens):
        # Assuming the above conversion was done,
        # convert `["edge", a, b]` to an edge incident to a and b.
        return edge(tokens[1], tokens[2])

    def to_graph(text, position, tokens):
        # `tokens` is now a list of integers and edges, in any order.
        nodes = set()
        edges = set()
        for token in tokens:
            # If the token is an integer, it's a node id.
            if isinstance(token, int):
                nodes.add(token)
            else:
                edges.add(token)
        return network(edges, nodes)

    # Do the conversions as soon as the respective tokens are parsed.
    integer.setParseAction(to_int)
    anEdge.setParseAction(to_edge)
    graph.setParseAction(to_graph)

    # Parse the text with the main grammar rule.
    # Return the result as a list, not as a pyparsing object.
    return list.parseString(text).asList()

# Check that converting to GML and back returns the original graph.
assert [EMPTY] == from_gml(to_gml(EMPTY))
assert [N1] == from_gml(to_gml(N1))
assert [LOOP] == from_gml(to_gml(LOOP))
assert [K3] == from_gml(to_gml(K3))
assert [EMPTY, LOOP] == from_gml('''
graph [directed 1]  # an empty graph that is (strangely enough) directed
graph [ node [id 1 label "ego"] edge [source 1 target 1 weight 0.5] ]
''')

# - Write helper functions that read from and write to a file.


# Projects
# --------
# These exercises require using or changing large parts of the library.
#
# - Write a program that imports SINPLE and reports various properties
# of the Arpanet network, e.g. the most highly connected nodes.
#
# - Write a program that imports SINPLE, reads a network from a GML file,
# e.g. provided by [Mark Newman](http://www-personal.umich.edu/~mejn/netdata/),
# and analyses it.
#
# - For example, download the Karate Club graph and check that it has
# two highly connected members, the student founder and
# the instructor, who are not friends of each other,
# and that most members are not friends of both.
# This conflict led to the split into two clubs.
#
# - Write a program that imports SINPLE, creates several random graphs
# with fixed `n` and increasing `p`, and reports their properties.
# Check if the mean degree is the theoretical expected value `p*(n-1)`.
# Check if a giant component emerges when `p > 1/n`.
# Check if the network becomes connected when `p > math.log(n) / n`.
#
# - Functions like `degree(node, graph)` assume that `node` is in `graph`.
# What happens if it isn't?
# Decide whether to return a sensible value or to raise an exception.
# If the former, add unit tests.
#
# - Extend SINPLE to support **weighted edges**,
# i.e. edges with an associated number called the edge's weight.
# An unweighted network can be seen
# as a weighted network where each edge has the default weight 1.
# Hints: Add a default parameter to `edge()`.
# Add a function to return the weight of an edge.
# Modify `geodesics()` to return the paths with the smallest total weight.
#
# - Extend SINPLE to store node and edge data, like node and edge labels,
# edge weight, etc. Change the GML functions to read and write such data.
# One possibility is to represent each node by a pair `(id, data)`
# and each edge by a pair `(endpoints, data)`,
# where `data` is a dictionary with a key-value pair for each datum.
# Another option is to put all data in a single dictionary,
# indexed by nodes and edges, of dictionaries, i.e. nodes and edges
# remain as now and the graph becomes a triple `(nodes, edges, data)`.
#
# - Make SINPLE work with Python 2.7.
#
# - Rewrite SINPLE using classes.
