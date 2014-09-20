"""
SImple Network Python Library for Education (SINPLE)

Michel Wermelinger

Many systems can be represented as networks, e.g. transport connections,
gene regulation, social interactions, the World Wide Web. Network science is an
important field to analyse and predict the behaviour of such systems.

SINPLE does not support analysis of very large networks;
there are various state-of-the-art libraries for that.
Instead SINPLE's aims are pedagogical:

- it illustrates basic Python constructs: loops, functions, tuples, sets, etc;
- it illustrates simple forms of
[unit testing](http://en.wikipedia.org/wiki/Unit_testing) and
[preconditions](http://en.wikipedia.org/wiki/Precondition) using assertions;
- it explains basic network concepts through executable definitions
(the functions) and concrete examples (the unit tests);
- it includes over 40 exercises;
- it is extensively documented in a simple [literate
programming](http://en.wikipedia.org/wiki/Literate_programming) style.

SINPLE is written for Python 3.

Running [Pycco](http://fitzgen.github.io/pycco/) on the source code
(`pycco -d . sinple.py`) generates `sinple.html`,
a side-by-side view of documentation and code.

The latest version of SINPLE can be found on http://tiny.cc/sinple.
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

# Network representation and creation
# -----------------------------------
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
# SINPLE represents a network by a pair `(nodes, edges)`,
# where `nodes` is a Python set of node names (e.g. integers or strings)
# and `edges` is a Python set of edges.
# Each edge is represented by a Python set of node names `{node1, node2}`,
# because the edge has no direction.
#
# The functions in this section encapsulate this representation
# from the rest of the library.


def edge(node1, node2):
    """Create an edge connecting the two nodes."""
    return frozenset({node1, node2})

assert edge(1, 2) == edge(2, 1)

# - Consult the Python documentation about sets.
#   Why is an edge represented by a frozen, i.e. immutable, set?
#   Hint: read again the definition of graph at the start of this section.


def endpoints(edge):
    """Return the set of the edge's **endpoints** (the nodes it connects)."""
    return set(edge)

assert endpoints(edge(2, 1)) == {1, 2}
assert endpoints(edge(1, 1)) == {1}
assert endpoints(edge("me", "my friend")) == {"my friend", "me"}


def edge_set(graph):
    """Return the set of the graph's edges."""
    (nodes, edges) = graph
    return edges


def node_set(graph):
    """Return the set of the graph's nodes."""
    (nodes, edges) = graph
    return nodes


def add_edges(edges, graph):
    """
    Create a new graph by adding the set of edges and their endpoints to graph.
    """
    nodes = set()
    for edge in edges:
        nodes = nodes | endpoints(edge)
    return (node_set(graph) | nodes, edge_set(graph) | edges)

# - Explain how `add_edges()` works. Why is `set()` used instead of `{}`?


def network(edges, nodes=set()):
    """Create a graph with the given nodes, edges, and their endpoints."""
    return add_edges(edges, (nodes, set()))

assert network({edge(1, 2)}) == network({edge(2, 1)}, {1, 2})


def delete_edges(edges, graph):
    """Create a new graph by removing the set of edges from graph."""
    return (node_set(graph), edge_set(graph) - edges)

assert delete_edges({edge(1, 2)}, network({edge(2, 1)})) \
    == network(set(), {1, 2})

# - Write a function to remove a set of nodes from a graph.


# ### Example networks
#
# These networks are used in further unit tests.


def null(n):
    """Create a **null graph** with n nodes numbered 1 to n and no edges."""
    return network(set(), {x for x in range(1, n+1)})

# The **empty graph** has no nodes and hence no edges.
EMPTY = null(0)
N1 = null(1)
N2 = null(2)

# The smallest non-simple graph has a single loop.
LOOP = network({edge(1, 1)})

# A **complete graph** K*n* connects each node to all other *n-1* nodes.
K2 = network({edge(1, 2)})
K3 = network({edge(1, 2), edge(1, 3), edge(2, 3)})

# - Define K1.
# - Write a function to create a complete graph, given the number of nodes.
#   Compare the output to the 'manually' created K1, K2, K3.

# A **cycle graph** C*n* has *n* nodes connected in a 'round-robin' fashion.
# Cycle graphs can be drawn as regular shapes: triangle, square, pentagon, etc.
C3 = network({edge(1, 2), edge(2, 3), edge(3, 1)})
C4 = add_edges({edge(3, 4), edge(4, 1)}, delete_edges({edge(3, 1)}, C3))
C5 = add_edges({edge(4, 5), edge(5, 1)}, delete_edges({edge(4, 1)}, C4))

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


# Edge properties
# ---------------


def incident(edge, node):
    """
    Check if edge and node are **incident** to each other (edge connects node).
    """
    return node in endpoints(edge)

assert incident(edge(1, 2), 1)
assert incident(edge(1, 2), 2)
assert not incident(edge(1, 2), 3)


def is_loop(edge):
    """Check if edge is a **loop** (connects a node to itself)."""
    return len(endpoints(edge)) == 1

assert is_loop(edge("A", "A"))
assert not is_loop(edge("A", "B"))


# Node properties
# ---------------


def adjacent(node1, node2, graph):
    """Check if the nodes are **adjacent** (an edge connects them)."""
    return edge(node1, node2) in edge_set(graph)

assert adjacent(1, 2, K2)
assert adjacent(1, 4, C4)
assert adjacent(1, 1, LOOP)
assert not adjacent(1, 1, K2)
assert not adjacent(1, 3, C4)


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


def neighbours(node, graph):
    """Return the set of the node's **neighbours** (adjacent nodes)."""
    return {node2 for node2 in node_set(graph) if adjacent(node, node2, graph)}

assert neighbours(1, LOOP) == {1}
assert neighbours(1, PETERSEN) == {2, 5, 6}

# - Write a function to check if a node is **isolated**
#   (not connected to any node). Implement 4 versions,
#   using `incident()`, `adjacent()`, `degree()`, and `neighbours()`.
# - Which version is the most efficient?


# Further graph construction
# --------------------------


def subgraph(nodes, graph):
    """
    Return the graph's **subgraph induced** by the set of nodes.

    It includes those edges of graph that connect the given nodes.
    """
    edges = {edge(n1, n2) for n1 in nodes for n2 in nodes}
    return network(edges & edge_set(graph), nodes)

assert subgraph({1}, N2) == N1
assert subgraph({1, 2, 3, 4, 5}, PETERSEN) == C5

# - Explain how `subgraph()` works.
# - Make it more efficient.
# - Write a function to check if *g1* is a **subgraph** of *g2*,
#   i.e. if the nodes and edges of *g1* are included in those of *g2*.


def neighbourhood(node, graph):
    """
    Return the node's **open neighbourhood** (subgraph induced by neighbours).
    """
    return subgraph(neighbours(node, graph), graph)

assert neighbourhood(1, N2) == EMPTY
assert neighbourhood(3, K3) == K2

# - Write a function that returns a node's **closed neighbourhood**,
#   which includes the node itself.


def renumber(graph, n=1):
    """Create a new graph by renaming the graph's nodes as n, n+1, n+2, etc."""
    number = n
    map = {}
    for node in node_set(graph):
        map[node] = number
        number = number + 1
    new_graph = network(set(), set(range(n, number)))
    for e in edge_set(graph):
        nodes = endpoints(e)
        node1 = map[nodes.pop()]
        node2 = map[nodes.pop()] if nodes else node1
        new_graph = add_edges({edge(node1, node2)}, new_graph)
    return new_graph

assert renumber(EMPTY, -1) == EMPTY
assert renumber(N2, 5) == network(set(), {5, 6})
assert renumber(network({edge("A", "B")})) == network({edge(1, 2)})


def gml(graph):
    """Return a string representing the graph in GML format.

    The nodes of the graph must be represented by integers.
    """
    # Concatenating many strings is slow; join a list of strings instead.
    lines = ["graph [", "  directed 0"]  # in GML, 0 means false
    for node in node_set(graph):
        lines.append("  node [id {}]".format(node))
    for edge in edge_set(graph):
        nodes = endpoints(edge)
        # Take one of the nodes as the edge's source.
        source = nodes.pop()
        # Take the other as target, unless set is now empty (edge is a loop).
        target = nodes.pop() if nodes else source
        lines.append("  edge [source {} target {}]".format(source, target))
    lines.append("]")
    # Join the lines, separated by newlines.
    return "\n".join(lines)

assert gml(LOOP) == """\
graph [
  directed 0
  node [id 1]
  edge [source 1 target 1]
]"""
assert gml(N2) == """\
graph [
  directed 0
  node [id 1]
  node [id 2]
]"""

# - Write a function that given a file name and a graph,
#   writes the graph in GML format to that file.


# Graph properties
# ----------------


def size(graph):
    """Return the graph's **size** (number of edges)."""
    return len(edge_set(graph))

assert size(EMPTY) == 0
assert size(N1) == 0
assert size(C3) == 3
assert size(PETERSEN) == 15


def order(graph):
    """Return the graph's **order** (number of nodes)."""
    return len(node_set(graph))

assert order(EMPTY) == 0
assert order(N1) == 1
assert order(C3) == 3
assert order(PETERSEN) == 10


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

# - Write a function to check if a graph is **bipartite**, i.e.
#   its nodes can be separated into two parts so that each edge
#   connects nodes in different parts.
#   Include a unit test showing that the utility graph is bipartite.


def density(graph):
    """
    Return the graph's **density** (ratio of actual to potential edges).

    It's the size divided by the size of a complete graph of the same order.
    The density of the empty graph is set to zero.
    """
    o = order(graph)
    return size(graph) / (o * (o - 1) / 2) if o > 0 else 0

assert density(EMPTY) == 0
# Every null graph has density 0.
assert density(N2) == 0
# Every complete graph has density 1.
assert density(K3) == 1
# A graph with loops can have density > 1.
assert density(add_edges({edge(1, 1)}, K3)) > 1


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


def degree_distribution(graph):
    """Return a list of length order(graph)+2 of a node's degree probabilities.

    The n-th element is the probability of a node having degree n.
    """
    o = order(graph)
    max_degree = o + 1
    if o == 0:
        return [0, 0]
    distribution = [0 for d in range(max_degree + 1)]
    for node in node_set(graph):
        distribution[degree(node, graph)] += 1
    for d in range(max_degree + 1):
        distribution[d] /= o
    return distribution

assert degree_distribution(EMPTY) == [0, 0]
assert degree_distribution(N2) == [1, 0, 0, 0]
assert degree_distribution(LOOP) == [0, 0, 1]

# - Explain how `degree_distribution()` works.
#   Why is the list never longer than `order(graph) + 2`?
# - Change `degree_distribution()` so that the list has no trailing zeros.
# - Write a function to compute the degree distribution from a degree sequence.


def mean_degree(graph):
    """Return the graph's **average degree** (mean of all node degrees)."""
    degrees = degree_sequence(graph)
    return sum(degrees) / len(degrees) if degrees else 0

assert mean_degree(EMPTY) == 0
assert mean_degree(N2) == 0
assert mean_degree(C3) == 2

# - Rewrite `mean_degree()` to make it much more efficient.
#   Hint: each edge contributes 2 to the total degree count.


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

# - Write a function to check if a graph is **regular**,
#   i.e. if all its nodes have the same degree.
# - Write a function to check if a degree sequence forms a regular graph.
# - Write a function to check if a graph is a null graph.
# - Write a function to check if a graph is **cubic**, i.e. 3-regular.
#   Include a unit test showing that the Petersen graph is cubic.
# - Write a function to check if a graph is a complete graph.
#   Include a unit test with a 2-regular graph of order 3 that is not K3.
# - Write a function to check if a graph is a cycle graph.
#   Include a unit test with a 2-regular graph of order 2 that is not C2.


def local_clustering_coefficient(node, graph):
    """
    Return node's **local clustering coefficient** (density of neighbourhood).
    """
    return density(neighbourhood(node, graph))

assert local_clustering_coefficient(1, K3) == 1
assert local_clustering_coefficient(3, PETERSEN) == 0


def shortest_path(node1, node2, graph):
    """
    Return a shortest path [node1, ..., node2] between the two nodes.

    A **path** is a node sequence n_1, n_2, ... with n_i adjacent to n_i+1.
    Return [] if there is no path. Return [node1] if both nodes are the same.
    """
    paths = [[node1]]
    visited = []
    while paths != []:
        path = paths.pop(0)
        last = path[-1]
        visited.append(last)
        if last == node2:
            return path
        for neighbour in neighbours(last, graph):
            if neighbour not in visited:
                paths.append(path + [neighbour])
    return []

assert shortest_path(1, 2, N2) == []
assert shortest_path(1, 1, LOOP) == [1]
assert shortest_path(3, 1, K3) == [3, 1]
assert shortest_path(1, 4, C5) == [1, 5, 4]

# - Explain how `shortest_path` works.
# - Write a function to return *all* shortest paths between two nodes.


def distance(node1, node2, graph):
    """
    Return the **distance** (length of shortest path) between the nodes.

    The length of a path is the number of its edges.
    If there's no path, the distance is infinite.
    """
    path = shortest_path(node1, node2, graph)
    return len(path) - 1 if path else float("infinity")

assert distance(1, 1, C3) == 0
assert distance(2, 4, C4) == 2
assert distance(3, 2, K3) == 1
assert distance(1, 2, N2) == float("infinity")


def diameter(graph):
    """Return the graph's **diameter** (longest of all pairwise distances)."""
    nodes = node_set(graph)
    distances = [distance(n1, n2, graph) for n1 in nodes for n2 in nodes]
    return max(distances) if distances else float("infinity")

assert diameter(EMPTY) == float("infinity")
assert diameter(N1) == 0
assert diameter(N2) == float("infinity")
assert diameter(K3) == 1
assert diameter(C5) == 2

# - Explain how `diameter()` works.
# - Rewrite `diameter()` to make it twice as fast.
#   Hint: the distance from A to B is the distance from B to A.


def connected(graph):
    """
    Check if graph is **connected** (there's a path between any pair of nodes).
    """
    for node1 in node_set(graph):
        for node2 in node_set(graph):
            if shortest_path(node1, node2, graph) == []:
                return False
    return True

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

assert components(K3) == [K3]
assert [order(c) for c in components(N2)] == [1, 1]
assert [size(c) for c in components(N2)] == [0, 0]

# - Explain how `components()` works. Explain its unit tests.
# - Rewrite `connected()` using `components()`.


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

assert giant_component(C3) == C3
assert giant_component(N2) == null(0)


# Import a random number generator with uniform distribution
from random import uniform


def random(n, p):
    """
    Return a simple **random graph** with n nodes and edge probability p.

    The two numbers can't be negative, and p can't be larger than 1.
    """
    graph = null(n)
    for node1 in range(1, n+1):
        for node2 in range(node1+1, n+1):
            if uniform(0, 1) < p:
                graph = add_edges({edge(node1, node2)}, graph)
    return graph

assert order(random(4, 0)) == 4
assert size(random(4, 0)) == 0
assert size(random(4, 1)) == 6
assert simple(random(4, 0.5))

# - Write a function to create a simple random graph of given order and size.


# Projects
# --------
#
# These are exercises that require using or changing the whole library.
#
# - Functions like `degree(node, graph)` assume that `node` is in `graph`.
# What happens if it isn't?
# Decide whether to return some value in those cases or to raise an exception.
# If the former, add unit tests.
#
# - Write a program that imports SINPLE and reports various properties
# of the Arpanet network, e.g. the most highly connected nodes.
#
# - Write a program that imports SINPLE, creates several random graphs
# with fixed `n` but increasing `p`, and reports their properties.
# Check if the mean degree is the theoretical expected value `p*(n-1)`.
# Check if a giant component emerges when `p > 1/n`.
# Check if the network becomes connected when `p > math.log(n) / n`.
#
# - Extend SINPLE to support **weighted edges**,
# i.e. edges that have an associated number, called the edge's weight.
# An unweighted network can be seen
# as a weighted network where each edge has the default weight 1.
# Add a default parameter to `edge()`.
# Add a function to return the weight of an edge.
# Modify `shortest_path()` to return the path with the smallest total weight.
#
# - Make SINPLE work with Python 2.7.
#
# - Rewrite SINPLE using classes.
