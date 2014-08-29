"""
SImple Network Python Library for Education (SINPLE)

Michel Wermelinger

Many systems can be represented as networks, e.g. transport connections,
gene regulation, social interactions, the World Wide Web.
[Network science](http://en.wikipedia.org/wiki/Network_science) is an
important field to analyse and predict the behaviour of such systems.

SINPLE does *not* enable analysis of large networks;
there are various state-of-the-art libraries for that.
Instead SINPLE's aims are pedagogical:

- it illustrates basic Python constructs: loops, functions, tuples, sets, etc;
- it illustrates simple forms of
[unit testing](http://en.wikipedia.org/wiki/Unit_testing) and
[preconditions](http://en.wikipedia.org/wiki/Precondition) using assertions;
- it explains basic network concepts through executable definitions
(the functions) and concrete examples (the unit tests);
- it includes small exercises and larger projects;
- it is extensively documented in a simple [literate
programming](http://en.wikipedia.org/wiki/Literate_programming) style.

Running [Pycco](http://fitzgen.github.io/pycco/) on the source code
(`pycco -d . sinple.py`) generates `sinple.html`,
a side-by-side view of documentation and code.

The latest version of SINPLE can be found on http://tiny.cc/sinple.
"""

# Coding conventions
# ------------------
# Functions named with adjectives or `is_`noun return booleans.
# Functions named with nouns return network entities or metrics.
# The code passes the [PEP8 style checker](https://github.com/jcrocholl/pep8).

# Network representation and creation
# -----------------------------------
# The mathematical basis of network science is graph theory.
# Graph terminology varies among authors.
# We will mostly follow the Wikipedia [glossary of graph
# theory](http://en.wikipedia.org/wiki/Glossary_of_graph_theory).
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


def edges(graph):
    """Return the set of the graph's edges"""
    (nodes, edges) = graph
    return edges


def nodes(graph):
    """Return the set of the graph's nodes."""
    (nodes, edges) = graph
    return nodes


def edge(node1, node2):
    """Create an edge connecting the two nodes."""
    return frozenset({node1, node2})

assert edge(1, 2) == edge(2, 1)

# **EXERCISES:**
#
# - Consult the Python documentation about sets.
#   Why is an edge represented by a frozen (i.e. immutable) set?


def endpoints(edge):
    """Return the edge's **endpoints**, i.e. the set of nodes it connects."""
    return edge

assert endpoints(edge(1, 2)) == {1, 2}


def null(n):
    """Return the **null graph** with nodes numbered 1 to *n* and no edges."""
    assert n >= 0
    return ({x for x in range(1, n+1)}, set())

assert nodes(null(0)) == set()
assert edges(null(2)) == set()
assert nodes(null(2)) == {1, 2}

# **EXERCISES:**
#
# - Explain how `null()` works. Why is `set()` used instead of `{}`?


def add_edges(edgeset, graph):
    """
    Add a set of edges and their endpoints to a graph to create a new graph.
    """
    nodeset = set()
    for edge in edgeset:
        nodeset = nodeset | endpoints(edge)
    return (nodes(graph) | nodeset, edges(graph) | edgeset)

assert add_edges({edge(1, 2)}, null(0)) == add_edges({edge(2, 1)}, null(2))

# **EXERCISES:**
#
# - Explain how `add_edges()` works. Why isn't `edgseset` just named `edges`?


def delete_edges(edgeset, graph):
    """Remove a set of edges from a graph to create a new graph."""
    return (nodes(graph), edges(graph) - edgeset)

assert delete_edges({edge(1, 2)}, add_edges({edge(2, 1)}, null(0))) == null(2)

# **EXERCISES:**
#
# - Write a function to remove a set of nodes from a graph.
# - Write a function `subgraph(nodeset, graph)` that returns the **subgraph**
#   with the nodes of `graph` that are in `nodeset` and the edges of `graph`
#   that connect those nodes.


def network(edgeset):
    """Create a network from the given set of edges."""
    return add_edges(edgeset, null(0))

# ### Example networks
#
# These examples are used in further unit tests.

# The **empty graph** has no nodes and hence no edges.
EMPTY = null(0)

# The **trivial graph** has a single node and no edges.
N1 = null(1)
N2 = null(2)

# The smallest non-simple graph has a single loop.
LOOP = network({edge(1, 1)})

# A **complete graph** K*n* connects each node to all other *n-1* nodes.
K2 = network({edge(1, 2)})
K3 = network({edge(1, 2), edge(1, 3), edge(2, 3)})

# **EXERCISES:**
#
# - Define K1.
# - Write a function to create a complete graph, given the number of nodes.
#   Don't forget to write a precondition and unit tests.
#   Hint: compare the output to the 'manually' created K1, K2, K3.

# A **cycle graph** C*n* has *n* nodes connected in a 'round-robin' fashion.
# Cycle graphs can be drawn as regular shapes: triangle, square, pentagon, etc.
C3 = network({edge(1, 2), edge(2, 3), edge(3, 1)})
C4 = add_edges({edge(3, 4), edge(4, 1)}, delete_edges({edge(3, 1)}, C3))
C5 = add_edges({edge(4, 5), edge(5, 1)}, delete_edges({edge(4, 1)}, C4))

# **EXERCISES:**
#
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

# **EXERCISES:**
#
# - Draw the Petersen graph on paper and label the nodes as in `PETERSEN`.

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
    """An edge is **incident** to the nodes it connects, and vice-versa."""
    return node in endpoints(edge)

assert incident(edge(1, 2), 1)
assert incident(edge(1, 2), 2)
assert not incident(edge(1, 2), 3)


def is_loop(edge):
    """An edge is a **loop** if it connects a node to itself."""
    return len(endpoints(edge)) == 1

assert is_loop(edge("A", "A"))
assert not is_loop(edge("A", "B"))


# Node properties
# ---------------


def degree(node, graph):
    """A node's **degree** is the number of incident edge tips.

    Loops count twice, so that each edge contributes 2 to the degree sum.
    """
    degree = 0
    for edge in edges(graph):
        if incident(edge, node):
            degree = degree + 1
            if is_loop(edge):
                degree = degree + 1
    return degree

assert degree(1, N1) == 0
assert degree(1, LOOP) == 2
assert degree(1, C4) == 2
assert degree(6, PETERSEN) == 3

# **EXERCISES**:
#
# - Explain why all nodes in a cycle graph have degree 2.
# - What is the degree of each node in K*n*? Write some unit tests to confirm.


def adjacent(node1, node2, graph):
    """Two nodes are **adjacent** if an edge connects them."""
    return edge(node1, node2) in edges(graph)

assert adjacent(1, 2, K2)
assert adjacent(1, 4, C4)
assert adjacent(1, 1, LOOP)
assert not adjacent(1, 1, K2)
assert not adjacent(1, 3, C4)


def isolated(node, graph):
    """A node is **isolated** if it is not incident to any edge."""
    for edge in edges(graph):
        if incident(edge, node):
            return False
    return True

assert isolated(1, N1)
assert not isolated(1, K2)

# **EXERCISES:**
#
# - Rewrite `isolated()` using `adjacent()`.
# - Rewrite `isolated()` using `degree()`.
# - Which of the 3 versions is more efficient?


# Graph properties
# ----------------


def size(graph):
    """The graph's **size** is the number of its edges."""
    return len(edges(graph))

assert size(EMPTY) == 0
assert size(N1) == 0
assert size(C3) == 3
assert size(PETERSEN) == 15


def order(graph):
    """The graph's **order** is the number of its nodes."""
    return len(nodes(graph))

assert order(EMPTY) == 0
assert order(N1) == 1
assert order(C3) == 3
assert order(PETERSEN) == 10


def simple(graph):
    """A graph is **simple** if it has no loops."""
    for edge in edges(graph):
        if is_loop(edge):
            return False
    return True

assert simple(EMPTY)
assert simple(N1)
assert simple(PETERSEN)
assert not simple(LOOP)


def degree_sequence(graph):
    """A graph's **degree sequence** is a descending list of nodes' degrees."""
    degrees = [degree(node, graph) for node in nodes(graph)]
    degrees.sort(reverse=True)
    return degrees

assert degree_sequence(EMPTY) == []
assert degree_sequence(N2) == [0, 0]
assert degree_sequence(K3) == [2, 2, 2]
assert degree_sequence(ARPANET) == [4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2]

# **EXERCISES:**
#
# - Write two functions to compute respectively the maximum and minimum degrees
#   of a degree sequence.
# - Write two functions to compute respectively the maximum and minimum degrees
#   of a graph. Hint: it's not always necessary to loop over all nodes.


def degree_distribution(graph):
    """Return a list of length order(graph)+2 of a node's degree probabilities.

    The n-th element is the probability of a node having degree n.
    """
    o = order(graph)
    max_degree = o + 1
    if o == 0:
        return [0, 0]
    distribution = [0 for d in range(max_degree + 1)]
    for node in nodes(graph):
        distribution[degree(node, graph)] += 1
    for d in range(max_degree + 1):
        distribution[d] /= o
    return distribution

assert degree_distribution(EMPTY) == [0, 0]
assert degree_distribution(N2) == [1, 0, 0, 0]
assert degree_distribution(LOOP) == [0, 0, 1]

# **EXERCISES:**
#
# - Explain how `degree_distribution()` works.
#   Why is the list never longer than `order(graph) + 2`?
# - Change `degree_distribution` so that trailing zeros are removed.
# - Write a function to compute the degree distribution from a degree sequence.


def mean_degree(graph):
    """Return a graph's average degree."""
    assert order(graph) > 0
    total = 0
    for node in nodes(graph):
        total = total + degree(node, graph)
    return total / order(graph)

assert mean_degree(N2) == 0
assert mean_degree(C3) == 2

# **EXERCISES:**
#
# - Rewrite `mean_degree()` to make it much more efficient.
#   Hint: each edge contributes 2 to the total degree count.


def k_regular(graph):
    """A graph is **k-regular** if all its nodes have degree k.

    Return k if the graph is k-regular, otherwise return -1.
    """
    k = -1
    for node in nodes(graph):
        if k == -1:
            k = degree(node, graph)
        elif k != degree(node, graph):
            return -1
    return k

assert k_regular(EMPTY) == -1
assert k_regular(N2) == 0
assert k_regular(LOOP) == 2
assert k_regular(C3) == 2
assert k_regular(PETERSEN) == 3

# **EXERCISES:**
#
# - Write a function to check whether a graph is **regular**,
#   i.e. whether all its nodes have the same degree.
# - Write a function to check whether a degree sequence forms a regular graph.
# - Write a function to check whether a graph is a null graph.
# - Write a function to check whether a graph is a complete graph.
#   Hint: find a 2-regular graph with 3 nodes that is not K3, to see that
#   just using `k_regular` is not enough.
# - Write a function to check whether a graph is a cycle graph.
# - Write a function to check whether a graph is **cubic**, i.e. 3-regular.


def neighbourhood(node, graph):
    """The **neighbourhood** of a node is the set of its adjacent nodes."""
    return {node2 for node2 in nodes(graph) if adjacent(node, node2, graph)}

assert neighbourhood(1, LOOP) == {1}
assert neighbourhood(1, PETERSEN) == {2, 5, 6}


def shortest_path(node1, node2, graph):
    """Return the first shortest path found between the two nodes.

    A **path** is a sequence of nodes n_1, n_2, ... with n_i adjacent to n_i+1.
    If there is no path, return `[]`, otherwise return `[node1, ..., node2]`.
    Return `[node1]` if node1 and node2 are the same.
    """
    paths = [[node1]]
    visited = []
    while paths != []:
        path = paths.pop(0)
        last = path[-1]
        visited.append(last)
        if last == node2:
            return path
        for neighbour in neighbourhood(last, graph):
            if neighbour not in visited:
                paths.append(path + [neighbour])
    return []

assert shortest_path(1, 2, N2) == []
assert shortest_path(1, 1, LOOP) == [1]
assert shortest_path(3, 1, K3) == [3, 1]
assert shortest_path(1, 4, C5) == [1, 5, 4]

# **EXERCISES:**
#
# - Write a function to return *all* shortest paths between two nodes.


def distance(node1, node2, graph):
    """The **distance** is the length of the shortest path from node1 to node2.

    The length of a path is the number of its edges.
    If there's no path, return -1 to represent infinite distance.
    """
    return len(shortest_path(node1, node2, graph)) - 1

assert distance(1, 1, LOOP) == 0
assert distance(2, 4, C4) == 2
assert distance(3, 2, K3) == 1
assert distance(1, 5, K3) == -1


def diameter(graph):
    """A graph's **diameter** is the longest of all pairwise distances."""
    diameter = -1
    for node1 in nodes(graph):
        for node2 in nodes(graph):
            d = distance(node1, node2, graph)
            if d == -1:
                return d
            elif d > diameter:
                diameter = d
    return diameter

assert diameter(EMPTY) == -1
assert diameter(N2) == -1
assert diameter(K3) == 1
assert diameter(C5) == 2

# **EXERCISES:**
#
# - Rewrite `diameter()` to make it twice as fast.
#   Hint: the distance from A to B is the distance from B to A.


def connected(graph):
    """A network is **connected** if every node is reachable from any node."""
    for node1 in nodes(graph):
        for node2 in nodes(graph):
            if distance(node1, node2, graph) == -1:
                return False
    return True


# Import a random number generator with uniform distribution
from random import uniform


def random(n, p):
    """Return a simple **random graph** with n nodes and edge probability p."""
    assert n >= 0
    assert 0 <= p and p <= 1
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

# **EXERCISES:**
#
# - Write a program that imports SINPLE, creates several random graphs,
#   collects their mean degrees and diameters, and checks whether the
#   empirical values match the expected theoretical values.
#   The expected mean degree is `n/p`.
# - Write a function to create a simple random graph of given order and size.

# ### Exercises
#
# 1. What is the mean degree of the empty graph?
# 1. Write a program that imports SINPLE and reports various properties
#    of the Arpanet network, e.g. the most connected nodes.

# Projects
#
# Extend SINPLE to support **weighted edges**,
# i.e. edges that have an associated number, called the edge's weight.
# By default the edge weight is 1, i.e. an unweighted network can be seen
# as a weighted network where each edge has weight 1.
# Learn about default parameter values in Python and change `edge()`.
# Add a function to return the weight of an edge.
# Modify `shortest_path()` to return the path with the smallest total weight.
#
# Check that ARPANET is robust to any communication link failing, i.e.
# the network remains connected if a single edge is removed, i.e.
# each edge belongs to a cycle.
#
# Check that KARATE has two highly connected members, the student founder and
# the instructor, who are not friends of each other,
# and that most members are not friends of both.
# This conflict led to the split into two rival clubs.
