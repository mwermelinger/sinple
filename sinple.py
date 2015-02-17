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


def edge(node1, node2):
    """Create an edge connecting the two nodes."""
    return frozenset({node1, node2})

assert edge(1, 2) == edge(2, 1)

# - Consult the Python documentation about sets.
#   Why is an edge represented by a frozen, i.e. immutable, set?
#   Hint: read again the definition of graph at the start of this section.


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
        nodes.update(endpoints(edge))
    return (node_set(graph) | nodes, edge_set(graph) | edges)

# - Explain how `add_edges()` works. Why is `set()` used instead of `{}`?


def network(edges, nodes=set()):
    """Create a graph with the given nodes, edges, and their endpoints."""
    return add_edges(edges, (nodes, set()))

assert network({edge(1, 2)}) == network({edge(2, 1)}, {1, 2})

# - Change `network()` to also accept lists of nodes and edges.


def delete_edges(edges, graph):
    """Create a new graph by removing the set of edges from graph."""
    return (node_set(graph), edge_set(graph) - edges)

assert delete_edges({edge(1, 2)}, network({edge(2, 1)})) \
    == network(set(), {1, 2})

# - Write a function to remove a set of nodes from a graph.


# ### Example networks
#
# These networks are used in further unit tests.

# The **empty graph** has no nodes and hence no edges.
EMPTY = network(set())
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


# Basics
# ------


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

    The induced subgraph has the edges of `graph` that connect the given nodes.
    """
    edges = {edge(n1, n2) for n1 in nodes for n2 in nodes}
    return network(edges & edge_set(graph), nodes)

assert subgraph({1}, N2) == N1
assert subgraph({1, 2, 3, 4, 5}, PETERSEN) == C5

# - Explain how `subgraph()` works.
# - Write a function to check if *g1* is a **subgraph** of *g2*,
#   i.e. if the nodes and edges of *g1* are included in those of *g2*.


# Import a random number generator with uniform distribution
from random import uniform


def random(n, p):
    """
    Return a simple **random graph** with n nodes and edge probability p.

    The two numbers can't be negative, and p can't be larger than 1.
    """
    edgeset = set()
    for node1 in range(1, n+1):
        for node2 in range(node1+1, n+1):
            if uniform(0, 1) < p:
                edgeset.add(edge(node1, node2))
    return add_edges(edgeset, null(n))

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


# Independence
# ------------
# These functions relate to non-adjacency.

# - Write a function to check if a node is **isolated**
#   (not connected to any node). Implement 4 versions,
#   using `incident()`, `adjacent()`, `degree()`, and `neighbours()`.
# - Which version is the most efficient?
# - Write a function to check if a graph is **bipartite**, i.e.
#   if its nodes can be separated into two parts so that each edge
#   connects nodes in different parts.
#   Include a unit test showing that the utility graph is bipartite.


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
assert sorted(geodesics(5, 1, network(
    {edge(1, 2), edge(2, 3), edge(2, 4), edge(3, 5), edge(4, 5)}))) \
    == [[5, 3, 2, 1], [5, 4, 2, 1]]

# - Explain how `geodesics` works.
# - Write a simplified `geodesic` (singular!) function that returns
# the first shortest path found, or `[]` if there is none.


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
