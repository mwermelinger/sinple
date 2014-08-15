"""
SImple Network Python Library for Education (SINPLE), version 0.1

Michel Wermelinger, 15 August 2014

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

# Network representation
# ----------------------
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
# Each edge is represented by a pair of node names `(node1, node2)`.
# The order of the two nodes does not matter.
#
# The functions in this and the next section encapsulate this representation
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
    return (node1, node2)


def endpoints(edge):
    """Return the **endpoints** of the edge, i.e. the nodes it connects."""
    (node1, node2) = edge
    return {node1, node2}


def incident(edge, node):
    """An edge is **incident** to the nodes it connects, and vice-versa."""
    return node in endpoints(edge)

assert(incident(edge(1, 2), 1))
assert(incident(edge(1, 2), 2))
assert(not incident(edge(1, 2), 3))


def is_loop(edge):
    """An edge is a **loop** if it connects a node to itself."""
    return len(endpoints(edge)) == 1

assert(is_loop(edge("A", "A")))
assert(not is_loop(edge("A", "B")))

# ### Network creation and example networks
#
# The following functions construct new graphs.
# The example graphs are used in further unit tests.


def null(n):
    """Return the **null graph** with nodes numbered 1 to *n* and no edges."""
    return ({x for x in range(1, n+1)}, set())

# The **empty graph** has no nodes and hence no edges.
EMPTY = null(0)

# The **trivial graph** has a single node and no edges.
N1 = null(1)
N2 = null(2)


def add_edges(edgeset, graph):
    """Add the edge set and their incident nodes to the graph."""
    nodeset = set()
    for edge in edgeset:
        nodeset = nodeset | endpoints(edge)
    return (nodes(graph) | nodeset, edges(graph) | edgeset)


def network(edges):
    """Create a network from the given edges."""
    return add_edges(edges, EMPTY)

# The smallest non-simple graph has a single loop.
LOOP = network({edge(1, 1)})

# A **complete graph** K*n* connects each node to all other *n-1* nodes.
K2 = network({edge(1, 2)})
K3 = network({edge(1, 2), edge(1, 3), edge(2, 3)})


def delete_edges(edgeset, graph):
    """Remove the edge set, but not their incident nodes, from the graph."""
    return (nodes(graph), edges(graph) - edgeset)

# A **cycle graph** C*n* has *n* nodes connected in a 'round-robin' fashion.
# Cycle graphs can be drawn as regular shapes: triangle, square, pentagon, etc.
C3 = network({edge(1, 2), edge(2, 3), edge(3, 1)})
C4 = add_edges({edge(3, 4), edge(4, 1)}, delete_edges({edge(3, 1)}, C3))
C5 = add_edges({edge(4, 5), edge(5, 1)}, delete_edges({edge(4, 1)}, C4))

# The **Petersen graph** is a (counter-)example of various graph concepts.
PETERSEN = add_edges({
    # It is usually drawn as an inner pentagram (C5)
    edge(6, 8), edge(8, 10), edge(10, 7), edge(7, 9), edge(9, 6),
    # connected to
    edge(1, 6), edge(2, 7), edge(3, 8), edge(4, 9), edge(5, 10)},
    # an outer pentagon (C5).
    C5)

# In December 1970, the Arpanet (the precursor of the Internet) had 13 nodes.
# Source: http://som.csudh.edu/cis/lpress/history/arpamaps.
ARPANET = network({
    edge("SRI", "Stanford"), edge("SRI", "UCSB"),
    edge("SRI", "UCLA"), edge("SRI", "Utah"),
    edge("Stanford", "UCLA"), edge("UCLA", "UCSB"), edge("UCLA", "RAND"),
    edge("RAND", "SDC"), edge("SDC", "Utah"), edge("UTAH", "MIT"),
    edge("RAND", "BBN"), edge("MIT", "BBN"), edge("BBN", "Harvard"),
    edge("Harvard", "Carnegie"), edge("Carnegie", "CASE"),
    edge("CASE", "Lincoln"), edge("Lincoln", "MIT")
    })

# ### Exercises
#
# 1. Consult the Python documentation about sets.
# Explain how `null()` works, and why `{}` wasn't used instead of `set()`.
# 1. Explain how `add_edges()` works. Why isn't `edgseset` just named `edges`?
# 1. Draw the Petersen graph on paper and label the nodes as in `PETERSEN`.
# 1. Explain how C4 is constructed from C3.
# 1. Define K1 and C2.

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

assert(degree(1, N1) == 0)
assert(degree(1, LOOP) == 2)
assert(degree(1, C4) == 2)
assert(degree(6, PETERSEN) == 3)


def adjacent(node1, node2, graph):
    """Two nodes are **adjacent** if an edge connects them."""
    return edge(node1, node2) in edges(graph) \
        or edge(node2, node1) in edges(graph)

assert(adjacent(1, 2, K2))
assert(adjacent(1, 4, C4))
assert(adjacent(1, 1, LOOP))
assert(not adjacent(1, 1, K2))
assert(not adjacent(1, 3, C4))


def isolated(node, graph):
    """A node is **isolated** if it is not adjacent to any other node."""
    for node2 in nodes(graph):
        if adjacent(node, node2, graph):
            return False
    return True

assert(isolated(1, N1))
assert(not isolated(1, K2))

# ### Exercises
#
# 1. Make `adjacent()` more efficient, searching the edge set once.
# 1. Is every node adjacent to itself?
#    Does that depend on the existence of a loop?
#    Look for answers on the internet, modify `adjacent()` if needed,
#    and create further unit tests.
# 1. Rewrite `isolated()` without using `adjacent()`.
#    Which version is more efficient?

# Graph properties
# ----------------


def size(graph):
    """The graph's **size** is the number of its edges."""
    return len(edges(graph))

assert(size(EMPTY) == 0)
assert(size(N1) == 0)
assert(size(C3) == 3)
assert(size(PETERSEN) == 15)


def order(graph):
    """The graph's **order** is the number of its nodes."""
    return len(nodes(graph))

assert(order(EMPTY) == 0)
assert(order(N1) == 1)
assert(order(C3) == 3)
assert(order(PETERSEN) == 10)


def simple(graph):
    """A graph is **simple** if it has no loops."""
    for edge in edges(graph):
        if is_loop(edge):
            return False
    return True

assert(simple(EMPTY))
assert(simple(N1))
assert(simple(PETERSEN))
assert(not simple(LOOP))
