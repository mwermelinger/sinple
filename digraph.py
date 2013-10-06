"""
This is a pedagogical implementation of directed graphs in Python.

It is pedagogical at various levels:

- It illustrates basic Python constructs: loops, functions, tuples, sets, etc.
- It illustrates a simple form of unit testing.
- It explains basic graph theoretical concepts through executable definitions (the functions) and concrete examples (the unit tests).
- It is extensively documented. You can run [Pycco](http://fitzgen.github.io/pycco/) on this file to generate a side-by-side view of documentation and code. 
- It includes exercises.

This implementation does *not* aim to provide efficient operations for large graphs: there are several state-of-the-art libraries for that. 

The terminology was taken from the Wikipedia pages on [digraphs](http://en.wikipedia.org/wiki/Digraph_(mathematics)).

The latest version of this file can be found on https://bitbucket.org/mwermelinger/digraph-py.

Digraph representation
----------------------
A digraph (short for 'directed graph') is represented by a pair `(nodes, arcs)`, where `nodes` is the *set* of unique node names (e.g. integers or strings) and `arcs` is the *set* of arcs. Each arc is a pair of node names `(tail, head)`, stating there is a directed arc from the `tail` node to the `head` node.
"""

# Example graphs
# --------------
# The following digraphs will be used as unit tests for functions. 

# The **empty graph** has no nodes and hence no arcs.
EMPTY = (set(), set())		# {} is the empty dictionary, not the empty set

# A **null graph** Nn has n nodes but no edges.
N1 = ({1}, set())			# also called the trivial graph
N2 = ({1,2}, set())

# The smallest non-simple graph.
LOOP = ({1}, {(1,1)})

# A **cycle graph** Cn has n nodes connected in a 'round-robin' fashion. 
# Cycle graphs can be drawn as regular shapes: triangle, square, pentagon, hexagon, etc.
C3 = ({1,2,3}, {(1,2), (2,3), (3,1)})               		
C4 = ({1,2,3,4}, {(1,2), (2,3), (3,4), (4,1)})      		
C5 = ({1,2,3,4,5},{(1,2), (2,3), (3,4), (4,5), (5,1)})     

# In a **complete graph** Kn, each of the n nodes is connected to all other n-1 nodes.
K2 = ({1,2}, {(1,2), (2,1)})
K3 = ({1,2,3}, {(1,2), (1,3), (2,1), (2,3), (3,1), (3,2)})

# The **Petersen graph** is a (counter-)example of various graph theory concepts.
# It is usually drawn as a pentagram connected to a pentagon.
PETERSEN = ({1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
            {(1,2), (2,3), (3,4), (4,5), (5,1),     # the outer C5 (pentagon)
             (6,8), (8, 10), (10, 7), (7,9), (9,6), # the inner C5 (pentagram)
             (1,6), (2,7), (3,8), (4,9), (5,10)     # connected to one another
            })

# ### Exercises
#
# 1. Define K1 and C2.

# Basic functions
# ---------------
# These functions make the remaining functions independent of the actual data structure.

def arcs(graph):
	"""Return the set of the graph's arcs"""
	(nodes, arcs) = graph
	return arcs
	
def nodes(graph):
	"""Return the set of the graph's nodes."""
	(nodes, arcs) = graph
	return nodes
	
def tail(arc):
	"""The **tail** of an arc is the node it points *from*."""
	(tail, head) = arc
	return tail
	
def head(arc):
	"""The **head** of an arc is the node it points *to*."""
	(tail, head) = arc
	return head
	
def make_arc(tail, head):
	"""Create an arc from tail node to head node."""
	return (tail, head)
	
# Arc properties
# --------------

def inverted(arc):
	"""The **inverted arc** has its head and tail swapped."""
	return make_arc(head(arc), tail(arc))
	
def loop(arc):
	"""An arc is a **loop** if it connects a node to itself."""
	return tail(arc) == head(arc)
	
# Node properties
# ---------------

def adjacent(node1, node2, graph):
	"""Two nodes are **adjacent** if an arc connects them in either direction."""
	return 	make_arc(node1, node2) in arcs(graph) or \
			make_arc(node2, node1) in arcs(graph)
		
# Unit tests for `adjacent()`
assert(adjacent(1, 2, K2))
assert(adjacent(1, 4, C4))
assert(not adjacent(1, 3, C4))

def isolated (node, graph):
	"""A node is **isolated** if it is not adjacent to any other node."""
	for node2 in nodes(graph):
		if adjacent(node, node2, graph):
			return False
	return True
	
# Unit tests for `isolated()`
assert(isolated(1, N1))
assert(not isolated(1, K2))

def indegree(node, graph):
	"""The **indegree** of a node is the number of incoming arcs."""
	degree = 0
	for arc in arcs(graph):
		if head(arc) == node:
			degree = degree + 1
	return degree
	
def outdegree(node, graph):
	"""The **outdegree** of a node is the number of outgoing arcs."""
	degree = 0
	for arc in arcs(graph):
		if tail(arc) == node:
			degree = degree + 1
	return degree

def source(node, graph):
	"""A node is a **source** if it has no incoming arcs."""
	return indegree(node, graph) == 0
	
# Unit tests for `source()`
assert(source(1, N1))
assert(not source(3, C4))

def sink(node, graph):
	"""A node is a **sink** if it has no outgoing arcs."""
	return outdegree(node, graph) == 0

# Unit tests for `sink()`
assert(sink(1, N1))
assert(not sink(3, C4))

# ### Exercises
#
# 1. Implement `outdegree()` and `sink()` and un-comment its unit tests.
# 1. Rewrite `adjacent()` to make it more efficient, by searching the arc set only once.
# 1. Is every node adjacent to itself? Does that depend on the existence of a loop? Look for answers on the internet, modify `adjacent()` if needed, and create further unit tests. 
# 1. Write another version of `isolated()` using the degree concept. Which is more efficient?

# Graph properties
# ----------------

def size(graph):
	"""The graph's **size** is the number of its arcs."""
	return len(arcs(graph))

# Unit tests for `size()`
assert(size(EMPTY) == 0)
assert(size(N1) == 0)
assert(size(C3) == 3)
assert(size(PETERSEN) == 15)
	
def order(graph):
	"""The graph's **order** is the number of its nodes."""
	return len(nodes(graph))
	
# Unit tests for `order()`
assert(order(EMPTY) == 0)
assert(order(N1) == 1)
assert(order(C3) == 3)
assert(order(PETERSEN) == 10)

def simple(graph):
	"""A graph is **simple** if it has no loops and no duplicate arcs."""
	for arc in arcs(graph):
		if loop(arc):
			return False
	return True
	
# Unit tests for `simple()`
assert(simple(EMPTY))
assert(simple(N1))
assert(simple(PETERSEN))
assert(not simple(LOOP))
	
def oriented(graph):
	"""An **oriented** graph is an undirected graph with a direction imposed on each edge."""
	for arc in arcs(graph):
		if not loop(arc) and inverted(arc) in arcs(graph):
			return False
	return True		
	
assert(oriented(EMPTY))
assert(oriented(N1))
assert(oriented(LOOP))
assert(not oriented(K2))

def symmetric(graph):
	"""A digraph is **symmetric** if for every arc it also has the inverted arc."""
	for arc in arcs(graph):
		if inverted(arc) not in arcs(graph):
			return False
	return True
	
# Unit tests for `symmetric()`
assert(symmetric(EMPTY))
assert(symmetric(N1))
assert(symmetric(LOOP))
assert(symmetric(K2))
assert(not symmetric(C3))

def balanced(graph):
	"""A digraph is **balanced** if every node has the same in and out degrees."""
	for node in nodes(graph):
		if indegree(node, graph) != outdegree(node, graph):
			return False
	return True

# Unit tests for `balanced()`
assert(balanced(EMPTY))
assert(balanced(N1))
assert(balanced(K2))
assert(balanced(C3))
assert(not balanced(PETERSEN))

def tournament(graph):
	"""A tournament is a digraph with exactly one arc for each pair of distinct nodes."""
	for node1 in nodes(graph):
		for node2 in nodes(graph):
			if node1 != node2:
				has_arc = make_arc(node1, node2) in arcs(graph)
				has_inverted = inverted(arc) in arcs(graph)
				if (not has_arc and not has_inverted) or (has_arc and has_inverted):
					return False
	return True
	
def acyclic(graph):
	"""A graph is **acyclic** if it has no cycles."""
	return False
	
# ### Exercises
#
# 1. Why doesn't `simple()` check for duplicate arcs?
# 1. Rewrite `simple()` to make it more efficient. Hint: usually the order is smaller than the size. 
