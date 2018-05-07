Many systems can be represented as networks, e.g. transport connections, 
gene regulation, social interactions, the World Wide Web. 
Network science is an increasingly important field to analyse 
and predict the behaviour of such systems. 

SINPLE is a SImple Network Python Library for Education. 
It doesn't support analysis of very large networks: 
there are various state-of-the-art libraries for that.
Instead SINPLE's aims are pedagogical:

- The code is written to be understandable: 
it follows closely the formal network definitions, 
passes the [PEP8 style checker](https://github.com/jcrocholl/pep8),
and uses assertions as a simple form of unit testing. 
- The unit tests provide concrete examples and 
counter-examples of the network concepts.
- Besides assertions, the code illustrates the use of sets, comprehensions,
default parameter values, and string formatting in Python.
- SINPLE includes almost 50 exercises. 
Please don't make your solutions publicly available.
- SINPLE can read and write networks in
[GML](http://en.wikipedia.org/wiki/Graph_Modelling_Language) format,
with some restrictions.
- SINPLE includes real-world networks.

SINPLE aims to be usable for courses introducing Python, 
data structures and algorithms, or network/graph theory.

## Using SINPLE

SINPLE is provided as a single Python 3 file.
Use the button above to download it.

You may wish to first read the 
[API documentation](https://mwermelinger.github.io/sinple/api.html), 
to get a sense of the functionality provided. 
The documentation was generated from the functions' docstrings, using
[pydoc](https://docs.python.org/3/library/pydoc.html).

SINPLE is written in a 
semi-[literate programming](https://en.wikipedia.org/wiki/Literate_programming) 
style, with abundant comments. You can read through a 
[side-by-side view](https://mwermelinger.github.io/sinple/sinple.html)
of code and comments, generated with [Pycco](http://fitzgen.github.io/pycco/). 

SINPLE is released under the Apache license 
(see the downloaded LICENSE.txt file).

<!--Use the buttons on the left side to download SINPLE and 
to report bugs and suggest enhancements on the issue tracker.-->
