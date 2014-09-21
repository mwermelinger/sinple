# SINPLE

Many systems can be represented as networks, e.g. transport connections, gene regulation, social interactions, the World Wide Web. Network science is an increasingly important field to analyse and predict the behaviour of such systems. 

SINPLE is a SImple Network Python Library for Education. It doesn't support analysis of very large networks: there are various state-of-the-art libraries for that. Instead SINPLE's aims are pedagogical.

SINPLE is written to be understandable, not necessarily highly efficiency. The code follows closely the formal network definitions, passes the [PEP8 style checker](https://github.com/jcrocholl/pep8), and illustrates a simple form of [unit testing](http://en.wikipedia.org/wiki/Unit_testing), using assertions. The unit tests provide concrete examples and counter-examples of the network concepts. 
Besides assertions, the code illustrates the use of sets, comprehensions, default parameter values, and string formatting.

The code is written in Python 3 in a semi-literate programming style to generate a side-by-side view of [documentation and code](http://mwermelinger.bitbucket.org/sinple/api.html) with [Pycco](http://fitzgen.github.io/pycco/). All functions have 'docstrings' to generate the [API documentation](http://mwermelinger.bitbucket.org/sinple/api.html) with pydoc.

SINPLE includes over 40 exercises. Please don't make your solutions publicly available.

SINPLE can write networks in [GML](http://en.wikipedia.org/wiki/Graph_Modelling_Language) format. Reading GML files is on the to-do list.

SINPLE is released under the Apache license (see the LICENSE.txt file).

Use the buttons on the left side to download SINPLE and to report bugs and suggest enhancements on the issue tracker.
