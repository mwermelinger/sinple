# SINPLE

Many systems can be represented as networks, e.g. transport connections, gene regulation, social interactions, the World Wide Web. [Network science](http://en.wikipedia.org/wiki/Network_science) is an increasingly important field to analyse and predict the behaviour of such systems. 

SINPLE is a SImple Network Python Library for Education. SINPLE does *not* enable analysis of large networks; there are various state-of-the-art libraries for that. Instead SINPLE's aims are pedagogical:

- it illustrates basic Python constructs: loops, functions, tuples, sets, etc;
- it illustrates simple forms of [unit testing](http://en.wikipedia.org/wiki/Unit_testing) and [preconditions](http://en.wikipedia.org/wiki/Precondition) using assertions;
- it explains basic network concepts through executable definitions (the functions) and concrete examples (the unit tests);
- it includes small exercises and larger projects;
- it is extensively documented in a simple [literate programming](http://en.wikipedia.org/wiki/Literate_programming) style. 

Running [Pycco](http://fitzgen.github.io/pycco/) on the source code (`pycco -d . sinple.py`) generates `sinple.html` with a side-by-side view of documentation and code.