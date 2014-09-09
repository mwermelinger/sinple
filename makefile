all:
	python sinple.py
	pep8 sinple.py
	pydoc sinple > sinple.txt
	pycco -d . sinple.py
