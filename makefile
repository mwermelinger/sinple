all:
	python3 sinple.py
	pep8 sinple.py
	pydoc3 sinple > sinple.txt
	pycco -d . sinple.py
