all:
	python3 sinple.py
	pydoc3 -w sinple
	sed 's|<font.*mw.*/font>||' sinple.html > docs/api.html
	pycco -d . sinple.py
	mv sinple.html docs/
	pep8 sinple.py
