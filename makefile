all:
	python3 sinple.py
	pydoc3 -w sinple
	sed 's|<font.*mw.*/font>||' sinple.html > api.html
	pycco -d . sinple.py
	cp -a pycco.css sinple.html api.html ../site/sinple/
	pep8 sinple.py
