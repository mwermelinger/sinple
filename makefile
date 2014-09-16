all:
	python3 sinple.py
	pydoc3 -w sinple
	mv sinple.html sinple-api.html
	pycco -d . sinple.py
	cp -a pycco.css sinple.html sinple-api.html ../site/sinple/
	pep8 sinple.py
