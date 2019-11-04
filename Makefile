.PHONY: install

install:
	pip install .

test:
	python -m unittest discover -s dfn/tests/
