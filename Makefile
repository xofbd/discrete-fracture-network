.PHONY: build

build:
	pip install .

test:
	python -m unittest discover -s dfn/tests/
