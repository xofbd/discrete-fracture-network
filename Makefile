SHELL := /bin/bash
POETRY_RUN := poetry run

.PHONY: all
all: clean install

poetry.lock: pyproject.toml
	poetry lock

.make.install.prod: poetry.lock
	poetry install --no-dev
	rm -f .make.install.prod
	touch $@

.make.install.dev: poetry.lock
	poetry install
	rm -f .make.install.dev
	touch $@

.PHONY: install
install:
	pip install .

.PHONY: tests
tests: test-lint test-unit

.PHONY: test-unit
test-unit: .make.install.dev
	$(POETRY_RUN) python3 -m unittest -v

.PHONY: test-lint
test-lint: .make.install.dev
	$(POETRY_RUN) flake8 dfn tests

.PHONY: clean
clean:
	rm -f .make.*
	find . | grep __pycache__ | xargs rm -rf
