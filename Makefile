SHELL := /bin/bash
POETRY_RUN := poetry run

.PHONY: all
all: clean install

poetry.lock: pyproject.toml
	poetry lock

.make.install.prod: poetry.lock
	poetry install --no-dev
	rm -f .make.install.dev
	touch $@

.make.install.dev: poetry.lock
	poetry install
	rm -f .make.install.prod
	touch $@

.PHONY: install
install:
	pip install .

.PHONY: tests
tests: test-lint test-unit

.PHONY: test-unit
test-unit: .make.install.dev
	$(POETRY_RUN) coverage run -m unittest -v
	$(POETRY_RUN) coverage report
	$(POETRY_RUN) coverage xml

.PHONY: test-lint
test-lint: .make.install.dev
	$(POETRY_RUN) flake8 dfn tests

.PHONY: tox
tox: .make.install.dev
	$(POETRY_RUN) tox

.PHONY: clean
clean:
	rm -f .coverage .make.* coverage.xml
	rm -rf .tox
	find . | grep __pycache__ | xargs rm -rf
