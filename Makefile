SHELL := /bin/bash
POETRY_RUN := poetry run

.PHONY: all
all: clean install

poetry.lock: pyproject.toml
	poetry lock

.make.install: poetry.lock
	poetry install
	touch $@

.PHONY: install
install: .make.install

.PHONY: tests
tests: test-lint test-unit

.PHONY: test-unit
test-unit: .make.install
	$(POETRY_RUN) python3 -m unittest -v

.PHONY: test-lint
test-lint: .make.install
	$(POETRY_RUN) flake8 dfn tests

.PHONY: clean
clean:
	rm -f .make.install
	find . | grep __pycache__ | xargs rm -rf
