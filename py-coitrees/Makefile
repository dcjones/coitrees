SHELL := /bin/bash

ts := $(shell date -u +"%Y-%m-%dT%H:%M:%SZ")

.PHONY: help
help: ## This help message
	@echo -e "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s :)"

.PHONY: build
build: dev-packages ## Builds Rust code and Python modules
	maturin build

.PHONY: build-release
build-release: dev-packages ## Build module in release mode
	maturin build --release

.PHONY: nightly
nightly: ## Set rust compiler to nightly version
	rustup override set nightly

.PHONY: install
install: dev-packages ## Install module into current virtualenv
	maturin develop --release

.PHONY: publish
publish: ## Publish crate on Pypi
	poetry run maturin publish

.PHONY: clean
clean: ## Clean up build artifacts
	cargo clean
	find .  \( -type f -name "*.py[co]" -o -type d -name "__pycache__" \) -delete && echo "Removed pycs and __pycache__"
	rm -f data/*log

.PHONY: dev-packages
dev-packages: ## Install Python development packages for project
	poetry install

.PHONY: cargo-test
cargo-test: ## Run cargo tests only
	cargo test

.PHONY: test
test: cargo-test dev-packages install quicktest ## Install rscannls module and run tests

.PHONY: quicktest
quicktest: ## Run tests on already installed hyperjson module
	poetry run pytest tests

.PHONY: bench
bench: ## Run benchmarks
	poetry run pytest benchmarks

.PHONY: bench-compare
bench-compare: nightly dev-packages install ## Run benchmarks and compare results with other JSON encoders
	poetry run pytest benchmarks --compare

.PHONY: build-profile
build-profile: ## Builds binary for profiling
	cd profiling && poetry run cargo build --release

# Setup instructions here:
# https://gist.github.com/dlaehnemann/df31787c41bd50c0fe223df07cf6eb89
.PHONY: profile
profile: OUTPUT_PATH = measurements/flame-$(ts).svg
profile: FLAGS=booleans --iterations 10000
profile: nightly build-profile ## Run perf-based profiling (only works on Linux!)
	perf record --call-graph dwarf,16384 -e cpu-clock -F 997 target/release/profiling $(FLAGS)
	time perf script | stackcollapse-perf.pl | c++filt | flamegraph.pl > $(OUTPUT_PATH)
	@echo "$(OUTPUT_PATH)"
