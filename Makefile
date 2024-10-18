SHELL := /bin/bash

.PHONY: setup
setup:
	@\
	python3 manage.py setup; \
	source .venv/bin/activate; \
	python3 manage.py install;

.PHONY: run
run:
	@\
	source ./.venv/bin/activate; \
	nohup python3 cabgen_pipeline_main.py & \

.PHONY: supervisord_run
supervisord_run:
	@\
	source ./.venv/bin/activate; \
	python3 cabgen_pipeline_main.py \
