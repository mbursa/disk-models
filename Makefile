default: all

all: models


.PHONY: models examples clean


# buid individual models
models:
	for model in `ls -d models/*/`; do \
	    echo "models/$${model%%/}/Makefile"; \
	    if [ -f "models/$${model%%/}/Makefile" ]; then \
	        echo "=== Building $${model%%/} model"; \
	        $(MAKE) -C models/$${model%%/}; \
	        echo "=== done"; \
	    fi; \
	done
#	@mv */*.so .


#diskmodel-test: diskmodel-test.c
#	@echo "=== Building diskmodel-test"
#	gcc $@.c -o $@ -ldl -rdynamic
#	@echo "=== done"



clean:
	@rm -rf diskmodel-test.o diskmodel-test
	@for model in `ls -d */`; do \
		$(MAKE) -s -C models/$${model%%/} clean; \
	done



