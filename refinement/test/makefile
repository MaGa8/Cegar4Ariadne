ROOTDIR = ../..
MODULES = $(REFINEMENT_DIR) $(TEST_RUNNER_DIR) $(TREE_DIR) $(GRAPH_DIR)
# TARGET = runRefinementTreeTest

-include $(ROOTDIR)/common.mk

INCFLAGS += -I $(ARIADNE_DIR)/ -I $(ARIADNE_DIR)/source/
LIBFLAGS = -L $(ARIADNE_DIR)/build/ -l ariadne -l ariadne-kernel
OPTFLAGS = -std=c++17 -g -w

.PHONY: reftree
reftree: check depend
	$(makeDependencies)
	make build

.PHONY: cegar
cegar: check depend
	$(makeDependencies)
	make build
	$(call makeExecutable,runCegarTest)

.PHONY: test
test:
	@echo $(LOCAL_DEPS)
