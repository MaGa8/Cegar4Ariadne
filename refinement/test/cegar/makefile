ROOTDIR = ../../..
TARGET = runCegarTest
MODULES = $(REFINEMENT_DIR) $(TEST_RUNNER_DIR) $(TREE_DIR) $(GRAPH_DIR)

# only set this variables if not in recursion (-> duplication)
INCFLAGS ?= -isystem $(ARIADNE_DIR)/ -isystem $(ARIADNE_DIR)/source/
LDFLAGS ?= -L $(ARIADNE_DIR)/build/
LDLIBS ?= -l ariadne

include $(ROOTDIR)/common.mk

.PHONY: dummy
dummy:
	@echo "got right mkfile"
