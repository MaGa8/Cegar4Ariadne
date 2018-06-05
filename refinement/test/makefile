ROOTDIR = ../..
MODULES = $(REFINEMENT_DIR) $(TEST_RUNNER_DIR) $(TREE_DIR) $(GRAPH_DIR)

# only set this variables if not in recursion (-> duplication)
INCFLAGS ?= -isystem $(ARIADNE_DIR)/ -isystem $(ARIADNE_DIR)/source/
LDFLAGS ?= -L $(ARIADNE_DIR)/build/
LDLIBS ?= -l ariadne

include $(ROOTDIR)/common.mk

.PHONY: reftree_debug
reftree_debug: TARGET = runRefinementTreeTest
reftree_debug: debug

.PHONY: cegar_debug
cegar_debug: TARGET = runCegarTest
cegar_debug: debug

.PHONY: reftree_release
reftree_release: TARGET = runRefinementTreeTest
reftree_release: release

.PHONY: cegar_release
cegar_release: TARGET = runCegarTest
cegar_release: release

# .PHONY: reftree
# reftree: check depend
# 	$(makeDependencies)
# 	make build
# 	$(call makeExecutable,runRefinementTreeTest)

# .PHONY: cegar
# cegar: check depend
# 	$(makeDependencies)
# 	make build
# 	$(call makeExecutable,runCegarTest)

.PHONY: test
test:
	@echo $(LOCAL_DEPS)
