ROOTDIR = ..
TESTDIRS = $(TREE_DIR) $(GRAPH_DIR)
MODULES = $(TEST_RUNNER_DIR) $(UTIL_DIR) $(REFINEMENT_DIR) $(TESTDIRS) $(addsuffix /test,$(TESTDIRS)) $(REFINEMENT_DIR)/test/refinementTree $(REFINEMENT_DIR)/test/cegar
TARGET = runAllTests

INCFLAGS ?= -isystem $(ARIADNE_DIR) -isystem $(ARIADNE_DIR)/source
LDFLAGS += -L $(ARIADNE_DIR)/build
LDLIBS += -l ariadne

-include $(ROOTDIR)/common.mk



