ROOTDIR = ..
TESTDIRS = $(TREE_DIR) $(GRAPH_DIR) $(REFINEMENT_DIR)
MODULES = $(TEST_RUNNER_DIR) $(TESTDIRS) $(addsuffix /test,$(TESTDIRS))
TARGET = runAllTests

INCFLAGS ?= -isystem $(ARIADNE_DIR) -isystem $(ARIADNE_DIR)/source
LDFLAGS += -L $(ARIADNE_DIR)/build
LDLIBS += -l ariadne

-include $(ROOTDIR)/common.mk



