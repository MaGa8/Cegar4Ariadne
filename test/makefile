ROOTDIR = ..
TESTDIRS = $(TREE_DIR) $(GRAPH_DIR) $(REFINEMENT_DIR)
MODULES = $(TEST_RUNNER_DIR) $(TESTDIRS) $(addsuffix /test,$(TESTDIRS))
TARGET = runAllTests

-include $(ROOTDIR)/common.mk

OPTFLAGS += -g -w -std=c++17
INCFLAGS += -I $(ARIADNE_DIR) -I $(ARIADNE_DIR)/source
LIBFLAGS += -L $(ARIADNE_DIR)/build -l ariadne -l ariadne-kernel
