ROOTDIR = ../..
MODULES = $(REFINEMENT_DIR) $(TREE_DIR) $(GRAPH_DIR) $(UTIL_DIR)
TARGET = cegarDemo

-include $(ROOTDIR)/common.mk

INCFLAGS += -isystem $(ARIADNE_DIR)/ -isystem $(ARIADNE_DIR)/source/
LDFLAGS += -L $(ARIADNE_DIR)/build/
LDLIBS += -l ariadne
