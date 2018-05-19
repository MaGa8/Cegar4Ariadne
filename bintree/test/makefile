# variables required to set
# ROOTDIR relative path to cegar root directory
# MODULES directories in which subcomponents reside that are required to be built first
# INCFLAGS include paths, other than modules specified and local include directory
# LIBFLAGS library paths and concrete libraries
#
# variables to set optionally
# TARGET target executable to make, leave unset if no executable is to be made
# OPTFLAGS optional flags passed to compiler

# shorthands for modules in this project listed below
TEST_RUNNER_DIR = $(ROOTDIR)/testRunner
TREE_DIR = $(ROOTDIR)/bintree
GRAPH_DIR = $(ROOTDIR)/digraph
REFINEMENT_DIR = $(ROOTDIR)/refinement
ARIADNE_DIR = $(ROOTDIR)/../ariadne

# added these
TARGET = runLinkedFixedBranchTreeTest
ROOTDIR = ../..
MODULES = $(TREE_DIR) $(TEST_RUNNER_DIR)
OPTFLAGS = -std=c++17 -Wall -g

SRCDIR = src
INCDIR = include
DEPDIR = depend
OBJDIR = obj
BINDIR = bin

vpath %.cpp $(SRCDIR)/
vpath %.hpp $(INCDIR)/
vpath %.o   $(OBJDIR)/
vpath %.d   $(DEPDIR)/

CXX = g++
INCFLAGS = -I $(INCDIR)/ $(foreach mod,$(MODULES),-I $(mod)/$(INCDIR)/ )
CXXFLAGS = $(INCFLAGS) $(OPTFLAGS)
LINK_FLAGS = $(CXXFLAGS) $(LIBFLAGS)

LOCAL_SOURCES = $(notdir $(wildcard $(SRCDIR)/*.cpp) )
LOCAL_OBJECTS = $(addprefix $(OBJDIR)/,$(LOCAL_SOURCES:%.cpp=%.o))
LOCAL_DEPS = $(addprefix $(DEPDIR)/,$(LOCAL_SOURCES:%.cpp=%.d))

GLOBAL_MAINS = $(wildcard $(SRCDIR)/run*.cpp) $(shell find $(ROOTDIR)/ -name "run*.cpp")
GLOBAL_SOURCES = $(filter-out $(GLOBAL_MAINS),$(wildcard $(SRCDIR)/*.cpp) $(foreach mod,$(MODULES),$(wildcard $(mod)/$(SRCDIR)/*.cpp)))
GLOBAL_OBJECTS = $(subst $(SRCDIR),$(OBJDIR),$(GLOBAL_SOURCES:%.cpp=%.o))

MAKEFLAGS += --no-builtin-rules

define makeDependencies =
$(foreach mod,$(MODULES), echo $(PWD) && cd $(mod) && make depend && make build && cd $(PWD) &&) echo "made dependencies"
endef

define makeExecutable =
$(CXX) $(LINK_FLAGS) -o $(BINDIR)/$1 $(filter-out $(OBJDIR)/$1.o,$(GLOBAL_OBJECTS)) $(OBJDIR)/$1.o
endef

.PHONY: all
all: 	check depend dependencies build executable

$(DEPDIR)/%.d:	%.cpp
		$(CXX) $(CXXFLAGS) -M -MF $@ $^

$(OBJDIR)/%.o:	%.cpp $(DEPDIR)/%.d
		$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BINDIR)/$(TARGET): $(LOCAL_OBJECTS) $(GLOBAL_OBJECTS)
	$(call makeExecutable,$(notdir $@))

# makes dependency files
.PHONY: depend
depend: $(LOCAL_DEPS)

# recurses into module directories
.PHONY: dependencies
dependencies:
	$(makeDependencies)

# makes object files
.PHONY: build
build: 	check $(LOCAL_OBJECTS)

.PHONY: executable
executable: $(BINDIR)/$(TARGET)

# makes sure root is defined
.PHONY: check
check: 	
ifndef ROOTDIR
	$(error "relative path to project root not defined")
endif

.PHONY: clean
clean:	check
	find $(OBJDIR) $(BINDIR) $(DEPDIR) -mindepth 1 -exec rm \{\} \;
	$(foreach mod,$(MODULES),echo $PWD && cd $(mod) && make clean && cd $(PWD) &&) echo "cleaned"

include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(LOCAL_SOURCES))))
