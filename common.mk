# variables required to set
# ROOTDIR relative path to cegar root directory
# MODULES directories in which subcomponents reside that are required to be built first
# INCFLAGS include paths, other than modules specified and local include directory
# LDFLAGS library paths
# LDLIBS libraries
#
# variables to set optionally
# TARGET target executable to make, leave unset if no executable is to be made
# OPTFLAGS optional flags passed to compiler

# shorthands for modules in this project listed below
TEST_RUNNER_DIR = $(ROOTDIR)/testRunner
TREE_DIR = $(ROOTDIR)/bintree
GRAPH_DIR = $(ROOTDIR)/digraph
REFINEMENT_DIR = $(ROOTDIR)/refinement
UTIL_DIR = $(ROOTDIR)/util
ARIADNE_DIR = $(ROOTDIR)/../ariadne

SRCDIR = src
INCDIR = include
DEPDIR = depend
OBJDIR = obj
BINDIR = bin
RELDIR = release
DBGDIR = debug
export BUILDDIR

export TARGET

vpath %.cpp $(SRCDIR)/
vpath %.hpp $(INCDIR)/

export CXX = g++
export OPTFLAGS = -std=c++17 -Wall
# includes should be locally defined
RELFLAGS += -O
DBGFLAGS += -g
INCFLAGS += -I $(INCDIR)/ $(foreach mod,$(MODULES),-I $(mod)/$(INCDIR)/ )
CXXFLAGS += $(INCFLAGS) $(OPTFLAGS) -g -fopenmp
# libraries should be locally defined

LOCAL_SOURCES = $(notdir $(wildcard $(SRCDIR)/*.cpp) )
LOCAL_OBJECTS = $(addprefix $(BUILDDIR)/$(OBJDIR)/,$(LOCAL_SOURCES:%.cpp=%.o))
LOCAL_DEPS = $(addprefix $(DEPDIR)/,$(LOCAL_SOURCES:%.cpp=%.d))

GLOBAL_MAINS = $(wildcard $(SRCDIR)/run*.cpp) $(shell find $(ROOTDIR)/ -name "run*.cpp")
GLOBAL_SOURCES = $(filter-out $(GLOBAL_MAINS),$(wildcard $(SRCDIR)/*.cpp) $(foreach mod,$(MODULES),$(wildcard $(mod)/$(SRCDIR)/*.cpp)))
GLOBAL_OBJECTS = $(subst $(SRCDIR),$(BUILDDIR)/$(OBJDIR),$(GLOBAL_SOURCES:%.cpp=%.o))

define makeDependencies =
$(foreach mod,$(MODULES), echo $(PWD) && cd $(mod) && make depend && make build && cd $(PWD) &&) echo "made dependencies"
endef

# parameters: target file name (without directory)
define makeExecutable =
$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) -o $(BUILDDIR)/$(BINDIR)/$1 $(filter-out $(BUILDDIR)/$(OBJDIR)/$1.o,$(GLOBAL_OBJECTS)) $(BUILDDIR)/$(OBJDIR)/$1.o
endef

define localObjects =
$(addprefix $(BUILDDIR)/$(OBJDIR)/,$(LOCAL_SOURCES:%.cpp=%.o))
endef

.PHONY: debug
debug: OPTFLAGS += $(DBGFLAGS)
debug: BUILDDIR := $(DBGDIR)
debug: 	setup depend dependencies
	$(MAKE) executable
# need recursion: exported builddir, so makefile can be reread and prerequisites adapted

.PHONY: release
release: OPTFLAGS += $(RELFLAGS)
release: BUILDDIR = $(RELDIR)
release: setup depend dependencies
	$(MAKE) executable

$(DEPDIR)/%.d:	%.cpp
		$(CXX) $(CXXFLAGS) -MM -MT '$$(BUILDDIR)'/$(OBJDIR)/$(basename $(notdir $<)).o -MF $@ $^
#		sed -i 's/\($$(BUILDDIR)\/$(OBJDIR)\/[a-Z]\+.o\)/$(subst .,\.,$(subst /,\/,$@)) \1/' $@
		echo -e '\t' '$$(COMPILE.cpp)' -o '$$(BUILDDIR)'/'$$(OBJDIR)'/$(basename $(notdir $@)).o >> $@ $<

$(BUILDDIR)/$(BINDIR)/$(TARGET): $(LOCAL_OBJECTS) $(GLOBAL_OBJECTS)
	@echo $(CXXFLAGS)
	@echo $(OPTFLAGS)
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
build: 	setup $(LOCAL_OBJECTS)
	@echo "build dir is " $(BUILDDIR)

.PHONY: executable
executable: $(BUILDDIR)/$(BINDIR)/$(TARGET)

# makes sure root is defined
.PHONY: check
check: 	
ifndef ROOTDIR
	$(error "relative path to project root not defined")
endif

.PHONY: setup
setup: check
	if [ ! -d $(BUILDDIR) ]; then mkdir $(BUILDDIR); fi
	if [ ! -d $(BUILDDIR)/$(OBJDIR) ]; then mkdir $(BUILDDIR)/$(OBJDIR); fi
	if [ ! -d $(DEPDIR) ]; then mkdir $(DEPDIR); fi
	if [ ! -d $(BUILDDIR)/$(BINDIR) ]; then mkdir $(BUILDDIR)/$(BINDIR); fi

.PHONY: clean
clean:
	find $(addprefix $(DBGDIR)/,$(OBJDIR) $(BINDIR)) $(addprefix $(RELDIR)/,$(OBJDIR) $(BINDIR)) $(DEPDIR) -mindepth 1 -exec rm \{\} \;
	$(foreach mod,$(MODULES),echo $(PWD) && cd $(mod) && make clean && cd $(PWD) &&) echo "cleaned"

-include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(LOCAL_SOURCES))))

.PHONY: test
test:
	@echo "local sources " $(LOCAL_SOURCES)
	@echo "global sources " $(GLOBAL_SOURCES)
	@echo "cxx flags " $(CXXFLAGS)
	@echo "lib flags " $(LDFLAGS) $(LDLIBS)
