########## defaut ##########

CC      = g++ --std=c++14 -fPIC -no-pie
DEBUG   = -g -Wall
CFLAGS  = -c
CFLAGS_NO_C  = -Wall $(DEBUG)
LFLAGS  = -Wall

#-I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include

################################################################################
# ROOT Libs and links
################################################################################
GLIB         = `pkg-config --cflags --libs glib-2.0`
ROOT_LINKS   = `root-config --cflags --glibs`
#ROOFIT       = -lRooFit -lHtml -lMinuit -lRooFitCore -lRooStats
ROOT_LIBS    = -L${ROOTSYS}/lib -lTreePlayer -lCore
ROOT_INCLUDE = -I${ROOTSYS}/include
BOOST_LIB    = -L/usr/include/boost/property_tree

################################################################################
# top level folder path
################################################################################
BASE_PATH = /home/federico/Scrivania/BetaAnalysis/analisi 
MAIN_DIR = ./main
BIN_DIR = ./bin
################################################################################
# General Dir & Src
################################################################################
GENERAL_INCLUDE_DIR = ./include
GENERAL_SRC_DIR = ./src
GENERAL_BUILD_DIR = ./src
GENERAL_SRC_FILES = $(wildcard $(GENERAL_SRC_DIR)/*.cpp)
GENERAL_BUILD_FILES = $(patsubst $(GENERAL_SRC_DIR)/%.cpp, $(GENERAL_BUILD_DIR)/%.so, $(GENERAL_SRC_FILES))

################################################################################
# Chameleon Dir & Src
################################################################################
CHAMELEON_INCLUDE_DIR = ./include
CHAMELEON_SRC_DIR = ./src
CHAMELEON_BUILD_DIR = ./src
CHAMELEON_SRC_FILES = $(wildcard $(CHAMELEON_SRC_DIR)/*.cpp)
CHAMELEON_BUILD_FILES = $(patsubst $(CHAMELEON_SRC_DIR)/%.cpp, $(CHAMELEON_BUILD_DIR)/%.so, $(CHAMELEON_SRC_FILES))

################################################################################
# ConfigFile parser Dir & Src
################################################################################
CONFIG_INCLUDE_DIR = ./include
CONFIG_SRC_DIR = ./src
CONFIG_BUILD_DIR = ./src
CONFIG_SRC_FILES = $(wildcard $(CONFIG_SRC_DIR)/*.cpp)
CONFIG_BUILD_FILES = $(patsubst $(CONFIG_SRC_DIR)/%.cpp, $(CONFIG_BUILD_DIR)/%.so, $(CONFIG_SRC_FILES))

################################################################################
# Analyzer Dir & Src
################################################################################
ANALYZER_INCLUDE_DIR = ./src
ANALYZER_SRC_DIR = ./src
ANALYZER_BUILD_DIR = ./src
ANALYZER_SRC_FILES = $(wildcard $(ANALYZER_SRC_DIR)/*.cpp)
ANALYZER_BUILD_FILES = $(patsubst $(ANALYZER_SRC_DIR)/%.cpp, $(ANALYZER_BUILD_DIR)/%.so, $(ANALYZER_SRC_FILES))

################################################################################
# Main Dir & Src
################################################################################
#ANALISI_INCLUDE_DIR = ./include
#ANALISI_SRC_DIR = ./src
#ANALISI_BUILD_DIR = ./src
#ANALISI_SRC_FILES = $(wildcard $(ANALYZER_SRC_DIR)/*.cpp)
#ANALISI_BUILD_FILES = $(patsubst $(ANALYZER_SRC_DIR)/%.cpp, $(ANALYZER_BUILD_DIR)/%.so, $(ANALYZER_SRC_FILES))


#===============================================================================
# Make exe

all:$(BIN_DIR)/analisi.exe

SIMPLE_BETA_RUN_OBJ = $(GENERAL_BUILD_FILES)
SIMPLE_BETA_RUN_OBJ += $(CHAMELEON_BUILD_FILES)
SIMPLE_BETA_RUN_OBJ += $(CONFIG_BUILD_FILES)
SIMPLE_BETA_RUN_OBJ += $(ANALYZER_BUILD_FILES)
$(BIN_DIR)/analisi.exe: $(MAIN_DIR)/main_analisi.cpp $(SIMPLE_BETA_RUN_OBJ)
	$(call PRINT_BUILD,$@)
	$(CC) $(DEBUG) $(LFLAGS) $^ -o $@ $(ROOT_LINKS) $(ROOT_LIBS) $(ROOT_INCLUDE) $(GLIB) $(BOOST_LIB)
	@printf "\n"


#===============================================================================


#*******************************************************************************
# Build General object files.
#*******************************************************************************
$(GENERAL_BUILD_DIR)/%.so: $(GENERAL_SRC_DIR)/%.cpp
	$(call PRINT_BUILD,$@)
	$(CC) $(DEBUG) -c $^ -o $@ $(ROOT_LINKS) $(ROOT_LIBS) $(ROOT_INCLUDE)
	@printf "\n"

#*******************************************************************************
# Build Chaneleon object files.
#*******************************************************************************
$(CHAMELEON_BUILD_DIR)/%.so: $(CHAMELEON_SRC_DIR)/%.cpp
	$(call PRINT_BUILD,$@)
	$(CC) $(DEBUG) -c $^ -o $@ $(ROOT_LINKS) $(ROOT_LIBS) $(ROOT_INCLUDE)
	@printf "\n"

#*******************************************************************************
# Build Config object files.
#*******************************************************************************
$(CONFIG_BUILD_DIR)/%.so: $(CONFIG_SRC_DIR)/%.cpp
	$(call PRINT_BUILD,$@)
	$(CC) $(DEBUG) -c $^ -o $@ $(ROOT_LINKS) $(ROOT_LIBS) $(ROOT_INCLUDE)
	@printf "\n"

#*******************************************************************************
# Build Analyzer object files.
#*******************************************************************************
#$(ANALYZER_BUILD_DIR)/%.so: $(ANALYZER_SRC_DIR)/%.cpp
#	$(call PRINT_BUILD,$@)
#	$(CC) $(DEBUG) -c $^ -o $@ $(ROOT_LINKS) $(ROOT_LIBS) $(ROOT_INCLUDE)
#	@printf "\n"


$(BASE_PATH)/$(ANALYZER_BUILD_DIR)/Analyzer.cpp:
	$(call PRINT_BUILD,$(BASE_PATH)/$(ANALYZER_BUILD_DIR)/Analyzer.cpp)
	rootcint -f $(BASE_PATH)/$(ANALYZER_BUILD_DIR)/Analyzer.cpp -c $(BASE_PATH)/$(ANALYZER_INCLUDE_DIR)/Analyzer.hpp $(BASE_PATH)/$(ANALYZER_INCLUDE_DIR)/LinkDef.h

$(BASE_PATH)/$(ANALYZER_BUILD_DIR)/Analyzer.so: $(BASE_PATH)/$(ANALYZER_SRC_DIR)/Analyzer.cpp
	$(CC) $(DEBUG) $(CFLAGS) $< -o $@ $(ROOT_LINKS) $(ROOT_LIBS) $(ROOT_INCLUDE)
	@printf "\n"

$(BASE_PATH)/$(BETA_SCOPE_BUILD_DIR)/TimeWindow_Dict.o: $(BASE_PATH)/$(BETA_SCOPE_BUILD_DIR)/TimeWindow_Dict.cxx
	$(CC) $(DEBUG) $(CFLAGS) $< -o $@ $(ROOT_LINKS) $(ROOT_LIBS) $(ROOT_INCLUDE)
	@printf "\n"

$(BASE_PATH)/$(BETA_SCOPE_BUILD_DIR)/libTimeWindow_Dict.so: $(BASE_PATH)/$(BETA_SCOPE_BUILD_DIR)/TimeWindow_Dict.o $(BASE_PATH)/$(BETA_SCOPE_BUILD_DIR)/TimeWindow.o
	$(CC) -shared $^ -o $@ $(ROOT_LINKS) $(ROOT_LIBS) $(ROOT_INCLUDE)
	@printf "\033[0;33m $@ is created \e[0m \n"


analisi_cpp := $(SRC)$(BETA_SCOPE_ANALYSIS)$(CPP)main_analisi.cpp
$(BIN)analisi: $(analisi_cpp) $(GLIB_CONFIG_DEP)
	$(CC) $(DEBUG) $(LFLAGS) $^ -o $@ $(ROOT_LINKS) $(ROOT_LIBS) $(ROOT_INCLUDE) $(GLIB)


#.PHONY: clean
#clean:
#	@echo "cleaning .o files in build and binary files in bin."
#	@for file in ./*/build/*; do rm -rf $$file; echo Delete Object file: $$file; done
#	@for file in ./bin/*; do rm -rf $$file; echo Delete binary file: $$file; done
#	@echo "finished"
