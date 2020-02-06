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
#MAIN_DIR = ./main
#BIN_DIR = ./bin
################################################################################
# General Dir & Src
################################################################################
GENERAL_INCLUDE_DIR = ./include
GENERAL_SRC_DIR = ./src
GENERAL_BUILD_DIR = ./build
GENERAL_SRC_FILES = $(wildcard $(WAVEFORM_ANALYSIS_SRC_DIR)/*.cpp)
GENERAL_BUILD_FILES = $(patsubst $(WAVEFORM_ANALYSIS_SRC_DIR)/%.cpp, $(WAVEFORM_ANALYSIS_BUILD_DIR)/%.o, $(WAVEFORM_ANALYSIS_SRC_FILES))

################################################################################
# Chameleon Dir & Src
################################################################################
CHAMELEON_INCLUDE_DIR = ./include
CHAMELEON_SRC_DIR = ./src
CHAMELEON_BUILD_DIR = BetaScope/build
CHAMELEON_SRC_FILES = $(wildcard $(BETA_SCOPE_SRC_DIR)/*.cpp)
CHAMELEON_BUILD_FILES = $(patsubst $(BETA_SCOPE_SRC_DIR)/%.cpp, $(BETA_SCOPE_BUILD_DIR)/%.o, $(BETA_SCOPE_SRC_FILES))

################################################################################
# ConfigFile parser Dir & Src
################################################################################
CONFIG_INCLUDE_DIR = ./include
CONFIG_SRC_DIR = ./src
CONFIG_BUILD_DIR = ./build
CONFIG_SRC_FILES = $(wildcard $(BETA_CONFIG_PARSER_SRC_DIR)/*.cpp)
CONFIG_BUILD_FILES = $(patsubst $(BETA_CONFIG_PARSER_SRC_DIR)/%.cpp, $(BETA_CONFIG_PARSER_BUILD_DIR)/%.o, $(BETA_CONFIG_PARSER_SRC_FILES))

################################################################################
# Ana Dir & Src
################################################################################
ANALYZER_INCLUDE_DIR = ./include
ANALYZER_SRC_DIR = ./src
ANALYZER_BUILD_DIR = ./build
ANALYZER_SRC_FILES = $(wildcard $(BETA_CONFIG_PARSER_SRC_DIR)/*.cpp)
ANALYZER_BUILD_FILES = $(patsubst $(BETA_CONFIG_PARSER_SRC_DIR)/%.cpp, $(BETA_CONFIG_PARSER_BUILD_DIR)/%.o, $(BETA_CONFIG_PARSER_SRC_FILES))