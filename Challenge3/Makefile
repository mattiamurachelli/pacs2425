# Compiler 
CXX = mpic++

# Flags
CXXFLAGS  ?= -std=c++17 -O3 -fopenmp -Wall -Wno-conversion-null -Wno-deprecated-declarations -MMD -MP
CPPFLAGS  ?= -Iinclude

# Executables
EXEC_SEQ  = main_sequential
EXEC_PAR  = main_parallel

# Directories
SRC_DIR   = src
OBJ_DIR   = obj
OBJ_SEQ   = $(OBJ_DIR)/seq
OBJ_PAR   = $(OBJ_DIR)/par

# Source file classification
COMMON_SRCS = $(wildcard $(SRC_DIR)/*.cpp)
COMMON_SRCS := $(filter-out $(SRC_DIR)/main_%.cpp, $(COMMON_SRCS))

SRCS_SEQ = $(COMMON_SRCS) $(SRC_DIR)/main_sequential.cpp
SRCS_PAR = $(COMMON_SRCS) $(SRC_DIR)/main_parallel.cpp

OBJS_SEQ = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_SEQ)/%.o, $(SRCS_SEQ))
OBJS_PAR = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_PAR)/%.o, $(SRCS_PAR))

DEPS_SEQ = $(OBJS_SEQ:.o=.d)
DEPS_PAR = $(OBJS_PAR:.o=.d)

# Build both by default
all: $(EXEC_SEQ) $(EXEC_PAR)

# Build rules
$(EXEC_SEQ): $(OBJS_SEQ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(EXEC_PAR): $(OBJS_PAR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

# Compile sources for sequential
$(OBJ_SEQ)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_SEQ)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Compile sources for parallel
$(OBJ_PAR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_PAR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Create build directories
$(OBJ_SEQ) $(OBJ_PAR):
	mkdir -p $@

# Clean
clean:
	$(RM) -r $(OBJ_DIR) test/Results/

distclean: clean
	$(RM) $(EXEC_SEQ) $(EXEC_PAR)

# Include auto-generated dependency files
-include $(DEPS_SEQ) $(DEPS_PAR)
