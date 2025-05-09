# Compiler and flags
CXX      ?= g++
CXXFLAGS ?= -std=c++20 -O3 -Wall -Wno-conversion-null -Wno-deprecated-declarations -MMD -MP
CPPFLAGS ?= -Iinclude

# Directories
SRC_DIR     = src
INC_DIR     = include
OBJ_DIR     = obj

# Executables
MAIN_BIN    = main
ROW_BIN     = test_rowmajor
COL_BIN     = test_colmajor
BUILDER_BIN = matrixBuilder

# Source files
MAIN_SRC    = $(SRC_DIR)/main.cpp
ROW_SRC     = $(SRC_DIR)/testing_rowmajor.cpp
COL_SRC     = $(SRC_DIR)/testing_colmajor.cpp
BUILDER_SRC = $(SRC_DIR)/RandomMatrixBuilder.cpp

# Object files
MAIN_OBJ    = $(OBJ_DIR)/main.o
ROW_OBJ     = $(OBJ_DIR)/testing_rowmajor.o
COL_OBJ     = $(OBJ_DIR)/testing_colmajor.o
BUILDER_OBJ = $(OBJ_DIR)/RandomMatrixBuilder.o

# Dependency files
DEPS        = $(wildcard $(OBJ_DIR)/*.d)

# Targets
.PHONY: all main testing builder clean distclean main_clean testing_clean builder_clean

all: $(MAIN_BIN) $(ROW_BIN) $(COL_BIN) $(BUILDER_BIN)

main_exec: $(MAIN_BIN)
testing_exec: $(ROW_BIN) $(COL_BIN)
builder_exec: $(BUILDER_BIN)

# Executable build rules
$(MAIN_BIN): $(MAIN_OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(ROW_BIN): $(ROW_OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(COL_BIN): $(COL_OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(BUILDER_BIN): $(BUILDER_OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

# Object build rules
$(OBJ_DIR)/main.o: $(MAIN_SRC) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJ_DIR)/testing_rowmajor.o: $(ROW_SRC) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJ_DIR)/testing_colmajor.o: $(COL_SRC) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJ_DIR)/RandomMatrixBuilder.o: $(BUILDER_SRC) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Ensure obj directory exists
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Clean only main-related files
main_clean:
	$(RM) $(MAIN_BIN) $(MAIN_OBJ) $(OBJ_DIR)/main.d

# Clean only testing-related files
testing_clean:
	$(RM) $(ROW_BIN) $(COL_BIN)
	$(RM) $(ROW_OBJ) $(OBJ_DIR)/testing_rowmajor.d
	$(RM) $(COL_OBJ) $(OBJ_DIR)/testing_colmajor.d

# Clean only builder-related files
builder_clean:
	$(RM) $(BUILDER_BIN) $(BUILDER_OBJ) $(OBJ_DIR)/RandomMatrixBuilder.d

# Shorthand for user-friendly syntax
run_main: $(MAIN_BIN)
run_testing: $(ROW_BIN) $(COL_BIN)
run_builder: $(BUILDER_BIN)

run_main_clean: main_clean
run_testing_clean: testing_clean
run_builder_clean: builder_clean

# Generic cleanup (leaves obj folder)
clean:
	$(RM) $(MAIN_OBJ) $(ROW_OBJ) $(COL_OBJ) $(BUILDER_OBJ)
	$(RM) $(OBJ_DIR)/*.d

# Full cleanup
distclean: clean
	$(RM) -r $(OBJ_DIR)
	$(RM) $(MAIN_BIN) $(ROW_BIN) $(COL_BIN) $(BUILDER_BIN)

# Include dependency files if they exist
-include $(DEPS)