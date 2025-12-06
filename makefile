CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra
TARGET = qc_simulator
TEST_DIR = tests
TEST_SOURCES = $(wildcard $(TEST_DIR)/*.cpp)
TEST_TARGETS = $(TEST_SOURCES:.cpp=)

all: $(TARGET)

$(TARGET): main.cpp
	$(CXX) $(CXXFLAGS) main.cpp -o $(TARGET)

# Build individual tests
$(TEST_DIR)/%: $(TEST_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

# Build all tests
tests: $(TEST_TARGETS)

# Run all tests
test: tests
	@echo "Running all tests..."
	@for test in $(TEST_TARGETS); do \
		echo "\n=== Running $$test ==="; \
		./$$test || exit 1; \
	done
	@echo "✓ All tests passed!"

run: all
	./$(TARGET)

clean:
	rm -f $(TARGET) $(TEST_TARGETS)

.PHONY: all tests run_tests run clean