CC=g++
CFLAGS=-Wall -Wextra -Werror -std=c++17 -c
LIBS=

SRC_DIR=./
OBJS_DIR=./build/
SOURCES=$(wildcard $(SRC_DIR)*.cpp)
OBJS=$(patsubst $(SRC_DIR)%.cpp, $(OBJS_DIR)%.o, $(SOURCES))
GCOV_FLAGS=-fprofile-arcs -ftest-coverage
TARGET=s21_matrix_oop.a

TESTS_SRC_DIR=./tests/
TESTS_OBJS_DIR=./build/
TESTS_SOURCES=$(wildcard $(TESTS_SRC_DIR)*.cpp)
TESTS_OBJS=$(patsubst $(TESTS_SRC_DIR)%.cpp, $(TESTS_OBJS_DIR)%.o, $(TESTS_SOURCES))

TESTS_TARGET=$(TESTS_OBJS_DIR)s21_matrix_oop_tests.exe

DEV_LIBS=-lgtest -lpthread

all: $(TARGET)

$(TARGET): $(OBJS)
	@printf "[BUILD] \e[92mLinking objects\e[0m\n"
	@ar rcs $@ $(OBJS) $(LIBS)
	@printf "[BUILD] \e[95mBuilt target $@\e[0m\n"

$(OBJS_DIR)%.o: $(SRC_DIR)%.cpp
	@printf "[BUILD] \e[32mBuilding $< -> $@\e[0m\n"
	@$(CC) $(CFLAGS) $< -o $@

$(TESTS_TARGET): $(TARGET) $(TESTS_OBJS)
	@printf "[BUILD] \e[92mLinking objects\e[0m\n"
	@$(CC) -o $@ $(TESTS_OBJS) $(TARGET) $(DEV_LIBS)
	@printf "[BUILD] \e[95mBuilt target $@\e[0m\n"


$(TESTS_OBJS_DIR)%.o: $(TESTS_SRC_DIR)%.cpp
	@printf "[BUILD] \e[32mBuilding $< -> $@\e[0m\n"
	@$(CC) $(CFLAGS) $< -o $@


test: clean $(TESTS_TARGET)
	@./$(TESTS_TARGET)

clean: clean_target clean_tests
	@printf "[CLEAN] \e[33mCleaned\e[0m\n"

rebuild: clean all

format:
	@clang-format -i $(SOURCES) s21_matrix_oop.h ./tests/*.cpp ./tests/*.h

check-format:
	@clang-format -n $(SOURCES) s21_matrix_oop.h ./tests/*.cpp ./tests/*.h

valgrind: $(TESTS_TARGET)
	@valgrind --tool=memcheck --leak-check=yes ./$(TESTS_TARGET)

clean_target:
	@rm -rf $(OBJS) $(TARGET)

clean_tests:
	@rm -rf $(TESTS_OBJS) $(TESTS_TARGET)
	



