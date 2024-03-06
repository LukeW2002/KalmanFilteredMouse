
all: build
build:
	mkdir -p build
# Compiler
CC := gcc

# Compiler flags
CFLAGS := -Wall -Wextra -Werror -Wpedantic -g 
#-fsanitize=address 

# Library
LIBS := -lm -lSDL2

# Source files
SRC := src/main.c src/linAlg.c

# Object files
OBJ := $(patsubst src/%.c, build/%.o, $(SRC))

# Executable name
EXEC := multivariate

# Main target
all: $(EXEC)

# Rule to build the executable
$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS) 

# Rule to build object files
build/%.o: src/%.c
	$(CC) $(CFLAGS) -c -o $@ $< 


# Clean target
clean:
	rm -rf $(EXEC) build  

.PHONY: all clean
