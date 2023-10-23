
TEST_SRC = $(wildcard tests/*.c)
TESTS = $(patsubst %.c, %, $(TEST_SRC))
.PHONY : all
all : newton $(TESTS)

newton: a3.c
	gcc -march=native -lpthread -lm -O2 -o $@ $<
newton_debug: a3.c
	gcc -march=native -lpthread -lm -O0 -g -o $@ $<
$(TESTS): %:%.c
	gcc -O2 -g -lm -o $@ $<

.PHONY: clean
clean: 
	rm -rf newton
