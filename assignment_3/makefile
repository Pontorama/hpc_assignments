# - the GNU make manual starts has its first section at
#   https://www.gnu.org/software/make/manual/make.html#Rule-Introduction
# - comments in makefiles are prefaced by #

# The form of a rule in make is
# TARGET : PREREQUISITS
# 	RECIPE
# The recipe must be indended by TABS as opposed to any other whitespace.

TEST_SRC = $(wildcard tests/*.c)
TESTS = $(patsubst %.c, %, $(TEST_SRC))
.PHONY : all
all : newton $(TESTS)

newton: src/assignment_3.c
	gcc -march=native -lm -g -O2 -o $@ $<
$(TESTS): %:%.c
	gcc -O3 -g -o $@ $<

.PHONY: clean
clean: 
	rm -rf newton
