same_file:
	gcc -O2 same_file.c -o same_file.out

different_files:
	gcc -O2 different_file.c different_file_mul.c -o different_files.out

inline:
	gcc -O2 inline.c -o inline.out

disk_operations:
	gcc -O2 disk_writing.c ../lib/benchmark_utils.c -o disk_operations.out 

all: same_file different_files inline
