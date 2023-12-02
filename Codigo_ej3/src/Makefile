CC = gcc
CFLAGS = -O3 -Wall -std=c99 -g -fopenmp

N = 
N_GEN = 
T_POB = 
M_RATE =
NHI =
NHF =

EXEC = sec
OUTFILE = ../output/out.txt
C_FILES = main.c io.c ga.c

sec: $(C_FILES)
	$(CC) $(CFLAGS) $(C_FILES) -o $(EXEC) -DTIME -DDEBUG

test_sec:
	./$(EXEC) $(N) $(N_GEN) $(T_POB) $(M_RATE) $(NHI) $(NHF) >> $(OUTFILE)

clean:
	rm -f $(EXEC)
