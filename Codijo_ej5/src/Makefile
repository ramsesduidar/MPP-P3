CC = mpicc
CFLAGS = -O3 -Wall -std=c99 -g

N = 
N_GEN = 
T_POB = 
M_RATE =
NGM =
NEM =

EXEC = sec
OUTFILE = ../output/out.txt
C_FILES = main.c io.c ga.c crear.c

sec: $(C_FILES)
	$(CC) $(CFLAGS) $(C_FILES) -o $(EXEC) -DTIME -DDEBUG

test_sec:
	mpirun -np 1 --oversubscribe ./$(EXEC) 2000 40 36 0.3 1 2 >> $(OUTFILE)
	mpirun -np 2 --oversubscribe ./$(EXEC) 2000 40 36 0.3 1 2 >> $(OUTFILE)
	mpirun -np 3 --oversubscribe ./$(EXEC) 2000 40 36 0.3 1 2 >> $(OUTFILE)
	mpirun -np 4 --oversubscribe ./$(EXEC) 2000 40 36 0.3 1 2 >> $(OUTFILE)
	mpirun -np 5 --oversubscribe ./$(EXEC) 2000 40 36 0.3 1 2 >> $(OUTFILE)
	mpirun -np 6 --oversubscribe ./$(EXEC) 2000 40 36 0.3 1 2 >> $(OUTFILE)
	mpirun -np 8 --oversubscribe ./$(EXEC) 2000 40 36 0.3 1 2 >> $(OUTFILE)
	mpirun -np 10 --oversubscribe ./$(EXEC) 2000 40 36 0.3 1 2 >> $(OUTFILE)
	mpirun -np 12 --oversubscribe ./$(EXEC) 2000 40 36 0.3 1 2 >> $(OUTFILE)
	mpirun -np 14 --oversubscribe ./$(EXEC) 2000 40 36 0.3 1 2 >> $(OUTFILE)
	mpirun -np 16 --oversubscribe ./$(EXEC) 2000 40 36 0.3 1 2 >> $(OUTFILE)

clean:
	rm -f $(EXEC)
