test: func.o test.o debug.o
	g++ -o test func.o test.o debug.o
	rm -f *.o

func.o: func.cpp
	g++ -c func.cpp -DBLOCK_SIZE=$(if $(BLOCK_SIZE),$(BLOCK_SIZE),1)

test.o: test.cpp
	g++ -c test.cpp -DBLOCK_SIZE=$(if $(BLOCK_SIZE),$(BLOCK_SIZE),1)

debug.o: debug.cpp
	g++ -c debug.cpp

clean:
	rm -f test
	rm -f *.o