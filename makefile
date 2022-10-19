test: func.o test.o debug.o
	g++ -o test func.o test.o debug.o -O2
	rm -f *.o

func.o: func.cpp
	g++ -c func.cpp -O2

test.o: test.cpp
	g++ -c test.cpp -O2

debug.o: debug.cpp
	g++ -c debug.cpp -O2

clean:
	rm -f test
	rm -f *.o