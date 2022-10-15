test: func.o test.o
	g++ -o test func.o test.o
	rm -f *.o

func.o: func.cpp
	g++ -c func.cpp

test.o: test.cpp
	g++ -c test.cpp

clean:
	rm -f test
	rm -f *.o