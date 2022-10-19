tool = $(if $(toolset),$(toolset),g++)

test: func.o test.o debug.o
	@$(tool) -w -o test func.o test.o debug.o -O2
	@rm -f *.o

func.o: func.cpp
	@$(tool) -w -c func.cpp -O2

test.o: test.cpp
	@$(tool) -w -c test.cpp -O2

debug.o: debug.cpp
	@$(tool) -w -c debug.cpp -O2

clean:
	@rm -f test
	@rm -f *.o