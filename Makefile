all: main

main:
	f2py -c -m fnl neighborlist.f90

clean:
	rm -f fnl.*so
	rm -f *.o
	rm -f *.mod
