all: perspectiveTransform.o

perspectiveTransform.o:  perspectiveTransform.cpp perspectiveTransform.hpp
	g++ -c perspectiveTransform.cpp -o perspectiveTransform.o

clean:
	rm *.o

