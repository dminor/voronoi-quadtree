INCS = -I../../include 
LIBS = 
CFLAGS = -g
LDFLAGS = -L../../bin 
OBJS = main.o
TARGET = ../../bin/voronoi-pts

all: $(OBJS)
	g++ $(LDFLAGS) $(LIBS) $(OBJS) -o $(TARGET) 

.cpp.o:
	g++ $(INCS) $(CFLAGS) -c $< -o $@

main.o: ../../include/voronoi_quadtree.h

clean:
	rm *.o $(TARGET) 


