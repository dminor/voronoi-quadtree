INCS = -I../../include 
LIBS = -lshp
CFLAGS = -g
LDFLAGS = -L../../bin 
OBJS = shapefile.o
TARGET = ../../bin/render-shapefile

all: $(OBJS)
	g++ $(LDFLAGS) $(LIBS) $(OBJS) -o $(TARGET) 

.cpp.o:
	g++ $(INCS) $(CFLAGS) -c $< -o $@

shapefile.o: ../../include/voronoi_quadtree.h

clean:
	rm *.o $(TARGET) 


