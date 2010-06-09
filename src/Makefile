INCS = -I../include
FLAGS = -std=gnu99 -fpic
OBJS = voronoi_quadtree.o
BINDIR = ../bin
TARGET = libvoronoiquadtree.so

all: $(OBJS)
	gcc $(LIBS) $(FLAGS) -shared -Wl,-soname,$(TARGET) $(OBJS) -o $(BINDIR)/$(TARGET)

.c.o:
	gcc $(INCS) $(FLAGS) -c $< -o $@

clean:
	rm *.o $(TARGET) 

