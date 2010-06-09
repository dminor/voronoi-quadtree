INCS = -I../../include 
LIBS = -lvoronoiquadtree -lshp
CFLAGS = -std=gnu99
LDFLAGS = -L../../bin 
OBJS = render_tree.o
TARGET = ../../bin/render-tree

all: $(OBJS)
	gcc $(LDFLAGS) $(LIBS) $(OBJS) -o $(TARGET) 

.c.o:
	gcc $(INCS) $(CFLAGS) -c $< -o $@

clean:
	rm *.o $(TARGET) 


