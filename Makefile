
DIRS = examples 

all:
	mkdir -p bin
	for dir in $(DIRS); do cd $$dir; make; cd ..; done

clean:
	for dir in $(DIRS); do cd $$dir; make clean; cd ..; done







