SUBDIRS = bmpio entities grid sim lists settings threadlocals affinity 
LDFLAGS = -lpthread
CFLAGS = -O2 -Wall 
CC = g++

export LDFLAGS CFLAGS CC

all: 
	for i in $(SUBDIRS) ; do \
	(cd $$i; make ) ; \
	done && \
	$(CC) $(CFLAGS) -c main.cpp ; \
	$(CC) $(CFLAGS) $(LDFLAGS) -o limmsim main.o affinity/*.o threadlocals/*.o bmpio/*.o entities/*.o grid/*.o settings/*.o sim/*.o


clean:
	rm -rf output/dumps/ ;\
	mkdir -p output/dumps ;\
	mkdir -p output/dumps/ag ;\
	mkdir -p output/dumps/b ;\
	rm -f output/data/*\

distclean: clean
	for i in $(SUBDIRS); do \
	(cd $$i; make clean ); \
	done; \
	rm -f limmsim *.o\

