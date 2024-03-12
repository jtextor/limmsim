OBJECTS = affinity.o 

all: $(OBJECTS)

%.o: %.cpp %.h 
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -f *.o