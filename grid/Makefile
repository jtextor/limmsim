OBJECTS = gridpoint.o grid.o

all: $(OBJECTS)

%.o: %.cpp %.h 
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -f $(OBJECTS)