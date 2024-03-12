OBJECTS = antigen.o antibody.o entity.o bcell.o tcell.o macrophage.o immunecomplex.o

all: $(OBJECTS)

%.o: %.cpp %.h 
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -f *.o