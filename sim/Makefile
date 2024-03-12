OBJECTS = sim.o diffusion.o reaction.o 

all: $(OBJECTS)

reaction.o: reaction.cpp phagocytosis.cpp costimulation.cpp proliferation.cpp response.cpp sim.h
	$(CC) $(CFLAGS) -c reaction.cpp

%.o: %.cpp sim.h
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -f *.o