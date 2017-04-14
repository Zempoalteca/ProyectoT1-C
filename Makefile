CC=gcc
TARGET=predrg
#CFLAGS= -O3 -fopenmp
CFLAGS=-O3 -fopenmp
LIBS=-lm
OBJECTS= Imgs_Gplot_Rmol.o pred_do_MC.o pmain_MC.o
HEADERS= Imgs_Gplot_Rmol.h pred_MC.h

all: $(TARGET)

$(TARGET): $(OBJECTS)
	${CC} $(CFLAGS) $(OBJECTS) -o $(TARGET) $(LIBS)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf $(TARGET) *.o 


