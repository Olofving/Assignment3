CFLAGS=-Wall -O3
INCLUDES=-I/opt/X11/include
LDFLAGS=-L/opt/X11/lib -lX11 -lm
 
galsim: galsim.o graphics.o
	gcc -o galsim galsim.o graphics.o $(LDFLAGS)

galsim.o: galsim.c graphics.c graphics.h
	gcc $(CFLAGS) $(INCLUDES) -c galsim.c $(LDFLAGS)

graphics.o: graphics.c graphics.h
	gcc $(CFLAGS) $(INCLUDES) -c graphics.c $(LDFLAGS)

clean:
	rm -f galsim.o galsim graphics.o result.gal

