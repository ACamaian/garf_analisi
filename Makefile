OBJS = vedaloss.o ecorr_veda.o theta_flow.o Classe_geo.o Classe_evento.o Classe_analisi.o perdenc.o faianalisi.o relkin.o
HEADER = Classe_geo.h Classe_evento.h Classe_analisi.h Classe_formule.h Analisi_Principale.h


faianalisi.out : $(OBJS) $(HEADER)
	g++  -O3 -Wall $(OBJS) -I. $(shell $(ROOTSYS)/bin/root-config --glibs --cflags ) -lgfortran -Wl,-Map=a.map -o faianalisi.out

.C.o:
	g++ -c -fno-aggressive-loop-optimizations -O3 -Wall $< -I$(ROOTSYS)/include -o $@

.c.o:
	gcc -c -fno-aggressive-loop-optimizations -O3 -Wall $< -o $@

%.o: %.cxx $(HEADER)
	g++ -c -fno-aggressive-loop-optimizations -O3 -Wall $< -I$(ROOTSYS)/include -o $@

.f.o:
	gfortran -O3 -fno-second-underscore -c $< -o $@

clean:
	rm *.out *.o
