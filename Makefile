
CC = g++ -Wall -fPIC -no-pie
all: analisi

#src/general.cpp src/Chameleon.cpp src/ConfigFile.cpp src/Analyzer.cpp

general.o: src/general.cpp
	$(CC) -c $^ -o $@ -I./include/


Chameleon.o: src/Chameleon.cpp
	$(CC) -c $^ -o $@ -I./include/


ConfigFile.o: src/ConfigFile.cpp
	$(CC) -c $^ -o $@ -I./include/


Analyzer.o: src/Analyzer.cpp
	$(CC) -c $^ -o $@ -I./include/ -lstdc++fs `root-config --cflags --glibs` -L${ROOTSYS}/lib -I${ROOTSYS}/include -lCore -lTreePlayer


analisi: analisi.C Analyzer.o general.o Chameleon.o ConfigFile.o
	$(CC) -I./include/ -I./ -I./src/ $^ -o $@  -lstdc++fs `root-config --cflags --glibs` -L${ROOTSYS}/lib -I${ROOTSYS}/include -lCore -lTreePlayer

.PHONE: clean

clean:
	rm general.o
	rm Chameleon.o
	rm ConfigFile.o
	rm Analyzer.o
	rm analisi
