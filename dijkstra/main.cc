#include <stdlib.h>

#include "dijkstra.hh"

int main(int argc, char *argv[]){
	int pSize = 500000;

	if(argc == 2){
		pSize = atoi(argv[1]);
	}

	Elf::dijkstra(pSize);

	return 0;
}
