#ifndef _GA
#define _GA
	
	typedef struct {
		int array_int[1000];
		double fitness;
	} Individuo;
	
	void cruzar(Individuo *, Individuo *, Individuo *, Individuo *, int, int);
	void fitness_reduction(const double *, Individuo *, int, int);
	void mutar(Individuo *, int, double);
	void mergeSort(Individuo **poblacion, int izq, int der, int n);
#endif
