#include "ampi_stack.h"

void AMPI_push(ampi_stack *s, double val) {
    if(AMPI_full(s))
	AMPI_expand(s);
    /*printf("size: %d, top: %d, val: %f\n", s->size, s->top, val);*/
    s->v[s->top] = val; 
    (s->top)=(s->top)+1;    
}

double AMPI_pop(ampi_stack *s) {
    /*if(empty(s))*/
    /*shrink(s);*/
    (s->top)=(s->top)-1;
    /*printf("size: %d, top: %d, val: %f\n", s->size, s->top, s->v[s->top]);*/
    return (s->v[s->top]);
}

void AMPI_stack_init(ampi_stack *s) {
    s->top=0;
    s->v = (double*) malloc(sizeof(double)*CHUNK_SIZE);
    s->size=CHUNK_SIZE;
}

void AMPI_destroy(ampi_stack *s) {
    s->top=0;
    free(s->v);
    s->size=CHUNK_SIZE;
}

int AMPI_full(ampi_stack *s) {
    return (s->top >= s->size);
}

void AMPI_expand(ampi_stack *s) {
    double *tmp;
    s->size=s->size+CHUNK_SIZE;
    tmp=(double*) realloc(s->v,s->size*sizeof(double));
    if(tmp != NULL) {
	s->v = tmp;
    }
}

void AMPI_shrink(ampi_stack *s) {
    double *tmp;
    s->size=s->size-CHUNK_SIZE;
    tmp=(double*) realloc(s->v,s->size*sizeof(double));
    if(tmp != NULL) {
	s->v = tmp;
    }
}

int AMPI_empty(ampi_stack *s) {
    return (s->top <= s->size-CHUNK_SIZE);
}
