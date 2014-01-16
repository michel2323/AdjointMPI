#include "ampi_stack.h"

void push(ampi_stack *s, double val) {
    if(full(s))
	expand(s);
    /*printf("size: %d, top: %d, val: %f\n", s->size, s->top, val);*/
    s->v[s->top] = val; 
    (s->top)=(s->top)+1;    
}

double pop(ampi_stack *s) {
    /*if(empty(s))*/
    /*shrink(s);*/
    (s->top)=(s->top)-1;
    /*printf("size: %d, top: %d, val: %f\n", s->size, s->top, s->v[s->top]);*/
    return (s->v[s->top]);
}

void stack_init(ampi_stack *s) {
    s->top=0;
    s->v = (double*) malloc(sizeof(double)*CHUNK_SIZE);
    s->size=CHUNK_SIZE;
}

void destroy(ampi_stack *s) {
    s->top=0;
    free(s->v);
    s->size=CHUNK_SIZE;
}

int full(ampi_stack *s) {
    return (s->top >= s->size);
}

void expand(ampi_stack *s) {
    double *tmp;
    s->size=s->size+CHUNK_SIZE;
    tmp=(double*) realloc(s->v,s->size*sizeof(double));
    if(tmp != NULL) {
	s->v = tmp;
    }
}

void shrink(ampi_stack *s) {
    double *tmp;
    s->size=s->size-CHUNK_SIZE;
    tmp=(double*) realloc(s->v,s->size*sizeof(double));
    if(tmp != NULL) {
	s->v = tmp;
    }
}

int empty(ampi_stack *s) {
    return (s->top <= s->size-CHUNK_SIZE);
}
