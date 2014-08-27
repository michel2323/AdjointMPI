#ifndef AMPI_STACK_H
#define AMPI_STACK_H
#include <stdlib.h>

#define CHUNK_SIZE 1000


typedef struct {
    double *v;
    int top;
    size_t size;
} ampi_stack;

void push(ampi_stack *s, double val);
double pop(ampi_stack *s);
void stack_init(ampi_stack *s);
void destroy(ampi_stack *s);
int full(ampi_stack *s);
void expand(ampi_stack *s);
void shrink(ampi_stack *s);
int empty(ampi_stack *s);
#endif
