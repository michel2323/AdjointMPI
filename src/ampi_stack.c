#include "ampi_stack.h"

void AMPI_push(ampi_stack *s, double val) {
    s->v[s->top] = val;
    (s->top)=(s->top)+1;
}

double AMPI_pop(ampi_stack *s) {
    (s->top)=(s->top)-1;
    return (s->v[s->top]);
}

ampi_stack* AMPI_stack_create(size_t size) {
    ampi_stack* s = (ampi_stack*) malloc(sizeof(ampi_stack));
    s->top=0;
    s->v = (double*) malloc(sizeof(double)*size);
    s->size=size;

    return s;
}

void AMPI_stack_delete(ampi_stack *s) {
    free(s->v);
    free(s);
}

void AMPI_stack_reset(ampi_stack *s) {
    /* We assume that the stack is always used to full size */
    s->top = s->size;
}
