#ifndef AMPI_STACK_H
#define AMPI_STACK_H
#include <stdlib.h>

#define CHUNK_SIZE 1000

/**
 * AMPI stack used to save values in particular for the reductions
 */
typedef struct {
    double *v;  /**< value */
    long int top; /**< top of the stack */
    long int size; /**< size of the stack */
} ampi_stack; /**< stack entry*/

/**
 * @brief Push value
 *
 * @param s stack
 * @param val value
 */
void AMPI_push(ampi_stack *s, double val);

/**
 * @brief Pop value 
 *
 * @param s stack
 *
 * @return value
 */
double AMPI_pop(ampi_stack *s);

/**
 * @brief Create stack
 *
 * @param s created stack
 */
void AMPI_stack_init(ampi_stack *s);

/**
 * @brief Destroy stack
 *
 * @param s stack to be destroyed
 */
void AMPI_destroy(ampi_stack *s);

/**
 * @brief Check whether stack is full. Only used by internally. 
 *
 * @param s stack to be checked
 *
 * @return non zero if stack is full
 */
int AMPI_full(ampi_stack *s);

/**
 * @brief Expand stack. Only used internally. 
 *
 * @param s stack to expanded
 */
void AMPI_expand(ampi_stack *s);

/**
 * @brief Shrink stack. Only used internally.
 *
 * @param s stack to be shrunk
 */
void AMPI_shrink(ampi_stack *s);

/**
 * @brief Check whether stack is empty. Only used internally.
 *
 * @param s stack to be checked
 *
 * @return non zero if stack is empty
 */
int AMPI_empty(ampi_stack *s);
#endif
