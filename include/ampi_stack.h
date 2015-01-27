#ifndef AMPI_STACK_H
#define AMPI_STACK_H
/** \file
 * \brief The AMPI stack is used for tracing the operations in the AMPI_Reduce()
 * and AMPI_Allreduce().
 */

#include <stdlib.h>

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
 * @param size Size of the stack. The stack should use the full size.
 */
ampi_stack* AMPI_stack_create(size_t size);

/**
 * @brief Destroy stack
 *
 * @param s stack to be destroyed
 */
void AMPI_stack_delete(ampi_stack *s);

/**
* @brief Reset the stack to the state after the pushing was finished
*
* @param s stack to be reseted
*/
void AMPI_stack_reset(ampi_stack *s);
#endif
