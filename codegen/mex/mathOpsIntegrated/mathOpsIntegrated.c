/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mathOpsIntegrated.c
 *
 * Code generation for function 'mathOpsIntegrated'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mathOpsIntegrated.h"
#include "adder.h"

/* Function Definitions */
real_T mathOpsIntegrated(const emlrtStack *sp, real_T in1, real_T in2, real_T
  len)
{
  (void)sp;

  /*  for code generation, preinitialize the output variable */
  /*  data type, size, and complexity  */
  /*  generate an include in the C code */
  /*  evaluate the C function */
  return adder(in1, in2, len);
}

/* End of code generation (mathOpsIntegrated.c) */
