/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mathOpsIntegrated_terminate.c
 *
 * Code generation for function 'mathOpsIntegrated_terminate'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mathOpsIntegrated.h"
#include "mathOpsIntegrated_terminate.h"
#include "_coder_mathOpsIntegrated_mex.h"
#include "mathOpsIntegrated_data.h"

/* Function Definitions */
void mathOpsIntegrated_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void mathOpsIntegrated_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (mathOpsIntegrated_terminate.c) */
