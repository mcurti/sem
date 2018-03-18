/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * BarycentricWeights.c
 *
 * Code generation for function 'BarycentricWeights'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "BarycentricWeights.h"
#include "BarycentricWeights_emxutil.h"

/* Function Definitions */
void BarycentricWeights(const emxArray_real_T *x, emxArray_real_T *w)
{
  int u0;
  int u1;
  int b_w;

  /* BarycentricWeights Is returning the weights for Lagrange interpolation */
  /* algorithm 30 */
  /*    For a set of points "x_i" on x axis the functions returns the weights */
  /*    for Lagrange interpolation */
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    u1 = 0;
  } else {
    u0 = x->size[0];
    u1 = x->size[1];
    if (u0 > u1) {
      u1 = u0;
    }
  }

  b_w = w->size[0] * w->size[1];
  w->size[0] = 1;
  w->size[1] = u1;
  emxEnsureCapacity_real_T(w, b_w);
  for (b_w = 0; b_w < u1; b_w++) {
    w->data[b_w] = 1.0;
  }

  for (u0 = 1; u0 - 1 <= u1 - 2; u0++) {
    for (b_w = 0; b_w < u0; b_w++) {
      w->data[b_w] *= x->data[b_w] - x->data[u0];
      w->data[u0] *= x->data[u0] - x->data[b_w];
    }
  }

  b_w = w->size[0] * w->size[1];
  w->size[0] = 1;
  emxEnsureCapacity_real_T(w, b_w);
  u0 = w->size[0];
  b_w = w->size[1];
  u0 *= b_w;
  for (b_w = 0; b_w < u0; b_w++) {
    w->data[b_w] = 1.0 / w->data[b_w];
  }
}

/* End of code generation (BarycentricWeights.c) */
