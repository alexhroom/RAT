/*
 * Non-Degree Granting Education License -- for use at non-degree
 * granting, nonprofit, education, and research organizations only. Not
 * for commercial or industrial use.
 *
 * strcmp.c
 *
 * Code generation for function 'strcmp'
 *
 */

/* Include files */
#include "strcmp.h"
#include "reflectivity_calculation_data.h"
#include "reflectivity_calculation_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo hb_emlrtRSI = {
    91,       /* lineNo */
    "strcmp", /* fcnName */
    "/usr/local/MATLAB/R2021b/toolbox/eml/eml/+coder/+internal/strcmp.m" /* pathName
                                                                          */
};

static emlrtRSInfo ib_emlrtRSI = {
    167,          /* lineNo */
    "loc_strcmp", /* fcnName */
    "/usr/local/MATLAB/R2021b/toolbox/eml/eml/+coder/+internal/strcmp.m" /* pathName
                                                                          */
};

static emlrtRSInfo jb_emlrtRSI = {
    240,       /* lineNo */
    "charcmp", /* fcnName */
    "/usr/local/MATLAB/R2021b/toolbox/eml/eml/+coder/+internal/strcmp.m" /* pathName
                                                                          */
};

/* Function Definitions */
boolean_T b_strcmp(const emxArray_char_T *a)
{
  static const char_T b_cv[10] = {'s', 't', 'a', 'n', 'd',
                                  'a', 'r', 'd', 'T', 'F'};
  int32_T exitg1;
  int32_T kstr;
  const char_T *a_data;
  boolean_T b_bool;
  a_data = a->data;
  b_bool = false;
  if (a->size[1] == 10) {
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr < 10) {
        if (a_data[kstr] != b_cv[kstr]) {
          exitg1 = 1;
        } else {
          kstr++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
  return b_bool;
}

boolean_T c_strcmp(const emlrtStack *sp, const emxArray_char_T *a)
{
  static const char_T b_cv[13] = {'a', 'i', 'r', '/', 's', 'u', 'b',
                                  's', 't', 'r', 'a', 't', 'e'};
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  int32_T exitg1;
  int32_T i;
  int32_T kstr;
  const char_T *a_data;
  boolean_T b_bool;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  a_data = a->data;
  st.site = &hb_emlrtRSI;
  b_st.site = &ib_emlrtRSI;
  b_bool = false;
  if (a->size[1] == 13) {
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr < 13) {
        c_st.site = &jb_emlrtRSI;
        i = (uint8_T)a_data[kstr];
        if (i > 127) {
          emlrtErrorWithMessageIdR2018a(
              &c_st, &s_emlrtRTEI, "Coder:toolbox:unsupportedString",
              "Coder:toolbox:unsupportedString", 2, 12, 127);
        }
        if (cv[i] != cv[(int32_T)b_cv[kstr]]) {
          exitg1 = 1;
        } else {
          kstr++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
  return b_bool;
}

/* End of code generation (strcmp.c) */