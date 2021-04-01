//
//  Non-Degree Granting Education License -- for use at non-degree
//  granting, nonprofit, educational organizations only. Not for
//  government, commercial, or other organizational use.
//
//  power.cpp
//
//  Code generation for function 'power'
//


// Include files
#include "power.h"
#include "eml_int_forloop_overflow_check.h"
#include "mwmathutil.h"
#include "reflectivity_calculation.h"
#include "reflectivity_calculation_data.h"
#include "rt_nonfinite.h"
#include "sum.h"

// Variable Definitions
static emlrtRSInfo xb_emlrtRSI = { 70, // lineNo
  "power",                             // fcnName
  "/usr/local/MATLAB/R2020a/toolbox/eml/lib/matlab/ops/power.m"// pathName
};

static emlrtRSInfo yb_emlrtRSI = { 79, // lineNo
  "fltpower",                          // fcnName
  "/usr/local/MATLAB/R2020a/toolbox/eml/lib/matlab/ops/power.m"// pathName
};

static emlrtRSInfo ke_emlrtRSI = { 66, // lineNo
  "applyBinaryScalarFunction",         // fcnName
  "/usr/local/MATLAB/R2020a/toolbox/eml/eml/+coder/+internal/applyBinaryScalarFunction.m"// pathName 
};

static emlrtRSInfo le_emlrtRSI = { 188,// lineNo
  "flatIter",                          // fcnName
  "/usr/local/MATLAB/R2020a/toolbox/eml/eml/+coder/+internal/applyBinaryScalarFunction.m"// pathName 
};

static emlrtRTEInfo eh_emlrtRTEI = { 79,// lineNo
  5,                                   // colNo
  "power",                             // fName
  "/usr/local/MATLAB/R2020a/toolbox/eml/lib/matlab/ops/power.m"// pName
};

// Function Definitions
void power(const emlrtStack *sp, const coder::array<real_T, 1U> &a, coder::array<
           real_T, 1U> &y)
{
  int32_T nx;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &xb_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  b_st.site = &yb_emlrtRSI;
  y.set_size((&eh_emlrtRTEI), (&b_st), a.size(0));
  c_st.site = &ke_emlrtRSI;
  nx = a.size(0);
  d_st.site = &le_emlrtRSI;
  if ((1 <= a.size(0)) && (a.size(0) > 2147483646)) {
    e_st.site = &vb_emlrtRSI;
    check_forloop_overflow_error(&e_st);
  }

  for (int32_T k = 0; k < nx; k++) {
    y[k] = muDoubleScalarPower(a[k], 2.0);
  }
}

// End of code generation (power.cpp)