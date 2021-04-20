//
// Non-Degree Granting Education License -- for use at non-degree
// granting, nonprofit, educational organizations only. Not for
// government, commercial, or other organizational use.
// File: reflectivity_calculation_terminate.h
//
// MATLAB Coder version            : 5.0
// C/C++ source code generated on  : 24-Feb-2021 09:21:20
//
#ifndef REFLECTIVITY_CALCULATION_TERMINATE_H
#define REFLECTIVITY_CALCULATION_TERMINATE_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "reflectivity_calculation_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
REFLECTIVITY_CALCULATION_DLL_EXPORT extern void
  reflectivity_calculation_terminate();

#endif

//
// File trailer for reflectivity_calculation_terminate.h
//
// [EOF]
//