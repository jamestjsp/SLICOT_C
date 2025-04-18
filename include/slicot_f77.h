/**
 * @file slicot_f77.h
 * @brief Fortran 77 name mangling and utilities for SLICOT C wrappers.
 *
 * Defines the F77_FUNC macro for handling Fortran name mangling differences
 * across platforms (Windows/Intel vs others) and may include other
 * definitions needed by the SLICOT C wrapper library.
 */
#ifndef SLICOT_F77_H
#define SLICOT_F77_H

/*
 * Define macros for Fortran 77 name mangling.
 *
 * F77_FUNC(name, NAME)
 * name: Lowercase Fortran function name
 * NAME: Uppercase Fortran function name
 *
 * The macro expands depending on the platform and compiler:
 * - Windows + Intel Fortran (__INTEL_COMPILER defined): Expands to NAME (uppercase)
 * - Windows + Other Compilers (e.g., MinGW/gfortran, MSVC): Expands to name##_ (lowercase + underscore)
 * - Non-Windows (Linux, macOS, etc.): Expands to name##_ (lowercase + underscore)
 */

#ifdef _WIN32 /* Check if compiling on Windows */

  /* Check specifically for the Intel Fortran compiler on Windows */
  #if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER) /* Include check for newer Intel compilers */
    /* Windows and Intel Fortran: Use uppercase names */
    #define F77_FUNC(name, NAME) NAME
  #else
    /* Windows but not Intel Fortran (e.g., MinGW/gfortran, MSVC): Use lowercase + underscore */
    /* This assumes gfortran-style name mangling for non-Intel compilers on Windows. */
    #define F77_FUNC(name, NAME) name##_
  #endif /* defined(__INTEL_COMPILER) */

#else /* Not _WIN32 */

  /* Non-Windows systems (Linux, macOS, etc.): Use lowercase + underscore */
  /* This is common for gfortran, adjust if using a different compiler */
  #define F77_FUNC(name, NAME) name##_

#endif /* _WIN32 */

#endif /* SLICOT_F77_H */
