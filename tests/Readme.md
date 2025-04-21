# Parsing SLICOT Example Data into Standard Matrix Format

This document explains how to read the numerical data presented in the **"Program Data"** and **"Program Results"** sections of SLICOT HTML documentation examples and format it into standard mathematical matrices (row-major format).

---

## Background: Fortran vs. Mathematical Notation

SLICOT routines are written in Fortran, which stores multi-dimensional arrays in **column-major order**. This means elements of a column are stored contiguously in memory. However, standard mathematical notation and many other programming languages (like C/C++) typically use **row-major order**, where elements of a row are contiguous.

While the underlying Fortran code uses column-major, the presentation of data in the HTML documentation examples often appears **row-wise for readability**.

---

## Key Observation in Provided Examples

Looking at the `AB01MD.html`, `AB01ND.html`, `AB01OD.html`, `AB04MD.html`, and `AB05MD.html` files, the numerical data for matrices in the "Program Data" and "Program Results" sections is generally listed **row by row**.

---

## How to Parse (Row-wise Presentation)

When the data is presented row-wise, parsing it into a standard mathematical (row-major) matrix is straightforward:

1. **Identify Dimensions**: Find the dimensions of the matrix (e.g., N, M, P). These are usually specified just before the numerical data or can be inferred from the context (e.g., A is NxN, B is NxM).

2. **Read Sequentially**: Read the numbers in the order they appear in the HTML section.

3. **Fill Row by Row**: Fill your matrix row by row, from left to right.

---

## Examples

### Example 1: AB01ND (`AB01ND.html`)

**Program Data:**
```
 AB01ND EXAMPLE PROGRAM DATA
   3     2     0.0     I  --> N=3, M=2
  -1.0   0.0   0.0      --> A, Row 1
  -2.0  -2.0  -2.0      --> A, Row 2
  -1.0   0.0  -3.0      --> A, Row 3
   1.0   0.0   0.0      --> B, Row 1 (Data seems transposed here, see note below)
   0.0   2.0   1.0      --> B, Row 2 (Data seems transposed here)
```

**Parsing A (3x3):**
```
A = [ -1.0   0.0   0.0 ]
    [ -2.0  -2.0  -2.0 ]
    [ -1.0   0.0  -3.0 ]
```

**Parsing B (3x2):**

> Note: The data listing `1.0 0.0 0.0` then `0.0 2.0 1.0` seems inconsistent for a 3x2 matrix B. Assuming it should represent the columns read row-wise, or perhaps it's listed transposed in the HTML data section. Based on standard notation and the Fortran code reading `((B(I,J), J=1,M), I=1,N)`, it should be read column-wise. If we read the numbers sequentially `1.0, 0.0, 0.0, 0.0, 2.0, 1.0` and assume it's column-major data for B(3,2), the row-major form is:

```
B = [ 1.0   0.0 ]
    [ 0.0   2.0 ]
    [ 0.0   1.0 ]
```

**Program Results:**
```
 The transformed state dynamics matrix ... is
  -3.0000   2.2361  --> Acont, Row 1
   0.0000  -1.0000  --> Acont, Row 2
 ...
 The transformed input/state matrix B ... is
   0.0000  -2.2361  --> Bcont, Row 1
   1.0000   0.0000  --> Bcont, Row 2
 ...
 The similarity transformation matrix Z is
   0.0000   1.0000   0.0000  --> Z, Row 1
  -0.8944   0.0000  -0.4472  --> Z, Row 2
  -0.4472   0.0000   0.8944  --> Z, Row 3
```

**Parsing Acont (2x2):**
```
Acont = [ -3.0000   2.2361 ]
        [  0.0000  -1.0000 ]
```

**Parsing Bcont (2x2):**
```
Bcont = [  0.0000  -2.2361 ]
        [  1.0000   0.0000 ]
```

**Parsing Z (3x3):**
```
Z = [  0.0000   1.0000   0.0000 ]
    [ -0.8944   0.0000  -0.4472 ]
    [ -0.4472   0.0000   0.8944 ]
```

---

### Example 2: AB05MD (`AB05MD.html`)

**Program Data:** Defines dimensions N1=3, M1=2, P1=2, N2=3, P2=2. Then lists matrices A1, B1, C1, D1, A2, B2, C2, D2. Each matrix's data is presented row by row.

```
   1.0   0.0  -1.0  --> A1, Row 1
   0.0  -1.0   1.0  --> A1, Row 2
   1.0   1.0   2.0  --> A1, Row 3
   1.0   1.0   0.0  --> B1, Row 1 (Data seems transposed here, see note below)
   2.0   0.0   1.0  --> B1, Row 2 (Data seems transposed here)
   3.0  -2.0   1.0  --> C1, Row 1
   0.0   1.0   0.0  --> C1, Row 2
   1.0   0.0         --> D1, Row 1
   0.0   1.0         --> D1, Row 2
  -3.0   0.0   0.0  --> A2, Row 1
   1.0   0.0   1.0  --> A2, Row 2
   0.0  -1.0   2.0  --> A2, Row 3
   0.0  -1.0   0.0  --> B2, Row 1 (Data seems transposed here)
   1.0   0.0   2.0  --> B2, Row 2 (Data seems transposed here)
   1.0   1.0   0.0  --> C2, Row 1
   1.0   1.0  -1.0  --> C2, Row 2
   1.0   1.0         --> D2, Row 1
   0.0   1.0         --> D2, Row 2
```

**Parsing A1 (3x3):** Read the first 3 lines directly.
```
A1 = [  1.0   0.0  -1.0 ]
     [  0.0  -1.0   1.0 ]
     [  1.0   1.0   2.0 ]
```

**Parsing B1 (3x2):** Similar note as AB01ND - the data layout `1.0 1.0 0.0`, `2.0 0.0 1.0` doesn't fit 3x2 directly. Assuming it represents column-major data read sequentially: `1.0, 2.0, 1.0, 1.0, 0.0, 0.0`. The row-major form is:
```
B1 = [ 1.0   1.0 ]
     [ 2.0   0.0 ]
     [ 1.0   0.0 ]
```

**Parsing C1 (2x3):** Read the next two lines directly.
```
C1 = [ 3.0  -2.0   1.0 ]
     [ 0.0   1.0   0.0 ]
```
...and so on for the other matrices.

**Program Results:** The output matrices A, B, C, D are clearly printed row by row. You read them sequentially to form the standard row-major matrices.

---

## Handling Column-wise Presentation (Less Common in these HTMLs)

If you encounter data explicitly stated or formatted as column-wise (e.g., listing all elements of column 1, then all elements of column 2), you would need to perform a transpose operation mentally or programmatically to get the standard row-major mathematical representation.

**Example:** If 3x2 matrix B was listed column-wise as `1.0, 0.0, 0.0, 0.0, 2.0, 1.0`, you would parse it as:

- Column 1 = `[1.0, 0.0, 0.0]'`
- Column 2 = `[0.0, 2.0, 1.0]'`

Resulting in the column-major matrix:
```
[ 1 0
  0 2
  0 1 ]
```
Transposing this gives the row-major form:
```
B = [ 1.0   0.0 ]
    [ 0.0   2.0 ]
    [ 0.0   1.0 ]
```

---

## Conclusion

For the specific SLICOT HTML examples provided, the matrix data in the "Program Data" and "Program Results" sections is generally presented in a row-wise format. This makes parsing straightforward: identify the dimensions and read the numbers sequentially to fill the matrix row by row, left to right. Be mindful of potential inconsistencies in how data is presented in the "Program Data" section and cross-reference with the Fortran READ statements in the example code if available and if precision is critical.