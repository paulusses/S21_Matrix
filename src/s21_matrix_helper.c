#include "s21_matrix.h"

int ExistMatrix(matrix_t* A) {
  int fail = OK;
  if (A->matrix == NULL || A == NULL || A->rows < 1 || A->columns < 1 ||
      !A->matrix) {
    fail = INCORRECT_MATRIX;
  }
  return fail;
}

int CheckTypeMatrix(matrix_t* A, matrix_t* B) {
  int fail = OK;
  if (A->columns != B->columns || A->rows != B->rows) {
    fail = CALC_ERROR;
  }
  return fail;
}

int SubAndSumMatrix(matrix_t* A, matrix_t* B, matrix_t* result, int sign) {
  int fail = s21_create_matrix(A->rows, A->columns, result);
  if (ExistMatrix(A) || ExistMatrix(B)) {
    fail = INCORRECT_MATRIX;
  } else if (CheckTypeMatrix(A, B) == 2) {
    fail = CALC_ERROR;
  } else {
    if (fail == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + (sign * (B->matrix[i][j]));
        }
      }
    }
  }
  return fail;
}

int MinorMatrix(matrix_t* A, matrix_t* minor_mat, int n, int m) {
  int fail = OK;
  fail = s21_create_matrix(A->rows - 1, A->columns - 1, minor_mat);
  if (ExistMatrix(A) != 0) {
    fail = INCORRECT_MATRIX;
  } else if (A->rows == 1 && A->columns == 1) {
    fail = CALC_ERROR;
  } else {
    int k = 0;
    int q = 0;
    if (fail == OK) {
      for (int i = 0; i < A->rows; i++) {
        if (i != n) {
          q = 0;
          for (int j = 0; j < A->columns; j++) {
            if (j != m) {
              minor_mat->matrix[k][q] = A->matrix[i][j];
              q++;
            }
          }
          k++;
        }
      }
    }
  }
  return fail;
}

int MinorCalc(matrix_t* A, matrix_t* result) {
  int fail = OK;
  result->matrix[0][0] = 1;
  if (A->rows != 1) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        double det = 0.0;
        matrix_t minor;
        fail = MinorMatrix(A, &minor, i, j);
        if (fail == OK) {
          fail = s21_determinant(&minor, &det);
          if (fail == OK) {
            result->matrix[i][j] = pow((-1), i + j) * det;
          }
        }
        s21_remove_matrix(&minor);
      }
    }
  } else {
    fail = INCORRECT_MATRIX;
  }
  return fail;
}

double Determinant(matrix_t* A) {
  double res = 0.0;
  if (A->rows == 2) {
    res = A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
  } else {
    for (int i = 0; i < A->columns; i++) {
      matrix_t minor_mat;
      MinorMatrix(A, &minor_mat, 0, i);
      res += pow(-1.0, i) * Determinant(&minor_mat) * A->matrix[0][i];
      s21_remove_matrix(&minor_mat);
    }
  }
  return res;
}
