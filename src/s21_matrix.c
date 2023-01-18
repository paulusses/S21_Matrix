#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int fail = OK;
  result->matrix = NULL;
  result->rows = 0;
  result->columns = 0;
  if (rows > 0 && columns > 0) {
    result->matrix = calloc((size_t)rows, sizeof(double *));
    if (result->matrix) {
      result->rows = rows;
      result->columns = columns;
      for (int i = 0; i < result->rows; i++) {
        result->matrix[i] = calloc((size_t)columns, sizeof(double));
        if (result->matrix[i] == NULL) {
          s21_remove_matrix(result);
          fail = INCORRECT_MATRIX;
          break;
        }
      }
    } else {
      fail = INCORRECT_MATRIX;
    }
  } else {
    fail = INCORRECT_MATRIX;
  }
  return fail;
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix) {
    for (int i = 0; i < A->rows; i++) {
      free(A->matrix[i]);
    }
    free(A->matrix);
    A->matrix = NULL;
  }
  if (A->rows) {
    A->rows = 0;
  }
  if (A->columns) {
    A->columns = 0;
  }
  A = NULL;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int fail = SUCCESS;
  if (ExistMatrix(A) == OK && ExistMatrix(B) == OK) {
    if (CheckTypeMatrix(A, B) != OK) {
      fail = FAILURE;
    } else {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-7) {
            fail = FAILURE;
            break;
          }
        }
      }
    }
  } else {
    fail = FAILURE;
  }
  return fail;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  return SubAndSumMatrix(A, B, result, 1);
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  return SubAndSumMatrix(A, B, result, -1);
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int fail = OK;
  fail = s21_create_matrix(A->rows, A->columns, result);
  if (ExistMatrix(A) == OK && fail == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  } else {
    fail = INCORRECT_MATRIX;
  }
  return fail;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int fail = s21_create_matrix(A->rows, B->columns, result);
  if (ExistMatrix(A) != OK || ExistMatrix(B) != OK || fail != OK) {
    fail = INCORRECT_MATRIX;
  } else if (A->columns != B->rows) {
    fail = CALC_ERROR;
  } else {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < B->columns; j++) {
        result->matrix[i][j] = 0;
        for (int k = 0; k < B->rows; k++) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    }
  }
  return fail;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int fail = s21_create_matrix(A->columns, A->rows, result);
  if (ExistMatrix(A)) {
    fail = INCORRECT_MATRIX;
  } else {
    for (int i = 0; i < result->rows; i++) {
      for (int j = 0; j < result->columns; j++) {
        result->matrix[i][j] = A->matrix[j][i];
      }
    }
  }
  return fail;
}

int s21_determinant(matrix_t *A, double *result) {
  int fail = OK;
  if (ExistMatrix(A) != OK) {
    fail = INCORRECT_MATRIX;
  } else if (A->columns != A->rows) {
    fail = CALC_ERROR;
  } else {
    if (A->rows == 1) {
      *result = A->matrix[0][0];
    } else {
      *result = Determinant(A);
    }
  }
  return fail;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int fail = s21_create_matrix(A->rows, A->columns, result);
  if (ExistMatrix(A) != OK) {
    fail = INCORRECT_MATRIX;
  } else if (A->rows != A->columns || A->rows == 1) {
    fail = CALC_ERROR;
  } else {
    if (fail == OK) {
      fail = MinorCalc(A, result);
    } else {
      fail = INCORRECT_MATRIX;
    }
  }
  return fail;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int fail = s21_create_matrix(A->rows, A->columns, result);
  double det = 0.0;
  if (ExistMatrix(A) != 0) {
    fail = INCORRECT_MATRIX;
  } else if (A->columns != A->rows) {
    fail = CALC_ERROR;
  } else {
    fail = s21_determinant(A, &det);
    if (fail == OK) {
      if (A->columns == 1 && A->rows == 1) {
        result->matrix[0][0] = 1 / det;
      } else {
        matrix_t algeb_dop;
        matrix_t transpose;
        fail = s21_calc_complements(A, &algeb_dop);
        if (fail == OK && det != 0) {
          fail = s21_transpose(&algeb_dop, &transpose);
          if (fail == OK) {
            s21_remove_matrix(result);
            s21_mult_number(&transpose, 1 / det, result);
          }
          s21_remove_matrix(&transpose);
        } else {
          fail = CALC_ERROR;
        }
        s21_remove_matrix(&algeb_dop);
      }
    } else {
      fail = CALC_ERROR;
    }
  }
  return fail;
}