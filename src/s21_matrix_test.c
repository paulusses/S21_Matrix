
#include "s21_matrix.h"

#include <check.h>

START_TEST(create_1) {
  matrix_t m1;
  int check = s21_create_matrix(2, 4, &m1);
  int check_origin = 0;
  s21_remove_matrix(&m1);
  ck_assert_int_eq(check, check_origin);
}
END_TEST

START_TEST(create_2) {
  matrix_t m1;
  int check = s21_create_matrix(-2, 4, &m1);
  int check_origin = 1;
  s21_remove_matrix(&m1);
  ck_assert_int_eq(check, check_origin);
}
END_TEST

START_TEST(create_3) {
  matrix_t m1;
  int check = s21_create_matrix(2, -4, &m1);
  int check_origin = 1;
  ck_assert_int_eq(check, check_origin);
}
END_TEST

START_TEST(create_4) {
  matrix_t m1;
  int check = s21_create_matrix(0, 15, &m1);
  int check_origin = 1;
  s21_remove_matrix(&m1);
  ck_assert_int_eq(check, check_origin);
}
END_TEST

START_TEST(create_5) {
  matrix_t m1;
  int check = s21_create_matrix(1, 1, &m1);
  int check_origin = 0;
  s21_remove_matrix(&m1);
  ck_assert_int_eq(check, check_origin);
}
END_TEST

START_TEST(eq_1) {
  matrix_t m1;
  matrix_t m2;
  int check_origin = 0;
  int check_1 = s21_create_matrix(3, 3, &m1);
  int check_2 = s21_create_matrix(3, 3, &m2);
  m1.matrix[0][0] = 0.0;
  m1.matrix[0][1] = 1.0;
  m1.matrix[0][2] = 2.0;
  m1.matrix[1][0] = 3.0;
  m1.matrix[1][1] = 4.0;
  m1.matrix[1][2] = 5.0;
  m1.matrix[2][0] = 6.0;
  m1.matrix[2][1] = 7.0;
  m1.matrix[2][2] = 8.0;
  m2.matrix[0][0] = 0.0;
  m2.matrix[0][1] = 1.0;
  m2.matrix[0][2] = 2.0;
  m2.matrix[1][0] = 3.0;
  m2.matrix[1][1] = 4.0;
  m2.matrix[1][2] = 5.0;
  m2.matrix[2][0] = 6.0;
  m2.matrix[2][1] = 7.0;
  m2.matrix[2][2] = 8.0;
  int result = s21_eq_matrix(&m1, &m2);
  int result_origin = SUCCESS;
  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(result, result_origin);
}
END_TEST

START_TEST(eq_2) {
  matrix_t m1;
  matrix_t m2;
  int check_origin = 0;
  int check_1 = s21_create_matrix(3, 3, &m1);
  int check_2 = s21_create_matrix(3, 3, &m2);
  m1.matrix[0][0] = 987654321.123456789;
  m1.matrix[0][1] = 1111111111111.1111111111111;
  m1.matrix[0][2] = 2.0;
  m1.matrix[1][0] = 3.0000001;
  m1.matrix[1][1] = 4.111111111;
  m1.matrix[1][2] = 5123454321.5123454321;
  m1.matrix[2][0] = 666666666666.0;
  m1.matrix[2][1] = 7.0;
  m1.matrix[2][2] = 987654321.9123456789;
  m2.matrix[0][0] = 987654321.123456789;
  m2.matrix[0][1] = 1111111111111.1111111111111;
  m2.matrix[0][2] = 2.0;
  m2.matrix[1][0] = 3.0000001;
  m2.matrix[1][1] = 4.111111111;
  m2.matrix[1][2] = 5123454321.5123454321;
  m2.matrix[2][0] = 666666666666.0;
  m2.matrix[2][1] = 7.0;
  m2.matrix[2][2] = 987654321.9123456789;
  int result = s21_eq_matrix(&m1, &m2);
  int result_origin = SUCCESS;
  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(result, result_origin);
}
END_TEST

START_TEST(eq_3) {
  matrix_t m1;
  matrix_t m2;
  int check_origin = 0;
  int check_1 = s21_create_matrix(3, 3, &m1);
  int check_2 = s21_create_matrix(3, 3, &m2);
  m1.matrix[0][0] = 987654321.1234567;
  m1.matrix[0][1] = 1111111111111.1111111;
  m1.matrix[0][2] = 2.0000000;
  m1.matrix[1][0] = 3.0000001;
  m1.matrix[1][1] = 4.1111111;
  m1.matrix[1][2] = 5123454321.5123454;
  m1.matrix[2][0] = 666666666666.0000000;
  m1.matrix[2][1] = 7.0000000;
  m1.matrix[2][2] = 987654321.9123456;
  m2.matrix[0][0] = 987654321.1234567;
  m2.matrix[0][1] = 1111111111111.1111111;
  m2.matrix[0][2] = 2.0000000;
  m2.matrix[1][0] = 3.0000001;
  m2.matrix[1][1] = 4.1111111;
  m2.matrix[1][2] = 5123454321.5123454;
  m2.matrix[2][0] = 666666666666.0000000;
  m2.matrix[2][1] = 7.0000000;
  m2.matrix[2][2] = 987654321.9123456;
  int result = s21_eq_matrix(&m1, &m2);
  int result_origin = SUCCESS;
  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(result, result_origin);
}
END_TEST

START_TEST(eq_4) {
  matrix_t m1;
  matrix_t m2;
  int check_origin = 0;
  int check_1 = s21_create_matrix(3, 3, &m1);
  int check_2 = s21_create_matrix(3, 3, &m2);
  m1.matrix[0][0] = 987654321.1234567;
  m1.matrix[0][1] = 1111111111111.111111;
  m1.matrix[0][2] = 2.0000005;
  m1.matrix[1][0] = 3.0000001;
  m1.matrix[1][1] = 4.1111111;
  m1.matrix[1][2] = 5123454321.5123454;
  m1.matrix[2][0] = 666666666666.0;
  m1.matrix[2][1] = 7.0;
  m1.matrix[2][2] = 987654321.9123456;
  m2.matrix[0][0] = 987654321.1234567;
  m2.matrix[0][1] = 1111111111111.111111;
  m2.matrix[0][2] = 2.0;
  m2.matrix[1][0] = 3.0000001;
  m2.matrix[1][1] = 4.1111111;
  m2.matrix[1][2] = 5123454321.5123454;
  m2.matrix[2][0] = 666666666666.0;
  m2.matrix[2][1] = 7.0;
  m2.matrix[2][2] = 987654321.9123456;
  int result = s21_eq_matrix(&m1, &m2);
  int result_origin = FAILURE;
  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(result, result_origin);
}
END_TEST

START_TEST(s21_eq_matrix_1) {
  matrix_t A;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 3.3333335;
  A.matrix[0][1] = -5;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = -3;
  matrix_t B;
  s21_create_matrix(2, 2, &B);
  B.matrix[0][0] = 3.33333356;
  B.matrix[0][1] = -5;
  B.matrix[1][0] = 6;
  B.matrix[1][1] = -3;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(result, 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_eq_matrix_2) {
  matrix_t A;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 3.33333357;
  A.matrix[0][1] = -5;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = -3;
  matrix_t B;
  s21_create_matrix(2, 2, &B);
  B.matrix[0][0] = 3.33333356;
  B.matrix[0][1] = -5;
  B.matrix[1][0] = 6;
  B.matrix[1][1] = -3;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(result, 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_eq_matrix_3) {
  matrix_t A;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 3.0;
  A.matrix[0][1] = -5;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = -3;
  matrix_t B;
  s21_create_matrix(2, 2, &B);
  B.matrix[0][0] = 2.0;
  B.matrix[0][1] = -5;
  B.matrix[1][0] = 6;
  B.matrix[1][1] = -3;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(result, 0);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_eq_matrix_4) {
  matrix_t A;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 3.3333333;
  A.matrix[0][1] = -5;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3.3333334;
  matrix_t B;
  s21_create_matrix(2, 2, &B);
  B.matrix[0][0] = 3.3333333;
  B.matrix[0][1] = -5;
  B.matrix[1][0] = 6;
  B.matrix[1][1] = 3.3333335;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(result, 0);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_eq_matrix_6) {
  matrix_t A;
  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = 1;
  matrix_t B;
  s21_create_matrix(1, 1, &B);
  B.matrix[0][0] = 1;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(result, 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_eq_matrix_7) {
  matrix_t A;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 1;
  A.matrix[1][0] = 1;
  A.matrix[1][1] = 1;
  matrix_t B;
  s21_create_matrix(1, 1, &B);
  B.matrix[0][0] = 1;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(result, 0);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(sum_1) {
  matrix_t m1;
  matrix_t m2;
  matrix_t origin;
  matrix_t r;
  int check_origin = 0;
  int check_1 = s21_create_matrix(3, 3, &m1);
  int check_2 = s21_create_matrix(3, 3, &m2);
  int check_3 = s21_create_matrix(3, 3, &origin);
  m1.matrix[0][0] = 0.1;
  m1.matrix[0][1] = 0.2;
  m1.matrix[0][2] = 0.3;
  m1.matrix[1][0] = 0.4;
  m1.matrix[1][1] = 0.5;
  m1.matrix[1][2] = 0.6;
  m1.matrix[2][0] = 0.7;
  m1.matrix[2][1] = 0.8;
  m1.matrix[2][2] = 0.9;
  m2.matrix[0][0] = 1.0;
  m2.matrix[0][1] = 2.0;
  m2.matrix[0][2] = 3.0;
  m2.matrix[1][0] = 4.0;
  m2.matrix[1][1] = 5.0;
  m2.matrix[1][2] = 6.0;
  m2.matrix[2][0] = 7.0;
  m2.matrix[2][1] = 8.0;
  m2.matrix[2][2] = 9.0;
  origin.matrix[0][0] = 1.1;
  origin.matrix[0][1] = 2.2;
  origin.matrix[0][2] = 3.3;
  origin.matrix[1][0] = 4.4;
  origin.matrix[1][1] = 5.5;
  origin.matrix[1][2] = 6.6;
  origin.matrix[2][0] = 7.7;
  origin.matrix[2][1] = 8.8;
  origin.matrix[2][2] = 9.9;
  int result = s21_sum_matrix(&m1, &m2, &r);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(check_3, check_origin);
  ck_assert_int_eq(result, check_origin);
  for (int i = 0; i < r.rows; i++) {
    for (int j = 0; j < r.columns; j++) {
      ck_assert_double_eq_tol(origin.matrix[i][j], r.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
  s21_remove_matrix(&r);
  s21_remove_matrix(&origin);
}
END_TEST

START_TEST(sum_2) {
  matrix_t m1;
  matrix_t m2;
  matrix_t origin;
  matrix_t r;
  int check_origin = 0;
  int check_1 = s21_create_matrix(3, 3, &m1);
  int check_2 = s21_create_matrix(3, 3, &m2);
  int check_3 = s21_create_matrix(3, 3, &origin);
  m1.matrix[0][0] = 0.1;
  m1.matrix[0][1] = 0.2;
  m1.matrix[0][2] = 0.3;
  m1.matrix[1][0] = 0.4;
  m1.matrix[1][1] = 0.5;
  m1.matrix[1][2] = 0.6;
  m1.matrix[2][0] = 0.7;
  m1.matrix[2][1] = 0.8;
  m1.matrix[2][2] = 0.9;
  m2.matrix[0][0] = 0.0;
  m2.matrix[0][1] = 0.0;
  m2.matrix[0][2] = 0.0;
  m2.matrix[1][0] = 0.0;
  m2.matrix[1][1] = 0.0;
  m2.matrix[1][2] = 0.0;
  m2.matrix[2][0] = 0.0;
  m2.matrix[2][1] = 0.0;
  m2.matrix[2][2] = 0.0;
  origin.matrix[0][0] = 0.1;
  origin.matrix[0][1] = 0.2;
  origin.matrix[0][2] = 0.3;
  origin.matrix[1][0] = 0.4;
  origin.matrix[1][1] = 0.5;
  origin.matrix[1][2] = 0.6;
  origin.matrix[2][0] = 0.7;
  origin.matrix[2][1] = 0.8;
  origin.matrix[2][2] = 0.9;
  int result = s21_sum_matrix(&m1, &m2, &r);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(check_3, check_origin);
  ck_assert_int_eq(result, check_origin);
  for (int i = 0; i < r.rows; i++) {
    for (int j = 0; j < r.columns; j++) {
      ck_assert_double_eq_tol(origin.matrix[i][j], r.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
  s21_remove_matrix(&r);
  s21_remove_matrix(&origin);
}
END_TEST

START_TEST(sum_3) {
  matrix_t m1;
  matrix_t m2;
  matrix_t origin;
  matrix_t r;
  int check_origin = 0;
  int check_1 = s21_create_matrix(3, 3, &m1);
  int check_2 = s21_create_matrix(3, 3, &m2);
  int check_3 = s21_create_matrix(3, 3, &origin);
  m1.matrix[0][0] = 0.1;
  m1.matrix[0][1] = 0.2;
  m1.matrix[0][2] = 0.3;
  m1.matrix[1][0] = 0.4;
  m1.matrix[1][1] = 0.5;
  m1.matrix[1][2] = 0.6;
  m1.matrix[2][0] = 0.7;
  m1.matrix[2][1] = 0.8;
  m1.matrix[2][2] = 0.9;
  m2.matrix[0][0] = -1.0;
  m2.matrix[0][1] = -2.0;
  m2.matrix[0][2] = -3.0;
  m2.matrix[1][0] = -4.0;
  m2.matrix[1][1] = -5.0;
  m2.matrix[1][2] = -6.0;
  m2.matrix[2][0] = -7.0;
  m2.matrix[2][1] = -8.0;
  m2.matrix[2][2] = -9.0;
  origin.matrix[0][0] = -0.9;
  origin.matrix[0][1] = -1.8;
  origin.matrix[0][2] = -2.7;
  origin.matrix[1][0] = -3.6;
  origin.matrix[1][1] = -4.5;
  origin.matrix[1][2] = -5.4;
  origin.matrix[2][0] = -6.3;
  origin.matrix[2][1] = -7.2;
  origin.matrix[2][2] = -8.1;
  int result = s21_sum_matrix(&m1, &m2, &r);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(check_3, check_origin);
  ck_assert_int_eq(result, check_origin);
  for (int i = 0; i < r.rows; i++) {
    for (int j = 0; j < r.columns; j++) {
      ck_assert_double_eq_tol(origin.matrix[i][j], r.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
  s21_remove_matrix(&r);
  s21_remove_matrix(&origin);
}
END_TEST

START_TEST(s21_sum_matrix_1) {
  matrix_t result, B, A;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 1;
  A.matrix[1][0] = 1;
  A.matrix[1][1] = 1;
  B.matrix[0][0] = 3;
  B.matrix[0][1] = 3;
  B.matrix[1][0] = 3;
  B.matrix[1][1] = 3;
  int check = s21_sum_matrix(&A, &B, &result);
  int check_origin = 0;
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < A.rows; i++) {
    for (int j = 0; j < B.columns; j++) {
      ck_assert_int_eq(result.matrix[i][j], 4);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(s21_sum_matrix_2) {
  matrix_t result, B, A;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 0;
  A.matrix[0][1] = 0;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 0;
  B.matrix[0][0] = 0;
  B.matrix[0][1] = 0;
  B.matrix[1][0] = 0;
  B.matrix[1][1] = 0;
  int check = s21_sum_matrix(&A, &B, &result);
  int check_origin = 0;
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < A.rows; i++) {
    for (int j = 0; j < B.columns; j++) {
      ck_assert_int_eq(result.matrix[i][j], 0);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(s21_sum_matrix_3) {
  matrix_t result, B, A;
  s21_create_matrix(1, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 1;
  B.matrix[0][0] = 3;
  B.matrix[0][1] = 3;
  B.matrix[1][0] = 3;
  B.matrix[1][1] = 3;
  int check = s21_sum_matrix(&A, &B, &result);
  int check_origin = 2;
  ck_assert_int_eq(check, check_origin);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(s21_sum_matrix_4) {
  matrix_t result, res_origin, B, A;
  s21_create_matrix(4, 4, &A);
  s21_create_matrix(4, 4, &B);
  s21_create_matrix(4, 4, &res_origin);
  A.matrix[0][0] = 2;
  A.matrix[0][1] = 5;
  A.matrix[0][2] = 7;
  A.matrix[0][3] = -4;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;
  A.matrix[1][3] = 14;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;
  A.matrix[2][3] = -1;
  A.matrix[3][0] = 10;
  A.matrix[3][1] = 3;
  A.matrix[3][2] = 9;
  A.matrix[3][3] = -3;
  B.matrix[0][0] = 8;
  B.matrix[0][1] = 8;
  B.matrix[0][2] = 3;
  B.matrix[0][3] = 5;
  B.matrix[1][0] = -7;
  B.matrix[1][1] = 7;
  B.matrix[1][2] = 1;
  B.matrix[1][3] = 1;
  B.matrix[2][0] = 7;
  B.matrix[2][1] = 3;
  B.matrix[2][2] = 4;
  B.matrix[2][3] = -4;
  B.matrix[3][0] = 10;
  B.matrix[3][1] = 3;
  B.matrix[3][2] = 9;
  B.matrix[3][3] = -3;
  res_origin.matrix[0][0] = 10;
  res_origin.matrix[0][1] = 13;
  res_origin.matrix[0][2] = 10;
  res_origin.matrix[0][3] = 1;
  res_origin.matrix[1][0] = -1;
  res_origin.matrix[1][1] = 10;
  res_origin.matrix[1][2] = 5;
  res_origin.matrix[1][3] = 15;
  res_origin.matrix[2][0] = 12;
  res_origin.matrix[2][1] = 1;
  res_origin.matrix[2][2] = 1;
  res_origin.matrix[2][3] = -5;
  res_origin.matrix[3][0] = 20;
  res_origin.matrix[3][1] = 6;
  res_origin.matrix[3][2] = 18;
  res_origin.matrix[3][3] = -6;
  int check = s21_sum_matrix(&A, &B, &result);
  int check_origin = 0;
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_uint_eq(res_origin.matrix[i][j], result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_sum_matrix_5) {
  matrix_t result, res_origin, B, A;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  s21_create_matrix(2, 2, &res_origin);
  A.matrix[0][0] = 2.987654;
  A.matrix[0][1] = -5.3456;
  A.matrix[1][0] = 6.23;
  A.matrix[1][1] = -3.9292;
  B.matrix[0][0] = 4.838;
  B.matrix[0][1] = 8.0987;
  B.matrix[1][0] = -8;
  B.matrix[1][1] = -0.921;
  res_origin.matrix[0][0] = 7.825654;
  res_origin.matrix[0][1] = 2.7531;
  res_origin.matrix[1][0] = -1.77;
  res_origin.matrix[1][1] = -4.8502;
  int check = s21_sum_matrix(&A, &B, &result);
  int check_origin = 0;
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_uint_eq(res_origin.matrix[i][j], result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_sum_matrix_7) {
  matrix_t result, res_origin, B, A;
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(5, 5, &B);
  s21_create_matrix(5, 5, &res_origin);
  for (int i = 0; i < A.rows; i++) {
    for (int j = 0; j < A.columns; j++) {
      A.matrix[i][j] = i;
    }
  }
  for (int i = 0; i < B.rows; i++) {
    for (int j = 0; j < B.columns; j++) {
      B.matrix[i][j] = i;
    }
  }
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      res_origin.matrix[i][j] = i * 2;
    }
  }
  int check = s21_sum_matrix(&A, &B, &result);
  int check_origin = 0;
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_uint_eq(res_origin.matrix[i][j], result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_sum_matrix_8) {
  matrix_t result, res_origin, B, A;
  s21_create_matrix(5, 4, &A);
  s21_create_matrix(5, 4, &B);
  s21_create_matrix(5, 4, &res_origin);
  for (int i = 0; i < A.rows; i++) {
    for (int j = 0; j < A.columns; j++) {
      A.matrix[i][j] = i;
      res_origin.matrix[i][j] = i;
    }
  }
  int check = s21_sum_matrix(&A, &B, &result);
  int check_origin = 0;
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_uint_eq(res_origin.matrix[i][j], result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_sum_matrix_9) {
  matrix_t result, B, A;
  s21_create_matrix(-5, 5, &A);
  s21_create_matrix(5, 4, &B);
  int check = s21_sum_matrix(&A, &B, &result);
  int check_origin = 1;
  ck_assert_int_eq(check, check_origin);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(sub_1) {
  matrix_t m1;
  matrix_t m2;
  matrix_t origin;
  matrix_t r;
  int check_origin = 0;
  int check_1 = s21_create_matrix(3, 3, &m1);
  int check_2 = s21_create_matrix(3, 3, &m2);
  int check_3 = s21_create_matrix(3, 3, &origin);
  m1.matrix[0][0] = 0.1;
  m1.matrix[0][1] = 0.2;
  m1.matrix[0][2] = 0.3;
  m1.matrix[1][0] = 0.4;
  m1.matrix[1][1] = 0.5;
  m1.matrix[1][2] = 0.6;
  m1.matrix[2][0] = 0.7;
  m1.matrix[2][1] = 0.8;
  m1.matrix[2][2] = 0.9;
  m2.matrix[0][0] = 1.0;
  m2.matrix[0][1] = 2.0;
  m2.matrix[0][2] = 3.0;
  m2.matrix[1][0] = 4.0;
  m2.matrix[1][1] = 5.0;
  m2.matrix[1][2] = 6.0;
  m2.matrix[2][0] = 7.0;
  m2.matrix[2][1] = 8.0;
  m2.matrix[2][2] = 9.0;
  origin.matrix[0][0] = -0.9;
  origin.matrix[0][1] = -1.8;
  origin.matrix[0][2] = -2.7;
  origin.matrix[1][0] = -3.6;
  origin.matrix[1][1] = -4.5;
  origin.matrix[1][2] = -5.4;
  origin.matrix[2][0] = -6.3;
  origin.matrix[2][1] = -7.2;
  origin.matrix[2][2] = -8.1;
  int result = s21_sub_matrix(&m1, &m2, &r);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(check_3, check_origin);
  ck_assert_int_eq(result, check_origin);
  for (int i = 0; i < r.rows; i++) {
    for (int j = 0; j < r.columns; j++) {
      ck_assert_double_eq_tol(origin.matrix[i][j], r.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
  s21_remove_matrix(&r);
  s21_remove_matrix(&origin);
}
END_TEST

START_TEST(sub_2) {
  matrix_t m1;
  matrix_t m2;
  matrix_t origin;
  matrix_t r;
  int check_origin = 0;
  int check_1 = s21_create_matrix(3, 3, &m1);
  int check_2 = s21_create_matrix(3, 3, &m2);
  int check_3 = s21_create_matrix(3, 3, &origin);
  m1.matrix[0][0] = 0.1;
  m1.matrix[0][1] = 0.2;
  m1.matrix[0][2] = 0.3;
  m1.matrix[1][0] = 0.4;
  m1.matrix[1][1] = 0.5;
  m1.matrix[1][2] = 0.6;
  m1.matrix[2][0] = 0.7;
  m1.matrix[2][1] = 0.8;
  m1.matrix[2][2] = 0.9;
  m2.matrix[0][0] = 0.0;
  m2.matrix[0][1] = 0.0;
  m2.matrix[0][2] = 0.0;
  m2.matrix[1][0] = 0.0;
  m2.matrix[1][1] = 0.0;
  m2.matrix[1][2] = 0.0;
  m2.matrix[2][0] = 0.0;
  m2.matrix[2][1] = 0.0;
  m2.matrix[2][2] = 0.0;
  origin.matrix[0][0] = 0.1;
  origin.matrix[0][1] = 0.2;
  origin.matrix[0][2] = 0.3;
  origin.matrix[1][0] = 0.4;
  origin.matrix[1][1] = 0.5;
  origin.matrix[1][2] = 0.6;
  origin.matrix[2][0] = 0.7;
  origin.matrix[2][1] = 0.8;
  origin.matrix[2][2] = 0.9;
  int result = s21_sub_matrix(&m1, &m2, &r);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(check_3, check_origin);
  ck_assert_int_eq(result, check_origin);
  for (int i = 0; i < r.rows; i++) {
    for (int j = 0; j < r.columns; j++) {
      ck_assert_double_eq_tol(origin.matrix[i][j], r.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
  s21_remove_matrix(&r);
  s21_remove_matrix(&origin);
}
END_TEST

START_TEST(sub_3) {
  matrix_t m1;
  matrix_t m2;
  matrix_t origin;
  matrix_t r;
  int check_origin = 0;
  int check_1 = s21_create_matrix(3, 3, &m1);
  int check_2 = s21_create_matrix(3, 3, &m2);
  int check_3 = s21_create_matrix(3, 3, &origin);
  m1.matrix[0][0] = 0.1;
  m1.matrix[0][1] = 0.2;
  m1.matrix[0][2] = 0.3;
  m1.matrix[1][0] = 0.4;
  m1.matrix[1][1] = 0.5;
  m1.matrix[1][2] = 0.6;
  m1.matrix[2][0] = 0.7;
  m1.matrix[2][1] = 0.8;
  m1.matrix[2][2] = 0.9;
  m2.matrix[0][0] = 1.0;
  m2.matrix[0][1] = 2.0;
  m2.matrix[0][2] = 3.0;
  m2.matrix[1][0] = 4.0;
  m2.matrix[1][1] = 5.0;
  m2.matrix[1][2] = 6.0;
  m2.matrix[2][0] = 7.0;
  m2.matrix[2][1] = 8.0;
  m2.matrix[2][2] = 9.0;
  origin.matrix[0][0] = -0.9;
  origin.matrix[0][1] = -1.8;
  origin.matrix[0][2] = -2.7;
  origin.matrix[1][0] = -3.6;
  origin.matrix[1][1] = -4.5;
  origin.matrix[1][2] = -5.4;
  origin.matrix[2][0] = -6.3;
  origin.matrix[2][1] = -7.2;
  origin.matrix[2][2] = -8.1;
  int result = s21_sub_matrix(&m1, &m2, &r);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(check_3, check_origin);
  ck_assert_int_eq(result, check_origin);
  for (int i = 0; i < r.rows; i++) {
    for (int j = 0; j < r.columns; j++) {
      ck_assert_double_eq_tol(origin.matrix[i][j], r.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
  s21_remove_matrix(&r);
  s21_remove_matrix(&origin);
}
END_TEST

START_TEST(sub_4) {
  matrix_t m1;
  matrix_t m2;
  matrix_t origin;
  matrix_t r;
  int check_origin = 0;
  int check_1 = s21_create_matrix(3, 3, &m1);
  int check_2 = s21_create_matrix(3, 3, &m2);
  int check_3 = s21_create_matrix(3, 3, &origin);
  m1.matrix[0][0] = 0.1;
  m1.matrix[0][1] = 0.2;
  m1.matrix[0][2] = 0.3;
  m1.matrix[1][0] = 0.4;
  m1.matrix[1][1] = 0.5;
  m1.matrix[1][2] = 0.6;
  m1.matrix[2][0] = 0.7;
  m1.matrix[2][1] = 0.8;
  m1.matrix[2][2] = 0.9;
  m2.matrix[0][0] = 0.0;
  m2.matrix[0][1] = 0.0;
  m2.matrix[0][2] = 0.0;
  m2.matrix[1][0] = 0.0;
  m2.matrix[1][1] = 0.0;
  m2.matrix[1][2] = 0.0;
  m2.matrix[2][0] = 0.0;
  m2.matrix[2][1] = 0.0;
  m2.matrix[2][2] = 0.0;
  origin.matrix[0][0] = 0.1;
  origin.matrix[0][1] = 0.2;
  origin.matrix[0][2] = 0.3;
  origin.matrix[1][0] = 0.4;
  origin.matrix[1][1] = 0.5;
  origin.matrix[1][2] = 0.6;
  origin.matrix[2][0] = 0.7;
  origin.matrix[2][1] = 0.8;
  origin.matrix[2][2] = 0.9;
  int result = s21_sub_matrix(&m2, &m1, &r);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(check_3, check_origin);
  ck_assert_int_eq(result, check_origin);
  for (int i = 0; i < r.rows; i++) {
    for (int j = 0; j < r.columns; j++) {
      ck_assert_double_eq_tol(-(origin.matrix[i][j]), r.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&m1);
  s21_remove_matrix(&m2);
  s21_remove_matrix(&r);
  s21_remove_matrix(&origin);
}
END_TEST

START_TEST(s21_sub_matrix_1) {
  matrix_t result, B, A;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 1;
  A.matrix[1][0] = 1;
  A.matrix[1][1] = 1;
  B.matrix[0][0] = 3;
  B.matrix[0][1] = 3;
  B.matrix[1][0] = 3;
  B.matrix[1][1] = 3;
  int check = s21_sub_matrix(&A, &B, &result);
  int check_origin = 0;
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < A.rows; i++) {
    for (int j = 0; j < B.columns; j++) {
      ck_assert_int_eq(result.matrix[i][j], -2);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(s21_sub_matrix_2) {
  matrix_t result, B, A;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 0;
  A.matrix[0][1] = 0;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 0;
  B.matrix[0][0] = 0;
  B.matrix[0][1] = 0;
  B.matrix[1][0] = 0;
  B.matrix[1][1] = 0;
  int check = s21_sub_matrix(&A, &B, &result);
  int check_origin = 0;
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < A.rows; i++) {
    for (int j = 0; j < B.columns; j++) {
      ck_assert_int_eq(result.matrix[i][j], 0);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(s21_sub_matrix_3) {
  matrix_t result, B, A;
  s21_create_matrix(1, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 1;
  B.matrix[0][0] = 3;
  B.matrix[0][1] = 3;
  B.matrix[1][0] = 3;
  B.matrix[1][1] = 3;
  int check = s21_sub_matrix(&A, &B, &result);
  int check_origin = 2;
  ck_assert_int_eq(check, check_origin);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(s21_sub_matrix_9) {
  matrix_t result, B, A;
  s21_create_matrix(-2, 2, &A);
  s21_create_matrix(2, 2, &B);
  int check = s21_sub_matrix(&A, &B, &result);
  int check_origin = 1;
  ck_assert_int_eq(check, check_origin);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(s21_sub_matrix_10) {
  matrix_t result, B, A;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(-2, 2, &B);
  int check = s21_sub_matrix(&A, &B, &result);
  int check_origin = 1;
  ck_assert_int_eq(check, check_origin);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(s21_sub_matrix_4) {
  matrix_t result, res_origin, B, A;
  s21_create_matrix(4, 4, &A);
  s21_create_matrix(4, 4, &B);
  s21_create_matrix(4, 4, &res_origin);
  A.matrix[0][0] = 2;
  A.matrix[0][1] = 5;
  A.matrix[0][2] = 7;
  A.matrix[0][3] = -4;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;
  A.matrix[1][3] = 14;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;
  A.matrix[2][3] = -1;
  A.matrix[3][0] = 10;
  A.matrix[3][1] = 3;
  A.matrix[3][2] = 9;
  A.matrix[3][3] = -3;
  B.matrix[0][0] = 8;
  B.matrix[0][1] = 8;
  B.matrix[0][2] = 3;
  B.matrix[0][3] = 5;
  B.matrix[1][0] = -7;
  B.matrix[1][1] = 7;
  B.matrix[1][2] = 1;
  B.matrix[1][3] = 1;
  B.matrix[2][0] = 7;
  B.matrix[2][1] = 3;
  B.matrix[2][2] = 4;
  B.matrix[2][3] = -4;
  B.matrix[3][0] = 10;
  B.matrix[3][1] = 3;
  B.matrix[3][2] = 9;
  B.matrix[3][3] = -3;
  res_origin.matrix[0][0] = -6;
  res_origin.matrix[0][1] = -3;
  res_origin.matrix[0][2] = 4;
  res_origin.matrix[0][3] = -9;
  res_origin.matrix[1][0] = 13;
  res_origin.matrix[1][1] = -4;
  res_origin.matrix[1][2] = 3;
  res_origin.matrix[1][3] = 13;
  res_origin.matrix[2][0] = -2;
  res_origin.matrix[2][1] = -5;
  res_origin.matrix[2][2] = -7;
  res_origin.matrix[2][3] = 3;
  res_origin.matrix[3][0] = 0;
  res_origin.matrix[3][1] = 0;
  res_origin.matrix[3][2] = 0;
  res_origin.matrix[3][3] = 0;
  int check = s21_sub_matrix(&A, &B, &result);
  int check_origin = 0;
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_uint_eq(res_origin.matrix[i][j], result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_sub_matrix_5) {
  matrix_t result, res_origin, B, A;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  s21_create_matrix(2, 2, &res_origin);
  A.matrix[0][0] = 2.987654;
  A.matrix[0][1] = -5.3456;
  A.matrix[1][0] = 6.23;
  A.matrix[1][1] = -3.9292;
  B.matrix[0][0] = 4.838;
  B.matrix[0][1] = 8.0987;
  B.matrix[1][0] = -8;
  B.matrix[1][1] = -0.921;
  res_origin.matrix[0][0] = -1.850346;
  res_origin.matrix[0][1] = -13.4443;
  res_origin.matrix[1][0] = 14.23;
  res_origin.matrix[1][1] = -3.0082;
  int check = s21_sub_matrix(&A, &B, &result);
  int check_origin = 0;
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_uint_eq(res_origin.matrix[i][j], result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_sub_matrix_7) {
  matrix_t result, res_origin, B, A;
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(5, 5, &B);
  s21_create_matrix(5, 5, &res_origin);
  for (int i = 0; i < A.rows; i++) {
    for (int j = 0; j < A.columns; j++) {
      A.matrix[i][j] = 4 * i;
    }
  }
  for (int i = 0; i < B.rows; i++) {
    for (int j = 0; j < B.columns; j++) {
      B.matrix[i][j] = 2 * i;
    }
  }
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      res_origin.matrix[i][j] = i * 2;
    }
  }
  int check = s21_sub_matrix(&A, &B, &result);
  int check_origin = 0;
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_uint_eq(res_origin.matrix[i][j], result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_sub_matrix_8) {
  matrix_t result, res_origin, B, A;
  s21_create_matrix(5, 4, &A);
  s21_create_matrix(5, 4, &B);
  s21_create_matrix(5, 4, &res_origin);
  for (int i = 0; i < A.rows; i++) {
    for (int j = 0; j < A.columns; j++) {
      A.matrix[i][j] = i;
      res_origin.matrix[i][j] = i;
    }
  }
  int check = s21_sub_matrix(&A, &B, &result);
  int check_origin = 0;
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_uint_eq(res_origin.matrix[i][j], result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(mul_num_1) {
  matrix_t m1;
  matrix_t origin;
  matrix_t r;
  double num = 2.0;
  int check_origin = 0;
  int check_1 = s21_create_matrix(3, 3, &m1);
  int check_2 = s21_create_matrix(3, 3, &origin);
  m1.matrix[0][0] = 0.1;
  m1.matrix[0][1] = 0.2;
  m1.matrix[0][2] = 0.3;
  m1.matrix[1][0] = 0.4;
  m1.matrix[1][1] = 0.5;
  m1.matrix[1][2] = 0.6;
  m1.matrix[2][0] = 0.7;
  m1.matrix[2][1] = 0.8;
  m1.matrix[2][2] = 0.9;
  origin.matrix[0][0] = 0.2;
  origin.matrix[0][1] = 0.4;
  origin.matrix[0][2] = 0.6;
  origin.matrix[1][0] = 0.8;
  origin.matrix[1][1] = 1.0;
  origin.matrix[1][2] = 1.2;
  origin.matrix[2][0] = 1.4;
  origin.matrix[2][1] = 1.6;
  origin.matrix[2][2] = 1.8;
  int result = s21_mult_number(&m1, num, &r);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(result, check_origin);
  for (int i = 0; i < r.rows; i++) {
    for (int j = 0; j < r.columns; j++) {
      ck_assert_double_eq_tol(origin.matrix[i][j], r.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&m1);
  s21_remove_matrix(&r);
  s21_remove_matrix(&origin);
}
END_TEST

START_TEST(mul_num_2) {
  matrix_t m1;
  matrix_t r;
  double num = 0.0;
  int check_origin = 0;
  int check_1 = s21_create_matrix(3, 3, &m1);
  m1.matrix[0][0] = 0.1;
  m1.matrix[0][1] = 0.2;
  m1.matrix[0][2] = 0.3;
  m1.matrix[1][0] = 0.4;
  m1.matrix[1][1] = 0.5;
  m1.matrix[1][2] = 0.6;
  m1.matrix[2][0] = 0.7;
  m1.matrix[2][1] = 0.8;
  m1.matrix[2][2] = 0.9;
  int result = s21_mult_number(&m1, num, &r);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(result, check_origin);
  for (int i = 0; i < r.rows; i++) {
    for (int j = 0; j < r.columns; j++) {
      ck_assert_double_eq_tol(0.0, r.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&m1);
  s21_remove_matrix(&r);
}
END_TEST

START_TEST(mul_num_3) {
  matrix_t m1;
  matrix_t r;
  double num = 1.0;
  int check_origin = 0;
  int check_1 = s21_create_matrix(3, 3, &m1);
  m1.matrix[0][0] = 0.1;
  m1.matrix[0][1] = 0.2;
  m1.matrix[0][2] = 0.3;
  m1.matrix[1][0] = 0.4;
  m1.matrix[1][1] = 0.5;
  m1.matrix[1][2] = 0.6;
  m1.matrix[2][0] = 0.7;
  m1.matrix[2][1] = 0.8;
  m1.matrix[2][2] = 0.9;
  int result = s21_mult_number(&m1, num, &r);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(result, check_origin);
  for (int i = 0; i < r.rows; i++) {
    for (int j = 0; j < r.columns; j++) {
      ck_assert_double_eq_tol(m1.matrix[i][j], r.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&m1);
  s21_remove_matrix(&r);
}
END_TEST

START_TEST(s21_mult_number_1) {
  matrix_t A, result, res_origin;
  double num = 2.0;
  int check_origin = 0;
  int check_1 = s21_create_matrix(5, 4, &A);
  int check_3 = s21_create_matrix(5, 4, &res_origin);
  A.matrix[0][0] = 2;
  A.matrix[0][1] = 5;
  A.matrix[0][2] = 7;
  A.matrix[0][3] = -4;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;
  A.matrix[1][3] = 14;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;
  A.matrix[2][3] = -1;
  A.matrix[3][0] = 10;
  A.matrix[3][1] = 3;
  A.matrix[3][2] = 9;
  A.matrix[3][3] = -3;
  A.matrix[4][0] = -3;
  A.matrix[4][1] = 3;
  A.matrix[4][2] = 1;
  A.matrix[4][3] = 5;
  res_origin.matrix[0][0] = 4;
  res_origin.matrix[0][1] = 10;
  res_origin.matrix[0][2] = 14;
  res_origin.matrix[0][3] = -8;
  res_origin.matrix[1][0] = 12;
  res_origin.matrix[1][1] = 6;
  res_origin.matrix[1][2] = 8;
  res_origin.matrix[1][3] = 28;
  res_origin.matrix[2][0] = 10;
  res_origin.matrix[2][1] = -4;
  res_origin.matrix[2][2] = -6;
  res_origin.matrix[2][3] = -2;
  res_origin.matrix[3][0] = 20;
  res_origin.matrix[3][1] = 6;
  res_origin.matrix[3][2] = 18;
  res_origin.matrix[3][3] = -6;
  res_origin.matrix[4][0] = -6;
  res_origin.matrix[4][1] = 6;
  res_origin.matrix[4][2] = 2;
  res_origin.matrix[4][3] = 10;
  int check = s21_mult_number(&A, num, &result);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_3, check_origin);
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_uint_eq(res_origin.matrix[i][j], result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_mult_number_2) {
  matrix_t A, result, res_origin;
  double num = 4.25;
  int check_origin = 0;
  int check_1 = s21_create_matrix(5, 4, &A);
  int check_3 = s21_create_matrix(5, 4, &res_origin);
  A.matrix[0][0] = 2;
  A.matrix[0][1] = 5;
  A.matrix[0][2] = 7;
  A.matrix[0][3] = -4;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;
  A.matrix[1][3] = 14;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;
  A.matrix[2][3] = -1;
  A.matrix[3][0] = 10;
  A.matrix[3][1] = 3;
  A.matrix[3][2] = 9;
  A.matrix[3][3] = -3;
  A.matrix[4][0] = -3;
  A.matrix[4][1] = 3;
  A.matrix[4][2] = 1;
  A.matrix[4][3] = 5;
  res_origin.matrix[0][0] = 8.5;
  res_origin.matrix[0][1] = 21.25;
  res_origin.matrix[0][2] = 29.75;
  res_origin.matrix[0][3] = -17;
  res_origin.matrix[1][0] = 25.5;
  res_origin.matrix[1][1] = 12.75;
  res_origin.matrix[1][2] = 17;
  res_origin.matrix[1][3] = 59.5;
  res_origin.matrix[2][0] = 21.25;
  res_origin.matrix[2][1] = -8.5;
  res_origin.matrix[2][2] = -12.75;
  res_origin.matrix[2][3] = -4.25;
  res_origin.matrix[3][0] = 42.5;
  res_origin.matrix[3][1] = 12.75;
  res_origin.matrix[3][2] = 38.25;
  res_origin.matrix[3][3] = -12.75;
  res_origin.matrix[4][0] = -12.75;
  res_origin.matrix[4][1] = 12.75;
  res_origin.matrix[4][2] = 4.25;
  res_origin.matrix[4][3] = 21.25;
  int check = s21_mult_number(&A, num, &result);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_3, check_origin);
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_uint_eq(res_origin.matrix[i][j], result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_mult_number_3) {
  matrix_t A, result, res_origin;
  double num = 4.25;
  int check_origin = 1;
  s21_create_matrix(-5, 4, &A);
  s21_create_matrix(5, 4, &res_origin);
  int check = s21_mult_number(&A, num, &result);
  ck_assert_int_eq(check, check_origin);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(transpose_1) {
  matrix_t m1;
  matrix_t r;
  matrix_t origin;
  int check_origin = 0;
  int check_1 = s21_create_matrix(3, 2, &m1);
  int check_2 = s21_create_matrix(2, 3, &origin);
  m1.matrix[0][0] = 1;
  m1.matrix[0][1] = 4;
  m1.matrix[1][0] = 2;
  m1.matrix[1][1] = 5;
  m1.matrix[2][0] = 3;
  m1.matrix[2][1] = 6;
  origin.matrix[0][0] = 1;
  origin.matrix[0][1] = 2;
  origin.matrix[0][2] = 3;
  origin.matrix[1][0] = 4;
  origin.matrix[1][1] = 5;
  origin.matrix[1][2] = 6;
  int result = s21_transpose(&m1, &r);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(result, check_origin);
  for (int i = 0; i < r.rows; i++) {
    for (int j = 0; j < r.columns; j++) {
      ck_assert_double_eq_tol(origin.matrix[i][j], r.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&m1);
  s21_remove_matrix(&r);
  s21_remove_matrix(&origin);
}
END_TEST

START_TEST(transpose_2) {
  // 1 4 7    1 2 3
  // 2 5 8    4 5 6
  // 3 6 9    7 8 9
  matrix_t m1;
  matrix_t r;
  matrix_t origin;
  int check_origin = 0;
  int check_1 = s21_create_matrix(3, 3, &m1);
  int check_2 = s21_create_matrix(3, 3, &origin);
  m1.matrix[0][0] = 1;
  m1.matrix[0][1] = 4;
  m1.matrix[0][2] = 7;
  m1.matrix[1][0] = 2;
  m1.matrix[1][1] = 5;
  m1.matrix[1][2] = 8;
  m1.matrix[2][0] = 3;
  m1.matrix[2][1] = 6;
  m1.matrix[2][2] = 9;
  origin.matrix[0][0] = 1;
  origin.matrix[0][1] = 2;
  origin.matrix[0][2] = 3;
  origin.matrix[1][0] = 4;
  origin.matrix[1][1] = 5;
  origin.matrix[1][2] = 6;
  origin.matrix[2][0] = 7;
  origin.matrix[2][1] = 8;
  origin.matrix[2][2] = 9;
  int result = s21_transpose(&m1, &r);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(result, check_origin);
  for (int i = 0; i < r.rows; i++) {
    for (int j = 0; j < r.columns; j++) {
      ck_assert_double_eq_tol(origin.matrix[i][j], r.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&m1);
  s21_remove_matrix(&r);
  s21_remove_matrix(&origin);
}
END_TEST

START_TEST(transpose_3) {
  matrix_t m1;
  matrix_t r;
  matrix_t origin;
  int check_origin = 0;
  int check_1 = s21_create_matrix(1, 1, &m1);
  int check_2 = s21_create_matrix(1, 1, &origin);
  m1.matrix[0][0] = 1.1;
  origin.matrix[0][0] = 1.1;
  int result = s21_transpose(&m1, &r);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(result, check_origin);
  for (int i = 0; i < r.rows; i++) {
    for (int j = 0; j < r.columns; j++) {
      ck_assert_double_eq_tol(origin.matrix[i][j], r.matrix[i][j], 1e-7);
    }
  }
  s21_remove_matrix(&m1);
  s21_remove_matrix(&r);
  s21_remove_matrix(&origin);
}
END_TEST

START_TEST(s21_transpose_4) {
  matrix_t A, result, result_origin;
  s21_create_matrix(5, 4, &A);
  s21_create_matrix(4, 5, &result_origin);
  A.matrix[0][0] = 2;
  A.matrix[0][1] = 5;
  A.matrix[0][2] = 7;
  A.matrix[0][3] = -4;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;
  A.matrix[1][3] = 14;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;
  A.matrix[2][3] = -1;
  A.matrix[3][0] = 10;
  A.matrix[3][1] = 3;
  A.matrix[3][2] = 9;
  A.matrix[3][3] = -3;
  A.matrix[4][0] = -3;
  A.matrix[4][1] = 3;
  A.matrix[4][2] = 1;
  A.matrix[4][3] = 5;
  result_origin.matrix[0][0] = 2;
  result_origin.matrix[0][1] = 6;
  result_origin.matrix[0][2] = 5;
  result_origin.matrix[0][3] = 10;
  result_origin.matrix[0][4] = -3;
  result_origin.matrix[1][0] = 5;
  result_origin.matrix[1][1] = 3;
  result_origin.matrix[1][2] = -2;
  result_origin.matrix[1][3] = 3;
  result_origin.matrix[1][4] = 3;
  result_origin.matrix[2][0] = 7;
  result_origin.matrix[2][1] = 4;
  result_origin.matrix[2][2] = -3;
  result_origin.matrix[2][3] = 9;
  result_origin.matrix[2][4] = 1;
  result_origin.matrix[3][0] = -4;
  result_origin.matrix[3][1] = 14;
  result_origin.matrix[3][2] = -1;
  result_origin.matrix[3][3] = -3;
  result_origin.matrix[3][4] = 5;
  s21_transpose(&A, &result);
  for (int i = 0; i < result_origin.rows; i++) {
    for (int j = 0; j < result_origin.columns; j++) {
      ck_assert_uint_eq(result_origin.matrix[i][j], result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&result_origin);
}
END_TEST

START_TEST(s21_transpose_5) {
  matrix_t A, result, res_origin;
  s21_create_matrix(4, 7, &A);
  s21_create_matrix(7, 4, &res_origin);
  for (int i = 0; i < A.rows; i++) {
    for (int j = 0; j < A.columns; j++) {
      A.matrix[i][j] = i;
    }
  }
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      res_origin.matrix[i][j] = j;
    }
  }
  s21_transpose(&A, &result);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_uint_eq(res_origin.matrix[i][j], result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_transpose_6) {
  matrix_t A, result;
  s21_create_matrix(4, -7, &A);
  int check = s21_transpose(&A, &result);
  int check_origin = 1;
  ck_assert_uint_eq(check, check_origin);
  s21_remove_matrix(&result);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_mult_matrix_1) {
  matrix_t A, B, result, res_origin;
  int check_origin = 0;
  int check_1 = s21_create_matrix(5, 4, &A);
  int check_2 = s21_create_matrix(4, 4, &B);
  int check_4 = s21_create_matrix(5, 4, &res_origin);
  A.matrix[0][0] = 2;
  A.matrix[0][1] = 5;
  A.matrix[0][2] = 7;
  A.matrix[0][3] = -4;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;
  A.matrix[1][3] = 14;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;
  A.matrix[2][3] = -1;
  A.matrix[3][0] = 10;
  A.matrix[3][1] = 3;
  A.matrix[3][2] = 9;
  A.matrix[3][3] = -3;
  A.matrix[4][0] = -3;
  A.matrix[4][1] = 3;
  A.matrix[4][2] = 1;
  A.matrix[4][3] = 5;
  B.matrix[0][0] = 8;
  B.matrix[0][1] = 8;
  B.matrix[0][2] = 3;
  B.matrix[0][3] = 5;
  B.matrix[1][0] = -7;
  B.matrix[1][1] = 7;
  B.matrix[1][2] = 1;
  B.matrix[1][3] = 1;
  B.matrix[2][0] = 7;
  B.matrix[2][1] = 3;
  B.matrix[2][2] = 4;
  B.matrix[2][3] = -4;
  B.matrix[3][0] = 10;
  B.matrix[3][1] = 3;
  B.matrix[3][2] = 9;
  B.matrix[3][3] = -3;
  res_origin.matrix[0][0] = -10;
  res_origin.matrix[0][1] = 60;
  res_origin.matrix[0][2] = 3;
  res_origin.matrix[0][3] = -1;
  res_origin.matrix[1][0] = 195;
  res_origin.matrix[1][1] = 123;
  res_origin.matrix[1][2] = 163;
  res_origin.matrix[1][3] = -25;
  res_origin.matrix[2][0] = 23;
  res_origin.matrix[2][1] = 14;
  res_origin.matrix[2][2] = -8;
  res_origin.matrix[2][3] = 38;
  res_origin.matrix[3][0] = 92;
  res_origin.matrix[3][1] = 119;
  res_origin.matrix[3][2] = 42;
  res_origin.matrix[3][3] = 26;
  res_origin.matrix[4][0] = 12;
  res_origin.matrix[4][1] = 15;
  res_origin.matrix[4][2] = 43;
  res_origin.matrix[4][3] = -31;
  int check = s21_mult_matrix(&A, &B, &result);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(check_4, check_origin);
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_uint_eq(res_origin.matrix[i][j], result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_mult_matrix_2) {
  matrix_t A, B, result, res_origin;
  int check_origin = 0;
  int check_1 = s21_create_matrix(2, 2, &A);
  int check_2 = s21_create_matrix(2, 4, &B);
  int check_4 = s21_create_matrix(2, 4, &res_origin);
  A.matrix[0][0] = 2.987654;
  A.matrix[0][1] = -5.3456;
  A.matrix[1][0] = 6.23;
  A.matrix[1][1] = -3.9292;
  B.matrix[0][0] = 4.838;
  B.matrix[0][1] = 8.0987;
  B.matrix[0][2] = 8.09;
  B.matrix[0][3] = -987;
  B.matrix[1][0] = -8;
  B.matrix[1][1] = -0.921;
  B.matrix[1][2] = -0.91;
  B.matrix[1][3] = 21;
  res_origin.matrix[0][0] = 57.219070052;
  res_origin.matrix[0][1] = 29.1194110498;
  res_origin.matrix[0][2] = 29.03461686;
  res_origin.matrix[0][3] = -3061.072098;
  res_origin.matrix[1][0] = 61.57434;
  res_origin.matrix[1][1] = 54.0736942;
  res_origin.matrix[1][2] = 53.976272;
  res_origin.matrix[1][3] = -6231.5232;
  int check = s21_mult_matrix(&A, &B, &result);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(check_4, check_origin);
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_uint_eq(res_origin.matrix[i][j], result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_mult_matrix_3) {
  matrix_t A, B, result, res_origin;
  int check_origin = 0;
  int check_1 = s21_create_matrix(3, 3, &A);
  int check_2 = s21_create_matrix(3, 3, &B);
  int check_4 = s21_create_matrix(3, 3, &res_origin);
  A.matrix[0][0] = 0;
  A.matrix[0][1] = 9;
  A.matrix[0][2] = 1;
  A.matrix[1][0] = 1;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 7;
  A.matrix[2][0] = 2;
  A.matrix[2][1] = 7;
  A.matrix[2][2] = 0;
  B.matrix[0][0] = 0;
  B.matrix[0][1] = 9;
  B.matrix[0][2] = 8;
  B.matrix[1][0] = 0;
  B.matrix[1][1] = 9;
  B.matrix[1][2] = 9;
  B.matrix[2][0] = 1;
  B.matrix[2][1] = 9;
  B.matrix[2][2] = 6;
  res_origin.matrix[0][0] = 1;
  res_origin.matrix[0][1] = 90;
  res_origin.matrix[0][2] = 87;
  res_origin.matrix[1][0] = 7;
  res_origin.matrix[1][1] = 99;
  res_origin.matrix[1][2] = 77;
  res_origin.matrix[2][0] = 0;
  res_origin.matrix[2][1] = 81;
  res_origin.matrix[2][2] = 79;
  int check = s21_mult_matrix(&A, &B, &result);
  ck_assert_int_eq(check_1, check_origin);
  ck_assert_int_eq(check_2, check_origin);
  ck_assert_int_eq(check_4, check_origin);
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_uint_eq(res_origin.matrix[i][j], result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_mult_matrix_4) {
  matrix_t A, B, result;
  int check_origin = 2;
  s21_create_matrix(2, 3, &A);
  s21_create_matrix(4, 2, &B);
  int check = s21_mult_matrix(&A, &B, &result);
  ck_assert_int_eq(check, check_origin);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(s21_mult_matrix_5) {
  matrix_t A, B, result;
  int check_origin = 0;
  s21_create_matrix(2, 3, &A);
  s21_create_matrix(3, 12, &B);
  for (int i = 0; i < B.rows; i++) {
    for (int j = 0; j < B.columns; j++) {
      B.matrix[i][j] = i + j;
    }
  }
  int check = s21_mult_matrix(&A, &B, &result);
  ck_assert_int_eq(check, check_origin);
  for (int i = 0; i < result.rows; i++) {
    for (int j = 0; j < result.columns; j++) {
      ck_assert_int_eq(result.matrix[i][j], 0.0);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

// START_TEST(s21_mult_matrix_6) {
//   matrix_t A, B, result;
//   int check_origin = 1;
//   s21_create_matrix(-12, 3, &A);
//   s21_create_matrix(3, 12, &B);
//   int check = s21_mult_matrix(&A, &B, &result);
//   ck_assert_int_eq(check, check_origin);
//   s21_remove_matrix(&B);
//   s21_remove_matrix(&A);
//   s21_remove_matrix(&result);
// }
// END_TEST

START_TEST(s21_determinant_1) {
  matrix_t A;
  double result;
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 5;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = -5;
  A.matrix[2][1] = 3;
  A.matrix[2][2] = 3;
  int check = s21_determinant(&A, &result);
  int check_origin = 0;
  double result_origin = 43;
  ck_assert_uint_eq(check, check_origin);
  ck_assert_uint_eq(result, result_origin);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_2) {
  matrix_t A;
  double result;
  s21_create_matrix(4, 4, &A);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[0][3] = 3;
  A.matrix[1][0] = 5;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 2;
  A.matrix[1][3] = 2;
  A.matrix[2][0] = -5;
  A.matrix[2][1] = 3;
  A.matrix[2][2] = 3;
  A.matrix[2][3] = 3;
  A.matrix[3][0] = -12.12;
  A.matrix[3][1] = 12.12;
  A.matrix[3][2] = 2.12;
  A.matrix[3][3] = -2.12;
  int check = s21_determinant(&A, &result);
  int check_origin = 0;
  double result_origin = -182.32;
  ck_assert_uint_eq(check, check_origin);
  ck_assert_uint_eq(result, result_origin);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_3) {
  matrix_t A;
  double result;
  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = 1;
  int check = s21_determinant(&A, &result);
  int check_origin = 0;
  double result_origin = 1;
  ck_assert_uint_eq(check, check_origin);
  ck_assert_uint_eq(result, result_origin);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_4) {
  matrix_t A;
  double result;
  s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 55.5;
  A.matrix[0][1] = 14.4;
  A.matrix[1][0] = -0.752;
  A.matrix[1][1] = 156.2;
  int check = s21_determinant(&A, &result);
  int check_origin = 0;
  double result_origin = 8679.9288;
  ck_assert_uint_eq(check, check_origin);
  ck_assert_uint_eq(result, result_origin);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_5) {
  matrix_t A;
  double result;
  s21_create_matrix(4, 4, &A);
  A.matrix[0][0] = 11.1;
  A.matrix[0][1] = 2.2;
  A.matrix[0][2] = 3.3;
  A.matrix[0][3] = 4.4;
  A.matrix[1][0] = 8.8;
  A.matrix[1][1] = -9.9;
  A.matrix[1][2] = 10.1;
  A.matrix[1][3] = 11.11;
  A.matrix[2][0] = 15.15;
  A.matrix[2][1] = 16.16;
  A.matrix[2][2] = 17.17;
  A.matrix[2][3] = 18.18;
  A.matrix[3][0] = -22.22;
  A.matrix[3][1] = 23.23;
  A.matrix[3][2] = 24.24;
  A.matrix[3][3] = 25.25;
  int check = s21_determinant(&A, &result);
  int check_origin = 0;
  double result_origin = -11770.2443124;
  ck_assert_uint_eq(check, check_origin);
  ck_assert_double_eq_tol(result, result_origin, 1e-7);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_6) {
  matrix_t A;
  double result;
  s21_create_matrix(4, 4, &A);
  int check = s21_determinant(&A, &result);
  int check_origin = 0;
  double result_origin = 0;
  ck_assert_uint_eq(check, check_origin);
  ck_assert_double_eq_tol(result, result_origin, 1e-7);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_7) {
  matrix_t A;
  double result;
  s21_create_matrix(5, 4, &A);
  int check = s21_determinant(&A, &result);
  int check_origin = 2;
  ck_assert_uint_eq(check, check_origin);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_8) {
  matrix_t A;
  double result;
  s21_create_matrix(-4, 4, &A);
  int check = s21_determinant(&A, &result);
  int check_origin = 1;
  ck_assert_uint_eq(check, check_origin);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_9) {
  matrix_t A;
  double result;
  s21_create_matrix(4, 4, &A);
  A.matrix[0][0] = 0;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[0][3] = 5;
  A.matrix[1][0] = 3;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 0;
  A.matrix[1][3] = 4;
  A.matrix[2][0] = 8;
  A.matrix[2][1] = 6;
  A.matrix[2][2] = 0;
  A.matrix[2][3] = 4;
  A.matrix[3][0] = 6;
  A.matrix[3][1] = 8;
  A.matrix[3][2] = 9;
  A.matrix[3][3] = 4;
  int check = s21_determinant(&A, &result);
  ck_assert_int_eq(check, 0);
  ck_assert_double_eq_tol(result, 438.0, 1e-7);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_10) {
  matrix_t A;
  double result;
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 0;
  A.matrix[1][0] = 0;
  A.matrix[2][0] = 0;
  A.matrix[0][1] = 3;
  A.matrix[1][1] = 5;
  A.matrix[2][1] = 9;
  A.matrix[0][2] = 1;
  A.matrix[1][2] = 56.23456;
  A.matrix[2][2] = 12345432.343;
  int check = s21_determinant(&A, &result);
  ck_assert_int_eq(check, 0);
  ck_assert_double_eq_tol(result, 0.0, 1e-7);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_11) {
  matrix_t A;
  double result;
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 123.456;
  A.matrix[1][0] = 123.455;
  A.matrix[2][0] = 123.454;
  A.matrix[0][1] = 123.453;
  A.matrix[1][1] = 123.452;
  A.matrix[2][1] = 123.451;
  A.matrix[0][2] = 123.450;
  A.matrix[1][2] = 1113.449;
  A.matrix[2][2] = 15.448;
  int check = s21_determinant(&A, &result);
  ck_assert_int_eq(check, 0);
  ck_assert_double_eq_tol(result, 0.006264, 1e-7 * 10.0);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_calc_complements_1) {
  matrix_t A, result, res_origin;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &res_origin);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = 2;
  A.matrix[2][2] = 1;
  res_origin.matrix[0][0] = 0;
  res_origin.matrix[0][1] = 10;
  res_origin.matrix[0][2] = -20;
  res_origin.matrix[1][0] = 4;
  res_origin.matrix[1][1] = -14;
  res_origin.matrix[1][2] = 8;
  res_origin.matrix[2][0] = -8;
  res_origin.matrix[2][1] = -2;
  res_origin.matrix[2][2] = 4;
  s21_calc_complements(&A, &result);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_uint_eq(res_origin.matrix[i][j], result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_calc_complements_2) {
  matrix_t A, result, res_origin;
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &res_origin);
  A.matrix[0][0] = 2;
  A.matrix[0][1] = 4;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  res_origin.matrix[0][0] = 3;
  res_origin.matrix[0][1] = -6;
  res_origin.matrix[1][0] = -4;
  res_origin.matrix[1][1] = 2;
  s21_calc_complements(&A, &result);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_uint_eq(res_origin.matrix[i][j], result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_calc_complements_3) {
  matrix_t A, result, res_origin;
  s21_create_matrix(-2, 2, &A);
  s21_create_matrix(2, 2, &res_origin);
  int check = s21_calc_complements(&A, &result);
  int check_origin = 1;
  ck_assert_uint_eq(check, check_origin);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_calc_complements_4) {
  matrix_t A, result, res_origin;
  s21_create_matrix(2, 4, &A);
  s21_create_matrix(2, 2, &res_origin);
  int check = s21_calc_complements(&A, &result);
  int check_origin = 2;
  ck_assert_uint_eq(check, check_origin);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_inverse_matrix_1) {
  matrix_t A, result, res_origin;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &res_origin);
  A.matrix[0][0] = 2;
  A.matrix[0][1] = 5;
  A.matrix[0][2] = 7;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;
  res_origin.matrix[0][0] = 1;
  res_origin.matrix[0][1] = -1;
  res_origin.matrix[0][2] = 1;
  res_origin.matrix[1][0] = -38;
  res_origin.matrix[1][1] = 41;
  res_origin.matrix[1][2] = -34;
  res_origin.matrix[2][0] = 27;
  res_origin.matrix[2][1] = -29;
  res_origin.matrix[2][2] = 24;
  s21_inverse_matrix(&A, &result);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_double_eq_tol(result.matrix[i][j], res_origin.matrix[i][j],
                              1e-7);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_inverse_matrix_2) {
  matrix_t A, result, res_origin;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &res_origin);
  A.matrix[0][0] = 2;
  A.matrix[0][1] = 1.3;
  A.matrix[0][2] = 0.7;
  A.matrix[1][0] = -1;
  A.matrix[1][1] = 0;
  A.matrix[1][2] = 3;
  A.matrix[2][0] = 0.3;
  A.matrix[2][1] = -0.3;
  A.matrix[2][2] = 2;
  res_origin.matrix[0][0] = 45.0 / 289.0;
  res_origin.matrix[0][1] = -281.0 / 578.0;
  res_origin.matrix[0][2] = 195.0 / 289.0;
  res_origin.matrix[1][0] = 145.0 / 289.0;
  res_origin.matrix[1][1] = 379.0 / 578.0;
  res_origin.matrix[1][2] = -335.0 / 289.0;
  res_origin.matrix[2][0] = 15.0 / 289.0;
  res_origin.matrix[2][1] = 99.0 / 578.0;
  res_origin.matrix[2][2] = 65.0 / 289.0;
  s21_inverse_matrix(&A, &result);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_double_eq_tol(result.matrix[i][j], res_origin.matrix[i][j],
                              1e-7);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_inverse_matrix_3) {
  matrix_t A, result, res_origin;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &res_origin);
  A.matrix[0][0] = 2.8;
  A.matrix[0][1] = 1.3;
  A.matrix[0][2] = 7.01;
  A.matrix[1][0] = -1.03;
  A.matrix[1][1] = -2.3;
  A.matrix[1][2] = 3.01;
  A.matrix[2][0] = 0;
  A.matrix[2][1] = -3;
  A.matrix[2][2] = 2;
  res_origin.matrix[0][0] = 44300.0 / 367429.0;
  res_origin.matrix[0][1] = -236300.0 / 367429.0;
  res_origin.matrix[0][2] = 200360.0 / 367429.0;
  res_origin.matrix[1][0] = 20600.0 / 367429.0;
  res_origin.matrix[1][1] = 56000.0 / 367429.0;
  res_origin.matrix[1][2] = -156483.0 / 367429.0;
  res_origin.matrix[2][0] = 30900.0 / 367429.0;
  res_origin.matrix[2][1] = 84000.0 / 367429.0;
  res_origin.matrix[2][2] = -51010.0 / 367429.0;
  s21_inverse_matrix(&A, &result);
  for (int i = 0; i < res_origin.rows; i++) {
    for (int j = 0; j < res_origin.columns; j++) {
      ck_assert_double_eq_tol(result.matrix[i][j], res_origin.matrix[i][j],
                              1e-7);
    }
  }
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_inverse_matrix_4) {
  matrix_t A, result, res_origin;
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &res_origin);
  A.matrix[0][0] = 777;
  res_origin.matrix[0][0] = 1.0 / 777.0;
  int check = s21_inverse_matrix(&A, &result);
  int check_origin = 0;
  ck_assert_double_eq_tol(res_origin.matrix[0][0], result.matrix[0][0], 1e-07);
  ck_assert_uint_eq(check, check_origin);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&res_origin);
}
END_TEST

START_TEST(s21_inverse_matrix_5) {
  matrix_t A, result;
  s21_create_matrix(-5, 5, &A);
  int check = s21_inverse_matrix(&A, &result);
  int check_origin = 1;
  ck_assert_uint_eq(check, check_origin);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(s21_inverse_matrix_6) {
  matrix_t A, result;
  s21_create_matrix(1, 2, &A);
  A.matrix[0][0] = 777;
  A.matrix[0][1] = 777;
  int check = s21_inverse_matrix(&A, &result);
  int check_origin = 2;
  ck_assert_uint_eq(check, check_origin);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(s21_inverse_matrix_7) {
  matrix_t A, result;
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 0;
  A.matrix[0][1] = 0;
  A.matrix[0][2] = 0;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 0;
  A.matrix[1][2] = 0;
  A.matrix[2][0] = 0;
  A.matrix[2][1] = 0;
  A.matrix[2][2] = 0;
  int check = s21_inverse_matrix(&A, &result);
  int check_origin = 2;
  ck_assert_uint_eq(check, check_origin);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(s21_remove_matrix_1) {
  matrix_t A;
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 75301501;
  s21_remove_matrix(&A);
  int check = (A.rows == 3) ? 0 : 1;
  int check_origin = 1;
  ck_assert_uint_eq(check, check_origin);
}
END_TEST

int main(void) {
  Suite *s1 = suite_create("Core");
  TCase *tc1_1 = tcase_create("Core");
  SRunner *sr = srunner_create(s1);
  int nf;

  suite_add_tcase(s1, tc1_1);
  tcase_add_test(tc1_1, create_1);
  tcase_add_test(tc1_1, create_2);
  tcase_add_test(tc1_1, create_3);
  tcase_add_test(tc1_1, create_4);
  tcase_add_test(tc1_1, create_5);
  tcase_add_test(tc1_1, eq_1);
  tcase_add_test(tc1_1, eq_2);
  tcase_add_test(tc1_1, eq_3);
  tcase_add_test(tc1_1, eq_4);
  tcase_add_test(tc1_1, s21_eq_matrix_1);
  tcase_add_test(tc1_1, s21_eq_matrix_2);
  tcase_add_test(tc1_1, s21_eq_matrix_3);
  tcase_add_test(tc1_1, s21_eq_matrix_4);
  tcase_add_test(tc1_1, s21_eq_matrix_6);
  tcase_add_test(tc1_1, s21_eq_matrix_7);
  tcase_add_test(tc1_1, sum_1);
  tcase_add_test(tc1_1, sum_2);
  tcase_add_test(tc1_1, sum_3);
  tcase_add_test(tc1_1, s21_sum_matrix_1);
  tcase_add_test(tc1_1, s21_sum_matrix_2);
  tcase_add_test(tc1_1, s21_sum_matrix_3);
  tcase_add_test(tc1_1, s21_sum_matrix_4);
  tcase_add_test(tc1_1, s21_sum_matrix_5);
  tcase_add_test(tc1_1, s21_sum_matrix_7);
  tcase_add_test(tc1_1, s21_sum_matrix_8);
  tcase_add_test(tc1_1, s21_sum_matrix_9);
  tcase_add_test(tc1_1, sub_1);
  tcase_add_test(tc1_1, sub_2);
  tcase_add_test(tc1_1, sub_3);
  tcase_add_test(tc1_1, sub_4);
  tcase_add_test(tc1_1, s21_sub_matrix_1);
  tcase_add_test(tc1_1, s21_sub_matrix_2);
  tcase_add_test(tc1_1, s21_sub_matrix_3);
  tcase_add_test(tc1_1, s21_sub_matrix_9);
  tcase_add_test(tc1_1, s21_sub_matrix_10);
  tcase_add_test(tc1_1, s21_sub_matrix_4);
  tcase_add_test(tc1_1, s21_sub_matrix_5);
  tcase_add_test(tc1_1, s21_sub_matrix_7);
  tcase_add_test(tc1_1, s21_sub_matrix_8);
  tcase_add_test(tc1_1, mul_num_1);
  tcase_add_test(tc1_1, mul_num_2);
  tcase_add_test(tc1_1, mul_num_3);
  tcase_add_test(tc1_1, s21_mult_number_1);
  tcase_add_test(tc1_1, s21_mult_number_2);
  tcase_add_test(tc1_1, s21_mult_number_3);
  tcase_add_test(tc1_1, transpose_1);
  tcase_add_test(tc1_1, transpose_2);
  tcase_add_test(tc1_1, transpose_3);
  tcase_add_test(tc1_1, s21_transpose_4);
  tcase_add_test(tc1_1, s21_transpose_5);
  tcase_add_test(tc1_1, s21_transpose_6);
  tcase_add_test(tc1_1, s21_mult_matrix_1);
  tcase_add_test(tc1_1, s21_mult_matrix_2);
  tcase_add_test(tc1_1, s21_mult_matrix_3);
  tcase_add_test(tc1_1, s21_mult_matrix_4);
  tcase_add_test(tc1_1, s21_mult_matrix_5);
  tcase_add_test(tc1_1, s21_determinant_1);
  tcase_add_test(tc1_1, s21_determinant_2);
  tcase_add_test(tc1_1, s21_determinant_3);
  tcase_add_test(tc1_1, s21_determinant_4);
  tcase_add_test(tc1_1, s21_determinant_5);
  tcase_add_test(tc1_1, s21_determinant_6);
  tcase_add_test(tc1_1, s21_determinant_7);
  tcase_add_test(tc1_1, s21_determinant_8);
  tcase_add_test(tc1_1, s21_determinant_9);
  tcase_add_test(tc1_1, s21_determinant_10);
  tcase_add_test(tc1_1, s21_determinant_11);
  tcase_add_test(tc1_1, s21_calc_complements_1);
  tcase_add_test(tc1_1, s21_calc_complements_2);
  tcase_add_test(tc1_1, s21_calc_complements_3);
  tcase_add_test(tc1_1, s21_calc_complements_4);
  tcase_add_test(tc1_1, s21_inverse_matrix_1);
  tcase_add_test(tc1_1, s21_inverse_matrix_2);
  tcase_add_test(tc1_1, s21_inverse_matrix_3);
  tcase_add_test(tc1_1, s21_inverse_matrix_4);
  tcase_add_test(tc1_1, s21_inverse_matrix_5);
  tcase_add_test(tc1_1, s21_inverse_matrix_6);
  tcase_add_test(tc1_1, s21_inverse_matrix_7);
  tcase_add_test(tc1_1, s21_remove_matrix_1);

  // srunner_set_fork_status(sr, CK_NOFORK);
  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}
