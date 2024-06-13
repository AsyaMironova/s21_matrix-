#ifndef __S21_MATRIX_TEST_H__
#define __S21_MATRIX_TEST_H__

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>

#include "s21_matrix_oop.h"

class S21MatrixTest : public testing::Test {
 protected:
  S21Matrix *matrix_1x1;
  S21Matrix *matrix_3x2;
  S21Matrix *matrix_3x3;
  S21Matrix *matrix_5x5;

  void SetUp();
  void TearDown();
  void FillRandom(S21Matrix &matrix);
};

#endif