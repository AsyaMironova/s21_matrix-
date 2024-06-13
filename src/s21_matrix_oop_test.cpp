#include "s21_matrix_oop_test.h"

void S21MatrixTest::SetUp() {
  matrix_1x1 = new S21Matrix(1, 1);
  (*matrix_1x1)(0, 0) = -2;

  matrix_3x2 = new S21Matrix(3, 2);
  (*matrix_3x2)(0, 0) = 1;
  (*matrix_3x2)(0, 1) = -2;
  (*matrix_3x2)(1, 0) = -3;
  (*matrix_3x2)(1, 1) = 4;
  (*matrix_3x2)(2, 0) = 5;
  (*matrix_3x2)(2, 1) = -6;

  matrix_3x3 = new S21Matrix(3, 3);
  (*matrix_3x3)(0, 0) = 2;
  (*matrix_3x3)(0, 1) = 5;
  (*matrix_3x3)(0, 2) = 7;
  (*matrix_3x3)(1, 0) = 6;
  (*matrix_3x3)(1, 1) = 3;
  (*matrix_3x3)(1, 2) = 4;
  (*matrix_3x3)(2, 0) = 5;
  (*matrix_3x3)(2, 1) = -2;
  (*matrix_3x3)(2, 2) = -3;

  matrix_5x5 = new S21Matrix(5, 5);
  FillRandom(*matrix_5x5);
}

void S21MatrixTest::TearDown() {
  delete matrix_1x1;
  delete matrix_3x2;
  delete matrix_3x3;
  delete matrix_5x5;
}

void S21MatrixTest::FillRandom(S21Matrix &matrix) {
  srand(time(nullptr));
  int rows = matrix.GetRows();
  int cols = matrix.GetCols();
  for (int row = 0; row < rows; row++)
    for (int col = 0; col < cols; col++) matrix(row, col) = rand() % 10;
}

// Constructors and destructor
TEST_F(S21MatrixTest, ConstructorNoParameters) {
  S21Matrix matrix;
  EXPECT_EQ(3, matrix.GetRows());
  EXPECT_EQ(3, matrix.GetCols());
}

TEST_F(S21MatrixTest, ConstructorWithParameters1) {
  S21Matrix matrix(3, 1);
  EXPECT_EQ(3, matrix.GetRows());
  EXPECT_EQ(1, matrix.GetCols());
}

TEST_F(S21MatrixTest, ConstructorWithParameters2) {
  EXPECT_THROW(S21Matrix matrix(0, 3), std::invalid_argument);
}

TEST_F(S21MatrixTest, ConstructorWithParameters3) {
  EXPECT_THROW(S21Matrix matrix(0, -1), std::invalid_argument);
}

TEST_F(S21MatrixTest, ConstructorCopy) {
  S21Matrix matrix(*matrix_5x5);
  EXPECT_EQ(5, matrix.GetRows());
  EXPECT_EQ(5, matrix.GetCols());
  EXPECT_EQ(true, matrix == *matrix_5x5);
}

TEST_F(S21MatrixTest, ConstructorMove) {
  S21Matrix src(*matrix_5x5);
  S21Matrix dest(std::move(src));
  EXPECT_EQ(0, src.GetRows());
  EXPECT_EQ(0, src.GetCols());
  EXPECT_EQ(5, dest.GetRows());
  EXPECT_EQ(5, dest.GetCols());
  EXPECT_EQ(true, dest == *matrix_5x5);
}

// Getters and setters
TEST_F(S21MatrixTest, GetRows) {
  S21Matrix matrix;
  EXPECT_EQ(3, matrix.GetRows());
}

TEST_F(S21MatrixTest, GetCols) {
  S21Matrix matrix;
  EXPECT_EQ(3, matrix.GetCols());
}

TEST_F(S21MatrixTest, SetRows1) {
  S21Matrix matrix;
  matrix.SetRows(5);
  EXPECT_EQ(5, matrix.GetRows());
}

TEST_F(S21MatrixTest, SetRows2) {
  S21Matrix matrix;
  EXPECT_THROW(matrix.SetRows(-1), std::out_of_range);
}

TEST_F(S21MatrixTest, SetCols1) {
  S21Matrix matrix;
  matrix.SetCols(5);
  EXPECT_EQ(5, matrix.GetCols());
}

TEST_F(S21MatrixTest, SetCols2) {
  S21Matrix matrix;
  EXPECT_THROW(matrix.SetCols(-1), std::out_of_range);
}

// Operators
TEST_F(S21MatrixTest, OperatorSum1) {
  S21Matrix matrixA(*matrix_3x3), expected;
  expected(0, 0) = 4;
  expected(0, 1) = 10;
  expected(0, 2) = 14;
  expected(1, 0) = 12;
  expected(1, 1) = 6;
  expected(1, 2) = 8;
  expected(2, 0) = 10;
  expected(2, 1) = -4;
  expected(2, 2) = -6;
  EXPECT_EQ(true, (matrixA + matrixA) == matrixA * 2);
}

TEST_F(S21MatrixTest, OperatorSum2) {
  S21Matrix matrixA(*matrix_3x3), matrixB(*matrix_3x2);
  EXPECT_THROW(matrixA + matrixB, std::out_of_range);
}

TEST_F(S21MatrixTest, OperatorSub1) {
  S21Matrix matrixA(*matrix_3x3), expected;
  EXPECT_EQ(true, (matrixA - matrixA) == expected);
}

TEST_F(S21MatrixTest, OperatorSub2) {
  S21Matrix matrixA(*matrix_3x3), matrixB(*matrix_3x2);
  EXPECT_THROW(matrixA - matrixB, std::out_of_range);
}

TEST_F(S21MatrixTest, OperatorMulMatrix1) {
  S21Matrix matrixA(*matrix_3x3), expected;
  expected(0, 0) = 69;
  expected(0, 1) = 11;
  expected(0, 2) = 13;
  expected(1, 0) = 50;
  expected(1, 1) = 31;
  expected(1, 2) = 42;
  expected(2, 0) = -17;
  expected(2, 1) = 25;
  expected(2, 2) = 36;
  // std::cout << matrixA * matrixA << '\n' << expected;
  EXPECT_EQ(true, (matrixA * matrixA) == expected);
}

TEST_F(S21MatrixTest, OperatorMulMatrix2) {
  S21Matrix matrixA(*matrix_3x3), matrixB(*matrix_3x2);
  EXPECT_THROW(matrixB * matrixA, std::out_of_range);
}

TEST_F(S21MatrixTest, OperatorMulNumber1) {
  S21Matrix matrixA(*matrix_3x3), expected;
  expected(0, 0) = -4;
  expected(0, 1) = -10;
  expected(0, 2) = -14;
  expected(1, 0) = -12;
  expected(1, 1) = -6;
  expected(1, 2) = -8;
  expected(2, 0) = -10;
  expected(2, 1) = 4;
  expected(2, 2) = 6;
  EXPECT_EQ(true, (matrixA * (-2)) == expected);
}

TEST_F(S21MatrixTest, OperatorMulNumber2) {
  S21Matrix matrixA(*matrix_3x3), expected;
  expected(0, 0) = -4;
  expected(0, 1) = -10;
  expected(0, 2) = -14;
  expected(1, 0) = -12;
  expected(1, 1) = -6;
  expected(1, 2) = -8;
  expected(2, 0) = -10;
  expected(2, 1) = 4;
  expected(2, 2) = 6;
  EXPECT_EQ(true, ((-2.) * matrixA) == expected);
}

TEST_F(S21MatrixTest, OperatorAssign) {
  S21Matrix matrixA(*matrix_3x3), matrixB;
  matrixB = matrixA;
  EXPECT_EQ(true, matrixA == matrixB);
}

TEST_F(S21MatrixTest, OperatorAssignSum) {
  S21Matrix matrixA(*matrix_3x3), expected;
  expected(0, 0) = 4;
  expected(0, 1) = 10;
  expected(0, 2) = 14;
  expected(1, 0) = 12;
  expected(1, 1) = 6;
  expected(1, 2) = 8;
  expected(2, 0) = 10;
  expected(2, 1) = -4;
  expected(2, 2) = -6;
  matrixA += matrixA;
  EXPECT_EQ(true, matrixA == expected);
}

TEST_F(S21MatrixTest, OperatorAssignSub) {
  S21Matrix matrixA(*matrix_3x3), expected;
  S21Matrix matrixB(*matrix_3x3);
  matrixA -= matrixB;
  EXPECT_EQ(true, matrixA == expected);
}

TEST_F(S21MatrixTest, OperatorAssignMulMatrix) {
  S21Matrix matrixA(*matrix_3x3), expected;
  expected(0, 0) = 69;
  expected(0, 1) = 11;
  expected(0, 2) = 13;
  expected(1, 0) = 50;
  expected(1, 1) = 31;
  expected(1, 2) = 42;
  expected(2, 0) = -17;
  expected(2, 1) = 25;
  expected(2, 2) = 36;
  matrixA *= matrixA;
  EXPECT_EQ(true, matrixA == expected);
}

TEST_F(S21MatrixTest, OperatorAssignMulNumber) {
  S21Matrix matrixA(*matrix_3x3), expected;
  expected(0, 0) = -4;
  expected(0, 1) = -10;
  expected(0, 2) = -14;
  expected(1, 0) = -12;
  expected(1, 1) = -6;
  expected(1, 2) = -8;
  expected(2, 0) = -10;
  expected(2, 1) = 4;
  expected(2, 2) = 6;
  matrixA *= -2;
  EXPECT_EQ(true, matrixA == expected);
}

TEST_F(S21MatrixTest, OperatorAssignIndexes1) {
  S21Matrix matrixA;
  EXPECT_DOUBLE_EQ(0, matrixA(1, 1));
  matrixA(1, 1) = -4;
  EXPECT_DOUBLE_EQ(-4, matrixA(1, 1));
}

TEST_F(S21MatrixTest, OperatorAssignIndexes2) {
  S21Matrix matrixA;
  EXPECT_THROW(matrixA(-1, 1), std::out_of_range);
}

TEST_F(S21MatrixTest, OperatorAssignIndexes3) {
  S21Matrix matrixA;
  EXPECT_THROW(matrixA(4, 1), std::out_of_range);
}

TEST_F(S21MatrixTest, OperatorAssignIndexes4) {
  S21Matrix matrixA;
  EXPECT_THROW(matrixA(1, -1), std::out_of_range);
}

TEST_F(S21MatrixTest, OperatorAssignIndexes5) {
  S21Matrix matrixA;
  EXPECT_THROW(matrixA(1, 4), std::out_of_range);
}

// Methods
TEST_F(S21MatrixTest, EqMatrix1) {
  S21Matrix matrixA, matrixB;
  EXPECT_EQ(true, matrixA.EqMatrix(matrixB));
}

TEST_F(S21MatrixTest, EqMatrix2) {
  S21Matrix matrixA, matrixB;
  matrixB(1, 1) = 1;
  EXPECT_EQ(false, matrixA.EqMatrix(matrixB));
}

TEST_F(S21MatrixTest, EqMatrix3) {
  S21Matrix matrixA, matrixB;
  matrixB.SetCols(2);
  EXPECT_EQ(false, matrixA.EqMatrix(matrixB));
}

TEST_F(S21MatrixTest, EqMatrix4) {
  S21Matrix matrixA, matrixB;
  matrixB.SetRows(2);
  EXPECT_EQ(false, matrixA.EqMatrix(matrixB));
}

TEST_F(S21MatrixTest, SumMatrix1) {
  S21Matrix matrixA(*matrix_3x3), expected;
  expected(0, 0) = 4;
  expected(0, 1) = 10;
  expected(0, 2) = 14;
  expected(1, 0) = 12;
  expected(1, 1) = 6;
  expected(1, 2) = 8;
  expected(2, 0) = 10;
  expected(2, 1) = -4;
  expected(2, 2) = -6;
  matrixA.SumMatrix(matrixA);
  EXPECT_EQ(true, matrixA == expected);
}

TEST_F(S21MatrixTest, SumMatrix2) {
  S21Matrix matrixA(*matrix_3x3), matrixB(*matrix_3x2);
  EXPECT_THROW(matrixA.SumMatrix(matrixB), std::out_of_range);
}

TEST_F(S21MatrixTest, SubMatrix1) {
  S21Matrix matrixA(*matrix_3x3), expected;
  matrixA.SubMatrix(matrixA);
  EXPECT_EQ(true, matrixA == expected);
}

TEST_F(S21MatrixTest, SubMatrix2) {
  S21Matrix matrixA(*matrix_3x3), matrixB(*matrix_3x2);
  EXPECT_THROW(matrixA.SubMatrix(matrixB), std::out_of_range);
}

TEST_F(S21MatrixTest, MulMatrix1) {
  S21Matrix matrixA(*matrix_3x3), expected;
  expected(0, 0) = 69;
  expected(0, 1) = 11;
  expected(0, 2) = 13;
  expected(1, 0) = 50;
  expected(1, 1) = 31;
  expected(1, 2) = 42;
  expected(2, 0) = -17;
  expected(2, 1) = 25;
  expected(2, 2) = 36;
  matrixA.MulMatrix(matrixA);
  EXPECT_EQ(true, matrixA == expected);
}

TEST_F(S21MatrixTest, MulMatrix2) {
  S21Matrix matrixA(*matrix_3x3), matrixB(*matrix_3x2);
  EXPECT_THROW(matrixB.MulMatrix(matrixA), std::out_of_range);
}

TEST_F(S21MatrixTest, MulNumber) {
  S21Matrix matrixA(*matrix_3x3), expected;
  expected(0, 0) = -4;
  expected(0, 1) = -10;
  expected(0, 2) = -14;
  expected(1, 0) = -12;
  expected(1, 1) = -6;
  expected(1, 2) = -8;
  expected(2, 0) = -10;
  expected(2, 1) = 4;
  expected(2, 2) = 6;
  matrixA.MulNumber(-2);
  EXPECT_EQ(true, matrixA == expected);
}

TEST_F(S21MatrixTest, Transpose) {
  S21Matrix matrixA(1, 3), expected(3, 1);
  EXPECT_EQ(true, matrixA.Transpose() == expected);
}

TEST_F(S21MatrixTest, Determinant1) {
  S21Matrix matrixA;
  EXPECT_DOUBLE_EQ(0, matrixA.Determinant());
}

TEST_F(S21MatrixTest, Determinant2) {
  S21Matrix matrixA(*matrix_3x3);
  EXPECT_DOUBLE_EQ(-1, matrixA.Determinant());
}

TEST_F(S21MatrixTest, Determinant3) {
  S21Matrix matrixA(*matrix_3x2);
  EXPECT_THROW(matrixA.Determinant(), std::out_of_range);
}

TEST_F(S21MatrixTest, CalcComplements1) {
  S21Matrix matrixA(*matrix_3x3), expected;
  expected(0, 0) = -1;
  expected(0, 1) = 38;
  expected(0, 2) = -27;
  expected(1, 0) = 1;
  expected(1, 1) = -41;
  expected(1, 2) = 29;
  expected(2, 0) = -1;
  expected(2, 1) = 34;
  expected(2, 2) = -24;
  EXPECT_EQ(true, matrixA.CalcComplements() == expected);
}

TEST_F(S21MatrixTest, CalcComplements2) {
  S21Matrix matrixA(*matrix_3x2);
  EXPECT_THROW(matrixA.CalcComplements(), std::out_of_range);
}

TEST_F(S21MatrixTest, Inverse1) {
  S21Matrix matrixA(*matrix_3x3), expected;
  expected(0, 0) = 1;
  expected(0, 1) = -1;
  expected(0, 2) = 1;
  expected(1, 0) = -38;
  expected(1, 1) = 41;
  expected(1, 2) = -34;
  expected(2, 0) = 27;
  expected(2, 1) = -29;
  expected(2, 2) = 24;
  EXPECT_EQ(true, matrixA.InverseMatrix() == expected);
}

TEST_F(S21MatrixTest, Inverse2) {
  S21Matrix matrixA(*matrix_3x2);
  EXPECT_THROW(matrixA.InverseMatrix(), std::out_of_range);
}

TEST_F(S21MatrixTest, Inverse3) {
  S21Matrix matrixA;
  EXPECT_THROW(matrixA.InverseMatrix(), std::out_of_range);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}