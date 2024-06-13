#include "s21_matrix_oop.h"

// база
S21Matrix::S21Matrix() : S21Matrix(3, 3) {}

// параметры
S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows < 1 || cols < 1)
    throw std::invalid_argument(
        "CreationError: Matrix sizes should be equal 1 or greater");  // exeption
  else {
    CreateMatrix(rows_, cols_);
  }
}

// копирование
S21Matrix::S21Matrix(const S21Matrix& other)
    : S21Matrix(other.rows_, other.cols_) {
  for (int i = 0; i < other.rows_; ++i) {
    for (int k = 0; k < other.cols_; k++) {
      matrix_[i][k] = other.matrix_[i][k];
    }
  }
}

// движение
S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

// деструктор
S21Matrix::~S21Matrix() { ClearMatrix(); }

// получить количество строк
int S21Matrix::GetRows() const { return rows_; }

// получить количество столбцов
int S21Matrix::GetCols() const { return cols_; }

// Сеттер для установления количество строк
void S21Matrix::SetCols(int cols) {
  if (cols < 1)
    throw std::out_of_range("Cols of matrix should be more or equal 1");
  if (cols != cols_) {
    S21Matrix result(rows_, cols);
    for (int i = 0; i < rows_; ++i)
      for (int k = 0; k < (cols_ < cols ? cols_ : cols); ++k)
        result.matrix_[i][k] = matrix_[i][k];
    *this = std::move(result);
  }
}

// Сеттер для установления количество столбцов
void S21Matrix::SetRows(int rows) {
  if (rows < 1)
    throw std::out_of_range("Rows of matrix should be more or equal 1");
  if (rows != rows_) {
    S21Matrix result(rows, cols_);
    for (int i = 0; i < (rows_ < rows ? rows_ : rows); ++i)
      for (int k = 0; k < cols_; ++k) result.matrix_[i][k] = matrix_[i][k];
    *this = std::move(result);
  }
}

// фукции

bool S21Matrix::EqMatrix(const S21Matrix& other) const noexcept {
  if (rows_ <= 0 || cols_ <= 0 || other.rows_ <= 0 || other.cols_ <= 0) {
    return false;
  }
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    return false;
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (std::abs(matrix_[i][j] - other.matrix_[i][j]) > EPS) {
        return false;
      }
    }
  }
  return true;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  else {
    for (int i = 0; i < rows_; ++i)
      for (int k = 0; k < cols_; ++k) matrix_[i][k] += other.matrix_[i][k];
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  else {
    for (int i = 0; i < rows_; ++i)
      for (int k = 0; k < cols_; ++k) matrix_[i][k] -= other.matrix_[i][k];
  }
}

void S21Matrix::MulNumber(const double num) noexcept {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] * num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_)
    throw std::out_of_range(
        "Columns of matrix 1 should be the same size as rows of matrix 2");
  S21Matrix result(rows_, other.cols_);
  for (int i = 0; i < result.rows_; ++i)
    for (int k = 0; k < result.cols_; ++k)
      for (int p = 0; p < cols_; ++p)
        result.matrix_[i][k] += matrix_[i][p] * other.matrix_[p][k];
  *this = result;
}

S21Matrix S21Matrix::Transpose() const {
  if (rows_ <= 0 || cols_ <= 0) {
    throw std::invalid_argument("Matrix must have positive dimensions");
  }

  S21Matrix transposed(cols_, rows_);

  for (int i = 0; i < cols_; i++) {
    for (int j = 0; j < rows_; j++) {
      transposed.matrix_[i][j] = matrix_[j][i];
    }
  }

  return transposed;
}

// operators
double& S21Matrix::operator()(int i, int j) {
  if (i < 0 || i >= rows_ || j < 0 || j >= cols_)
    throw std::out_of_range(
        "Arguments i(j) should be more or equal then 0 and less then"
        "rows_(cols_)");
  else
    return matrix_[i][j];
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this != &other) {
    ClearMatrix();
    CreateMatrix(other.rows_, other.cols_);
    for (int i = 0; i < other.rows_; ++i)
      for (int k = 0; k < other.cols_; ++k) matrix_[i][k] = other.matrix_[i][k];
  }
  return *this;
}

bool S21Matrix::operator==(const S21Matrix& other) { return EqMatrix(other); }

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double num) {
  MulNumber(num);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix result(*this);
  result.SubMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix result(*this);
  result.MulNumber(num);
  return result;
}

S21Matrix operator*(double num, S21Matrix& matrix) { return matrix * num; }

double S21Matrix::Determinant() const {
  if (cols_ != rows_)
    throw std::out_of_range("Exception: The matrix is not square");
  double det = 0;
  if (rows_ == 1) {
    det = matrix_[0][0];
  } else if (rows_ == 2) {
    det = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  } else {
    for (int i = 0; i < rows_; i++) {
      S21Matrix M(Minor_(1, i + 1));
      det += pow((-1), i) * matrix_[0][i] * M.Determinant();
    }
  }
  return det;
}

S21Matrix S21Matrix::CalcComplements() const {
  if (rows_ != cols_)
    throw std::out_of_range("CalcComplementsError: The matrix must be square");

  S21Matrix result(rows_, cols_);
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) {
      S21Matrix M(Minor_(i + 1, j + 1));
      result.matrix_[i][j] = M.Determinant() * pow(-1, i + j);
    }

  return result;
}

S21Matrix S21Matrix::InverseMatrix() const {
  double det = Determinant();
  if (det == 0) throw std::out_of_range("InverseError: Determinant is 0");

  S21Matrix result(rows_, cols_);
  result = CalcComplements();
  result = result.Transpose();
  return result * (1.0 / det);
}

void S21Matrix::ClearMatrix() noexcept {
  if (matrix_) {
    for (int i = 0; i < rows_; ++i) delete[] matrix_[i];  // clear memory
  }
  delete[] matrix_;
}

void S21Matrix::CreateMatrix(int rows, int cols) {
  rows_ = rows;
  cols_ = cols;
  matrix_ = new double*[rows_];
  for (int i = 0; i < rows_; ++i) {
    matrix_[i] = new double[cols_];
    for (int j = 0; j < cols_; ++j) matrix_[i][j] = 0;
  }
}

S21Matrix S21Matrix::Minor_(int rows, int cols) const {
  S21Matrix M(rows_ - 1, cols_ - 1);
  int m = 0, n = 0;

  for (int i = 0; i < rows_; i++) {
    m = i;
    if (i > rows - 1) {
      m--;
    }
    for (int j = 0; j < cols_; j++) {
      n = j;
      if (j > cols - 1) {
        n--;
      }
      if (i != rows - 1 && j != cols - 1) {
        M.matrix_[m][n] = matrix_[i][j];
      }
    }
  }
  return M;
}