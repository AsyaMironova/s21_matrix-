#ifndef S21_MATRIX_OOP_H_
#define S21_MATRIX_OOP_H_
#define EPS 1e-7

#include <cmath>
#include <iostream>

class S21Matrix {
 private:
  // атрибуты
  int rows_, cols_;
  double** matrix_;
  void ClearMatrix() noexcept;
  void CreateMatrix(int rows, int cols);
  S21Matrix Minor_(int rows, int cols) const;

 public:
  // конструктор
  S21Matrix();                        // база
  S21Matrix(int rows, int cols);      // параметры
  S21Matrix(const S21Matrix& other);  // копирование
  S21Matrix(S21Matrix&& other);       // движение
  // деструктор
  ~S21Matrix();
  // фукции
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num) noexcept;
  void MulMatrix(const S21Matrix& other);
  bool EqMatrix(const S21Matrix& other) const noexcept;
  S21Matrix Transpose() const;
  S21Matrix CalcComplements() const;
  double Determinant() const;
  S21Matrix InverseMatrix() const;
  // getters
  int GetRows() const;
  int GetCols() const;
  // setters
  void SetRows(int rows);
  void SetCols(int cols);
  // oоператоры
  double& operator()(int i, int j);
  S21Matrix& operator=(const S21Matrix& other);
  bool operator==(const S21Matrix& other);
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix operator*(const double num);
  friend S21Matrix operator*(double num, S21Matrix& matrix);
  S21Matrix& operator*=(const double num);
  S21Matrix& operator*=(const S21Matrix& other);
};

#endif  // S21_MATRIX_OOP_H_