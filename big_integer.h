#pragma once

#include <functional>
#include <iosfwd>
#include <string>
#include <vector>

using uint = unsigned int;

struct big_integer {
private:
  std::vector<uint> n;

  unsigned int back() const;
  void expand_size(size_t size);
  void pop_back_unused();

  big_integer&
  apply_bitwise_operator(big_integer const& rhs,
                         const std::function<uint(uint, uint)>& op);

  bool negative() const;
  void negate();

  bool is_zero() const;

  void short_multiply(uint rhs);

  uint short_divide(uint rhs);
  big_integer divide(big_integer const& rhs);

  void short_add(uint rhs);
  void short_add_segment(size_t from, uint rhs, uint add);
  uint add_segment(size_t from, size_t to, const big_integer& rhs);

public:
  big_integer();
  big_integer(big_integer const& other);
  big_integer(int a);
  big_integer(uint a);
  big_integer(long a);
  big_integer(unsigned long a);
  big_integer(long long a);
  big_integer(unsigned long long a);

  explicit big_integer(std::string const& str);
  ~big_integer();

  void swap(big_integer& other);

  big_integer& operator=(big_integer const& other);

  big_integer& operator+=(big_integer const& rhs);
  big_integer& operator-=(big_integer const& rhs);
  big_integer& operator*=(big_integer const& rhs);
  big_integer& operator/=(big_integer const& rhs);
  big_integer& operator%=(big_integer const& rhs);

  big_integer& operator&=(big_integer const& rhs);
  big_integer& operator|=(big_integer const& rhs);
  big_integer& operator^=(big_integer const& rhs);

  big_integer& operator<<=(int shift);
  big_integer& operator>>=(int shift);

  big_integer operator+() const;
  big_integer operator-() const;
  big_integer operator~() const;

  big_integer& operator++();
  big_integer operator++(int);

  big_integer& operator--();
  big_integer operator--(int);

  friend bool operator==(big_integer const& a, big_integer const& b);
  friend bool operator!=(big_integer const& a, big_integer const& b);
  friend bool operator<(big_integer const& a, big_integer const& b);
  friend bool operator>(big_integer const& a, big_integer const& b);
  friend bool operator<=(big_integer const& a, big_integer const& b);
  friend bool operator>=(big_integer const& a, big_integer const& b);

  std::string to_bin_string() const;
  friend std::string to_string(big_integer const& a);
};

big_integer operator+(big_integer a, big_integer const& b);
big_integer operator-(big_integer a, big_integer const& b);
big_integer operator*(big_integer a, big_integer const& b);
big_integer operator/(big_integer a, big_integer const& b);
big_integer operator%(big_integer a, big_integer const& b);

big_integer operator&(big_integer a, big_integer const& b);
big_integer operator|(big_integer a, big_integer const& b);
big_integer operator^(big_integer a, big_integer const& b);

big_integer operator<<(big_integer a, int b);
big_integer operator>>(big_integer a, int b);

bool operator==(big_integer const& a, big_integer const& b);
bool operator!=(big_integer const& a, big_integer const& b);
bool operator<(big_integer const& a, big_integer const& b);
bool operator>(big_integer const& a, big_integer const& b);
bool operator<=(big_integer const& a, big_integer const& b);
bool operator>=(big_integer const& a, big_integer const& b);

std::string to_string(big_integer const& a);
std::ostream& operator<<(std::ostream& s, big_integer const& a);
