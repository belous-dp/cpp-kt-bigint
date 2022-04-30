#include "big_integer.h"
#include <bitset>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <ostream>
#include <stdexcept>

using uint = unsigned int;

const big_integer ZERO = big_integer();
const big_integer ONE = big_integer(1);

big_integer::big_integer() = default;

big_integer::big_integer(big_integer const& other) = default;

big_integer::big_integer(int a) : n(1, a), sign(a < 0) {}

//big_integer::big_integer(uint a) : n(1, a) {}

std::invalid_argument err(std::string const& str) {
  return std::invalid_argument("Invalid 10-based integer: " + str);
}

big_integer::big_integer(std::string const& str) : sign(str[0] == '-') {
  auto cnt = std::count(str.begin(), str.end(), '-');
  if ((str.size() - sign) == 0 || cnt > 1 || (cnt == 1 && str.find('-') != 0)) {
    throw err(str);
  }
  const size_t CHUNK = 9;
  const big_integer MUL = (uint)1e9;
  size_t shift = (str.size() - sign) % CHUNK;
  if (shift > 0) {
    size_t pos = 0;
    n.push_back(std::stoi(str.substr(sign, shift), &pos));
    if (pos != shift) {
      throw err(str);
    }
  } else {
    n.push_back(0);
  }
  for (size_t i = sign + shift; i < str.size(); i += CHUNK) {
    *this *= MUL;
    size_t pos = 0;
    *this += std::stoi(str.substr(i, CHUNK), &pos);
    if (pos != CHUNK) {
      throw err(str);
    }
  }
  if (sign) {
    negate();
  }
}

big_integer::~big_integer() = default;

big_integer& big_integer::operator=(big_integer const& other) {
  if (this != &other) {
    big_integer copy(other);
    swap(copy);
  }
  return *this;
}

void big_integer::swap(big_integer& other) {
  std::swap(n, other.n);
}

std::string big_integer::to_bin_string() const {
  if (n.empty()) {
    return "0";
  }
  std::string res;
  for (uint i : n) {
    std::string str = std::bitset<std::numeric_limits<uint>::digits>(i).to_string();
    std::reverse(str.begin(), str.end());
    res += str;
  }
  return res;
}

uint msb(uint x) {
  return x >> (std::numeric_limits<uint>::digits - 1);
}

uint sext(uint x) {   // sign extension
  return -1 * msb(x); // perhaps UINT_MAX * msb(x) or
                      // std::numeric_limits<uint>::max() * msb(x)
}

uint big_integer::back() const {
  return n.empty() ? 0 : n.back();
}

big_integer& big_integer::operator+=(big_integer const& rhs) {
  expand_size(rhs.n.size());
  add_from(0, n.size(), rhs);
  pop_back_unused();
  return *this;
}

uint big_integer::add_from(size_t from, size_t to, const big_integer& rhs) {
  uint carry = 0;
  for (size_t i = from, j = 0; i < to; ++i, ++j) {
    n[i] += carry;
    carry = (n[i] < carry);
    uint add = j < rhs.n.size() ? rhs.n[j] : sext(rhs.back());
    n[i] += add;
    carry += n[i] < add;
  }
  return carry;
}

void big_integer::expand_size(size_t size) {
  n.resize(std::max(n.size() + 1, size), sext(back()));
}

void big_integer::pop_back_unused() {
  while (n.size() >= 2 && sext(n[n.size() - 2]) == n[n.size() - 1]) {
    n.pop_back();
  }
}

big_integer& big_integer::operator-=(big_integer const& rhs) {
  return *this += -rhs;
}

bool big_integer::negative() const {
  return msb(back());
}

void big_integer::negate() {
  if (n.empty())
    return;
  for (uint& i : n) {
    i = ~i;
  }
  ++*this;
}

uint lo_word(unsigned long long x) {
  uint mask = -1;
  return x & mask;
}

uint hi_word(unsigned long long x) {
  return x >> std::numeric_limits<uint>::digits;
}

big_integer& big_integer::operator*=(big_integer const& rhs) {
  bool sign1 = negative(), sign2 = rhs.negative();
  if (sign1) {
//    std::cerr << "this is negative" << std::endl;
    negate();
  }
  big_integer const* prhs = &rhs;
  if (sign2) {
//    std::cerr << "rhs is negative" << std::endl;
    auto* tmp = new big_integer(rhs);
    tmp->negate();
    prhs = tmp;
  }
//  std::cerr << "lhs:\n" << to_bin_string() << std::endl;
//  std::cerr << "rhs:\n" << prhs->to_bin_string() << std::endl;
  big_integer res = 0;
  for (size_t i = 0; i < prhs->n.size(); ++i) {
    big_integer add;
    add.n.resize(i, 0);
    uint carry = 0;
    for (size_t j = 0; j < n.size(); ++j) {
      unsigned long long mult = n[j] * 1ULL * prhs->n[i];
      mult += carry;
      add.n.push_back(lo_word(mult));
      carry = hi_word(mult);
    }
    add.n.push_back(carry);
    //pop_back_unused(); // ???
    res += add;
  }
//  std::cerr << "intermediate res:\n" << res.to_bin_string() << std::endl;
  if (sign1 ^ sign2) {
    res.negate();
  }
  if (sign2) { //todo smart pointers
    delete prhs;
  }
  *this = res;
  return *this;
}

/// divides this by integer rhs, assuming that this > rhs
uint big_integer::short_divide(uint rhs) {
  return 0;
}

/// divides this by big_integer rhs
// algorithm described in The Art of Computer Programming (vol.2, 3rd ed.),
// chapter 4.3.1, p. 272
big_integer big_integer::divide(big_integer const& rhs) {
//  if (rhs == ZERO) {
//    throw std::invalid_argument("division by zero");
//  }
//  // todo отрицательные числа, случай при rhs.back() == 0
//  if (rhs.n.size() == 1) {
//    return short_divide(rhs.back());
//  }
//  if (*this < rhs) {
//    big_integer tmp = *this;
//    *this = ZERO;
//    return tmp;
//  } else if (*this == rhs) {
//    *this = ONE;
//    return ZERO;
//  }
//  bool sign1 = negative(), sign2 = rhs.negative();
//  if (sign1) {
//    //    std::cerr << "this is negative" << std::endl;
//    negate();
//  }
//  big_integer const* v = &rhs;
//  if (sign2) {
//    //    std::cerr << "rhs is negative" << std::endl;
//    auto* tmp = new big_integer(rhs);
//    tmp->negate();
//    v = tmp;
//  }
//  if (n.back() == 0) {
//    n.pop_back();
//  }
//  if (v->n.back() == 0) {
//    v->n.pop_back();
//  }
//
//  const size_t sn = v->n.size();
//  const size_t sm = n.size() - sn;
//  const uint b = std::numeric_limits<uint>::max();
//  //D1
//  uint d = b / (v->n[sn - 1] + 1);
//  *this *= d;
//  assert(n.size() == sm + sn + 1);
//  std::vector<uint> quot(sm);
//  //D2-D7
//  for (size_t j = sm; j >= 0; --j) {
//    //D3
//    unsigned long long q0 = (n[j + sn] * b + n[j + sn - 1]) / v->n[sn - 1];
//    unsigned long long r0 = (n[j + sn] * b + n[j + sn - 1]) % v->n[sn - 1];
//    while ((q0 == b || q0 * v->n[sn - 2] > b * r0 + n[j + sn - 2]) && r0 < b) {
//      q0--;
//      r0 += v->n[sn - 1];
//    }
//    //D4-D5
//    // at this step q <= b
//    uint q = lo_word(q0);
//    big_integer subt = rhs * q;
//    subt.negate();
//    uint neg = add_from(j, j + sn, subt);
//    if (neg) {
//      //D6
//      q--;
//      add_from(j, j + sn, rhs);
//    }
//    quot[j] = q;
//  }
//  short_divide(d);
//  big_integer rem = *this;
//  n = quot;
//  if (sign2) {
//    delete v;
//  }
//  return rem;
  return *this;
}

big_integer& big_integer::operator/=(big_integer const& rhs) {
  divide(rhs);
  return *this;
}

big_integer& big_integer::operator%=(big_integer const& rhs) {
  *this = divide(rhs);
  return *this;
}

big_integer& big_integer::apply_bitwise_operator(
    big_integer const& rhs, const std::function<uint(uint, uint)>& op) {
  // вопрос: можно ли здесь как-то вместо functional
  // использовать указатели на функции?
  expand_size(rhs.n.size());
  for (size_t i = 0; i < n.size(); ++i) {
    uint other = i < rhs.n.size() ? rhs.n[i] : sext(rhs.back());
    n[i] = op(n[i], other);
  }
  pop_back_unused();
  return *this;
}

big_integer& big_integer::operator&=(big_integer const& rhs) {
  return apply_bitwise_operator(rhs, std::bit_and<>());
}

big_integer& big_integer::operator|=(big_integer const& rhs) {
  return apply_bitwise_operator(rhs, std::bit_or<>());
}

big_integer& big_integer::operator^=(big_integer const& rhs) {
  return apply_bitwise_operator(rhs, std::bit_xor<>());
}

big_integer& big_integer::operator<<=(int shift) {
  assert(shift >= 0);
  int d = shift / std::numeric_limits<uint>::digits;
  int r = shift % std::numeric_limits<uint>::digits;
  if (d > 0) {
    std::reverse(n.begin(), n.end());
    while (d--) {
      n.push_back(0);
    }
    std::reverse(n.begin(), n.end());
  }
  if (r > 0) {
    expand_size(n.size() + 1);
    for (size_t i = n.size(); i > 1; --i) {
      n[i - 1] <<= r;
      n[i - 1] |= n[i - 2] >> (std::numeric_limits<uint>::digits - r);
    }
    if (!n.empty()) {
      n[0] <<= r;
    }
  }
  pop_back_unused();
  return *this;
}


big_integer& big_integer::operator>>=(int shift) {
  assert(shift >= 0);
  int d = shift / std::numeric_limits<uint>::digits;
  int r = shift % std::numeric_limits<uint>::digits;
  if (d > 0) {
    std::reverse(n.begin(), n.end());
    while (!n.empty() && d--) {
      n.pop_back();
    }
    std::reverse(n.begin(), n.end());
  }
  if (r > 0) {
    for (size_t i = 1; i < n.size(); ++i) {
      n[i - 1] >>= r;
      n[i - 1] |= n[i] << (std::numeric_limits<uint>::digits - r);
    }
    if (!n.empty()) {
      uint tail = sext(n.back());
      n.back() >>= r;
      tail >>= (std::numeric_limits<uint>::digits - r);
      tail <<= (std::numeric_limits<uint>::digits - r);
      n.back() |= tail;
    }
  }
  pop_back_unused();
  return *this;
}

big_integer big_integer::operator+() const {
  return *this;
}

big_integer big_integer::operator-() const {
  if (n.empty())
    return *this;
  return ++~*this;
}

big_integer big_integer::operator~() const {
  big_integer res(*this);
  for (uint& i : res.n) {
    i = ~i;
  }
  return res;
}

big_integer& big_integer::operator++() {
  return *this += ONE;
}

big_integer big_integer::operator++(int) {
  big_integer copy = *this;
  ++*this;
  return copy;
}

big_integer& big_integer::operator--() {
  return *this -= ONE;
}

big_integer big_integer::operator--(int) {
  big_integer copy = *this;
  --*this;
  return copy;
}

big_integer operator+(big_integer a, big_integer const& b) {
  return a += b;
}

big_integer operator-(big_integer a, big_integer const& b) {
  return a -= b;
}

big_integer operator*(big_integer a, big_integer const& b) {
  return a *= b;
}

big_integer operator/(big_integer a, big_integer const& b) {
  return a /= b;
}

big_integer operator%(big_integer a, big_integer const& b) {
  return a %= b;
}

big_integer operator&(big_integer a, big_integer const& b) {
  return a &= b;
}

big_integer operator|(big_integer a, big_integer const& b) {
  return a |= b;
}

big_integer operator^(big_integer a, big_integer const& b) {
  return a ^= b;
}

big_integer operator<<(big_integer a, int b) {
  return a <<= b;
}

big_integer operator>>(big_integer a, int b) {
  return a >>= b;
}

bool operator==(big_integer const& a, big_integer const& b) {
//  std::cerr << "a:" << std::endl;
//  std::cerr << a.to_bin_string() << std::endl;
//  std::cerr << "b:" << std::endl;
//  std::cerr << b.to_bin_string() << std::endl;
  return a.n == b.n || (a.n.size() <= 1 && b.n.size() <= 1 && a.back() == b.back());
}

bool operator!=(big_integer const& a, big_integer const& b) {
  return !(a == b);
}

bool operator<(big_integer const& a, big_integer const& b) {
  if (a.negative() ^ b.negative()) {
    return a.negative() && !b.negative();
  }
  bool inv = a.negative();
  if (a.n.size() != b.n.size()) {
    return (a.n.size() < b.n.size()) ^ inv;
  }
  for (size_t i = a.n.size(); i > 0; --i) {
    if (a.n[i - 1] != b.n[i - 1]) {
      return (a.n[i - 1] < b.n[i - 1]) ^ inv;
    }
  }
  return false;
}

bool operator>(big_integer const& a, big_integer const& b) {
  return !(a < b || a == b);
}

bool operator<=(big_integer const& a, big_integer const& b) {
  return !(a > b);
}

bool operator>=(big_integer const& a, big_integer const& b) {
  return !(a < b);
}

std::string to_string(big_integer const& a) { // todo
  return "";
}

std::ostream& operator<<(std::ostream& s, big_integer const& a) {
  return s << to_string(a);
}
