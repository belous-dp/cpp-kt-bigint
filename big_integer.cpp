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

//todo пересмотреть код, потому что конструктор по-умолчанию не создает вектор
// но должен считаться как 0

const big_integer ZERO = big_integer(0);
const big_integer ONE = big_integer(1);

big_integer::big_integer() = default;

big_integer::big_integer(big_integer const& other) : n(other.n) {}

big_integer::big_integer(int a) : n(1, a) {}

big_integer::big_integer(std::string const& str) {} // todo

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
  uint carry = 0;
  expand_size(rhs.n.size());
  for (size_t i = 0; i < n.size(); ++i) {
    n[i] += carry;
    carry = (n[i] < carry);
    uint add = i < rhs.n.size() ? rhs.n[i] : sext(rhs.back());
    n[i] += add;
    carry += n[i] < add;
  }
  pop_back_unused();
  return *this;
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
    std::cerr << "this is negative" << std::endl;
    negate();
  }
  big_integer const* prhs = &rhs;
  if (sign2) {
    std::cerr << "rhs is negative" << std::endl;
    auto* tmp = new big_integer(rhs);
    tmp->negate();
    prhs = tmp;
  }
  std::cerr << "prhs:\n" << prhs->to_bin_string() << std::endl;
  big_integer res = 0;
  for (size_t i = 0; i < prhs->n.size(); ++i) {
    big_integer add;
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
  std::cerr << "intermediate res:\n" << res.to_bin_string() << std::endl;
  if (sign1 ^ sign2) {
    res.negate();
  }
  if (sign2) {
    delete prhs;
  }
  *this = res;
  return *this;
}

big_integer& big_integer::operator/=(big_integer const& rhs) { // todo
  return *this;
}

big_integer& big_integer::operator%=(big_integer const& rhs) { // todo
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

big_integer& big_integer::operator<<=(int rhs) { // todo
  return *this;
}

big_integer& big_integer::operator>>=(int rhs) { // todo
  return *this;
}

big_integer big_integer::operator+() const {
  return *this;
}

big_integer big_integer::operator-() const {
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
  return *this + ONE;
}

big_integer& big_integer::operator--() {
  return *this -= ONE;
}

big_integer big_integer::operator--(int) {
  return *this - ONE;
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
  return a.n == b.n;
}

bool operator!=(big_integer const& a, big_integer const& b) {
  return !(a == b);
}

bool operator<(big_integer const& a, big_integer const& b) { // todo
  return true;
}

bool operator>(big_integer const& a, big_integer const& b) { // todo
  return true;
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
