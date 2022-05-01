/**
 * Числа хранятся следующим образом:
 * вектор из 32-битных беззнаковых разрядов;
 * младший разряд лежит в n[0];
 * дополнение до двух, но немного модифицированное
 * если число отрицательное, его запись совпадает с обычным дополнением до двух
 * если число положительное, возможны две ситуации
 * 1. старший 32-битный разряд числа помещается в int.
 *    то есть старший разряд выглядит как *******************************0
 *      в этом случае записывается его обычная двоичная запись (само число)
 * 2. старший 32-битный разряд числа НЕ помещается в int.
 *    то есть старший разряд выглядит как *******************************1
 *      в этом случае записывается его обычная двоичная запись, но
 *      также добавляется лидирующий разряд, в котором записан 0
 * это сделано для того, чтобы отличать отрицательные числа
 * от больших (> MAX_INT, но < MAX_UINT) положительных
 *
*/

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

uint lo_word(unsigned long long x) {
  uint mask = -1;
  return x & mask;
}

uint hi_word(unsigned long long x) {
  return x >> std::numeric_limits<uint>::digits;
}

big_integer::big_integer(int a) : n(1, a) { }

big_integer::big_integer(uint a) : n({a, 0}) {
  pop_back_unused();
}

big_integer::big_integer(long a) : n({lo_word(a), hi_word(a)}) {
  pop_back_unused();
}

big_integer::big_integer(unsigned long a) : n({lo_word(a), hi_word(a), 0}) {
  pop_back_unused();
}

big_integer::big_integer(long long a) : n({lo_word(a), hi_word(a)}) {
  pop_back_unused();
}

big_integer::big_integer(unsigned long long a) : n({lo_word(a), hi_word(a), 0}) {
  pop_back_unused();
}



std::invalid_argument err(std::string const& str) {
  return std::invalid_argument("Invalid 10-based integer: " + str);
}

big_integer::big_integer(std::string const& str) {
  bool sign = str[0] == '-';
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
    carry += (n[i] < add);
  }
  return carry;
}

void big_integer::expand_size(size_t size) {
  n.resize(std::max(n.size() + 1, size), sext(back()));
}

void big_integer::pop_back_unused() {
  // если последний разряд все 1, и предпоследний точно отрицательный,
  // либо если последний разряд все 0 и предпоследний точно положительный
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

// обычное умножение в столбик. Умножаем this на последний разряд,
// складываем, учитывая смещение на B^шаг
// здесь и далее B или b это размер разряда, т.е. MAX_UINT
big_integer& big_integer::operator*=(big_integer const& rhs) {
  bool sign1 = negative(), sign2 = rhs.negative();
  if (sign1) {
    negate();
  }
  big_integer const* prhs = &rhs;
  if (sign2) {
    auto* tmp = new big_integer(rhs);
    tmp->negate();
    prhs = tmp;
  }
  big_integer res = 0;
  for (size_t i = 0; i < prhs->n.size(); ++i) {
    big_integer add;
    add.n.resize(i, 0);
    uint carry = 0;
    for (unsigned int j : n) {
      unsigned long long mult = j * 1ULL * prhs->n[i];
      mult += carry;
      add.n.push_back(lo_word(mult));
      carry = hi_word(mult);
    }
    add.n.push_back(carry);
    res += add;
  }
  if (sign1 ^ sign2) {
    res.negate();
  }
  if (sign2) {
    delete prhs;
  }
  *this = res;
  return *this;
}

/// divides this by integer rhs, assuming that this > rhs > 0
/// returns reminder
uint big_integer::short_divide(uint rhs) {
  uint r = 0;
  n.push_back(0);
  const unsigned long long b = std::numeric_limits<uint>::max() + 1ULL;
  for (size_t j = n.size(); j > 0; --j) {
    uint tmp = n[j - 1];
    n[j - 1] = (r * b + tmp) / rhs; // точно помещается в uint
    r = (r * b + tmp) % rhs; // точно помещается в uint
  }
  pop_back_unused();
  return r;
}

size_t div_size(std::vector<uint> const& v) {
  return v.size() - (!v.empty() && v.back() == 0);
}

/// divides this by big_integer rhs
/// returns reminder
// algorithm described in The Art of Computer Programming (vol.2, 3rd ed.),
// chapter 4.3.1, p. 272
big_integer big_integer::divide(big_integer const& rhs) {
  if (rhs == ZERO) {
    throw std::invalid_argument("division by zero");
  }
  if (*this == rhs) {
    *this = ONE;
    return ZERO;
  }

  bool sign1 = negative(), sign2 = rhs.negative();
  bool sign = sign1 ^ sign2;
  if (sign1) {
    negate();
  }
  big_integer v = rhs;
  if (sign2) {
    v.negate();
  }

  if (div_size(v.n) == 1) {
    big_integer rem = short_divide(v.n[0]);
    if (sign) {
      negate();
    }
    return (sign ^ sign2) ? -rem : rem;
  }
  if (div_size(n) < div_size(v.n)) {
    big_integer rem = (sign ^ sign2) ? -*this : *this;
    *this = ZERO;
    return rem;
  }


  const size_t sn = div_size(v.n);
  const size_t sm = div_size(n) - sn;
  const unsigned long long b = std::numeric_limits<uint>::max() + 1ULL;
  // D1
  uint d = b / (v.n[sn - 1] + 1);
  *this *= d;
  v *= d;
  if (n.size() < sm + sn + 1) {
    n.push_back(0);
  }
  assert(n.size() == sm + sn + 1);
  // quotient stored in the greatest sm + 1 blocks of n
  // reminder stored in the lowest sn blocks of n
  // D2-D7
  for (size_t j = sm; j >= 0; --j) {
    // D3
    unsigned long long q0 = (n[j + sn] * b + n[j + sn - 1]) / v.n[sn - 1];
    unsigned long long r0 = (n[j + sn] * b + n[j + sn - 1]) % v.n[sn - 1];
    while ((q0 == b || q0 * v.n[sn - 2] > b * r0 + n[j + sn - 2]) && r0 < b) {
      q0--;
      r0 += v.n[sn - 1];
    }
    // D4-D5
    //  at this step q <= b
    uint q = lo_word(q0);

    big_integer subt = v * q;
    subt.negate();
    uint carry = add_from(j, j + sn + 1, subt);
    if (!carry && subt != 0) {
      // D6
      // нам нужно как-то понять, происходил заём или нет
      // вместо вычитания из u v*q сделаем u += v*q (в соответствующих разрядах)
      // как я рассуждал:
      // пусть в последнем разряде произошёл перенос
      // это значит, что в следующем разряде после последнего мы сделаем
      // u[x] + 1 + (v*q)[y],
      // но (v*q)[y] = -1, т.к. у нас дополнение до двух, а мы прибавляем отрицательное число.
      // Тогда u[x] не изменится, значит при вычитании мы из него ничего не занимали
      // единственное исключение из этого правила -- когда вычитаем 0.
      q--;
      add_from(j, j + sn + 1, v);
    }
    n[j + sn] = q;
    if (j == 0) {
      break;
    }
  }
  big_integer rem;
  std::reverse(n.begin(), n.end());
  while (n.size() > sm + 1) {
    rem.n.push_back(n.back());
    n.pop_back();
  }
  rem.short_divide(d);
  std::reverse(n.begin(), n.end());
  pop_back_unused();
  if (sign) {
    negate();
  }
  return (sign ^ sign2) ? -rem : rem;
}

big_integer& big_integer::operator/=(big_integer const& rhs) {
  divide(rhs);
  return *this;
}

big_integer& big_integer::operator%=(big_integer const& rhs) {
  *this = divide(rhs);
  return *this;
}

big_integer&
big_integer::apply_bitwise_operator(big_integer const& rhs,
                                    const std::function<uint(uint, uint)>& op) {
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
  return a.n == b.n ||
         (a.n.size() <= 1 && b.n.size() <= 1 && a.back() == b.back());
}

bool operator!=(big_integer const& a, big_integer const& b) {
  return !(a == b);
}

bool operator<(big_integer const& a, big_integer const& b) {
  if (a.negative() ^ b.negative()) {
    return a.negative() && !b.negative();
  }
  if (a.n.size() != b.n.size()) {
    return (a.n.size() < b.n.size()) ^ a.negative();
  }
  for (size_t i = a.n.size(); i > 0; --i) {
    if (a.n[i - 1] != b.n[i - 1]) {
      return a.n[i - 1] < b.n[i - 1];
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

std::string big_integer::to_bin_string() const {
  if (n.empty()) {
    return "0";
  }
  std::string res;
  for (uint i : n) {
    std::string str =
        std::bitset<std::numeric_limits<uint>::digits>(i).to_string();
    std::reverse(str.begin(), str.end());
    res += str;
  }
  return res;
}

std::string to_string(big_integer const& a) {
  big_integer b = a;
  bool sign = b.negative();
  if (sign) {
    b.negate();
  }
  const size_t CHUNK = 9;
  const big_integer DIV = (uint)1e9;
  std::vector<uint> v;
  while (b != ZERO) {
    v.push_back((b % DIV).n[0]);
    b /= DIV;
  }
  std::string res = sign ? "-" : "";
  for (size_t i = v.size(); i > 0; --i) {
    if (v[i - 1] == 0 && res.empty()) {
      continue;
    }
    std::string cur = std::to_string(v[i - 1]);
    if (res.size() > sign && cur.size() < CHUNK) {
      std::reverse(cur.begin(), cur.end());
      while (!res.empty() && cur.size() < CHUNK) {
        cur.push_back('0');
      }
      std::reverse(cur.begin(), cur.end());
    }
    res += cur;
  }
  if (res.empty()) {
    res += '0';
  }
  return res;
}

std::ostream& operator<<(std::ostream& s, big_integer const& a) {
  return s << to_string(a);
}
