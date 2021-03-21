#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <complex>

using std::cin;
using std::cout;
using std::string;
using std::vector;
using std::istream;
using std::ostream;
using std::complex;
using std::max;

typedef complex<double> base;

constexpr double pi = 3.1415;

class BigInteger {
public:
    BigInteger();

    ~BigInteger() = default;

    BigInteger(const string &value);

    BigInteger(int value);

    explicit operator bool() const;

    string toString() const;

    BigInteger operator-() const;

    BigInteger &operator++();

    BigInteger &operator--();

    BigInteger operator++(int);

    BigInteger operator--(int);

    friend BigInteger operator+(const BigInteger &left, const BigInteger &right);

    friend BigInteger operator-(const BigInteger &left, const BigInteger &right);

    friend BigInteger operator*(const BigInteger &left, const BigInteger &right);

    friend BigInteger operator/(const BigInteger &left, const BigInteger &right);

    friend BigInteger operator%(const BigInteger &left, const BigInteger &right);

    BigInteger &operator+=(const BigInteger &other);

    BigInteger &operator-=(const BigInteger &other);

    BigInteger &operator*=(const BigInteger &other);

    BigInteger &operator/=(const BigInteger &other);

    BigInteger &operator%=(const BigInteger &other);

    bool operator==(const BigInteger &other) const;

    bool operator!=(const BigInteger &other) const;

    bool operator<(const BigInteger &other) const;

    bool operator<=(const BigInteger &other) const;

    bool operator>(const BigInteger &other) const;

    bool operator>=(const BigInteger &other) const;

    friend ostream &operator<<(ostream &stream, const BigInteger &value);

    friend istream &operator>>(istream &stream, BigInteger &value);

//private:
    vector<int> num_;

    int sign_;

    friend BigInteger PositiveSum(const BigInteger &left, const BigInteger &right);

    friend BigInteger PositiveSub(const BigInteger &left, const BigInteger &right);

    friend BigInteger Abs(const BigInteger &number);

    void RemoveZeroes();
};

inline int ascii_letter_to_number(char letter) {
    return letter - '0';
}

inline char number_to_ascii_letter(int number) {
    return char(number + 48);
}

BigInteger::BigInteger() : sign_(0) {
    num_.push_back(0);
}

BigInteger::BigInteger(const string &value) {
    if (value == "-0") {
        sign_ = 0;
    } else {
        if (value.length() == 1 && value[0] == '0') {
            sign_ = 0;
        } else {
            sign_ = (value[0] == '-') ? -1 : 1;
        }
    }

    for (int i = static_cast<int>(value.length()) - 1; i >= 0; --i) {
        if (value[i] != '-') {
            num_.push_back(ascii_letter_to_number(value[i]));
        }
    }
}

BigInteger::BigInteger(int value) {
    sign_ = (value == 0) ? 0 : 1;
    if (value < 0) {
        sign_ = -1;
        value *= -1;
    }

    if (value == 0) {
        num_.push_back(0);
    } else {
        while (value != 0) {
            num_.push_back(value % 10);
            value /= 10;
        }
    }
}

BigInteger::operator bool() const {
    return (this->sign_ != 0);
}

string BigInteger::toString() const {
    string s;
    if (sign_ == 0) {
        return "0";
    }

    if (sign_ < 0) {
        s += '-';
    }
    for (int i = static_cast<int>(num_.size()) - 1; i >= 0; --i) {
        s += number_to_ascii_letter(num_[i]);
    }

    return s;
}

void BigInteger::RemoveZeroes() {
    int significant_position = static_cast<int>(num_.size()) - 1;
    while (num_[significant_position] == 0 && significant_position > 0) {
        significant_position--;
    }
    num_.resize(significant_position + 1);
}

BigInteger PositiveSum(const BigInteger &left, const BigInteger &right) {
    if (left.sign_ == 0) {
        return right;
    }
    if (right.sign_ == 0) {
        return left;
    }

    BigInteger result;
    result.num_.resize(0);
    result.sign_ = 1;

    size_t max_discharges = 0;
    if (left.num_.size() > right.num_.size()) {
        max_discharges = left.num_.size();
    } else {
        max_discharges = right.num_.size();
    }

    int carry = 0;
    for (size_t i = 0; i < max_discharges; ++i) {
        int right_num = (i >= right.num_.size()) ? 0 : right.num_[i];
        int left_num = (i >= left.num_.size()) ? 0 : left.num_[i];
        result.num_.push_back((carry + right_num + left_num) % 10);
        carry = (right_num + left_num + carry > 9) ? 1 : 0;
    }
    if (carry > 0) {
        result.num_.push_back(1);
    }

    return result;
}

BigInteger PositiveSub(const BigInteger &left, const BigInteger &right) {
    int carry = 0;
    BigInteger result;
    result.num_.resize(0);
    for (size_t i = 0; i < left.num_.size(); i++) {
        int left_num = left.num_[i];
        int right_num = (size_t(i) >= right.num_.size()) ? 0 : right.num_[i];
        carry += left_num - right_num + 10;
        result.num_.push_back(carry % 10);
        carry = (carry < 10) ? -1 : 0;
    }

    result.RemoveZeroes();

    if (!(result.num_.size() == 1 && result.num_[0] == 0)) {
        result.sign_ = 1;
    }

    return result;
}

BigInteger BigInteger::operator-() const {
    BigInteger result(*this);
    result.sign_ = -1 * this->sign_;
    return result;
}

BigInteger &BigInteger::operator++() {
    return *this += 1;
}

BigInteger &BigInteger::operator--() {
    return *this -= 1;
}

BigInteger BigInteger::operator++(int) {
    BigInteger old(*this);
    ++*this;
    return old;
}

BigInteger BigInteger::operator--(int) {
    BigInteger old(*this);
    --*this;
    return old;
}

BigInteger Abs(const BigInteger &number) {
    if (number.sign_ >= 0) {
        return number;
    } else {
        return -number;
    }
}

BigInteger PositiveNegativeSum(const BigInteger &negative, const BigInteger &positive) {
    if (positive > Abs(negative)) {
        return PositiveSub(positive, negative);
    } else {
        BigInteger result = PositiveSub(negative, positive);
        result.sign_ = (result.sign_ == 0) ? 0 : -1;
        return result;
    }
}

BigInteger operator+(const BigInteger &left, const BigInteger &right) {
    BigInteger result;
    if (left.sign_ >= 0 && right.sign_ >= 0) {
        return PositiveSum(left, right);
    }
    if (left.sign_ < 0 && right.sign_ < 0) {
        result = PositiveSum(left, right);
        result.sign_ = -result.sign_;
        return result;
    }
    if (left.sign_ < 0) {
        return PositiveNegativeSum(left, right);
    } else {
        return PositiveNegativeSum(right, left);
    }
}

BigInteger operator-(const BigInteger &left, const BigInteger &right) {
    return left + (-right);
}

void FFT(vector<base> &coefs, bool invert) {
    int coefs_size = static_cast<int>(coefs.size());
    if (coefs_size == 1) {
        return;
    }
    vector<base> even_coefs(coefs_size / 2);
    vector<base> odd_coefs(coefs_size / 2);
    for (int i = 0, j = 0; i < coefs_size; i += 2, ++j) {
        even_coefs[j] = coefs[i];
        odd_coefs[j] = coefs[i + 1];
    }

    FFT(even_coefs, invert);
    FFT(odd_coefs, invert);

    double angle = 2 * pi / coefs_size * (invert ? -1 : 1);
    base w(1), wn(cos(angle), sin(angle));
    for (int i = 0; i < coefs_size / 2; ++i) {
        coefs[i] = even_coefs[i] + w * odd_coefs[i];
        coefs[i + coefs_size / 2] = even_coefs[i] - w * odd_coefs[i];
        if (invert) {
            coefs[i] /= 2, coefs[i + coefs_size / 2] /= 2;
        }
        w *= wn;
    }
}

BigInteger operator*(const BigInteger &left, const BigInteger &right) {
    BigInteger result(0);
    result.sign_ = 1;

    vector<base> left_copy(left.num_.begin(), left.num_.end());
    vector<base> right_copy(right.num_.begin(), right.num_.end());
    int final_size = 1;
    while (final_size < max(static_cast<int>(left.num_.size()), static_cast<int>(right.num_.size()))) {
        final_size <<= 1;
    }
    final_size <<= 1;
    left_copy.resize(final_size), right_copy.resize(final_size);

    FFT(left_copy, false);
    FFT(right_copy, false);

    for (int i = 0; i < final_size; ++i) {
        left_copy[i] *= right_copy[i];
    }
    FFT(left_copy, true);

    result.num_.resize(final_size);
    for (int i = 0; i < final_size; ++i)
        result.num_[i] = static_cast<int>(left_copy[i].real() + 0.5);

    int carry = 0;
    for (int i = 0; i < final_size; ++i) {
        result.num_[i] += carry;
        carry = result.num_[i] / 10;
        result.num_[i] %= 10;
    }

    result.sign_ = left.sign_ * right.sign_;
    result.RemoveZeroes();
    return result;
}

BigInteger operator/(const BigInteger &left, const BigInteger &right) {
    if (right.sign_ == 0) {
        return left;
    }

    if (Abs(left) < Abs(right)) {
        return BigInteger(0);
    }

    BigInteger copy_left(Abs(left));
    BigInteger copy_right(Abs(right));

    BigInteger result(0);

    while (copy_left >= copy_right) {
        BigInteger tmp;
        tmp.num_.resize(0);
        tmp.sign_ = 1;

        for (size_t i = 0; i < copy_right.num_.size(); i++) {
            tmp.num_.push_back(copy_left.num_[i + (copy_left.num_.size() - copy_right.num_.size())]);
        }
        if (tmp < copy_right) {
            tmp.num_.insert(tmp.num_.begin(), copy_left.num_[(copy_left.num_.size() - copy_right.num_.size()) - 1]);
        }

        bool good_num = false;
        for (int i = 9; i >= 1; --i) {
            BigInteger big_i(i);
            if (copy_right * big_i <= tmp && !good_num) {
                BigInteger to_ans(big_i);
                for (int j = 0; j < static_cast<int>(copy_left.num_.size() - tmp.num_.size()); ++j) {
                    to_ans.num_.insert(to_ans.num_.begin(), 0);
                }
                result = result + to_ans;
                copy_left = copy_left - to_ans * copy_right;
                good_num = true;
            }
        }
    }

    result.sign_ = left.sign_ * right.sign_;

    return result;
}

BigInteger operator%(const BigInteger &left, const BigInteger &right) {
    return left - (left / Abs(right)) * Abs(right);
}

BigInteger &BigInteger::operator+=(const BigInteger &other) {
    *this = *this + other;
    return *this;
}

BigInteger &BigInteger::operator-=(const BigInteger &other) {
    *this = *this - other;
    return *this;
}

BigInteger &BigInteger::operator*=(const BigInteger &other) {
    *this = *this * other;
    return *this;
}

BigInteger &BigInteger::operator/=(const BigInteger &other) {
    *this = *this / other;
    return *this;
}

BigInteger &BigInteger::operator%=(const BigInteger &other) {
    *this = *this % other;
    return *this;
}

bool BigInteger::operator==(const BigInteger &other) const {
    if (other.sign_ == 0 && this->sign_ == 0) {
        return true;
    }

    if (sign_ != other.sign_) {
        return false;
    }

    if (num_ != other.num_) {
        return false;
    }

    return true;
}

bool BigInteger::operator!=(const BigInteger &other) const {
    return !(*this == other);
}

bool BigInteger::operator<(const BigInteger &other) const {
    if (sign_ < other.sign_) {
        return true;
    }
    if (sign_ > other.sign_) {
        return false;
    }

    if (sign_ == 0 && other.sign_ == 0) {
        return false;
    }

    if (num_.size() < other.num_.size()) {
        return sign_ >= 0;
    }

    if (num_.size() > other.num_.size()) {
        return sign_ < 0;
    }

    for (int i = static_cast<int>(num_.size()) - 1; i >= 0; --i) {
        if (num_[i] < other.num_[i]) {
            return (sign_ != -1);
        }
        if (num_[i] > other.num_[i]) {
            return (sign_ == -1);
        }
    }

    return false;
}

bool BigInteger::operator<=(const BigInteger &other) const {
    return (*this < other) || (*this == other);
}

bool BigInteger::operator>(const BigInteger &other) const {
    return !(*this <= other);
}

bool BigInteger::operator>=(const BigInteger &other) const {
    return !(*this < other);
}

istream &operator>>(istream &stream, BigInteger &value) {
    string s;
    stream >> s;
    BigInteger tmp(s);
    value = tmp;
    return stream;
}

ostream &operator<<(ostream &stream, const BigInteger &value) {
    stream << value.toString();
    return stream;
}