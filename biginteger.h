//
// Created by Timur on 21.03.2020.
//

#ifndef BIGINT_BIGINTEGER_H
#define BIGINT_BIGINTEGER_H

#include <iostream>
#include <vector>
#include <string>


const int BASE = 1000000000;

class BigInteger {
public:
    BigInteger();
    BigInteger(long long longnumber);
    explicit BigInteger(const std::string& str);
    BigInteger(const BigInteger& b);

    std::string toString() const;
    BigInteger& operator=(const BigInteger& number);
    BigInteger& operator+=(const BigInteger& number);
    BigInteger& operator-=(const BigInteger& number);
    BigInteger& operator*=(const BigInteger& number);
    BigInteger& operator/=(const BigInteger& number);
    BigInteger& operator%=(const BigInteger& number);

    BigInteger& operator++();
    BigInteger operator++(int);
    BigInteger& operator--();
    BigInteger operator--(int);

    friend bool operator<(const BigInteger& first, const BigInteger& second);
    friend bool operator>(const BigInteger& first, const BigInteger& second);
    friend bool operator<=(const BigInteger& first, const BigInteger& second);
    friend bool operator>=(const BigInteger& first, const BigInteger& second);
    friend bool operator==(const BigInteger& first, const BigInteger& second);
    friend bool operator!=(const BigInteger& first, const BigInteger& second);

    BigInteger operator-() const;

    bool ifNegative() const;
    explicit operator bool() const;
    static BigInteger abs(const BigInteger& number);

private:
    static int moduleCompare(const BigInteger& first, const BigInteger& second);
    static int compare(const BigInteger& first, const BigInteger& second);
    void split(BigInteger& lessPart, BigInteger& bigPart, int length) const;
    static BigInteger multiple(const BigInteger& first, const BigInteger& second);
    static void div(const BigInteger& first, const BigInteger& second, BigInteger& result, BigInteger& mod);

    BigInteger& shift(int distance);
    void deleteLeadingZeros();
    BigInteger multipleOnInt(int number);
    void plus(const BigInteger& number);
    void minusSmaller(const BigInteger& number);
    void minusBigger(const BigInteger& number);

    std::vector<int> digits;
    bool sign;
};

bool BigInteger::ifNegative() const {
    return !sign;
}

std::istream& operator>>(std::istream& in, BigInteger& b) {
    std::string s;
    in >> s;
    b = BigInteger(s);
    return in;
}

std::ostream& operator<<(std::ostream& out,const BigInteger& b) {
    return out << b.toString();
}

BigInteger::BigInteger() : sign(true), digits(1, 0) {}

BigInteger::BigInteger(long long longnumber) {
    long long module = 0;
    if (longnumber > 0) {
        module = longnumber;
    } else module = -longnumber;
    long long ten = 1;
    while (module / ten > 0) ten *= 10;
    ten /= 10;
    std::string str;
    while (ten > 0) {
        str.push_back(char('0' + ((module / ten) % 10)));
        ten /= 10;
    }
    if (str.empty()) str.push_back(char('0'));
    BigInteger b = BigInteger(str);
    for (int digit : b.digits) {
        digits.push_back(digit);
    }
    sign = (longnumber >= 0);
}

BigInteger::BigInteger(const std::string& str) {
    int i = static_cast<int>(str.size()) - 1;
    int firstDigit = 0;
    if (str[0] == '-') {
        sign = false;
        firstDigit = 1;
    } else sign = true;
    while (i >= firstDigit) {
        int digit = 0;
        int ten = 1;
        for (int k = 0; k < 9; ++k) {
            if (i >= firstDigit) digit += ten * ((str[i] - '0'));
            ten *= 10;
            --i;
        }
        digits.push_back(digit);
    }
    (*this).deleteLeadingZeros();
}

BigInteger::BigInteger(const BigInteger& b) {
    sign = b.sign;
    digits.clear();
    for (int digit : b.digits) {
        digits.push_back(digit);
    }
}

std::string BigInteger::toString() const {
    std::string str;
    str.clear();
    if (!sign) {
        str.push_back('-');
    }
    bool ifFindDigit = false;
    for (int k = BASE / 10; k > 0; k /= 10) {
        int d = (digits[digits.size() - 1] / k) % 10;
        if (d != 0) ifFindDigit = true;
        if (ifFindDigit) str.push_back((char('0' + d)));
    }
    for (int i = static_cast<int>(digits.size()) - 2; i >= 0; --i) {
        for (int k = BASE / 10; k > 0; k /= 10) {
            int d = (digits[i] / k) % 10;
            str.push_back((char('0' + d)));
        }
    }
    if (str.empty()) str.push_back(char('0'));
    return str;
}

void BigInteger::deleteLeadingZeros() {
    while (digits.size() > 1 && digits[digits.size() - 1] == 0) {
        digits.pop_back();
    }
    if (digits[digits.size() - 1] == 0) {
        sign = true;
    }
}

BigInteger::operator bool() const {
    return (*this) != 0;
}

int BigInteger::moduleCompare(const BigInteger& first, const BigInteger& second) {
    if (first.digits.size() > second.digits.size()) return 1;
    if (first.digits.size() < second.digits.size()) return -1;
    int i = static_cast<int>(first.digits.size()) - 1;
    while (i >= 0 && first.digits[i] == second.digits[i]) --i;
    if (i < 0) return 0;
    if (first.digits[i] > second.digits[i]) return 1;
    return -1;
}

int BigInteger::compare(const BigInteger& first, const BigInteger& second) {
    if (first.sign != second.sign) {
        if (first.sign) return 1;
        return -1;
    }
    if (first.sign) return moduleCompare(first, second);
    return (-1) * moduleCompare(first, second);
}

bool operator<(const BigInteger& first, const BigInteger& second) {
    return BigInteger::compare(first, second) == -1;
}

bool operator>(const BigInteger& first, const BigInteger& second) {
    return BigInteger::compare(first, second) == 1;
}

bool operator<=(const BigInteger& first, const BigInteger& second) {
    return !(first > second);
}

bool operator>=(const BigInteger& first, const BigInteger& second) {
    return !(first < second);
}

bool operator==(const BigInteger& first, const BigInteger& second) {
    return BigInteger::compare(first, second) == 0;
}

bool operator!=(const BigInteger& first, const BigInteger& second) {
    return !(first == second);
}

BigInteger BigInteger::operator-() const{
    BigInteger temp = (*this);
    temp.sign = !temp.sign;
    temp.deleteLeadingZeros();
    return temp;
}

void BigInteger::plus(const BigInteger& number) {
    int maxLength = std::max(digits.size(), number.digits.size());
    for (int i = 0; i < maxLength; ++i) {
        if (i >= static_cast<int>(digits.size())) {
            digits.push_back(0);
        }
        if (i < static_cast<int>(number.digits.size())) {
            digits[i] += number.digits[i];
        }
        if (digits[i] >= BASE) {
            if (i == static_cast<int>(digits.size()) - 1) digits.push_back(0);
            digits[i] -= BASE;
            digits[i + 1] += 1;
        }
    }
    deleteLeadingZeros();
}

void BigInteger::minusSmaller(const BigInteger& number) {
    for (size_t i = 0; i < digits.size(); ++i) {
        if (i < number.digits.size()) {
            digits[i] -= number.digits[i];
        }
        if (digits[i] < 0) {
            digits[i] += BASE;
            digits[i + 1] -= 1;
        }
    }
    deleteLeadingZeros();
}

void BigInteger::minusBigger(const BigInteger& number) {
    for (size_t i = 0; i < number.digits.size(); ++i) {
        if (i == digits.size()) digits.push_back(0);
        digits[i] = number.digits[i] - digits[i];
        if (digits[i] < 0) {
            digits[i] += BASE;
            if (i + 1 == digits.size()) digits.push_back(0);
            digits[i + 1] += 1;
        }
    }
    deleteLeadingZeros();
}

BigInteger& BigInteger::operator+=(const BigInteger& number) {
    if (sign == number.sign) {
        plus(number);
    } else {
        if (moduleCompare(*this, number) != -1) {
            minusSmaller(number);
        } else {
            sign = !sign;
            minusBigger(number);
        }
    }
    return *this;
}

BigInteger& BigInteger::operator-=(const BigInteger& number) {
    if (*this != number) {
        sign = !sign;
        *this += number;
        sign = !sign;
        deleteLeadingZeros();
    } else (*this) = BigInteger();
    return *this;
}

BigInteger& BigInteger::operator=(const BigInteger& number) {
    if (*this != number) {
        sign = number.sign;
        digits.clear();
        for (size_t i = 0; i < number.digits.size(); ++i) {
            digits.push_back(number.digits[i]);
        }
    }
    return *this;
}

BigInteger& BigInteger::operator++() {
    *this += 1;
    return *this;
}

BigInteger BigInteger::operator++(int) {
    BigInteger old = *this;
    *this += 1;
    return old;
}

BigInteger& BigInteger::operator--() {
    *this -= 1;
    return *this;
}

BigInteger BigInteger::operator--(int) {
    BigInteger old = *this;
    *this -= 1;
    return old;
}

BigInteger operator+(const BigInteger& first, const BigInteger& second) {
    BigInteger temp = first;
    temp += second;
    return temp;
}

BigInteger operator-(const BigInteger& first, const BigInteger& second) {
    BigInteger temp = first;
    temp -= second;
    return temp;
}

void BigInteger::split(BigInteger& lessPart, BigInteger& bigPart, int length) const{
    lessPart.digits.clear();
    bigPart.digits.clear();
    int i = static_cast<int>(digits.size()) - 1;
    if (i >= length) {
        for (int j = length; j < static_cast<int>(digits.size()); ++j) {
            bigPart.digits.push_back(digits[j]);
        }
    } else bigPart.digits.push_back(0);
    for (int j = 0; j < length && j < static_cast<int>(digits.size()); ++j) {
        lessPart.digits.push_back(digits[j]);
    }
    lessPart.deleteLeadingZeros();
}

BigInteger& BigInteger::shift(int distance) {
    for (int i = 0; i < distance; ++i) digits.emplace(digits.begin(), 0);
    return *this;
}

BigInteger BigInteger::multiple(const BigInteger& first, const BigInteger& second) {
    if (first.digits.size() == 1 && first.digits[0] == 0) return BigInteger();
    if (second.digits.size() == 1 && second.digits[0] == 0) return BigInteger();
    int maxLength = std::max(first.digits.size(), second.digits.size());
    if (maxLength != 1) {
        int length = (maxLength + 1) / 2;
        BigInteger firstLessPart, firstBigPart;
        first.split(firstLessPart, firstBigPart, length);
        BigInteger secondLessPart, secondBigPart;
        second.split(secondLessPart, secondBigPart, length);

        BigInteger firstMult = multiple(firstLessPart, secondLessPart);
        BigInteger thirdMult = multiple(firstBigPart, secondBigPart);
        BigInteger secondMult = multiple(firstBigPart + firstLessPart, secondBigPart + secondLessPart);
        secondMult -= firstMult;
        secondMult -= thirdMult;
        secondMult.shift(length);
        thirdMult.shift(2 * length);
        BigInteger multiple(firstMult + secondMult + thirdMult);
        multiple.sign = (first.sign == second.sign);
        return multiple;
    } else {
        long long mult = static_cast<long long>(first.digits[0]) * static_cast<long long>(second.digits[0]);
        if (mult < 0) mult = -mult;
        BigInteger multiple(mult);
        multiple.sign = (first.sign == second.sign);
        return multiple;
    }
}

BigInteger& BigInteger::operator*=(const BigInteger& number) {
    *this = BigInteger::multiple(*this, number);
    return *this;
}

BigInteger operator*(const BigInteger& first, const BigInteger& second) {
    BigInteger tmp = first;
    tmp *= second;
    return tmp;
}

BigInteger BigInteger::abs(const BigInteger& number) {
    BigInteger tmp = number;
    tmp.sign = true;
    return tmp;
}

BigInteger BigInteger::multipleOnInt(int number) {
    BigInteger answer;
    long long transition = 0;
    for (size_t i = 0; i < digits.size(); ++i) {
        long long mult = static_cast<long long>(digits[i]) * static_cast<long long>(number);
        answer.digits.push_back(0);
        answer.digits[i] = static_cast<int>((mult + transition) % BASE);
        transition = (mult + transition) / BASE;
    }
    if (transition > 0) {
        answer.digits[answer.digits.size() - 1] = static_cast<int>(transition);
    }
    answer.deleteLeadingZeros();
    return answer;
}

void BigInteger::div(const BigInteger& first, const BigInteger& second, BigInteger& result, BigInteger& mod) {
    BigInteger module = abs(second);
    std::vector<int> tmpResult;
    tmpResult.clear();
    for (int i = static_cast<int>(first.digits.size()) - 1; i >= 0; --i) {
        mod *= BASE;
        mod += first.digits[i];
        int left = 0;
        int right = BASE;
        while (right - left > 1) {
            int middle = (right + left) / 2;
            if (mod >= module.multipleOnInt(middle)) {
                left = middle;
            } else {
                right = middle;
            }
        }
        tmpResult.push_back(left);
        mod -= (module * left);
    }
    result.digits.clear();
    for (int i = static_cast<int>(tmpResult.size()) - 1; i >= 0; --i) {
        result.digits.push_back(tmpResult[i]);
    }
    result.sign = (first.sign == second.sign);
    result.deleteLeadingZeros();
}

BigInteger& BigInteger::operator/=(const BigInteger& number) {
    BigInteger result;
    BigInteger mod;
    BigInteger::div(*this, number, result, mod);
    (*this) = result;
    deleteLeadingZeros();
    return *this;
}

BigInteger& BigInteger::operator%=(const BigInteger& number) {
    BigInteger tmp = (*this);
    tmp /= number;
    (*this) -= tmp * number;
    return *this;
}

BigInteger operator/(const BigInteger& first, const BigInteger& second) {
    BigInteger tmp = first;
    tmp /= second;
    return tmp;
}

BigInteger operator%(const BigInteger& first, const BigInteger& second) {
    BigInteger tmp = first;
    tmp %= second;
    return tmp;
}


class Rational {
public:
    Rational(const BigInteger& n, const BigInteger& d = BigInteger(1));
    Rational(int n, int d = 1);
    Rational() : numerator(0), denominator(1) {}
    Rational(const Rational& fraction) : numerator(fraction.numerator), denominator(fraction.denominator) {}

    Rational operator-() const;
    std::string toString() const;
    std::string asDecimal(size_t precision = 0) const;
    explicit operator double() const;

    Rational& operator=(const Rational& fraction) = default;
    Rational& operator+=(const Rational& fraction);
    Rational& operator-=(const Rational& fraction);
    Rational& operator/=(const Rational& fraction);
    Rational& operator*=(const Rational& fraction);

    friend bool operator<(const Rational& first, const Rational& second);
    friend bool operator>(const Rational& first, const Rational& second);
    friend bool operator<=(const Rational& first, const Rational& second);
    friend bool operator>=(const Rational& first, const Rational& second);
    friend bool operator==(const Rational& first, const Rational& second);
    friend bool operator!=(const Rational& first, const Rational& second);

private:
    static BigInteger nodRecursive(BigInteger& first, BigInteger& second);
    static BigInteger nod(const BigInteger& first, const BigInteger& second);
    void reduct();

    BigInteger numerator;
    BigInteger denominator;
};

Rational::Rational(const BigInteger &n, const BigInteger &d) : numerator(n), denominator(d) {
    if (denominator < 0) {
        denominator = -denominator;
        numerator = -numerator;
    }
    reduct();
}

Rational::Rational(int n, int d) : numerator(n), denominator(d) {
    if (denominator < 0) {
        denominator = -denominator;
        numerator = -numerator;
    }
    reduct();
}

BigInteger Rational::nodRecursive(BigInteger& first, BigInteger& second) {
    if (first == 0 && second == 0) return BigInteger(1);
    if (first == 0) return second;
    if (second == 0) return first;
    if (first <= second) {
        second %= first;
    } else {
        first %= second;
    }
    return nodRecursive(first, second);
}

BigInteger Rational::nod(const BigInteger& first, const BigInteger& second) {
    BigInteger tFirst = BigInteger::abs(first);
    BigInteger tSecond = BigInteger::abs(second);
    return nodRecursive(tFirst, tSecond);
}

void Rational::reduct() {
    BigInteger bi_nod = nod(numerator, denominator);
    numerator /= bi_nod;
    denominator /= bi_nod;
}

std::string Rational::toString() const {
    std::string str = numerator.toString();
    if (denominator != 1 && numerator != BigInteger(0)) {
        str.push_back('/');
        std::string tmp = denominator.toString();
        for (char ch : tmp) {
            str.push_back(ch);
        }
    }
    return str;
}

bool operator<(const Rational &first, const Rational &second) {
    if (first.numerator.ifNegative() != second.numerator.ifNegative()) {
        return first.numerator.ifNegative();
    }
    return (first.numerator * second.denominator) < (first.denominator * second.numerator);
}

bool operator>(const Rational &first, const Rational &second) {
    if (first.numerator.ifNegative() != second.numerator.ifNegative()) {
        return !first.numerator.ifNegative();
    }
    return (first.numerator * second.denominator) > (first.denominator * second.numerator);
}

bool operator<=(const Rational &first, const Rational &second) {
    return !(first > second);
}

bool operator>=(const Rational &first, const Rational &second) {
    return !(first < second);
}

bool operator==(const Rational &first, const Rational &second) {
    return (first.numerator == second.numerator) && (first.denominator == second.denominator);
}

bool operator!=(const Rational &first, const Rational &second) {
    return !(first == second);
}

Rational Rational::operator-() const {
    Rational tmp = (*this);
    tmp.numerator = -tmp.numerator;
    return tmp;
}

Rational& Rational::operator+=(const Rational& fraction) {
    if (*this == fraction) {
        *this *= 2;
    } else {
        numerator *= fraction.denominator;
        numerator += (denominator * fraction.numerator);
        denominator *= fraction.denominator;
        reduct();
    }
    return *this;
}

Rational& Rational::operator-=(const Rational &fraction) {
    if (*this == fraction) {
        *this = Rational();
    } else {
        numerator *= fraction.denominator;
        numerator -= (denominator * fraction.numerator);
        denominator *= fraction.denominator;
        reduct();
    }
    return *this;
}

Rational& Rational::operator/=(const Rational& fraction) {
    if (*this == fraction) {
        *this = Rational(1);
    } else {
        numerator *= fraction.denominator;
        denominator *= fraction.numerator;
        if (denominator < 0) {
            numerator = -numerator;
            denominator = -denominator;
        }
        reduct();
    }
    return *this;
}

Rational &Rational::operator*=(const Rational &fraction) {
    numerator *= fraction.numerator;
    denominator *= fraction.denominator;
    reduct();
    return *this;
}

Rational operator+(const Rational& first, const Rational& second) {
    Rational tmp = first;
    tmp += second;
    return tmp;
}

Rational operator-(const Rational& first, const Rational& second) {
    Rational tmp = first;
    tmp -= second;
    return tmp;
}

Rational operator/(const Rational& first, const Rational& second) {
    Rational tmp = first;
    tmp /= second;
    return tmp;
}

Rational operator*(const Rational& first, const Rational& second) {
    Rational tmp = first;
    tmp *= second;
    return tmp;
}

std::string Rational::asDecimal(size_t precision) const {
    BigInteger tmp = BigInteger(1);
    for (size_t i = 0; i < precision; ++i) {
        tmp *= BigInteger(10);
    }
    tmp *= numerator;
    tmp /= denominator;
    std::string str = tmp.toString();
    std::string tempStr;
    for (size_t i = 0; i < precision; ++i) {
        if (!(str.empty() || (str.size() == 1 && str[0] == '-'))) {
            tempStr.push_back(str[str.size() - 1]);
            str.pop_back();
        } else tempStr.push_back('0');
    }
    if (str.empty() || (str.size() == 1 && str[0] == '-')) str.push_back('0');
    if (precision != 0) str.push_back('.');
    for (size_t i = 0; i < precision; ++i) {
        str.push_back(tempStr[tempStr.size() - 1]);
        tempStr.pop_back();
    }
    return str;
}

Rational::operator double() const {
    std::string str = (*this).asDecimal(1000);
    return std::stod(str);
}

#endif //BIGINT_BIGINTEGER_H