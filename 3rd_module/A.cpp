#include <iostream>
#include <cmath>

using std::cin;
using std::cout;
using std::sqrt;
using std::istream;
using std::min;

constexpr long double precision = 1e-10;

class Point;

class Vector;

class Point {
public:
    Point() = default;

    Point(const Point &point) = default;

    Point(long double x, long double y, long double z) : x_(x), y_(y), z_(z) {}

    long double GetX() const;

    long double GetY() const;

    long double GetZ() const;

    Point operator+(const Vector &vector) const;

    friend istream &operator>>(istream &in, Point &point);

private:
    long double x_ = 0;
    long double y_ = 0;
    long double z_ = 0;
};


class Vector {
public:
    Vector() = default;

    Vector(const Vector &vector) = default;

    Vector(const Point &point1, const Point &point2) : first_point_(point1), second_point_(point2) {}

    Vector(long double x, long double y, long double z) {
        first_point_ = Point(0, 0, 0);
        second_point_ = Point(x, y, z);
    }

    long double GetX() const;

    long double GetY() const;

    long double GetZ() const;

    const Point &GetFirstPoint() const;

    const Point &GetSecondPoint() const;

    void SetFirstPoint(const Point &first_point);

    void SetSecondPoint(const Point &second_point);

    long double GetLength() const;

    Vector operator*(long double number);

    long double operator*(const Vector &vec);

private:
    Point first_point_;
    Point second_point_;
};


long double Point::GetX() const {
    return x_;
}

long double Point::GetY() const {
    return y_;
}

long double Point::GetZ() const {
    return z_;
}

Point Point::operator+(const Vector &vector) const {
    return Point(x_ + vector.GetX(), y_ + vector.GetY(), z_ + vector.GetZ());
}

istream &operator>>(istream &in, Point &point) {
    return in >> point.x_ >> point.y_ >> point.z_;
}


long double Vector::GetX() const {
    return second_point_.GetX() - first_point_.GetX();
}

long double Vector::GetY() const {
    return second_point_.GetY() - first_point_.GetY();
}

long double Vector::GetZ() const {
    return second_point_.GetZ() - first_point_.GetZ();
}

const Point &Vector::GetFirstPoint() const {
    return first_point_;
}

const Point &Vector::GetSecondPoint() const {
    return second_point_;
}

void Vector::SetFirstPoint(const Point &first_point) {
    first_point_ = first_point;
}

void Vector::SetSecondPoint(const Point &second_point) {
    second_point_ = second_point;
}

long double Vector::GetLength() const {
    return sqrt(GetX() * GetX() + GetY() * GetY() + GetZ() * GetZ());
}

Vector Vector::operator*(long double number) {
    long double tmp_x = GetX() * number;
    long double tmp_y = GetY() * number;
    long double tmp_z = GetZ() * number;

    return Vector(tmp_x, tmp_y, tmp_z);
}

long double Vector::operator*(const Vector &vector) {
    return GetX() * vector.GetX() + GetY() * vector.GetY() + GetZ() * vector.GetZ();
}


long double CalculateDistance(const Point &first_point, const Point &second_point) {
    return Vector(first_point, second_point).GetLength();
}

long double CalculateDistance(const Vector &vector, const Point &point) {
    Vector left_side = Vector(point, vector.GetFirstPoint());
    Vector right_side = Vector(point, vector.GetSecondPoint());

    if ((left_side * vector) * (right_side * vector) > 0) {
        return min(CalculateDistance(vector.GetFirstPoint(), point),
                   CalculateDistance(vector.GetSecondPoint(), point));
    }

    long double tmp_x =
            left_side.GetY() * right_side.GetZ() - left_side.GetZ() * right_side.GetY();
    long double tmp_y =
            left_side.GetZ() * right_side.GetX() - left_side.GetX() * right_side.GetZ();
    long double tmp_z =
            left_side.GetX() * right_side.GetY() - left_side.GetY() * right_side.GetX();
    long double square = Vector(tmp_x, tmp_y, tmp_z).GetLength();

    long double distance = square / vector.GetLength();

    return distance;
}

long double CalculateDistance(Vector &first_vector, Vector &second_vector) {
    if (first_vector.GetLength() == 0) {
        if (second_vector.GetLength() == 0) {
            return CalculateDistance(first_vector.GetFirstPoint(), second_vector.GetFirstPoint());
        } else {
            return CalculateDistance(second_vector, first_vector.GetFirstPoint());
        }
    }

    if (second_vector.GetLength() == 0) {
        return CalculateDistance(first_vector, second_vector.GetFirstPoint());
    }

    while (second_vector.GetLength() > precision) {
        Point left_border = second_vector.GetFirstPoint() + (second_vector * (1. / 3.));
        Point right_border = second_vector.GetFirstPoint() + (second_vector * (2. / 3.));

        if (CalculateDistance(first_vector, left_border) < CalculateDistance(first_vector, right_border)) {
            second_vector.SetSecondPoint(right_border);
        } else {
            second_vector.SetFirstPoint(left_border);
        }
    }

    return CalculateDistance(first_vector, second_vector.GetFirstPoint());
}


class Solver {
public:
    void operator()() {
        Point point1;
        Point point2;
        cin >> point1;
        cin >> point2;
        Vector first_vector(point1, point2);

        Point point3;
        Point point4;
        cin >> point3;
        cin >> point4;
        Vector second_vector(point3, point4);

        cout.precision(10);
        cout << CalculateDistance(first_vector, second_vector);
    }
};


int main() {
    Solver solver;
    solver();

    return 0;
}

