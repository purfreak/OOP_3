#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

using std::cin;
using std::cout;
using std::vector;
using std::istream;
using std::begin;
using std::end;

constexpr double pi = 3.1415;


class Point {
public:
    Point() = default;

    Point(double x, double y) : x_(x), y_(y) {}

    Point(const Point &point) = default;

    double GetX() const;

    double GetY() const;

    void SetX(double x);

    void SetY(double y);

    friend istream &operator>>(istream &in, Point &point);

private:
    double x_ = 0;
    double y_ = 0;
};


class Polygon {
public:
    Polygon() = default;

    Polygon(const Polygon &polygon) = default;

    void PushBack(const Point &point);

    size_t GetSize() const;

    Point GetElement(size_t index) const;

    void TransformPolygon();

    bool IsPointInPolygon(const Point &point) const;

private:
    vector<Point> polygon_;
};


class MinkowskiSumm {
public:
    MinkowskiSumm() = default;

    MinkowskiSumm(Polygon &first_polygon, Polygon &second_polygon) : first_polygon_(first_polygon),
                                                                     second_polygon_(second_polygon) {}

    Polygon operator()();

private:
    Polygon first_polygon_;
    Polygon second_polygon_;
};

bool operator<(const Point &point1, const Point &point2) {
    if (point1.GetY() < point2.GetY()) {
        return true;
    }
    if (point1.GetY() == point2.GetY()) {
        return point1.GetX() < point2.GetX();
    }
    return false;
}

istream &operator>>(istream &in, Point &point) {
    return in >> point.x_ >> point.y_;
}

Point operator+(const Point &point1, const Point &point2) {
    return Point(point1.GetX() + point2.GetX(), point1.GetY() + point2.GetY());
}

double Point::GetX() const {
    return x_;
}

double Point::GetY() const {
    return y_;
}

void Point::SetX(double x) {
    x_ = x;
}

void Point::SetY(double y) {
    y_ = y;
}

void Polygon::PushBack(const Point &point) {
    polygon_.push_back(point);
}

size_t Polygon::GetSize() const {
    return polygon_.size();
}

Point Polygon::GetElement(size_t index) const {
    return polygon_[index];
}

void Polygon::TransformPolygon() {
    reverse(polygon_.begin(), polygon_.end());

    auto min_point = min_element(polygon_.begin(), polygon_.end());
    int left_lower_index = distance(polygon_.begin(), min_point);

    rotate(polygon_.begin(), polygon_.begin() + left_lower_index, polygon_.end());

    PushBack(GetElement(0));
}

double CalculatePolarAngle(const Point &point1, const Point &point2) {
    double angle = atan2(point2.GetY() - point1.GetY(), point2.GetX() - point1.GetX());
    if (angle < 0) {
        angle += 2 * pi;
    }
    return angle;
}

Polygon MinkowskiSumm::operator()() {
    size_t fist_polygon_size = first_polygon_.GetSize() - 1;
    size_t second_polygon_size = second_polygon_.GetSize() - 1;

    Polygon result;

    int first_index = 0;
    int second_index = 0;
    while (first_index < fist_polygon_size && second_index < second_polygon_size) {
        result.PushBack(first_polygon_.GetElement(first_index) + second_polygon_.GetElement(second_index));

        double first_polygon_angle = CalculatePolarAngle(first_polygon_.GetElement(first_index),
                                                         first_polygon_.GetElement(first_index + 1));
        double second_polygon_angle = CalculatePolarAngle(second_polygon_.GetElement(second_index),
                                                          second_polygon_.GetElement(second_index + 1));
        if (first_polygon_angle < second_polygon_angle) {
            first_index++;
        } else {
            if (first_polygon_angle > second_polygon_angle) {
                second_index++;
            } else {
                first_index++;
                second_index++;
            }
        }
    }
    result.PushBack(first_polygon_.GetElement(first_index) + second_polygon_.GetElement(second_index));
    result.PushBack(result.GetElement(0));

    return result;
}

bool Polygon::IsPointInPolygon(const Point &point) const {
    for (size_t i = 0; i < GetSize() - 1; ++i) {
        Point a = GetElement(i);
        Point b = GetElement(i + 1);

        if ((b.GetX() - a.GetX()) * (point.GetY() - b.GetY()) - (point.GetX() - b.GetX()) * (b.GetY() - a.GetY()) < 0) {
            return false;
        }
    }
    return true;
}


class Solver {
public:
    void operator()() {
        Polygon polygon1;
        Polygon polygon2;
        size_t n = 0;
        cin >> n;

        for (size_t i = 0; i < n; ++i) {
            Point point;
            cin >> point;
            polygon1.PushBack(point);
        }

        size_t m = 0;
        cin >> m;
        for (size_t i = 0; i < m; ++i) {
            Point point;
            cin >> point;
            double x = point.GetX();
            double y = point.GetY();
            point.SetX(-x);
            point.SetY(-y);
            polygon2.PushBack(point);
        }

        polygon1.TransformPolygon();
        polygon2.TransformPolygon();

        MinkowskiSumm sum(polygon1, polygon2);

        Polygon result_polygon = sum();

        Point point(0.0, 0.0);
        cout << (result_polygon.IsPointInPolygon(point) ? "YES" : "NO");
    }
};

int main() {
    Solver solver;
    solver();

    return 0;
}