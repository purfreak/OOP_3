#include <algorithm>
#include <tuple>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <vector>

using std::cin;
using std::cout;
using std::vector;
using std::tuple;
using std::string;
using std::min;
using std::max;
using std::abs;
using std::swap;
using std::istream;
using std::ostream;
using std::numeric_limits;
using std::move;
using std::get;

constexpr int N = 3;


class Point {
public:
    double coordinates[N]{};

    Point() {
        memset(this->coordinates, 0, N * sizeof(double));
    }

    ~Point() = default;

    Point(const Point &point) {
        memcpy(this->coordinates, point.coordinates, N * sizeof(double));
    }

    Point(Point&& point) {
        memmove(this->coordinates, point.coordinates, N * sizeof(double));
    }

    explicit Point(double coords[]) {
        memcpy(this->coordinates, &coords[0], N * sizeof(double));
    }

    explicit Point(const vector<double> &coords) {
        memcpy(this->coordinates, &coords[0], N * sizeof(double));
    }

    Point(vector<double>&& coords) {
        memmove(this->coordinates, &coords[0], N * sizeof(double));
    }

    double operator[](int i) const {
        return coordinates[i];
    }

    bool operator==(const Point &point) const;

    bool operator!=(const Point &point) const;

    Point &operator=(const Point &point);

    Point operator-() const;

    Point operator+(const Point &point) const;

    Point operator-(const Point &point) const;

    Point operator*(double a) const;

    Point operator/(double a) const;

    Point &operator+=(const Point &point);

    Point &operator-=(const Point &point);

    Point &operator*=(double a);

    Point &operator/=(double a);
};


Point operator*(double a, const Point &base) {
    return base * a;
}

istream &operator>>(istream &in, Point &point) {
    for (auto &coord: point.coordinates) {
        in >> coord;
    }
    return in;
}

ostream &operator<<(ostream &out, const Point &point) {
    for (auto it = &point.coordinates[0]; it != &point.coordinates[N - 1]; ++it) {
        out << *it << " ";
    }
    out << point.coordinates[N - 1];

    return out;
}

Point &Point::operator=(const Point &point) {
    memcpy(this->coordinates, point.coordinates, N * sizeof(double));
    return *this;
}

bool are_equal(const double x, const double y) {
    return abs(x - y) <= numeric_limits<double>::epsilon();
}

bool Point::operator==(const Point &point) const {
    auto this_begin = &(this->coordinates[0]);

    for (const auto &x: point.coordinates) {
        if (!are_equal(*(this_begin++), x)) {
            return false;
        }
    }

    return true;
}

bool Point::operator!=(const Point &point) const {
    return !(*this == point);
}

Point Point::operator-() const {
    double neg_coords[N];

    auto neg_begin = &neg_coords[0];

    for (const auto &x: this->coordinates) {
        *(neg_begin++) = -x ?: 0;
    }

    return Point(neg_coords);
}

Point Point::operator+(const Point &point) const {
    auto result = *this;
    return result += point;
}

Point Point::operator-(const Point &point) const {
    auto result = *this;
    return result -= point;
}

Point Point::operator*(double a) const {
    auto result = *this;
    return result *= a;
}

Point Point::operator/(double a) const {
    auto result = *this;
    return result /= a;
}

Point &Point::operator+=(const Point &point) {
    auto this_begin = &(this->coordinates[0]);

    for (const auto &x: point.coordinates) {
        *(this_begin++) += x;
    }

    return *this;
}

Point &Point::operator-=(const Point &point) {
    return *this += (-point);
}

Point &Point::operator*=(double a) {
    auto end = &(this->coordinates[N]);
    ++end;
    for (auto this_begin = &(this->coordinates[0]); this_begin != end; ++this_begin) {
        *this_begin *= a;
    }

    return *this;
}

Point &Point::operator/=(double a) {
    auto end = &(this->coordinates[N]);
    ++end;

    for (auto this_begin = &(this->coordinates[0]); this_begin != end; ++this_begin) {
        *this_begin /= a;
    }

    return *this;
}


struct Vector : Point {
    Vector() = default;

    explicit Vector(const vector<double> &coordinates) : Point(coordinates) {}

    Vector(const Point &begin, const Point &end) : Point(end - begin) {}
};


Vector VectorProduct(const Vector &lhs, const Vector &rhs) {
    return Vector({lhs[1] * rhs[2] - lhs[2] * rhs[1],
                   -(lhs[0] * rhs[2] - lhs[2] * rhs[0]),
                   lhs[0] * rhs[1] - lhs[1] * rhs[0]});
}

using Shape = vector<Point>;

using Face = tuple<int, int, int>;

void MakeRightOrder(Face &face) {
    auto first = get<0>(face);
    auto second = get<1>(face);
    auto third = get<2>(face);

    if (second < first && second < third) {
        face = {second, third, first};
    } else if (third < first && third < second) {
        face = {third, first, second};
    }
}

ostream &operator<<(ostream &out, const Face &face) {
    return out << get<0>(face) << " " << get<1>(face) << " " << get<2>(face);
}


class ConvexHull {
private:
    struct NumberPoint : Point {
    public:
        int64_t ID = IDs_++;
        NumberPoint *prev = nullptr;
        NumberPoint *next = nullptr;

        static void ResetIds() {
            IDs_ = 0;
        }

        void Rotate(double angle);

        bool HoldMovie();

        NumberPoint() : Point() {}

        explicit NumberPoint(const Point &point) : Point(point) {}

    private:
        static void RotateAroundThirdDimension(double &x, double &y, const double &angle);

        static int IDs_;
    };

public:
    const double inf = 1e9;

    using Movie = NumberPoint;

    const vector<Face> &GetHull() {
        return hull_;
    }

    explicit ConvexHull(const vector<Point> &points);

private:
    vector<NumberPoint> points_;

    vector<Face> hull_;

    double Time(const NumberPoint *a, const NumberPoint *b, const NumberPoint *c);

    static bool IsAntiClockwise(const NumberPoint *a, const NumberPoint *b, const NumberPoint *c);

    static void BuildSupportingRib(NumberPoint *&u, NumberPoint *&v);

    vector<Movie *> MergeMovies(int64_t left_border, int64_t right_border);

    vector<Face> BuildConvexHull(bool negative);

    void BuildConvexHull();
};

int ConvexHull::NumberPoint::IDs_ = 0;

void ConvexHull::NumberPoint::RotateAroundThirdDimension(double &x, double &y, const double &angle) {
    const auto cos_angle = cos(angle);
    const auto sin_angle = sin(angle);

    double new_x, new_y;

    new_x = x * cos_angle + y * sin_angle;
    new_y = -x * sin_angle + y * cos_angle;

    x = new_x;
    y = new_y;
}

void ConvexHull::NumberPoint::Rotate(double angle) {
    RotateAroundThirdDimension(this->coordinates[1], this->coordinates[2], angle);

    RotateAroundThirdDimension(this->coordinates[0], this->coordinates[2], angle);

    RotateAroundThirdDimension(this->coordinates[0], this->coordinates[1], angle);
}

bool ConvexHull::NumberPoint::HoldMovie() {
    if (prev->next != this) {
        prev->next = next->prev = this;
        return true;
    } else {
        prev->next = next;
        next->prev = prev;
        return false;
    }
}

ConvexHull::ConvexHull(const vector<Point> &points) {
    points_.reserve(points.size());

    constexpr double angle_rad = 1e-2;

    for (auto &point : points) {
        points_.emplace_back(point);

        points_.back().Rotate(angle_rad);
    }

    BuildConvexHull();
    points_.back().ResetIds();
}

ostream &operator<<(ostream &out, ConvexHull &convex_hull) {
    out << convex_hull.GetHull().size() << "\n";
    for (const auto &face : convex_hull.GetHull()) {
        out << "3 " << face << "\n";
    }
    return out;
}

double
ConvexHull::Time(const NumberPoint *a, const NumberPoint *b, const NumberPoint *c) { // return time when abc is line
    if (a == nullptr || b == nullptr || c == nullptr) {
        return inf;
    }

    auto product = VectorProduct({*a, *b}, {*b, *c});
    return -product[1] / product[2];
}

bool ConvexHull::IsAntiClockwise(const NumberPoint *a, const NumberPoint *b, const NumberPoint *c) {
    if (a == nullptr || b == nullptr || c == nullptr) {
        return true;
    }

    return VectorProduct({*a, *b}, {*b, *c})[2] > 0;
}

void ConvexHull::BuildSupportingRib(NumberPoint *&u, NumberPoint *&v) {
    while (true) {
        if (!IsAntiClockwise(u->prev, u, v)) {
            u = u->prev;
        } else if (!IsAntiClockwise(u, v, v->next)) {
            v = v->next;
        } else {
            break;
        }
    }
}

vector<ConvexHull::Movie *> ConvexHull::MergeMovies(int64_t left_border, int64_t right_border) {
    if (right_border - left_border <= 1) {
        return vector<Movie *>();
    }

    const auto middle = (left_border + right_border) / 2;

    const vector<Movie *> left_hull = MergeMovies(left_border, middle);
    const vector<Movie *> right_hull = MergeMovies(middle, right_border);

    vector<Movie *> movies;

    auto u = &points_[middle - 1];
    auto v = &points_[middle];

    BuildSupportingRib(u, v);

    int64_t l_current = 0;
    int64_t r_current = 0;

    for (auto current_time = -inf; true;) {
        Movie *left = nullptr;
        Movie *right = nullptr;

        vector<double> times(6, inf);

        if (l_current < left_hull.size()) {
            left = left_hull[l_current];
            times[0] = Time(left->prev, left, left->next);
        }

        if (r_current < right_hull.size()) {
            right = right_hull[r_current];
            times[1] = Time(right->prev, right, right->next);
        }

        times[2] = Time(u->prev, u, v);
        times[3] = Time(u, u->next, v);
        times[4] = Time(u, v, v->next);
        times[5] = Time(u, v->prev, v);

        int64_t min_time_index = 0;
        auto min_time = times[min_time_index];

        for (int64_t i = 1; i < times.size(); ++i) {
            if (current_time < times[i] && times[i] < min_time) { // current_time?
                min_time_index = i;
                min_time = times[min_time_index];
            }
        }
        current_time = min_time;

        if (min_time >= inf) {
            break;
        }

        switch (min_time_index) {
            case 0:
                if (left->coordinates[0] < u->coordinates[0]) {
                    movies.emplace_back(left);
                }

                left->HoldMovie();
                ++l_current;
                break;

            case 1:
                if (right->coordinates[0] > v->coordinates[0]) {
                    movies.emplace_back(right);
                }

                right->HoldMovie();
                ++r_current;
                break;

            case 2:
                movies.emplace_back(u);
                u = u->prev;
                break;

            case 3:
                movies.emplace_back(u = u->next);
                break;

            case 4:
                movies.emplace_back(v);
                v = v->next;
                break;

            case 5:
                movies.emplace_back(v = v->prev);
                break;

            default:
                break;
        }
    }

    u->next = v;
    v->prev = u;

    for (auto it = movies.rbegin(); it != movies.rend(); ++it) {
        auto &current = *it;

        if (u->coordinates[0] < current->coordinates[0] && current->coordinates[0] < v->coordinates[0]) {
            u->next = v->prev = current;
            current->prev = u;
            current->next = v;

            if (current->coordinates[0] <= points_[middle - 1].coordinates[0]) {
                u = current;
            } else {
                v = current;
            }
        } else {
            current->HoldMovie();

            if (current == u) {
                u = u->prev;
            }

            if (current == v) {
                v = v->next;
            }
        }
    }

    return movies;
}

vector<Face> ConvexHull::BuildConvexHull(bool negative) {
    vector<Face> hull;

    vector<Movie *> movies = MergeMovies(0, points_.size());

    for (Movie *movie : movies) {
        auto face = Face({movie->prev->ID, movie->ID, movie->next->ID});

        if (movie->HoldMovie() == !negative) {
            swap(get<0>(face), get<1>(face));
        }

        hull.push_back(face);
    }

    return hull;
}

void ConvexHull::BuildConvexHull() {
    sort(points_.begin(), points_.end(),
         [](auto lhs, auto rhs) -> bool { return lhs.coordinates[0] < rhs.coordinates[0]; });

    hull_ = BuildConvexHull(true);

    for (NumberPoint &p : points_) {
        p.prev = nullptr;
        p.next = nullptr;
        p.coordinates[2] = -p.coordinates[2];
    }

    auto second_hull = BuildConvexHull(false);

    hull_.insert(hull_.end(), second_hull.begin(), second_hull.end());

    for (Face &triple : hull_) {
        MakeRightOrder(triple);
    }

    sort(hull_.begin(), hull_.end());
}


class Solver {
public:
    void operator()() {
        int m;
        cin >> m;

        for (int test = 0; test < m; ++test) {
            int n;
            cin >> n;

            vector<Point> points;
            points.reserve(n);

            for (int i = 0; i < n; ++i) {
                points.emplace_back();
                cin >> points.back();
            }

            ConvexHull convex_hull(points);
            cout << convex_hull << "\n";
        }
    }
};


int main() {
    Solver solver;
    solver();

    return 0;
}