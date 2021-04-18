//
// Created by Timur on 20.04.2020.
//

#ifndef GEOMETRY_GEOMETRY_H
#define GEOMETRY_GEOMETRY_H

#include <iostream>
#include <vector>
#include <cmath>

double EPS = 0.000000001;


class Line;

struct Point {
public:
    double x;
    double y;
    Point(double x, double y): x(x), y(y) {}
    Point(const Point &another) = default;

    void Rotate(const Point &center, double angle);
    void RotateCosForm(const Point &center, double cos);
    void Reflex(Point center);
    void Reflex(Line axis);
    void Scale(Point center, double coefficient);
};

bool operator==(const Point& first, const Point& second) {
    return (std::abs(first.x - second.x) < EPS) and (std::abs(first.y - second.y) < EPS);
}

bool operator!=(const Point& first, const Point& second) {
    return !(first==second);
}

void Point::Rotate(const Point &center, double angle) {
    x -= center.x;
    y -= center.y;
    double new_x = x * cos(angle) + y * sin(angle);
    double new_y = y * cos(angle) - x * sin(angle);
    x = new_x + center.x;
    y = new_y + center.y;
}

void Point::RotateCosForm(const Point &center, double cos) {
    x -= center.x;
    y -= center.y;
    double new_x = x * cos + y * std::sqrt(1.0 - cos * cos);
    double new_y = y * cos - x * std::sqrt(1.0 - cos * cos);
    x = new_x + center.x;
    y = new_y + center.y;
}

void Point::Reflex(Point center) {
    x = 2.0 * center.x - x;
    y = 2.0 * center.y - y;
}

void Point::Scale(Point center, double coefficient) {
    if (coefficient < 0) {
        coefficient = -coefficient;
        (*this).Reflex(center);
    }
    x -= center.x;
    y -= center.y;
    x *= coefficient;
    y *= coefficient;
    x += center.x;
    y += center.y;
}


class Line {
public:
    bool if_vertical;
    double coefficient;
    double shift;

    Line(double c, double s) : if_vertical(false), coefficient(c), shift(s) {}
    Line(const Point& p, double c) : if_vertical(false), coefficient(c), shift(p.y - p.x * c) {}
    Line(const Point& f, const Point& s);

    Line perpendicular(const Point& p) const;
};

Line::Line(const Point &f, const Point &s) {
    if (std::abs(f.x - s.x) < EPS) {
        if_vertical = true;
        coefficient = 1;
        shift = -f.x;
    } else {
        if_vertical = false;
        coefficient = (f.y - s.y) / (f.x - s.x);
        shift = f.y - f.x * coefficient;
    }
}

bool operator==(const Line& first, const Line& second) {
    if (first.if_vertical != second.if_vertical) return false;
    return (std::abs(first.coefficient - second.coefficient) < EPS) && (std::abs(first.shift - second.shift) < EPS);
}

bool operator!=(const Line& first, const Line& second) {
    return !(first == second);
}

Line Line::perpendicular(const Point &p) const {
    if (if_vertical) {
        return Line(0, p.y);
    } else {
        if (std::abs(coefficient) < EPS) {
            Line l(1, -p.x);
            l.if_vertical = true;
            return l;
        } else {
            return Line((-1) / coefficient, p.x / coefficient + p.y);
        }
    }
}

Point Intersection(const Line &first, const Line &second) {
    if (first.if_vertical)
        return Point(-first.shift, second.coefficient * (-first.shift) + second.shift);
    if (second.if_vertical)
        return Point(-second.shift, first.coefficient * (-second.shift) + first.shift);
    double x = (first.shift - second.shift) / (second.coefficient - first.coefficient);
    double y = first.coefficient * x + first.shift;
    return Point(x, y);
}

Point Projection(const Point &point, const Line &line) {
    if (line.if_vertical) return Point(line.shift, point.y);
    if (std::abs(line.coefficient) < EPS) return Point(point.x, line.shift);
    double x = (point.y + point.x / line.coefficient - line.shift) / (line.coefficient + 1 / line.coefficient);
    double y = line.coefficient * x + line.shift;
    return Point(x, y);
}

void Point::Reflex(Line axis) {
    if(axis.if_vertical) {
        x = 2.0 * (-axis.shift) - x;
    } else {
        if (std::abs(axis.coefficient) < EPS) {
            y = 2.0 * axis.shift - y;
        } else {
            double projection_x = (y + x / axis.coefficient - axis.shift) / (axis.coefficient + 1 / axis.coefficient);
            double projection_y = axis.coefficient * projection_x + axis.shift;
            (*this).Reflex(Point(projection_x, projection_y));
        }
    }
}


class Shape {
public:
    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool operator==(const Shape& another) const = 0;
    virtual bool operator!=(const Shape &another) const = 0;
    virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual bool isSimilarTo(const Shape& another) const = 0;

    virtual bool containsPoint(Point point) const = 0;

    virtual void rotate(Point center, double angle) = 0;
    virtual void reflex(Point center) = 0;
    virtual void reflex(Line axis) = 0;
    virtual void scale(Point center, double coefficient) = 0;

    virtual ~Shape() = default;
};


class Polygon : public Shape {
protected:
    std::vector<Point> vertexes;
    size_t vertexesCount;

public:
    int verticesCount() const;
    std::vector<Point> getVertices() const;

    Polygon(const std::vector<Point>& points);

    template<typename... Points>
    Polygon(Points... points);

    bool isConvex() const;

    double perimeter() const override;
    double area() const override;
    bool operator==(const Shape &another) const override;
    bool operator!=(const Shape &another) const override;
    bool isCongruentTo(const Shape& another) const override;
    bool isSimilarTo(const Shape& another) const override;
    bool containsPoint(Point point) const override;

    void rotate(Point center, double angle) override;
    void reflex(Point center) override;
    void reflex(Line axis) override;
    void scale(Point center, double coefficient) override;
};

Polygon::Polygon(const std::vector<Point> &points) : vertexesCount(points.size()) {
    for (Point point : points) {
        vertexes.push_back(point);
    }
}

int Polygon::verticesCount() const {
    return vertexesCount;
}

template<typename... Points>
Polygon::Polygon(Points... points) {
    for (auto p : std::initializer_list<Point>{points...}) {
        vertexes.push_back(p);
    }
    vertexesCount = vertexes.size();
}

std::vector<Point> Polygon::getVertices() const {
    return vertexes;
}

double PointsDist(const Point& first, const Point& second) {
    double length = (first.x - second.x) * (first.x - second.x);
    length += (first.y - second.y) * (first.y - second.y);
    return std::sqrt(length);
}

double Polygon::perimeter() const {
    double perimeter = 0;
    for (size_t i = 0; i < vertexes.size() - 1; ++i) {
        perimeter += PointsDist(vertexes[i], vertexes[i + 1]);
    }
    perimeter += PointsDist(vertexes[vertexes.size() - 1], vertexes[0]);
    return perimeter;
}

double Polygon::area() const {
    size_t n = vertexes.size();
    double area = 0;
    for (size_t i = 0; i < n - 1; ++i) {
        area += vertexes[i].x * vertexes[i + 1].y;
    }
    area += vertexes[n - 1].x * vertexes[0].y;
    for (size_t i = 0; i < n - 1; ++i) {
        area -= vertexes[i + 1].x * vertexes[i].y;
    }
    area -= vertexes[0].x * vertexes[n - 1].y;
    area = std::abs(area);
    area /= 2.0;
    return area;
}

bool Polygon::operator==(const Shape &another) const {
    const Polygon *p = dynamic_cast<const Polygon *>(&another);
    if (p == nullptr) return false;

    if (vertexesCount != p->vertexesCount) return false;
    bool if_found = false;
    size_t i = 0;
    while (!if_found && i < p->vertexesCount) {
        if (vertexes[0] == p->vertexes[i]) if_found = true;
        ++i;
    }
    --i;
    if (!if_found) return false;

    bool if_similar = true;
    for (size_t j = 0; j < vertexesCount; ++j) {
        if (i + j < p->vertexesCount) {
            if (vertexes[j] != p->vertexes[i + j]) if_similar = false;
        } else if (vertexes[j] != p->vertexes[i + j - p->vertexesCount]) if_similar = false;
    }
    if (if_similar) return true;

    if_similar = true;
    for (size_t j = 0; j < vertexesCount; ++j) {
        if (static_cast<int>(i) - static_cast<int>(j) >= 0) {
            if (vertexes[j] != p->vertexes[i - j]) if_similar = false;
        } else if (vertexes[j] != p->vertexes[i - j + p->vertexesCount]) if_similar = false;
    }
    return if_similar;
}

bool Polygon::operator!=(const Shape &another) const {
    return !((*this) == another);
}

bool Polygon::isSimilarTo(const Shape &another) const {
    const Polygon *p = dynamic_cast<const Polygon *>(&another);
    if (p == nullptr) return false;

    if (vertexesCount != p->vertexesCount) return false;

    double coef_of_similarity = std::sqrt(this->area() / p->area());

    bool if_found = false;
    for (size_t j = 0; j < vertexesCount; ++j) {
        bool if_current = true;
        for (size_t i = 1; i < vertexesCount; ++i) {
            size_t ind = j + i;
            if (ind >= vertexesCount) ind -= vertexesCount;

            double dist = PointsDist(vertexes[j], vertexes[ind]);
            if (std::abs(PointsDist(vertexes[0], vertexes[i]) - coef_of_similarity * dist) < EPS) {
                if_current = false;
                break;
            }
        }
        if (if_current) {
            if_found = true;
            break;
        }

        if_current = true;
        for (size_t i = 1; i < vertexesCount; ++i) {
            size_t ind = 0;
            if (static_cast<int>(j) - static_cast<int>(i) < 0)
                ind = j + vertexesCount - i;
            else ind = j - i;

            double dist = PointsDist(vertexes[j], vertexes[ind]);
            if (std::abs(PointsDist(vertexes[0], vertexes[i]) - coef_of_similarity * dist) < EPS) {
                if_current = false;
                break;
            }
        }
        if (if_current) {
            if_found = true;
            break;
        }
    }
    return if_found;
}

bool Polygon::isCongruentTo(const Shape &another) const {
    const Polygon *p = dynamic_cast<const Polygon *>(&another);
    if (p == nullptr) return false;

    if (vertexesCount != p->vertexesCount) return false;


    return (std::abs(this->area() - p->area()) < EPS) && this->isSimilarTo(*p);
}

bool direction(const Point &first, const Point &second, const Point &third) {
    double first_x = second.x - first.x;
    double first_y = second.y - first.y;
    double second_x = third.x - second.x;
    double second_y = third.y - second.y;
    return (first_x * second_y - first_y * second_x) >= 0;
}

bool Polygon::isConvex() const {
    size_t n = vertexesCount;

    bool dir = direction(vertexes[0], vertexes[1], vertexes[2]);
    for (size_t i = 0; i < n - 2; ++i) {
        if (dir != direction(vertexes[i], vertexes[i + 1], vertexes[i + 2])) return false;
    }
    if (dir != direction(vertexes[n - 2], vertexes[n - 1], vertexes[0])) return false;
    return (dir == direction(vertexes[n - 1], vertexes[0], vertexes[1]));
}

double getVectMultiple(const Point &start, const Point &first_end, const Point &second_end) {
    double first_x = first_end.x - start.x;
    double first_y = first_end.y - start.y;
    double second_x = second_end.x - start.x;
    double second_y = second_end.y - start.y;
    return (first_x * second_y - second_x * first_y) / 2.0;
}

bool Polygon::containsPoint(Point point) const {
    bool result = false;
    size_t j = vertexesCount - 1;
    for (size_t i = 0; i < vertexesCount; i++) {
        if (((vertexes[i].y < point.y && vertexes[j].y >= point.y) || (vertexes[j].y < point.y && vertexes[i].y >= point.y)) &&
             (vertexes[i].x + (point.y - vertexes[i].y) / (vertexes[j].y - vertexes[i].y) * (vertexes[j].x - vertexes[i].x) < point.x))
            result = !result;
        j = i;
    }
    return result;
}

void Polygon::rotate(Point center, double angle) {
    for (size_t i = 0; i < vertexesCount; ++i) {
        vertexes[i].Rotate(center, angle);
    }
}

void Polygon::reflex(Point center) {
    for (size_t i = 0; i < vertexesCount; ++i) {
        vertexes[i].Reflex(center);
    }
}

void Polygon::reflex(Line axis) {
    for (size_t i = 0; i < vertexesCount; ++i) {
        vertexes[i].Reflex(axis);
    }
}

void Polygon::scale(Point center, double coefficient) {
    for (size_t i = 0; i < vertexesCount; ++i) {
        vertexes[i].Scale(center, coefficient);
    }
}


class Rectangle : public Polygon {
protected:
    Point center;

public:
    Rectangle(const Point &first, const Point &second, double coef);
};

Rectangle::Rectangle(const Point &first, const Point &second, double coef) : Polygon(first), center(Point(0, 0)) {
    vertexesCount = 4;
    center = Point((first.x + second.x) / 2.0, (first.y + second.y) /2.0);
    if (coef > 1) coef = 1 / coef;

    double diagonal = PointsDist(first, second);
    double little_edge = diagonal * coef / std::sqrt(1 + coef * coef);
    double cos = 1 - little_edge * little_edge * 2.0 / (diagonal * diagonal);

    Point additional(first);
    additional.RotateCosForm(center, cos);
    vertexes.push_back(additional);

    vertexes.push_back(second);

    vertexes.emplace_back(2.0 * center.x - additional.x, 2.0 * center.y - additional.y);
}


class Ellipse : public Shape {
protected:
    std::pair<Point, Point> focus;
    double distances_sum;

    double focus_dist;
public:
    Ellipse(const Point& first, const Point& second, double d) : focus (std::make_pair(first, second)), distances_sum(d) {
        focus_dist = PointsDist(focus.first, focus.second) / 2.0;
    }

    std::pair<Point, Point> focuses();
    std::pair<Line, Line> directrices();
    double eccentricity();
    Point center();

    double perimeter() const override;
    double area() const override;
    bool operator==(const Shape &another) const override;
    bool operator!=(const Shape &another) const override;
    bool isCongruentTo(const Shape& another) const override;
    bool isSimilarTo(const Shape& another) const override;
    bool containsPoint(Point point) const override;

    void rotate(Point center, double angle) override;
    void reflex(Point center) override;
    void reflex(Line axis) override;
    void scale(Point center, double coefficient) override;
};

Point Ellipse::center() {
    double x = (focus.first.x + focus.second.x) / 2.0;
    double y = (focus.first.y + focus.second.y) / 2.0;
    return Point(x, y);
}

std::pair<Point, Point> Ellipse::focuses() {
    return focus;
}

double Ellipse::eccentricity() {
    return 2.0 * focus_dist / distances_sum;
}

std::pair<Line, Line> Ellipse::directrices() {
    Point center = this->center();
    double dist_coef = (distances_sum * distances_sum) / (4.0 * focus_dist);

    double x_shift = (focus.first.x - focus.second.x) / (2.0 * focus_dist);
    x_shift *= dist_coef;
    double y_shift = (focus.first.y - focus.second.y) / (2.0 * focus_dist);
    y_shift *= dist_coef;

    Point first(center.x + x_shift, center.y + y_shift);
    Point second(center.x - x_shift, center.y - y_shift);
    Line l(focus.first, focus.second);

    return std::make_pair(l.perpendicular(first), l.perpendicular(second));
}

double Ellipse::perimeter() const {
    double a = distances_sum / 2.0;
    double b = std::sqrt(a * a - focus_dist * focus_dist);
    return M_PI * (3 * (a + b) - std::sqrt((3 * a + b) * (a + 3 * b)));
}

double Ellipse::area() const {
    double a = distances_sum / 2.0;
    double b = std::sqrt(a * a - focus_dist * focus_dist);
    return M_PI * a * b;
}

bool Ellipse::operator==(const Shape &another) const {
    const Ellipse *el = dynamic_cast<const Ellipse *>(&another);
    if (el == nullptr) return false;

    bool if_focuses_similar = false;
    if (this->focus.first == el->focus.first && this->focus.second == el->focus.second) {
        if_focuses_similar = true;
    }
    if (this->focus.first == el->focus.second && this->focus.second == el->focus.first) {
        if_focuses_similar = true;
    }
    return if_focuses_similar && (std::abs(this->distances_sum - el->distances_sum) < EPS);
}

bool Ellipse::operator!=(const Shape &another) const {
    return !((*this) == another);
}

bool Ellipse::isCongruentTo(const Shape &another) const {
    const Ellipse *el = dynamic_cast<const Ellipse *>(&another);
    if (el == nullptr) return false;

    return (std::abs(this->focus_dist - el->focus_dist) < EPS) && (std::abs(this->distances_sum - el->distances_sum) < EPS);
}

bool Ellipse::isSimilarTo(const Shape &another) const {
    const Ellipse *el = dynamic_cast<const Ellipse *>(&another);
    if (el == nullptr) return false;

    if (this->focus_dist == 0) {
        return (el->focus_dist == 0);
    } else {
        return std::abs((el->focus_dist / this->focus_dist) - (el->distances_sum / this->distances_sum)) < EPS;
    }
}

bool Ellipse::containsPoint(Point point) const {
    double distance = 0;
    distance += PointsDist(point, focus.first);
    distance += PointsDist(point, focus.second);
    return (distance < distances_sum);
}

void Ellipse::rotate(Point center, double angle) {
    focus.first.Rotate(center, angle);
    focus.second.Rotate(center, angle);
}

void Ellipse::reflex(Point center) {
    focus.first.Reflex(center);
    focus.second.Reflex(center);
}

void Ellipse::reflex(Line axis) {
    focus.first.Reflex(axis);
    focus.second.Reflex(axis);
}

void Ellipse::scale(Point center, double coefficient) {
    focus.first.Scale(center, coefficient);
    focus.second.Scale(center, coefficient);
    distances_sum *= std::abs(coefficient);
    focus_dist *= std::abs(coefficient);
}


class Circle : public Ellipse {
public:
    Circle(const Point& center, double radius) : Ellipse(center, center, 2.0 * radius) {}

    double radius();
};

double Circle::radius() {
    return distances_sum / 2.0;
}


class Square : public Rectangle {
public:
    Square(const Point &first, const Point &second) : Rectangle(first, second, 1) {}

    Circle circumscribedCircle();
    Circle inscribedCircle();
};

Circle Square::circumscribedCircle() {
    return Circle(center, PointsDist(vertexes[0], vertexes[2]) / 2);
}

Circle Square::inscribedCircle() {
    return Circle(center, PointsDist(vertexes[0], vertexes[1]) / 2);
}


class Triangle : public Polygon {
public:
    Triangle(const Point &first, const Point &second, const Point &third) : Polygon(first, second, third) {}

    Circle circumscribedCircle();
    Circle inscribedCircle();
    Point centroid();
    Point orthocenter();
    Line EulerLine();
    Circle ninePointsCircle();
};

Circle Triangle::circumscribedCircle() {
    Point middle = Point((vertexes[0].x + vertexes[1].x) / 2.0, (vertexes[0].y + vertexes[1].y) / 2.0);
    Line first_perp(Line(vertexes[0], vertexes[1]).perpendicular(middle));

    middle = Point((vertexes[1].x + vertexes[2].x) / 2.0, (vertexes[1].y + vertexes[2].y) / 2.0);
    Line second_perp(Line(vertexes[1], vertexes[2]).perpendicular(middle));

    Point center = Intersection(first_perp, second_perp);
    return Circle(center, PointsDist(vertexes[0], center));
}

Point Triangle::centroid() {
    Point middle = Point((vertexes[0].x + vertexes[1].x) / 2, (vertexes[0].y + vertexes[1].y) / 2);
    Line first_median(vertexes[2], middle);
    middle = Point((vertexes[1].x + vertexes[2].x) / 2, (vertexes[1].y + vertexes[2].y) / 2);
    Line second_median(vertexes[0], middle);
    return Intersection(first_median, second_median);
}

Point Triangle::orthocenter() {
    Point first_proj = Projection(vertexes[0], Line(vertexes[1], vertexes[2]));
    Point second_proj = Projection(vertexes[1], Line(vertexes[0], vertexes[2]));
    return Intersection(Line(vertexes[0], first_proj), Line(vertexes[1], second_proj));
}

Line Triangle::EulerLine() {
    return Line((*this).centroid(), (*this).orthocenter());
}

Circle Triangle::ninePointsCircle() {
    Point f_middle = Point((vertexes[0].x + vertexes[1].x) / 2, (vertexes[0].y + vertexes[1].y) / 2);
    Point s_middle = Point((vertexes[1].x + vertexes[2].x) / 2, (vertexes[1].y + vertexes[2].y) / 2);
    Point t_middle = Point((vertexes[0].x + vertexes[2].x) / 2, (vertexes[0].y + vertexes[2].y) / 2);
    return Triangle(f_middle, s_middle, t_middle).circumscribedCircle();
}

Circle Triangle::inscribedCircle() {
    double f = PointsDist(vertexes[0], vertexes[1]);
    double s = PointsDist(vertexes[1], vertexes[2]);
    double t = PointsDist(vertexes[2], vertexes[0]);
    double x_cord = (f * vertexes[2].x + s * vertexes[0].x + t * vertexes[1].x) / (f + s + t);
    double y_cord = (f * vertexes[2].y + s * vertexes[0].y + t * vertexes[1].y) / (f + s + t);

    Point center(x_cord, y_cord);
    double radius = this->area() * 2.0 / (f + s + t);
    return Circle(center, radius);
}

#endif //GEOMETRY_GEOMETRY_H
