#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>


// Найти расстояние между двумя отрезками в трехмерном пространстве

static long double EPS = 0.000000001;

template<typename Type>
class Vector;

template<typename Type>
class Vertex {
public:
    Vertex(Type x, Type y, Type z) : x_cord(x), y_cord(y), z_cord(z){}
    Vertex(const Vertex<long long> &v) : x_cord(v.x_cord), y_cord(v.y_cord), z_cord(v.z_cord){}

    static long double Distance(const Vertex<Type>& first, const Vertex<Type>& second);
    static Type SquaredDistance(const Vertex<Type>& first, const Vertex<Type>& second);
    Vertex operator*(Type number) const;
    Vertex operator+(const Vector<Type>& vector) const;

    Type x_cord;
    Type y_cord;
    Type z_cord;
};

template<typename Type>
Type Vertex<Type>::SquaredDistance(const Vertex<Type> &first, const Vertex<Type> &second) {
    Type squared_dist = 0;
    squared_dist += (first.x_cord - second.x_cord) * (first.x_cord - second.x_cord);
    squared_dist += (first.y_cord - second.y_cord) * (first.y_cord - second.y_cord);
    squared_dist += (first.z_cord - second.z_cord) * (first.z_cord - second.z_cord);
    return squared_dist;
}

template<typename Type>
long double Vertex<Type>::Distance(const Vertex<Type> &first, const Vertex<Type> &second) {
    return std::sqrt(SquaredDistance(first, second));
}

template<typename Type>
Vertex<Type> Vertex<Type>::operator*(Type number) const {
    return Vertex(x_cord * number, y_cord * number, z_cord * number);
}

template<typename Type>
class Vector {
public:
    Vector(Type x, Type y, Type z) : x_cord(x), y_cord(y), z_cord(z){}
    Vector(const Vertex<Type> &ver) : x_cord(ver.x_cord), y_cord(ver.y_cord), z_cord(ver.z_cord) {}
    Vector(const Vertex<Type> &begin, const Vertex<Type> &end);
    Vector<Type>(const Vector<long long> &v) : x_cord(v.x_cord), y_cord(v.y_cord), z_cord(v.z_cord) {}

    static Vector<Type> VectorMultiple(const Vector<Type> &first, const Vector<Type> &second);
    static Type ScalarMultiple(const Vector<Type> &first, const Vector<Type> &second);
    bool IsZero() const;
    long double Length() const;
    Type SquaredLength() const;
    Vector<Type> operator*(Type number) const;

    template<typename NumberType>
    Vector<long double> operator/(NumberType number) const {
        return Vertex<long double>(x_cord / number, y_cord / number, z_cord / number);
    }

    Type x_cord;
    Type y_cord;
    Type z_cord;
};

template<typename Type>
Vector<Type>::Vector(const Vertex<Type> &begin, const Vertex<Type> &end) {
    x_cord = end.x_cord - begin.x_cord;
    y_cord = end.y_cord - begin.y_cord;
    z_cord = end.z_cord - begin.z_cord;
}

template<typename Type>
Vector<Type> Vector<Type>::VectorMultiple(const Vector<Type> &first, const Vector<Type> &second) {
    Type new_x_cord = (first.y_cord * second.z_cord - second.y_cord * first.z_cord);
    Type new_y_cord = (first.z_cord * second.x_cord - first.x_cord * second.z_cord);
    Type new_z_cord = (first.x_cord * second.y_cord - first.y_cord * second.x_cord);
    return Vector(new_x_cord, new_y_cord, new_z_cord);
}

template<typename Type>
Type Vector<Type>::ScalarMultiple(const Vector<Type> &first, const Vector<Type> &second) {
    return first.x_cord * second.x_cord + first.y_cord * second.y_cord + first.z_cord * second.z_cord;
}

template<typename Type>
bool Vector<Type>::IsZero() const {
    return (std::abs(x_cord) < EPS) && (std::abs(y_cord) < EPS) && (std::abs(z_cord) < EPS);
}

template<typename Type>
Type Vector<Type>::SquaredLength() const {
    return Vertex<Type>::SquaredDistance(Vertex<Type>(x_cord, y_cord, z_cord), Vertex<Type>(0, 0, 0));
}

template<typename Type>
long double Vector<Type>::Length() const {
    return std::sqrt(SquaredLength());
}

template<typename Type>
Vector<Type> Vector<Type>::operator*(Type number) const {
    return Vector<Type>(x_cord * number, y_cord * number, z_cord * number);
}

template<typename Type>
Vertex<Type> Vertex<Type>::operator+(const Vector<Type> &vector) const {
    Type x = x_cord + vector.x_cord;
    Type y = y_cord + vector.y_cord;
    Type z = z_cord + vector.z_cord;
    return Vertex<Type>(x, y, z);
}

template<typename Type>
struct LineSegment {
public:
    LineSegment(const Vertex<Type> &first, const Vertex<Type> &second) : first(first), second(second) {}
    LineSegment(const LineSegment<long long> &another) : first(another.first), second(another.second) {}

    template<typename VertType>
    bool ContainVertex(const Vertex<VertType> &vertex) const {
        bool is_contain = IsBetween(first.x_cord, second.x_cord, vertex.x_cord);
        is_contain = is_contain && IsBetween(first.y_cord, second.y_cord, vertex.y_cord);
        is_contain = is_contain && IsBetween(first.z_cord, second.z_cord, vertex.z_cord);
        return is_contain;
    }

    Vertex<Type> first;
    Vertex<Type> second;
private:
    template<typename NumberType>
    static bool IsBetween(Type f_border, Type s_border, NumberType number) {
        return (f_border <= number + EPS && number - EPS <= s_border)
               || (s_border <= number + EPS && number - EPS <= f_border);
    }
};


long long SquaredMinDistInEnds(const LineSegment<long long> &first_seg,
                               const LineSegment<long long> &second_seg) {
    long long min = Vertex<long long>::SquaredDistance(first_seg.first, second_seg.first);
    min = std::min(min, Vertex<long long>::SquaredDistance(first_seg.second, second_seg.first));
    min = std::min(min, Vertex<long long>::SquaredDistance(first_seg.first, second_seg.second));
    min = std::min(min, Vertex<long long>::SquaredDistance(first_seg.second, second_seg.second));
    return min;
}

long double SquaredDistOnLine(const LineSegment<long long> &first_seg,
                              const LineSegment<long long> &second_seg) {
    bool is_crossing = first_seg.ContainVertex(second_seg.first);
    is_crossing = is_crossing || first_seg.ContainVertex(second_seg.second);
    is_crossing = is_crossing || second_seg.ContainVertex(first_seg.first);
    is_crossing = is_crossing || second_seg.ContainVertex(first_seg.second);
    if (is_crossing) {
        return 0;
    } else {
        return SquaredMinDistInEnds(first_seg, second_seg);
    }
}

bool IsProjectionOnSegment(const LineSegment<long long> &first_seg, const Vertex<long long> &vertex) {
    Vector between_vector(first_seg.first, vertex);
    Vector direction_vector = Vector(first_seg.first, first_seg.second);
    long long squared_length = direction_vector.SquaredLength();

    long long alpha = Vector<long long>::ScalarMultiple(direction_vector, between_vector);
    Vertex<long double> new_vertex = Vertex<long double>(first_seg.first)
                                     + Vector<long double>(direction_vector) * alpha / squared_length;

    return first_seg.ContainVertex(new_vertex);
}

long double SquaredLengthOfProjectionOnLine(const LineSegment<long long> &first_seg,
                                            const Vertex<long long> &vertex) {
    Vector between_vector(first_seg.first, vertex);
    Vector direction_vector = Vector(first_seg.first, first_seg.second);
    long long squared_length = direction_vector.SquaredLength();
    long long dist = Vector<long long>::VectorMultiple(between_vector, direction_vector).SquaredLength();
    return (long double)dist / squared_length;
}

long double SquaredDistOnParallel(const LineSegment<long long> &first_seg,
                                  const LineSegment<long long> &second_seg) {
    if (IsProjectionOnSegment(first_seg, second_seg.first)
        || IsProjectionOnSegment(first_seg, second_seg.second)) {
        return SquaredLengthOfProjectionOnLine(first_seg, second_seg.first);
    } else if (IsProjectionOnSegment(second_seg, first_seg.first)
               || IsProjectionOnSegment(second_seg, first_seg.second)){
        return SquaredLengthOfProjectionOnLine(second_seg, first_seg.first);
    } else {
        return SquaredMinDistInEnds(first_seg, second_seg);
    }
}

template<typename Type>
Vertex<long double> CrossingVertex(const LineSegment<Type> &first_seg,
                                   const LineSegment<Type> &second_seg) {
    Vector<Type> first_vec(first_seg.first, first_seg.second);
    Vector<Type> second_vec(second_seg.first, second_seg.second);
    Vector<Type> btw_vector(second_seg.first, first_seg.first);

    Type module = first_vec.x_cord * second_vec.y_cord - first_vec.y_cord * second_vec.x_cord;
    Type coff = btw_vector.y_cord * second_vec.x_cord - btw_vector.x_cord * second_vec.y_cord;
    if (module == 0) {
        module = first_vec.z_cord * second_vec.y_cord - first_vec.y_cord * second_vec.z_cord;
        coff = btw_vector.y_cord * second_vec.z_cord - btw_vector.z_cord * second_vec.y_cord;
        if (module == 0) {
            module = first_vec.x_cord * second_vec.z_cord - first_vec.z_cord * second_vec.x_cord;
            coff = btw_vector.z_cord * second_vec.x_cord - btw_vector.x_cord * second_vec.z_cord;
        }
    }
    if (module < 0) {
        module = -module;
        coff = -coff;
    }
    Vertex<long double> vertex = Vertex<long double>(first_seg.first) + Vector<long double>(first_vec) * coff / module;
    return vertex;
}

long double SquaredDistIfNotCrossing(const LineSegment<long long> &first_seg,
                                     const LineSegment<long long> &second_seg) {
    auto min = (long double)SquaredMinDistInEnds(first_seg, second_seg);
    if (IsProjectionOnSegment(first_seg, second_seg.first)) {
        min = std::min(min, SquaredLengthOfProjectionOnLine(first_seg, second_seg.first));
    }
    if (IsProjectionOnSegment(first_seg, second_seg.second)) {
        min = std::min(min, SquaredLengthOfProjectionOnLine(first_seg, second_seg.second));
    }
    if (IsProjectionOnSegment(second_seg, first_seg.first)) {
        min = std::min(min, SquaredLengthOfProjectionOnLine(second_seg, first_seg.first));
    }
    if (IsProjectionOnSegment(second_seg, first_seg.second)) {
        min = std::min(min, SquaredLengthOfProjectionOnLine(second_seg, first_seg.second));
    }
    return min;
}

long double SquaredDistOnCrossing(const LineSegment<long long> &first_seg,
                                  const LineSegment<long long> &second_seg) {
    Vertex<long double> cross_vertex = CrossingVertex<long long>(first_seg, second_seg);

    if (first_seg.ContainVertex(cross_vertex) && second_seg.ContainVertex(cross_vertex)) {
        return 0;
    }
    return SquaredDistIfNotCrossing(first_seg, second_seg);
}

Vertex<long double> ProjectionOnPlane(const Vertex<long long> &vertex, const Vertex<long long> &vertex_on_plane,
                                      const Vector<long long> &normal_vector) {
    Vector between_vector(vertex, vertex_on_plane);
    long long alpha = Vector<long long>::ScalarMultiple(between_vector, normal_vector);
    return Vertex<long double>(vertex) + (Vector<long double>)normal_vector * alpha / normal_vector.SquaredLength();
}

long double DistOnOther(const LineSegment<long long> &first_seg,
                        const LineSegment<long long> &second_seg) {
    Vector first_vec(first_seg.first, first_seg.second);
    Vector second_vec(second_seg.first, second_seg.second);
    Vector normal_vector = Vector<long long>::VectorMultiple(first_vec, second_vec);

    Vector between_vector(first_seg.first, second_seg.first);
    long long scalar_mult = std::abs(Vector<long long>::ScalarMultiple(between_vector, normal_vector));
    long double dist_btw_planes = (long double)(scalar_mult) / normal_vector.Length();

    Vertex<long double> first_proj = ProjectionOnPlane(first_seg.first, second_seg.first, normal_vector);
    Vertex<long double> second_proj = ProjectionOnPlane(first_seg.second, second_seg.first, normal_vector);

    Vertex<long double> crossing_vertex = CrossingVertex<long double>
            (LineSegment<long double>(first_proj, second_proj), LineSegment<long double>(second_seg));
    if (LineSegment<long double>(first_proj, second_proj).ContainVertex(crossing_vertex)
        && second_seg.ContainVertex(crossing_vertex)) {
        return dist_btw_planes;
    }
    return std::sqrt(SquaredDistIfNotCrossing(first_seg, second_seg));
}

long double SquaredDistFromVertexToSegment(const LineSegment<long long> &segment,
                                           const Vertex<long long> &vertex) {
    if (IsProjectionOnSegment(segment, vertex)) {
        return SquaredLengthOfProjectionOnLine(segment, vertex);
    }
    return SquaredMinDistInEnds(LineSegment(vertex, vertex), segment);
}

long double FindDistBetweenSegments(const LineSegment<long long> &first_seg,
                                    const LineSegment<long long> &second_seg) {
    Vector first_vec(first_seg.first, first_seg.second);
    Vector second_vec(second_seg.first, second_seg.second);

    Vector normal_vector = Vector<long long>::VectorMultiple(first_vec, second_vec);
    Vector between_vector(second_seg.first, first_seg.first);
    if (normal_vector.IsZero()) {
        // один из векторов нулевой
        if (first_vec.IsZero() || second_vec.IsZero()) {
            if (!first_vec.IsZero()) {
                return std::sqrt(SquaredDistFromVertexToSegment(first_seg, second_seg.first));
            }
            if (!second_vec.IsZero()) {
                return std::sqrt(SquaredDistFromVertexToSegment(second_seg, first_seg.first));
            }
            return Vector(first_seg.first, second_seg.first).Length();
        }
        // отрезки лежат на одной // или // на параллельных прямых //
        if (Vector<long long>::VectorMultiple(between_vector, first_vec).IsZero() &&
            !first_vec.IsZero() && !second_vec.IsZero()) {
            return std::sqrt(SquaredDistOnLine(first_seg, second_seg));
        } else {
            return std::sqrt(SquaredDistOnParallel(first_seg, second_seg));
        }
    } else {
        // прямые отрезков пересекаются // или // скрещиваются //
        if (Vector<long long>::ScalarMultiple(normal_vector, between_vector) == 0) {
            return std::sqrt(SquaredDistOnCrossing(first_seg, second_seg));
        } else {
            return DistOnOther(first_seg, second_seg);
        }
    }
}

int main() {
    std::vector<Vertex<long long>> vec;
    long long x, y, z;
    for (size_t iter = 0; iter < 4; ++iter) {
        std::cin >> x >> y >> z;
        vec.emplace_back(x, y, z);
    }
    LineSegment first(vec[0], vec[1]);
    LineSegment second(vec[2], vec[3]);
    std::cout << std::fixed << std::setprecision(15) << FindDistBetweenSegments(first, second);
    return 0;
}