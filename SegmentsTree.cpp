#include <iostream>
#include <vector>


// Есть расскашенная полоса длины n.
// Каждым запросом требуется перекрасить данную часть полосы
// в другой цвет, и найти минимум на данном отрезке.

struct Vertex {
public:
    explicit Vertex(int v = -1, int m = -1, int f = -1, int l = -1)
            : value(v), min(m), first(f), last(l) {}

    int min;
    int value;
    int first;
    int last;
};

class SegmentsTree {
public:
    explicit SegmentsTree(int n, const std::vector<int> &start);
    int Minimum(int left, int right);
    void Insert(int left, int right, int number);

private:
    int vertex_count;
    std::vector<Vertex> vertexes;

    int Minimum(int left, int right, int v);
    int Insert(int left, int right, int number, int v);
    void Push(int v);
    void Init(int ind);
};

void SegmentsTree::Init(int ind) {
    if (ind > vertex_count - 1) return;
    if (ind * 2 + 1 > vertex_count - 1) return;
    Init(ind * 2 + 1);
    vertexes[ind].min = vertexes[ind * 2 + 1].min;
    vertexes[ind].first = vertexes[ind * 2 + 1].first;
    vertexes[ind].last = vertexes[ind * 2 + 1].last;
    if (ind * 2 + 2 < vertex_count) {
        Init(ind * 2 + 2);
        if (vertexes[ind].min > vertexes[ind * 2 + 2].min)
            if (vertexes[ind * 2 + 2].min != -1)
                vertexes[ind].min = vertexes[ind * 2 + 2].min;
        if (vertexes[ind * 2 + 2].last != -1)
            vertexes[ind].last = vertexes[ind * 2 + 2].last;
    }
}

SegmentsTree::SegmentsTree(int n, const std::vector<int> &start) {
    int two = 1;
    while (two < n) two *= 2;
    vertexes = std::vector<Vertex>(two - 1 + n);
    vertex_count = two - 1 + n;
    int ind = static_cast<int>(start.size()) - 1;
    for (int i = vertex_count - 1; ind >= 0; --i) {
        vertexes[i] = Vertex(start[ind], start[ind], ind, ind);
        --ind;
    }
    Init(0);
}

void SegmentsTree::Push(int v) {
    if (v >= vertex_count) return;
    if (vertexes[v].value == -1) return;
    if (v * 2 + 1 < vertex_count) {
        vertexes[v * 2 + 1].value = vertexes[v].value;
        vertexes[v * 2 + 1].min = vertexes[v].value;
    } else return;
    if (v * 2 + 2 < vertex_count) {
        vertexes[v * 2 + 2].value = vertexes[v].value;
        vertexes[v * 2 + 2].min = vertexes[v].value;
    }
    vertexes[v].value = -1;
}

int SegmentsTree::Minimum(int left, int right, int v) {
    if (vertexes[v].first == left && vertexes[v].last == right) {
        return vertexes[v].min;
    }
    if (vertexes[v].value != -1) return vertexes[v].min;
    if (left > vertexes[v * 2 + 1].last)
        return Minimum(left, right, v * 2 + 2);
    if (right <= vertexes[v * 2 + 1].last)
        return Minimum(left, right, v * 2 + 1);
    int left_min = Minimum(left, vertexes[v * 2 + 1].last, v * 2 + 1);
    int right_min = Minimum(vertexes[v * 2 + 2].first, right, v * 2 + 2);
    return std::min(left_min, right_min);
}

int SegmentsTree::Minimum(int left, int right) {
    return Minimum(left, right, 0);
}

int SegmentsTree::Insert(int left, int right, int number, int v) {
    if (vertexes[v].first == left && vertexes[v].last == right) {
        vertexes[v].value = number;
        vertexes[v].min = number;
        return vertexes[v].min;
    }
    if (vertexes[v].value != -1) Push(v);
    if (left > vertexes[v * 2 + 1].last) {
        int new_min = Insert(left, right, number, v * 2 + 2);
        vertexes[v].min = std::min(vertexes[v * 2 + 1].min, new_min);
        return vertexes[v].min;
    }
    if (right <= vertexes[v * 2 + 1].last) {
        int new_min = Insert(left, right, number, v * 2 + 1);
        if (v * 2 + 2 < vertex_count && vertexes[v * 2 + 2].min != -1)
            vertexes[v].min = std::min(vertexes[v * 2 + 2].min, new_min);
        else vertexes[v].min = new_min;
        return vertexes[v].min;
    }
    int left_min = Insert(left, vertexes[v * 2 + 1].last, number, v * 2 + 1);
    int right_min = Insert(vertexes[v * 2 + 2].first, right, number, v * 2 + 2);
    vertexes[v].min = std::min(left_min, right_min);
    return vertexes[v].min;
}

void SegmentsTree::Insert(int left, int right, int number) {
    Insert(left, right, number, 0);
}

int main() {
    int n;
    std::cin >> n;
    std::vector<int> start;
    int first, second, third;
    for (size_t i = 0; i < n; ++i) {
        std::cin >> first >> second >> third;
        start.Push_back(first + second + third);
    }
    SegmentsTree t(n, start);

    int m;
    std::cin >> m;
    int left, right;
    std::vector<int> answer(0);
    for (size_t i = 0; i < m; ++i) {
        std::cin >> left >> right;
        std::cin >> first >> second >> third;
        t.Insert(left, right, first + second + third);
        std::cin >> left >> right;
        answer.Push_back(t.Minimum(left, right));
    }

    for (int number : answer)
        std::cout << number << " ";
    return 0;
}