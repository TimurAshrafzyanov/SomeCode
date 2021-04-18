#include <iostream>
#include <vector>
#include <stack>


// Проверить, есть ли в графе цикл
// И если его нет, вывести отсортированные вершины
// так, что ребра идут из тех вершин что раньше
// только в те, что позже по порядку

enum Color {
    WHITE = 0,
    GREY = 1,
    BLACK = 2
};

class Graph {
public:
    Graph(int n) : vertexCount(n), edges(n){};

    void AddEdge(int begin, int end) {
        edges[begin].push_back(end);
    }

    bool HasEdge(int begin, int end) const;

    int VertexCount() const {
        return vertexCount;
    }

    const std::vector<int>& Neighbors(int vertex) const {
        return edges[vertex];
    }

private:
    int vertexCount;
    std::vector<std::vector<int>> edges;
};

bool Graph::HasEdge(int begin, int end) const {
    bool ifHere = false;
    for (int edge : edges[begin]) {
        if (edge == end) {
            ifHere = true;
        }
    }
    return ifHere;
}

bool PartTopSort(const Graph& g, std::vector<Color>& colors,  std::stack<int>& visited, int vertex) {
    colors[vertex] = GREY;
    const std::vector<int>& neighbors = g.Neighbors(vertex);
    for (int neighbor : neighbors) {
        if (colors[neighbor] == GREY) {
            return false;
        }
        if (colors[neighbor] == WHITE) {
            if (!PartTopSort(g, colors, visited, neighbor)) {
                return false;
            }
        }
    }
    colors[vertex] = BLACK;
    visited.push(vertex);
    return true;
}

bool TopSort(const Graph& g, std::vector<int>& sorted) {
    sorted.clear();

    std::vector<Color> colors(g.VertexCount(), WHITE);
    std::stack<int> visited;
    for (int i = 0; i < g.VertexCount(); i++) {
        if (colors[i] == WHITE) {
            if (!PartTopSort(g, colors, visited, i)) {
                return false;
            }
        }
    }

    while (!visited.empty()) {
        sorted.push_back(visited.top());
        visited.pop();
    }
    return true;
}

int main() {
    int n, m;
    std::cin >> n >> m;
    Graph g(n);
    int begin, end;
    for (int i = 0; i < m; i++) {
        std::cin >> begin >> end;
        g.AddEdge(begin, end);
    }

    std::vector<int> sorted;
    if (!TopSort(g, sorted)) {
        std::cout << "NO";
    } else {
        std::cout << "YES" << std::endl;
        for (int vertex : sorted){
            std::cout << vertex << " ";
        }
    }

    return 0;
}