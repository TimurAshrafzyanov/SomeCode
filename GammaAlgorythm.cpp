#include <iostream>
#include <vector>
#include <stack>


// Проверить, планарен ли граф
// Реализовано с помощью гамма-алгоритма

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
    for (const auto& e : edges[begin]) {
        if (e == end) {
            ifHere = true;
        }
    }
    return ifHere;
}

struct Edge {
    Edge(int begin, int end) : first(begin), second(end) {}

    int first;
    int second;
};

bool DetectCycleRecursively(const Graph& g, int vertex, int parent,
        std::vector<bool>& colours, std::stack<int>& visited) {
    colours[vertex] = true;
    const std::vector<int>& neighbors = g.Neighbors(vertex);
    for (int i = 0; i < neighbors.size(); ++i) {
        if (neighbors[i] != parent) {
            visited.push(neighbors[i]);
            if (colours[neighbors[i]]) {
                return true;
            }
            if (DetectCycleRecursively(g, neighbors[i], vertex, colours, visited)) {
                return true;
            }
            visited.pop();
        }
    }
    return false;
}

bool DetectCycle(const Graph& g, std::vector<int>& cycle) {
    std::vector<bool> colours(g.VertexCount(), false);
    std::stack<int> visited;
    bool ifFound = false;
    int i = 0;
    while (!ifFound && i < g.VertexCount()) {
        if (!colours[i]) {
            ifFound = DetectCycleRecursively(g, i, -1, colours, visited);
        }
        ++i;
    }
    if (!ifFound) return false;
    cycle.push_back(visited.top());
    visited.pop();
    while (!visited.empty() && cycle.front() != visited.top()) {
        cycle.push_back(visited.top());
        visited.pop();
    }
    return true;
}

void FindComponentsRecursively(const Graph& g, int vertex, int parent, std::vector<bool>& colours,
        std::vector<int>& tin, std::vector<int>& up, int& currentTime, std::stack<Edge>& edges,
        std::vector<std::vector<Edge>>& components, std::vector<std::vector<bool>>& ifEdgeDeleted) {
    const std::vector<int>& neighbors = g.Neighbors(vertex);
    colours[vertex] = true;
    tin[vertex] = currentTime;
    ++currentTime;
    up[vertex] = tin[vertex];

    for(int i = 0; i < neighbors.size(); ++i) {
        if (!ifEdgeDeleted[vertex][neighbors[i]] && neighbors[i] != parent) {
            edges.push(Edge(vertex, neighbors[i]));
        }
        if (!colours[neighbors[i]]) {
            FindComponentsRecursively(g, neighbors[i], vertex, colours, tin,
                    up, currentTime, edges, components, ifEdgeDeleted);
            up[vertex] = std::min(up[vertex], up[neighbors[i]]);

            if (up[neighbors[i]] >= tin[vertex]) {
                std::vector<Edge> currentComponent;
                while ((edges.top().first != vertex) || (edges.top().second != neighbors[i])) {
                    currentComponent.push_back(edges.top());
                    ifEdgeDeleted[edges.top().first][edges.top().second] = true;
                    ifEdgeDeleted[edges.top().second][edges.top().first] = true;
                    edges.pop();
                }
                currentComponent.push_back(edges.top());
                ifEdgeDeleted[edges.top().first][edges.top().second] = true;
                ifEdgeDeleted[edges.top().second][edges.top().first] = true;
                edges.pop();

                components.push_back(currentComponent);
            }

        } else {
            if (neighbors[i] != parent) {
                up[vertex] = std::min(up[vertex], tin[neighbors[i]]);
            }
        }
    }
}

void FindConnectivityComponents(const Graph& g, std::vector<std::vector<Edge>>& components) {
    int n = g.VertexCount();
    std::vector<bool> colours(n,false);
    std::vector<int> tin(n, 0);
    std::vector<int> up(n);
    std::stack<Edge> edges;
    std::vector<std::vector<bool>> ifEdgeDeleted(n, std::vector<bool>(n, false));

    for (int vertex = 0; vertex < n; ++vertex) {
        if(!colours[vertex]) {
            int currentTime = 0;
            FindComponentsRecursively(g,vertex, vertex, colours, tin, up,
                    currentTime,edges, components, ifEdgeDeleted);
        }
    }
}

void FindChain(const Graph& g, int vertex, int start, int parent,
        std::vector<bool>& colours, std::vector<std::vector<bool>>& deletedEdges,
        std::vector<int>& chain, bool& ifFound, std::vector<int>& contactVertexes) {
    colours[vertex] = true;
    chain.push_back(vertex);
    for (int i = 0; i < contactVertexes.size(); ++i) {
        if (contactVertexes[i] == vertex) {
            return;
        }
    }
    const std::vector<int>& neighbors = g.Neighbors(vertex);
    for (int i = 0; i < neighbors.size(); ++i) {
        if (!deletedEdges[vertex][neighbors[i]] && neighbors[i] != parent && neighbors[i] != start) {
            if (!ifFound) {
                bool ifContact = false;
                for (int j = 0; j < contactVertexes.size(); ++j) {
                    if (neighbors[i] == contactVertexes[j]) {
                        ifContact = true;
                    }
                }
                if (ifContact) {
                    chain.push_back(neighbors[i]);
                    ifFound = true;
                    break;
                }
            }
            if (!colours[neighbors[i]]) {
                if (!ifFound) {
                    FindChain(g, neighbors[i], start, vertex, colours,
                            deletedEdges, chain, ifFound, contactVertexes);
                }
            }
        }
    }
}

void MakeNewPlanes(std::vector<int>& oldPLane, std::vector<int>& chain,
                   std::vector<int>& firstNewPlane, std::vector<int>& secondNewPlane) {
    int i = 0;
    while (oldPLane[i] != chain[0] && oldPLane[i] != chain[chain.size() - 1]) {
        firstNewPlane.push_back(oldPLane[i]);
        ++i;
    }
    if (oldPLane[i] == chain[0]) {
        for (int j = 0; j < chain.size() - 1; ++j) {
            firstNewPlane.push_back(chain[j]);
        }
        while (oldPLane[i] != chain[chain.size() - 1]) {
            secondNewPlane.push_back(oldPLane[i]);
            ++i;
        }
        for (int j = i; j < oldPLane.size(); ++j) {
            firstNewPlane.push_back(oldPLane[j]);
        }
        for (int j = chain.size() - 1; j > 0; --j) {
            secondNewPlane.push_back(chain[j]);
        }
    } else {
        for (int j = chain.size() - 1; j > 0; --j) {
            firstNewPlane.push_back(chain[j]);
        }
        while (oldPLane[i] != chain[0]) {
            secondNewPlane.push_back(oldPLane[i]);
            ++i;
        }
        for (int j = i; j < oldPLane.size(); ++j) {
            firstNewPlane.push_back(oldPLane[j]);
        }
        for (int j = 0; j < chain.size() - 1; ++j) {
            secondNewPlane.push_back(chain[j]);
        }
    }
}

int GetNeighborOfContact(const Graph& g, int contVer, std::vector<int>& segment,
        std::vector<std::vector<bool>>& deletedEdges) {
    const std::vector<int>& neighbors = g.Neighbors(contVer);
    int i = 0;
    bool ifFound = false;
    while (!ifFound && i < neighbors.size()) {
        for (int j = 0; j < segment.size(); ++j) {
            if (!deletedEdges[contVer][neighbors[i]]) {
                if (segment[j] == neighbors[i]) {
                    ifFound = true;
                    break;
                }
            }
        }
        ++i;
    }
    return neighbors[i - 1];
}

void FindSegments(int& time, const Graph& g, int vertex, int start,
        std::vector<bool>& colours, std::vector<std::vector<bool>>& deletedEdges,
        std::vector<std::vector<int>>& segments, std::vector<int>& currentSegment,
        std::vector<int>& vertexesInGraph) {
    colours[vertex] = true;
    currentSegment.push_back(vertex);
    for (int i = 0; i < vertexesInGraph.size(); ++i) {
        if (vertexesInGraph[i] == vertex && time != 0) {
            return;
        }
    }
    ++time;
    const std::vector<int>& neighbors = g.Neighbors(vertex);
    for (int i = 0; i < neighbors.size(); ++i) {
        if (!colours[neighbors[i]] && !deletedEdges[vertex][neighbors[i]]) {
            FindSegments(time, g, neighbors[i], start, colours,
                         deletedEdges, segments, currentSegment, vertexesInGraph);

            if (vertex == start) {
                if (currentSegment.size() > 1) {
                    segments.push_back(currentSegment);
                }
                while (!currentSegment.empty()) {
                    currentSegment.pop_back();
                }
                currentSegment.push_back(vertex);
                for (int i = 0; i < vertexesInGraph.size(); ++i) {
                    if (vertexesInGraph[i] != vertex) {
                        colours[vertexesInGraph[i]] = false;
                    }
                }
            }
        }
    }
}

int FindSuitableSegment(int& minCount, std::vector<std::vector<int>>& segments,
        std::vector<std::vector<int>>& planes, std::vector<int>& oneOfSuitablePlanes,
        std::vector<std::vector<int>>& contactVertexes) {
    std::vector<int> countOfSuitablePlanes(segments.size());
    for (int i = 0; i < segments.size(); ++i) {
        if (segments[i].size() == 1) countOfSuitablePlanes[i] = planes.size() + 1;
        else {
            int count = 0;
            for (int j = 0; j < planes.size(); ++j) {
                bool ifSuitable = true;
                for (int k = 0; k < contactVertexes[i].size(); ++k) {
                    bool ifHere = false;
                    for (int t = 0; t < planes[j].size(); ++t) {
                        if (contactVertexes[i][k] == planes[j][t]) {
                            ifHere = true;
                            break;
                        }
                    }
                    if (!ifHere) {
                        ifSuitable = false;
                        break;
                    }
                }
                if (ifSuitable) {
                    oneOfSuitablePlanes[i] = j;
                    ++count;
                }
            }
            countOfSuitablePlanes[i] = count;
        }
    }

    minCount = static_cast<int>(planes.size()) + 1;
    int number = -1;
    for (int i = 0; i < countOfSuitablePlanes.size(); ++i) {
        if (countOfSuitablePlanes[i] < minCount) {
            minCount = countOfSuitablePlanes[i];
            number = i;
        }
    }
    return number;
}

bool InsertSegments(const Graph& g, std::vector<int>& vertexesInGraph,
        std::vector<std::vector<int>>& planes,
        std::vector<std::vector<bool>>& deletedEdges) {
    std::vector<bool> colours(g.VertexCount(), false);
    std::vector<std::vector<int>> segments;
    for (int i = 0; i < vertexesInGraph.size(); ++i) {
        std::vector<int> currentSegment;
        int time = 0;
        FindSegments(time, g, vertexesInGraph[i], vertexesInGraph[i],
                colours, deletedEdges, segments, currentSegment, vertexesInGraph);
    }
    std::vector<std::vector<int>> contactVertexes(segments.size());
    for (int i = 0; i < segments.size(); ++i) {
        for (int j = 0; j < segments[i].size(); ++j) {
            for (int k = 0; k < vertexesInGraph.size(); ++k) {
                if (segments[i][j] == vertexesInGraph[k]) {
                    contactVertexes[i].push_back(segments[i][j]);
                    break;
                }
            }
        }
    }

    std::vector<int> oneOfSuitablePlanes(segments.size());
    int minCount = 0;
    int number = FindSuitableSegment(minCount, segments,
            planes,oneOfSuitablePlanes, contactVertexes);
    if (number == -1) return true;
    if (minCount == 0) return false;

    std::vector<bool> newColours(g.VertexCount(), false);
    std::vector<int> chain;
    bool ifFound = false;
    int start = contactVertexes[number][0];
    chain.push_back(start);
    int neighborOfContact = GetNeighborOfContact(g, start, segments[number], deletedEdges);
    FindChain(g, neighborOfContact, start, start, newColours,deletedEdges,
            chain, ifFound, contactVertexes[number]);

    std::vector<int> firstNewPlane;
    std::vector<int> secondNewPlane;
    MakeNewPlanes(planes[oneOfSuitablePlanes[number]], chain,
            firstNewPlane, secondNewPlane);
    while (!planes[oneOfSuitablePlanes[number]].empty()) {
        planes[oneOfSuitablePlanes[number]].pop_back();
    }
    for (int i = 0; i < firstNewPlane.size(); ++i) {
        planes[oneOfSuitablePlanes[number]].push_back(firstNewPlane[i]);
    }
    planes.push_back(secondNewPlane);
    for (int i = 0; i < chain.size() - 1; ++i) {
        deletedEdges[chain[i]][chain[i + 1]] = true;
        deletedEdges[chain[i + 1]][chain[i]] = true;
        if (i != 0) vertexesInGraph.push_back(chain[i]);
    }

    return InsertSegments(g, vertexesInGraph, planes, deletedEdges);
}

bool IsComponentPlane(const Graph& g, std::vector<Edge> component) {
    Graph compGraph(g.VertexCount());
    for (int i = 0; i < component.size(); ++i) {
        if (!compGraph.HasEdge(component[i].first, component[i].second)) {
            compGraph.AddEdge(component[i].first, component[i].second);
            compGraph.AddEdge(component[i].second, component[i].first);
        }
    }
    std::vector<std::vector<int>> planes;
    std::vector<int> cycle;
    if (!DetectCycle(compGraph, cycle)) {
        return true;
    }
    planes.push_back(cycle);
    planes.push_back(cycle);

    std::vector<std::vector<bool>> deletedEdges(g.VertexCount(),
            std::vector<bool>(g.VertexCount(), false));
    std::vector<int> vertexesInGraph;
    for (int i = 1; i < cycle.size(); ++i) {
        deletedEdges[cycle[i]][cycle[i - 1]] = true;
        deletedEdges[cycle[i - 1]][cycle[i]] = true;
        vertexesInGraph.push_back(cycle[i]);
    }
    deletedEdges[cycle[cycle.size() - 1]][cycle[0]] = true;
    deletedEdges[cycle[0]][cycle[cycle.size() - 1]] = true;
    vertexesInGraph.push_back(cycle[0]);

    return InsertSegments(compGraph, vertexesInGraph, planes, deletedEdges);
}

bool IsPlaneGraph(const Graph& g) {
    std::vector<std::vector<Edge>> components;
    FindConnectivityComponents(g, components);

    for (int i = 0; i < components.size(); ++i) {
        if (!IsComponentPlane(g, components[i])) {
            return false;
        }
    }
    return true;
}

int main() {
    int n, m;
    std::cin >> n >> m;
    Graph g(n);
    int begin, end;
    for (int i = 0; i < m; ++i) {
        std::cin >> begin >> end;
        g.AddEdge(begin, end);
        g.AddEdge(end, begin);
    }

    if (IsPlaneGraph(g)) std::cout << "YES";
    else std::cout << "NO";

    return 0;
}