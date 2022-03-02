#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <vector>

typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int8_t   s8;
typedef int16_t  s16;
typedef int32_t  s32;
typedef int64_t  s64;

typedef float  f32;
typedef double f64;

struct GraphMatrix {
    s32 *cells;
    s32 dim;
};

void skip_lines(FILE *f, s32 count) {
    char c = 0;
    s32 lines_done = 0;
    while (lines_done < count) {
        s32 read_bytes = fread(&c, 1, 1, f);
        if (read_bytes != 1) break;
        if (c == '\n') lines_done++;
    }
}

struct Position2D {
    s32 x, y;
};

struct Instance {
    GraphMatrix graph;
    std::vector<Position2D> positions;
};

Instance parse_file(const char *filepath) {
    FILE *file = fopen(filepath, "r");
    if (!file) return {};

    // skip first 6 lines
    skip_lines(file, 6);

    std::vector<Position2D> positions;
    while (true) {
        s32 node_idx;
        Position2D pos;
        s32 read_elems = fscanf(file, " %d %d %d ", &node_idx, &pos.x, &pos.y);
        if (read_elems != 3) break;
        positions.push_back(pos);
    }

    Instance result;
    result.positions = positions;
    auto& graph = result.graph;
    graph.dim = positions.size();
    graph.cells = (s32 *)malloc(sizeof(s32) * graph.dim * graph.dim);

    for (s32 i = 0; i < graph.dim; i++) {
        Position2D a = positions[i];
        for (s32 j = 0; j < graph.dim; j++) {
            Position2D b = positions[j];
            s32 dx = a.x - b.x;
            s32 dy = a.y - b.y;
            f32 delta = sqrtf(dx*dx + dy*dy);
            s32 cost = (s32)(delta + 0.5f);
            graph.cells[graph.dim * i + j] = cost;
        }
    }

    return result;
}

struct Solution {
    std::vector<s32> loop_a;
    std::vector<s32> loop_b;
};

// returns the index into available vector
s32 closest_node(GraphMatrix g, s32 node, std::vector<s32>& available) {
    s32 result = 0;
    s32 result_dist = INT32_MAX;
    s32 base = g.dim * node;
    for (s32 i = 0; i < available.size(); i++) {
        s32 idx = available[i];
        s32 dist = g.cells[base + idx];
        if (dist < result_dist) {
            result_dist = dist;
            result = i;
        }
    }

    return result;
}

// returns the index into available vector
s32 furthest_node(GraphMatrix g, s32 node, std::vector<s32>& available) {
    s32 result = 0;
    s32 result_dist = 0;
    s32 base = g.dim * node;
    for (s32 i = 0; i < available.size(); i++) {
        s32 idx = available[i];
        s32 dist = g.cells[base + idx];
        if (dist > result_dist) {
            result_dist = dist;
            result = i;
        }
    }

    return result;
}

void vector_remove_idx(std::vector<s32>& v, s32 idx) {
    v[idx] = v[v.size()-1];
    v.pop_back();
}

// returns the node idx
s32 consume_closest(GraphMatrix g, s32 node, std::vector<s32>& available) {
    s32 i = closest_node(g, node, available);
    s32 result = available[i];
    vector_remove_idx(available, i);
    return result;
}

// returns the node idx
s32 consume_furthest(GraphMatrix g, s32 node, std::vector<s32>& available) {
    s32 i = furthest_node(g, node, available);
    s32 result = available[i];
    vector_remove_idx(available, i);
    return result;
}

Solution greedy_simple(GraphMatrix graph) {
    Solution result;
    auto& a = result.loop_a;
    auto& b = result.loop_b;

    s32 max_loop_a_size = graph.dim / 2 + (graph.dim & 1);
    s32 max_loop_b_size = graph.dim - max_loop_a_size;

    std::vector<s32> available(graph.dim);
    for (s32 i = 0; i < graph.dim; i++) {
        available[i] = i;
    }

    s32 first_idx = rand() % graph.dim;
    a.push_back(first_idx);
    vector_remove_idx(available, first_idx);

    s32 furthest = furthest_node(graph, first_idx, available);
    b.push_back(available[furthest]);
    vector_remove_idx(available, furthest);

    while (available.size() > 0) {
        s32 last_a = a[a.size()-1];
        s32 last_b = b[b.size()-1];
        s32 candidate_a = closest_node(graph, last_a, available);
        s32 candidate_b = closest_node(graph, last_b, available);

        s32 cost_a = graph.cells[last_a*graph.dim + available[candidate_a]];
        s32 cost_b = graph.cells[last_b*graph.dim + available[candidate_b]];

        bool can_expand_a = a.size() < max_loop_a_size;
        bool can_expand_b = b.size() < max_loop_b_size;
        bool should_expand_a = can_expand_a;
        if (can_expand_a && can_expand_b) {
            should_expand_a = cost_a <= cost_b;
        }
        if (should_expand_a) {
            a.push_back(available[candidate_a]);
            vector_remove_idx(available, candidate_a);
        } else {
            b.push_back(available[candidate_b]);
            vector_remove_idx(available, candidate_b);
        }
    }

    return result;
}

// returns index into the available vector
s32 best_intermediate_node(GraphMatrix graph, s32 a, s32 b, std::vector<s32>& available) {
    s32 best = 0;
    s32 best_cost = INT32_MAX;
    for (s32 i = 0; i < available.size(); i++) {
        s32 idx = available[i];
        s32 cost = graph.cells[idx*graph.dim + a] + graph.cells[idx*graph.dim + b];
        cost -= graph.cells[a*graph.dim + b];
        if (cost < best_cost) {
            best = i;
            best_cost = cost;
        }
    }
    return best;
}

Solution greedy_loop(GraphMatrix graph) {
    Solution result;
    auto& a = result.loop_a;
    auto& b = result.loop_b;

    s32 max_loop_a_size = graph.dim / 2 + (graph.dim & 1);
    s32 max_loop_b_size = graph.dim - max_loop_a_size;

    std::vector<s32> available(graph.dim);
    for (s32 i = 0; i < graph.dim; i++) {
        available[i] = i;
    }

    s32 loop_a_first = rand() % graph.dim;
    s32 first_node_idx = available[loop_a_first];
    a.push_back(first_node_idx);
    vector_remove_idx(available, loop_a_first);
    a.push_back(consume_closest(graph, first_node_idx, available));

    b.push_back(consume_furthest(graph, first_node_idx, available));
    b.push_back(consume_closest(graph, b[0], available));

    while (available.size() > 0) {
        s32 best_a = 0;
        s32 best_a_cost = INT32_MAX;
        s32 best_a_edge = 0;
        for (s32 edge = 0; edge < a.size(); edge++) {
            s32 node_a = a[edge];
            s32 node_b = a[(edge+1) % a.size()];
            s32 node_i = best_intermediate_node(graph, node_a, node_b, available);
            s32 node = available[node_i];
            s32 cost = graph.cells[node*graph.dim + node_a] + graph.cells[node*graph.dim + node_b];
            cost -= graph.cells[node_a*graph.dim + node_b];
            if (cost < best_a_cost) {
                best_a = node_i;
                best_a_cost = cost;
                best_a_edge = edge;
            }
        }

        s32 best_b = 0;
        s32 best_b_cost = INT32_MAX;
        s32 best_b_edge = 0;
        for (s32 edge = 0; edge < b.size(); edge++) {
            s32 node_a = b[edge];
            s32 node_b = b[(edge+1) % b.size()];
            s32 node_i = best_intermediate_node(graph, node_a, node_b, available);
            s32 node = available[node_i];
            s32 cost = graph.cells[node*graph.dim + node_a] + graph.cells[node*graph.dim + node_b];
            cost -= graph.cells[node_a*graph.dim + node_b];
            if (cost < best_b_cost) {
                best_b = node_i;
                best_b_cost = cost;
                best_b_edge = edge;
            }
        }

        bool can_expand_a = a.size() < max_loop_a_size;
        bool can_expand_b = b.size() < max_loop_b_size;
        bool should_expand_a = can_expand_a;
        if (can_expand_a && can_expand_b) {
            should_expand_a = best_a_cost <= best_b_cost;
        }
        if (should_expand_a) {
            s32 node = available[best_a];
            vector_remove_idx(available, best_a);
            a.insert(a.begin() + best_a_edge + 1, node);
        } else {
            s32 node = available[best_b];
            vector_remove_idx(available, best_b);
            b.insert(b.begin() + best_b_edge + 1, node);
        }
    }

    return result;
}

// returns true on success
bool write_entire_file(const char *filename, void *data, size_t size) {
    FILE *f = fopen(filename, "w");
    if (!f) return false;
    size_t written = fwrite(data, size, 1, f);
    fclose(f);
    if (written != 1) return false;
    return true;
}

template <typename T>
bool write_entire_file(const char *filename, std::vector<T>& v) {
    return write_entire_file(filename, &v[0], v.size()*sizeof(T));
}

s32 score(GraphMatrix graph, Solution s) {
    s32 cost = 0;
    for (s32 i = 0; i < s.loop_a.size(); i++) {
        s32 a = s.loop_a[i];
        s32 b = s.loop_a[(i+1) % s.loop_a.size()];
        cost += graph.cells[a*graph.dim + b];
    }
    for (s32 i = 0; i < s.loop_b.size(); i++) {
        s32 a = s.loop_b[i];
        s32 b = s.loop_b[(i+1) % s.loop_b.size()];
        cost += graph.cells[a*graph.dim + b];
    }
    return cost;
}

int main() {
    srand(time(0));

    auto parsed = parse_file("data/kroB100.tsp");
    auto solution1 = greedy_simple(parsed.graph);
    auto solution2 = greedy_loop(parsed.graph);

    printf("solution 1 cost = %d\n", score(parsed.graph, solution1));
    printf("solution 2 cost = %d\n", score(parsed.graph, solution2));

    write_entire_file("results/a.dat", solution2.loop_a);
    write_entire_file("results/b.dat", solution2.loop_b);
    write_entire_file("results/pos.dat", parsed.positions);

    return 0;
}
