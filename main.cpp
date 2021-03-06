#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <iterator>
#include <set>
#include <unordered_set>
#include <string>
#include <assert.h>
#define SOKOL_IMPL
#include "sokol_time.h"

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
    inline s32 get(s32 i, s32 j) {
        return cells[i*dim + j];
    };
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

        s32 cost_a = graph.get(last_a, available[candidate_a]);
        s32 cost_b = graph.get(last_b, available[candidate_b]);

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
        s32 cost = graph.get(idx, a) + graph.get(idx, b) - graph.get(a, b);
        if (cost < best_cost) {
            best = i;
            best_cost = cost;
        }
    }
    return best;
}

Solution greedy_loop(GraphMatrix graph, Solution *initial_solution = NULL) {
    Solution result;
    auto& a = result.loop_a;
    auto& b = result.loop_b;

    s32 max_loop_a_size = graph.dim / 2 + (graph.dim & 1);
    s32 max_loop_b_size = graph.dim - max_loop_a_size;

    std::vector<s32> available(graph.dim);
    for (s32 i = 0; i < graph.dim; i++) {
        available[i] = i;
    }

    if (initial_solution == NULL) {
        s32 loop_a_first = rand() % graph.dim;
        s32 first_node_idx = available[loop_a_first];
        a.push_back(first_node_idx);
        vector_remove_idx(available, loop_a_first);
        a.push_back(consume_closest(graph, first_node_idx, available));

        b.push_back(consume_furthest(graph, first_node_idx, available));
        b.push_back(consume_closest(graph, b[0], available));
    } else {
        a = initial_solution->loop_a;
        b = initial_solution->loop_b;

        for (s32 i = 0; i < a.size(); i++) {
            auto pos = std::find(available.begin(), available.end(), a[i]);
            if (pos != available.end()) available.erase(pos);
        }

        for (s32 i = 0; i < b.size(); i++) {
            auto pos = std::find(available.begin(), available.end(), b[i]);
            if (pos != available.end()) available.erase(pos);
        }
    }

    while (available.size() > 0) {
        s32 best_a = 0;
        s32 best_a_cost = INT32_MAX;
        s32 best_a_edge = 0;
        for (s32 edge = 0; edge < a.size(); edge++) {
            s32 node_a = a[edge];
            s32 node_b = a[(edge+1) % a.size()];
            s32 node_i = best_intermediate_node(graph, node_a, node_b, available);
            s32 node = available[node_i];
            s32 cost = graph.get(node, node_a) + graph.get(node, node_b) - graph.get(node_a, node_b);
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
            s32 cost = graph.get(node, node_a) + graph.get(node, node_b) - graph.get(node_a, node_b);
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

struct CycleEdgeCost {
    bool cycle = 0; // 0 - a, 1 -b
    s32 node_i = 0;
    s32 edge = 0;
    s32 cost = INT32_MAX;
};

Solution regret_loop(GraphMatrix graph) {
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
    a.push_back(consume_closest(graph, first_node_idx, available));

    b.push_back(consume_furthest(graph, first_node_idx, available));
    b.push_back(consume_closest(graph, b[0], available));
    b.push_back(consume_closest(graph, b[0], available));

    while (available.size() > 0) {
        CycleEdgeCost highest_regret_insertion;
        s32 highest_regret = 0;

        bool can_expand_a = a.size() < max_loop_a_size;
        bool can_expand_b = b.size() < max_loop_b_size;

        for (s32 node_i = 0; node_i < available.size(); node_i++) {
            s32 node = available[node_i];
            CycleEdgeCost best, second_best;

            if (can_expand_a) {
                for (s32 edge = 0; edge < a.size(); edge++) {
                    s32 node_a = a[edge];
                    s32 node_b = a[(edge+1) % a.size()];
                    s32 cost = graph.get(node, node_a) + graph.get(node, node_b) - graph.get(node_a, node_b);

                    if (cost <= best.cost) {
                        second_best = best;
                        best.cycle = 0;
                        best.node_i = node_i;
                        best.edge = edge;
                        best.cost = cost;
                    }
                    else if (cost <= second_best.cost) {
                        second_best.cycle = 0;
                        second_best.node_i = node_i;
                        second_best.edge = edge;
                        second_best.cost = cost;
                    }
                }
            }

            if (can_expand_b) {
                for (s32 edge = 0; edge < b.size(); edge++) {
                    s32 node_a = b[edge];
                    s32 node_b = b[(edge+1) % b.size()];
                    s32 cost = graph.get(node, node_a) + graph.get(node, node_b) - graph.get(node_a, node_b);

                    if (cost <= best.cost) {
                        second_best = best;
                        best.cycle = 1;
                        best.node_i = node_i;
                        best.edge = edge;
                        best.cost = cost;
                    }
                    else if (cost <= second_best.cost) {
                        second_best.cycle = 0;
                        second_best.node_i = node_i;
                        second_best.edge = edge;
                        second_best.cost = cost;
                    }
                }
            }

            second_best.cost = second_best.cost == INT32_MAX ? best.cost : second_best.cost;
            s32 regret = second_best.cost - best.cost;
            if (regret >= highest_regret) {
                highest_regret = regret;
                highest_regret_insertion = best;
            }
        }

        s32 node = available[highest_regret_insertion.node_i];
        vector_remove_idx(available, highest_regret_insertion.node_i);
        if (highest_regret_insertion.cycle == 0) {
            a.insert(a.begin() + highest_regret_insertion.edge + 1, node);
        }
        else {
            b.insert(b.begin() + highest_regret_insertion.edge + 1, node);
        }
    }

    return result;
}

// returns true on success
bool write_entire_file(const char *filename, void *data, size_t size) {
    FILE *f = fopen(filename, "wb");
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
        cost += graph.get(a, b);
    }
    for (s32 i = 0; i < s.loop_b.size(); i++) {
        s32 a = s.loop_b[i];
        s32 b = s.loop_b[(i+1) % s.loop_b.size()];
        cost += graph.get(a, b);
    }
    return cost;
}

void verify_solution(Instance instance, Solution s) {
    std::set<s32> a(s.loop_a.begin(), s.loop_a.end());
    std::set<s32> b(s.loop_b.begin(), s.loop_b.end());

    std::set<s32> u;
    std::set_union(a.begin(), a.end(),
                b.begin(), b.end(),
                std::inserter(u, u.begin()));

    if (u.size() != instance.graph.dim) {
        printf("Both loops do not contain all nodes of given instance\n");
        printf("Loop a size: %zu\n", a.size());
        printf("Loop b size: %zu\n", b.size());
        printf("Union of loops size: %zu\n", u.size());
        printf("Count of nodes that are not in any of the loops: %zu\n", instance.graph.dim-u.size());
    }
}

typedef Solution (*Solver)(GraphMatrix);

struct ExperimentResult {
    Solution best_solution;
    s32 min, max;
    f64 average;
    f64 average_time;
};

struct NodeExchange {
    s32 delta, i, j;
};

NodeExchange best_node_exchange(GraphMatrix graph, std::vector<s32>& loop_i, std::vector<s32>& loop_j, bool greedy) {
    NodeExchange result = {};
    s32 loop_i_count = loop_i.size();
    s32 loop_j_count = loop_j.size();
    s32 offset_i = rand() % loop_i_count;
    s32 offset_j = rand() % loop_j_count;
    for (s32 i = 0; i < loop_i.size(); i++) {
        s32 i_prev = loop_i[(i-1+loop_i_count+offset_i) % loop_i_count];
        s32 i_curr = loop_i[(i+offset_i) % loop_i_count];
        s32 i_next = loop_i[(i+1+offset_i) % loop_i_count];
        s32 i_cost = graph.get(i_curr, i_prev) + graph.get(i_curr, i_next);
        for (s32 j = 0; j < loop_j.size(); j++) {
            s32 j_prev = loop_j[(j-1+loop_j_count+offset_j) % loop_j_count];
            s32 j_curr = loop_j[(j+offset_j) % loop_j_count];
            s32 j_next = loop_j[(j+1+offset_j) % loop_j_count];
            s32 j_cost = graph.get(j_curr, j_prev) + graph.get(j_curr, j_next);
            s32 new_i_cost = graph.get(j_curr, i_prev) + graph.get(j_curr, i_next);
            s32 new_j_cost = graph.get(i_curr, j_prev) + graph.get(i_curr, j_next);
            s32 delta = new_i_cost + new_j_cost - i_cost - j_cost;
            if (delta < result.delta) {
                result.delta = delta;
                result.i = (i+offset_i) % loop_i_count;
                result.j = (j+offset_j) % loop_j_count;
                if (greedy) return result;
            }
        }
    }

    return result;
}

NodeExchange best_node_exchange(GraphMatrix graph, std::vector<s32>& loop, bool greedy) {
    NodeExchange result = {};
    s32 loop_count = loop.size();
    s32 offset = rand() % loop_count;
    for (s32 i = 0; i < loop_count; i++) {
        s32 i_prev = loop[(i-1+loop_count+offset) % loop_count];
        s32 i_curr = loop[(i+offset) % loop_count];
        s32 i_next = loop[(i+1+offset) % loop_count];
        s32 i_cost = graph.get(i_curr, i_prev) + graph.get(i_curr, i_next);

        for (s32 j = i+2; j < loop_count; j++) {
            s32 j_prev = loop[(j-1+offset) % loop_count];
            s32 j_curr = loop[(j+offset) % loop_count];
            s32 j_next = loop[(j+1+offset) % loop_count];
            s32 j_cost = graph.get(j_curr, j_prev) + graph.get(j_curr, j_next);
            s32 new_i_cost = graph.get(j_curr, i_prev) + graph.get(j_curr, i_next);
            s32 new_j_cost = graph.get(i_curr, j_prev) + graph.get(i_curr, j_next);
            s32 delta = new_i_cost + new_j_cost - i_cost - j_cost;
            if (j_curr == i_next) {
                s32 old_cost = graph.get(i_curr, i_prev) + graph.get(j_curr, j_next);
                s32 new_cost = graph.get(j_curr, i_prev) + graph.get(i_curr, j_next);
                delta = new_cost - old_cost;
            } else if (i_curr == j_next) {
                s32 old_cost = graph.get(j_curr, j_prev) + graph.get(i_curr, i_next);
                s32 new_cost = graph.get(i_curr, j_prev) + graph.get(j_curr, i_next);
                delta = new_cost - old_cost;
            }
            if (delta < result.delta) {
                result.delta = delta;
                result.i = (i+offset) % loop_count;
                result.j = (j+offset) % loop_count;
                if (greedy) return result;
            }
        }
    }
    return result;
}

Solution neighbour_search_node(GraphMatrix graph, Solution init, bool greedy) {
    Solution result = init;

    std::vector<s32> *loops[2];
    loops[0] = &result.loop_a;
    loops[1] = &result.loop_b;

    while (true) {
        NodeExchange best = {};
        std::vector<s32> *best_loop = NULL;
        NodeExchange change = {};

        s32 loop_offset = rand() & 1;
        for (s32 loop_idx = 0; loop_idx < 2; loop_idx++) {
            auto& loop = *loops[(loop_idx+loop_offset) & 1];
            change = best_node_exchange(graph, loop, greedy);
            if (change.delta < best.delta) {
                best = change;
                best_loop = &loop;
                if (greedy) goto apply_change;
            }
        }

        change = best_node_exchange(graph, result.loop_a, result.loop_b, greedy);
        if (change.delta < best.delta) {
            best = change;
            best_loop = NULL;
        }

apply_change:
        if (best.delta < 0) {
            if (best_loop) {
                auto& loop = *best_loop;
                s32 t = loop[best.i];
                loop[best.i] = loop[best.j];
                loop[best.j] = t;
            } else {
                s32 t = result.loop_a[best.i];
                result.loop_a[best.i] = result.loop_b[best.j];
                result.loop_b[best.j] = t;
            }
        } else {
            break;
        }
    }
    return result;
}

enum ExchangeType {
    EXCHANGE_LOOP_A  = 0,
    EXCHANGE_LOOP_B  = 1,
    EXCHANGE_LOOP_AB = 2,
};

struct Exchange {
    s32 delta;
    ExchangeType type;
    union {
        struct { // edge exchange
            s32 i_from, i_to;
            s32 j_from, j_to;
        };
        struct { // node exchange
            s32 i_prev, i_curr, i_next;
            s32 j_prev, j_curr, j_next;
        };
    };
};

bool operator <(const Exchange& lhs, const Exchange& rhs) {
    return lhs.delta < rhs.delta;
}

void print_move(Exchange move) {
    bool node_move = move.type == EXCHANGE_LOOP_AB;
    printf("%s move:\t", node_move ? "node" : "edge");
    printf("\tdelta = %d\t", move.delta);
    if (node_move) {
        printf("\ti = %d\t", move.i_curr);
        printf("\tj = %d\n", move.j_curr);
    } else {
        printf("\ti = (%d; %d)\t", move.i_from, move.i_to);
        printf("\tj = (%d; %d)\n", move.j_from, move.j_to);
    }
}

Solution neighbour_search_edge(GraphMatrix graph, Solution init, bool greedy = false) {
    Solution result = init;

    std::vector<s32> *loops[2];
    loops[0] = &result.loop_a;
    loops[1] = &result.loop_b;

    while (true) {
        NodeExchange best = {};
        NodeExchange best_cross_loop = {};
        std::vector<s32> *best_loop = NULL;

        s32 loop_offset = rand() & 1;
        for (s32 loop_i = 0; loop_i < 2; loop_i++) {
            auto& loop = *loops[(loop_i + loop_offset) & 1];
            s32 loop_count = loop.size();
            s32 offset = rand() % loop_count;
            for (s32 i = 0; i < loop_count; i++) {
                s32 i_from = loop[(i+offset) % loop_count];
                s32 i_to = loop[(i+1+offset) % loop_count];
                s32 i_cost = graph.get(i_from, i_to);
                for (s32 j = i+1; j < loop_count; j++) {
                    s32 j_from = loop[(j+offset) % loop_count];
                    s32 j_to = loop[(j+1+offset) % loop_count];
                    s32 j_cost = graph.get(j_from, j_to);

                    s32 new_i_cost = graph.get(j_from, i_from);
                    s32 new_j_cost = graph.get(j_to, i_to);
                    s32 delta = new_i_cost + new_j_cost - i_cost - j_cost;
                    if (delta < best.delta) {
                        best.delta = delta;
                        best.i = (i+offset) % loop_count;
                        best.j = (j+offset) % loop_count;
                        best_loop = &loop;
                        if (greedy) goto apply_change;
                    }
                }
            }
        }

        best_cross_loop = best_node_exchange(graph, result.loop_a, result.loop_b, greedy);
        if (best_cross_loop.delta < best.delta) {
            best = best_cross_loop;
            best_loop = NULL;
        }

apply_change:
        if (best.delta < 0) {
            if (best_loop) {
                s32 edge_i = (best.i < best.j) ? best.i : best.j;
                s32 edge_j = (best.i > best.j) ? best.i : best.j;
                auto& loop = *best_loop;
#if 0
                {
                    Exchange move;
                    move.type = EXCHANGE_LOOP_A;
                    move.i_from = loop[edge_i];
                    move.i_to = loop[(edge_i+1+loop.size()) % loop.size()];
                    move.j_from = loop[edge_j];
                    move.j_to = loop[(edge_j+1+loop.size()) % loop.size()];
                    move.delta = best.delta;
                    print_move(move);
                }
#endif
                s32 between = (edge_i+1 + edge_j) / 2;
                for (s32 i = edge_i+1; i <= between; i++) {
                    s32 j = edge_j + edge_i+1 - i;
                    s32 t = loop[i];
                    loop[i] = loop[j];
                    loop[j] = t;
                }
            } else {
                s32 t = result.loop_a[best.i];
                result.loop_a[best.i] = result.loop_b[best.j];
                result.loop_b[best.j] = t;
            }
        } else {
            break;
        }
    }
    return result;
}

void shuffle(std::vector<s32>& v) {
    for (s32 i = 0; i < v.size(); i++) {
        s32 j = i + rand() % (v.size() - i);
        s32 t = v[i];
        v[i] = v[j];
        v[j] = t;
    }
}

Solution random_solution(s32 node_count) {
    std::vector<s32> nodes(node_count);
    for (s32 i = 0; i < node_count; i++) {
        nodes[i] = i;
    }
    shuffle(nodes);
    Solution result;
    result.loop_a = std::vector<s32>(node_count / 2 + (node_count&1));
    result.loop_b = std::vector<s32>(node_count / 2);
    for (s32 i = 0; i < result.loop_a.size(); i++) {
        result.loop_a[i] = nodes[i];
    }
    s32 base = result.loop_a.size();
    for (s32 i = 0; i < result.loop_b.size(); i++) {
        result.loop_b[i] = nodes[base + i];
    }
    return result;
}

Solution random_walk(GraphMatrix graph, Solution init, s32 steps) {
    auto current_solution = init;

    std::vector<s32> *loops[2];
    loops[0] = &current_solution.loop_a;
    loops[1] = &current_solution.loop_b;

    Solution best_solution = init;
    auto best_score = score(graph, init);

    for (s32 i = 0; i < steps; i++) {
        s32 move_type = rand() % 3;
        switch(move_type) {
            case 0: {
                // switch two edges from a random loop
                auto& loop = *loops[rand()&1];
                s32 edge_i = rand() % loop.size();
                s32 edge_j = rand() % loop.size();
                while (abs(edge_i - edge_j) < 3)
                    edge_j = rand() % loop.size();
                if (edge_i > edge_j) {
                    s32 t = edge_i;
                    edge_i = edge_j;
                    edge_j = t;
                }
                s32 between = (edge_i+1 + edge_j) / 2;
                for (s32 i = edge_i+1; i <= between; i++) {
                    s32 j = edge_j + edge_i+1 - i;
                    s32 t = loop[i];
                    loop[i] = loop[j];
                    loop[j] = t;
                }
            } break;
            case 1: {
                // switch two nodes from a random loop
                auto& loop = *loops[rand()&1];
                s32 node_i = rand() % loop.size();
                s32 node_j = rand() % loop.size();
                while (node_i == node_j)
                    node_j = rand() % loop.size();
                s32 t = loop[node_i];
                loop[node_i] = loop[node_j];
                loop[node_j] = t;
            } break;
            case 2: {
                // switch two nodes from different loops
                s32 node_i = rand() % current_solution.loop_a.size();
                s32 node_j = rand() % current_solution.loop_b.size();
                s32 t = current_solution.loop_a[node_i];
                current_solution.loop_a[node_i] = current_solution.loop_b[node_j];
                current_solution.loop_b[node_j] = t;
            } break;

            auto new_score = score(graph, current_solution);
            if (new_score < best_score) {
                best_score = new_score;
                best_solution = current_solution;
            }
        }
    }
    return best_solution;
}

typedef Solution (*LS_Solver)(GraphMatrix graph, Solution init, bool greedy);

struct Node {
    s32 next;
    s32 prev;
};

// doubly linked list representing a single solution
struct SolutionList {
    std::vector<Node> nodes;
    std::vector<bool> loop; // loop[node] == 0 if node in loop_a, else 1

    Node& operator [](s32 idx) {
        return nodes[idx];
    }
};

s32 find(std::vector<bool>& v, bool value) {
    for (s32 i = 0; i < v.size(); i++) {
        if (v[i] == value) return i;
    }
    return v.size();
}

// convert a visit-order node list to node-order doubly linked list
SolutionList to_linked_list(Solution solution) {
    SolutionList result;
    s32 dim = solution.loop_a.size() + solution.loop_b.size();
    result.nodes.resize(dim);
    result.loop.resize(dim);
    s32 last = solution.loop_a.back();
    for (s32 next : solution.loop_a) {
        result[last].next = next;
        result[next].prev = last;
        result.loop[next] = 0;
        last = next;
    }
    last = solution.loop_b.back();
    for (s32 next : solution.loop_b) {
        result[last].next = next;
        result[next].prev = last;
        result.loop[next] = 1;
        last = next;
    }
    return result;
}

Solution to_visit_list(SolutionList list) {
    Solution result;
    s32 dim = list.nodes.size();

    result.loop_a.resize(dim/2 + (dim&1));
    s32 start = find(list.loop, 0);
    s32 node = start;
    s32 i = 0;
    do {
        node = list[node].next;
        result.loop_a[i] = node;
        i++;
    } while (node != start);

    result.loop_b.resize(dim/2);
    start = find(list.loop, 1);
    node = start;
    i = 0;
    do {
        node = list[node].next;
        result.loop_b[i] = node;
        i++;
    } while (node != start);

    assert(i == result.loop_b.size());

    return result;
}

inline s32 node_exchange_delta(GraphMatrix graph, Exchange move) {
    s32 new_i_next = (move.j_next != move.i_curr) ? move.j_next : move.j_curr;
    s32 new_j_next = (move.i_next != move.j_curr) ? move.i_next : move.i_curr;
    s32 i_cost = graph.get(move.i_curr, move.i_prev) +
                 graph.get(move.i_curr, move.i_next);
    s32 j_cost = graph.get(move.j_curr, move.j_prev) +
                 graph.get(move.j_curr, move.j_next);
    s32 new_i_cost = graph.get(move.j_curr, move.i_prev) +
                     graph.get(move.j_curr, new_j_next);
    s32 new_j_cost = graph.get(move.i_curr, move.j_prev) +
                     graph.get(move.i_curr, new_i_next);
    return new_i_cost + new_j_cost - i_cost - j_cost;
}

void print_loops(SolutionList& list) {
    printf("A: ");
    s32 start = find(list.loop, 0);
    s32 node = start;
    s32 count_a = 0;
    do {
        node = list[node].next;
        count_a++;
        printf("%d-", node);
    } while (node != start && count_a < list.nodes.size());

    printf("\nB: ");
    start = find(list.loop, 1);
    node = start;
    s32 count_b = 0;
    do {
        node = list[node].next;
        count_b++;
        printf("%d-", node);
    } while (node != start && count_b < list.nodes.size());
    printf("\nA count: %d\tB count: %d\n", count_a, count_b);
}

// returns the cost delta (should always be <= 0)
s32 neighbour_search_edge_cache(GraphMatrix graph, SolutionList& list) {
    // create a sorted list of all moves
    std::multiset<Exchange> moves;

    // determine some starting nodes for iteration in each loop
    s32 starting_nodes[2];
    starting_nodes[0] = 0;
    s32 loop_0 = list.loop[0];
    starting_nodes[1] = 1;
    while (list.loop[starting_nodes[1]] == loop_0) starting_nodes[1]++;

    // add all the edge exchange moves
    for (s32 loop_i = 0; loop_i < 2; loop_i++) {
        s32 start = starting_nodes[loop_i];
        Exchange move;
        move.type = (ExchangeType)(loop_i^loop_0);
        move.i_from = start;
        move.i_to = list[start].next;
        while (move.i_to != start) {
            s32 i_cost = graph.get(move.i_from, move.i_to);
            move.j_from = move.i_to;
            move.j_to = list[move.j_from].next;
            while (move.j_from != start) {
                s32 j_cost = graph.get(move.j_from, move.j_to);
                s32 new_i_cost = graph.get(move.j_from, move.i_from);
                s32 new_j_cost = graph.get(move.j_to, move.i_to);
                move.delta = new_i_cost + new_j_cost - i_cost - j_cost;
                if (move.delta < 0) moves.insert(move);

                // second orientation
                s32 t = move.j_from;
                move.j_from = move.j_to;
                move.j_to = t;
                new_i_cost = graph.get(move.j_from, move.i_from);
                new_j_cost = graph.get(move.j_to, move.i_to);
                move.delta = new_i_cost + new_j_cost - i_cost - j_cost;
                if (move.delta < 0) moves.insert(move);

                // j_to and j_from were already swapped above
                move.j_to = list[move.j_from].next;
            }
            move.i_from = move.i_to;
            move.i_to = list[move.i_from].next;
        }
    }

    // add all the node exchange moves
    {
        Exchange move;
        move.type = EXCHANGE_LOOP_AB;
        move.i_curr = 0;
        do {
            move.i_prev = list[move.i_curr].prev;
            move.i_next = list[move.i_curr].next;

            move.j_curr = starting_nodes[1];
            do {
                move.j_prev = list[move.j_curr].prev;
                move.j_next = list[move.j_curr].next;

                move.delta = node_exchange_delta(graph, move);
                if (move.delta < 0) moves.insert(move);

                move.j_curr = move.j_next;
            } while (move.j_curr != starting_nodes[1]);
            move.i_curr = move.i_next;
        } while (move.i_curr != 0);
    }

    s32 result_delta = 0;

    while (!moves.empty()) {
        s32 modified_nodes[6];
        s32 modified_nodes_count = 0;
        s32 new_edges[4];
        s32 new_edges_count = 0;
        bool move_applied = false;
        for (auto it = moves.begin(); it != moves.end();) {
            bool do_remove = false;
            bool do_break = false;
            auto move = *it;
            //print_move(move);
            switch (move.type) {
                case EXCHANGE_LOOP_A:
                case EXCHANGE_LOOP_B: {
                    bool same_loop = list.loop[move.i_from] == list.loop[move.j_from];
                    if (!same_loop) {
                        do_remove = true;
                        break;
                    }

                    bool i_correct = list[move.i_from].next == move.i_to;
                    bool i_reverse = list[move.i_from].prev == move.i_to;
                    bool j_correct = list[move.j_from].next == move.j_to;
                    bool j_reverse = list[move.j_from].prev == move.j_to;
                    if (!(i_correct || i_reverse) || !(j_correct || j_reverse)) {
                        do_remove = true;
                        break;
                    }
                    if (i_correct != j_correct) {
                        // the move may be applicable later but it is not now
                        // so we leave it in, but go to the next one
                        break; 
                    }

                    // if we got here, the move is applicable

                    //print_move(move);
                    //printf("i, j correct = (%d, %d)\n", i_correct, j_correct);
                    if (i_reverse) {
                        s32 t = move.i_to;
                        move.i_to = move.i_from;
                        move.i_from = t;
                        t = move.j_to;
                        move.j_to = move.j_from;
                        move.j_from = t;
                    }

                    list[move.i_from].next = move.j_from;
                    list[move.j_from].next = list[move.j_from].prev;
                    list[move.j_from].prev = move.i_from;

                    list[move.j_to].prev = move.i_to;
                    list[move.i_to].prev = list[move.i_to].next;
                    list[move.i_to].next = move.j_to;

                    s32 node = list[move.i_to].prev;
                    while (node != move.j_from) {
                        s32 t = list[node].next;
                        list[node].next = list[node].prev;
                        node = list[node].prev = t;
                    }
                    do_remove = true;
                    do_break = true;
                    modified_nodes[0] = move.i_from;
                    modified_nodes[1] = move.j_from;
                    modified_nodes[2] = move.i_to;
                    modified_nodes[3] = move.j_to;
                    modified_nodes_count = 4;
                    new_edges[0] = move.i_from;
                    new_edges[1] = move.i_to;
                    new_edges_count = 2;
                } break;
                case EXCHANGE_LOOP_AB: {
                    bool i_correct = (list[move.i_prev].next == move.i_curr) &&
                                     (list[move.i_curr].next == move.i_next);
                    bool i_reverse = (list[move.i_prev].prev == move.i_curr) &&
                                     (list[move.i_curr].prev == move.i_next);
                    bool j_correct = (list[move.j_prev].next == move.j_curr) &&
                                     (list[move.j_curr].next == move.j_next);
                    bool j_reverse = (list[move.j_prev].prev == move.j_curr) &&
                                     (list[move.j_curr].prev == move.j_next);
                    if (!(i_correct || i_reverse) || !(j_correct || j_reverse)) {
                        do_remove = true;
                        break;
                    }

                    // if we got here, the move is applicable

                    //print_move(move);

                    if (i_reverse) {
                        s32 t = move.i_next;
                        move.i_next = move.i_prev;
                        move.i_prev = t;
                    }
                    if (j_reverse) {
                        s32 t = move.j_next;
                        move.j_next = move.j_prev;
                        move.j_prev = t;
                    }

                    list.loop[move.i_curr] = !list.loop[move.i_curr];
                    list.loop[move.j_curr] = !list.loop[move.j_curr];

                    list[move.i_prev].next = move.j_curr;
                    list[move.j_curr].prev = move.i_prev;
                    list[move.j_curr].next = move.i_next;
                    list[move.i_next].prev = move.j_curr;

                    list[move.j_prev].next = move.i_curr;
                    list[move.i_curr].prev = move.j_prev;
                    list[move.i_curr].next = move.j_next;
                    list[move.j_next].prev = move.i_curr;

                    do_remove = true;
                    do_break = true;
                    modified_nodes[0] = move.i_prev;
                    modified_nodes[1] = move.i_curr;
                    modified_nodes[2] = move.i_next;
                    modified_nodes[3] = move.j_prev;
                    modified_nodes[4] = move.j_curr;
                    modified_nodes[5] = move.j_next;
                    modified_nodes_count = 6;
                    new_edges[0] = move.i_prev;
                    new_edges[1] = move.i_curr;
                    new_edges[2] = move.j_prev;
                    new_edges[3] = move.j_curr;
                    new_edges_count = 4;
                } break;
            }

            // C++ spec is not great
            auto it_copy = it;
            it++;
            if (do_remove) {
                //puts("removing");
                moves.erase(it_copy);
            }
            if (do_break) {
                move_applied = true;
                result_delta += move.delta;
                break;
            }
        }

        if (!move_applied) break;

        //print_loops(list);

        //puts("adding edge moves");
        Exchange move;
        for (s32 i = 0; i < new_edges_count; i++) {
            move.type = (ExchangeType)!!list.loop[new_edges[i]];
            move.i_from = new_edges[i];
            move.i_to = list[move.i_from].next;
            s32 i_cost = graph.get(move.i_from, move.i_to);
            move.j_from = move.i_to;
            while (move.j_from != move.i_from) {
                move.j_to = list[move.j_from].next;
                s32 j_cost = graph.get(move.j_from, move.j_to);
                s32 new_i_cost = graph.get(move.j_from, move.i_from);
                s32 new_j_cost = graph.get(move.j_to, move.i_to);
                move.delta = new_i_cost + new_j_cost - i_cost - j_cost;
                if (move.delta < 0) {
                    //puts("adding a move:");
                    //print_move(move);
                    moves.insert(move);
                }

                // second orientation
                s32 t = move.j_from;
                move.j_from = move.j_to;
                move.j_to = t;
                new_i_cost = graph.get(move.j_from, move.i_from);
                new_j_cost = graph.get(move.j_to, move.i_to);
                move.delta = new_i_cost + new_j_cost - i_cost - j_cost;
                if (move.delta < 0) moves.insert(move);
            }
        }

        //puts("adding node moves");
        move.type = EXCHANGE_LOOP_AB;
        for (s32 i = 0; i < modified_nodes_count; i++) {
            move.i_curr = modified_nodes[i];
            move.i_prev = list[move.i_curr].prev;
            move.i_next = list[move.i_curr].next;

            bool first_node_loop = list.loop[move.i_curr];
            s32 j_start = find(list.loop, !first_node_loop);
            move.j_curr = j_start;
            do {
                move.j_prev = list[move.j_curr].prev;
                move.j_next = list[move.j_curr].next;
                move.delta = node_exchange_delta(graph, move);
                if (move.delta < 0) {
                    //puts("adding a move:");
                    //print_move(move);
                    moves.insert(move);
                }
                move.j_curr = move.j_next;
            } while (move.j_curr != j_start);
        }
    }

    return result_delta;
}

Solution neighbour_search_edge_cache(GraphMatrix graph, Solution init) {
    auto list = to_linked_list(init);
    neighbour_search_edge_cache(graph, list);
    return to_visit_list(list);
}

struct NodeLoopPosition {
    s32 loop;
    s32 position;
};

NodeExchange best_node_exchange_candidates(GraphMatrix graph, std::vector<s32>& loop_i, std::vector<s32>& loop_j, std::vector<std::pair<s32, s32>>& best_neighbours, std::vector<NodeLoopPosition>& node_loop_positions, s32 max_neighbours) {
    NodeExchange result = {};
    s32 loop_i_count = loop_i.size();
    s32 loop_j_count = loop_j.size();
    s32 offset_i = rand() % loop_i_count;
    s32 offset_j = rand() % loop_j_count;
    for (s32 i = 0; i < loop_i.size(); i++) {
        s32 i_prev = loop_i[(i-1+loop_i_count+offset_i) % loop_i_count];
        s32 i_curr = loop_i[(i+offset_i) % loop_i_count];
        s32 i_next = loop_i[(i+1+offset_i) % loop_i_count];
        s32 i_cost = graph.get(i_curr, i_prev) + graph.get(i_curr, i_next);

        s32 neighbours_visited = 0;
        for (s32 k = 1; k < graph.dim; k++) {
            std::pair<s32, s32> candidate = best_neighbours[i_curr*graph.dim+k];
            s32 j_curr = candidate.first;
            if (node_loop_positions[j_curr].loop == 1) {
                neighbours_visited++;
                s32 j = node_loop_positions[j_curr].position;
                s32 j_prev = loop_j[(j-1 + loop_j_count) % loop_j_count];
                s32 j_next = loop_j[(j+1) % loop_j_count];
                s32 j_cost = graph.get(j_curr, j_prev) + graph.get(j_curr, j_next);
                s32 new_i_cost = graph.get(j_curr, i_prev) + graph.get(j_curr, i_next);
                s32 new_j_cost = graph.get(i_curr, j_prev) + graph.get(i_curr, j_next);
                s32 delta = new_i_cost + new_j_cost - i_cost - j_cost;
                if (delta < result.delta) {
                    result.delta = delta;
                    result.i = (i+offset_i) % loop_i_count;
                    result.j = j;
                }
            }
            if (neighbours_visited==max_neighbours) {
                break;
            }
        }
    }

    return result;
}

Solution neighbour_search_edge_candidates(GraphMatrix graph, Solution init, s32 max_neighbours) {
    Solution result = init;

    std::vector<s32> *loops[2];
    loops[0] = &result.loop_a;
    loops[1] = &result.loop_b;

    std::vector<NodeLoopPosition> node_loop_positions(graph.dim);

    for (s32 loop_i = 0; loop_i < 2; loop_i++) {
        auto& loop = *loops[loop_i];
        for (s32 i = 0; i < loop.size(); i++) {
            node_loop_positions[loop[i]] = NodeLoopPosition { .loop = loop_i, .position = i};
        }
    }

    std::vector<std::pair<s32, s32>> best_neighbours;
    best_neighbours.reserve(graph.dim*graph.dim);

    std::vector<std::pair<s32, s32>> node_neighbours(graph.dim);
    // prepare closest neighbours (node, cost)
    for (s32 node=0; node<graph.dim; node++) {
        for (s32 other_node=0; other_node<graph.dim; other_node++) {
            node_neighbours[other_node] = std::pair<s32, s32>(other_node, graph.get(node, other_node));
        }
        std::sort(node_neighbours.begin(), node_neighbours.end(), [](auto &left, auto &right) {
            return left.second < right.second;
        });
        std::copy(node_neighbours.begin(), node_neighbours.end(), std::back_inserter(best_neighbours));
    }

    while (true) {
        NodeExchange best = {};
        NodeExchange best_cross_loop = {};
        std::vector<s32> *best_loop = NULL;
        bool should_remove_successor = true;

        s32 loop_offset = rand() & 1;
        for (s32 loop_i = 0; loop_i < 2; loop_i++) {
            auto& loop = *loops[(loop_i + loop_offset) & 1];
            s32 loop_count = loop.size();
            s32 offset = rand() % loop_count;
            for (s32 i = 0; i < loop_count; i++) {

                s32 i_node = loop[(i+offset) % loop_count];
                s32 i_pred = loop[(i-1+offset) % loop_count];
                s32 i_succ = loop[(i+1+offset) % loop_count];

                s32 i_pred_cost = graph.get(i_pred, i_node);
                s32 i_succ_cost = graph.get(i_node, i_succ);

                s32 neighbors_visited = 0;
                for (s32 j = 1; j < graph.dim; j++) {
                    std::pair<s32, s32> candidate = best_neighbours[i_node*graph.dim+j];
                    s32 candidate_node = candidate.first;
                    if (node_loop_positions[candidate_node].loop == (loop_i + loop_offset) & 1) {
                        s32 candidate_position = node_loop_positions[candidate_node].position;
                        s32 candidate_pred = loop[(candidate_position-1 + loop_count) % loop_count];
                        s32 candidate_succ = loop[(candidate_position+1) % loop_count];
                        if (candidate_pred != i_node && candidate_succ != i_node) {
                            neighbors_visited++;
                            s32 candidate_pred_cost = graph.get(candidate_pred, candidate_node);
                            s32 candidate_succ_cost = graph.get(candidate_node, candidate_succ);
                            s32 candidate_cost = candidate.second;

                            s32 succ_to_succ_cost = graph.get(i_succ, candidate_succ);
                            s32 pred_to_pred_cost = graph.get(i_pred, candidate_pred);

                            s32 succ_removal_delta = candidate_cost + succ_to_succ_cost - i_succ_cost - candidate_succ_cost;
                            s32 pred_removal_delta = candidate_cost + pred_to_pred_cost - i_pred_cost - candidate_pred_cost;

                            s32 delta = succ_removal_delta < pred_removal_delta ? succ_removal_delta : pred_removal_delta;

                            if (delta < best.delta) {
                                best.delta = delta;
                                best.i = (i+offset) % loop_count;
                                best.j = candidate_position;
                                best_loop = &loop;
                                should_remove_successor = delta == succ_removal_delta;
                            }
                        
                            if (neighbors_visited == max_neighbours) {
                                break;
                            }
                        }
                    }           
                }
            }
        }

        best_cross_loop = best_node_exchange_candidates(graph, result.loop_a, result.loop_b, best_neighbours, node_loop_positions, max_neighbours);
        if (best_cross_loop.delta < best.delta) {
            best = best_cross_loop;
            best_loop = NULL;
        }

        if (best.delta < 0) {
            if (best_loop) {
                auto& loop = *best_loop;
                s32 edge_i = (best.i < best.j) ? best.i : best.j;
                s32 edge_j = (best.i > best.j) ? best.i : best.j;
                if (should_remove_successor) {    
                    s32 between = (edge_i+1 + edge_j) / 2;
                    for (s32 i = edge_i+1; i <= between; i++) {
                        s32 j = edge_j + edge_i+1 - i;
                        s32 t = loop[i];

                        node_loop_positions[loop[i]].position = j;
                        node_loop_positions[loop[j]].position = i;

                        loop[i] = loop[j];
                        loop[j] = t;
                    }
                } else {
                    s32 loop_count = loop.size();
                    s32 between = ((edge_i-edge_j+loop_count) % loop_count)/2;
                    for (s32 i = 0; i < between; i++) {
                        int loop_i = (edge_j+i) % loop_count;
                        int j = (edge_i-1-i+loop_count) % loop_count;
                        s32 t = loop[loop_i];

                        node_loop_positions[loop[loop_i]].position = j;
                        node_loop_positions[loop[j]].position = loop_i;

                        loop[loop_i] = loop[j];
                        loop[j] = t;
                    }
                }
            } else {
                s32 t = result.loop_a[best.i];

                node_loop_positions[t].loop=1;
                node_loop_positions[t].position=best.j;
                node_loop_positions[result.loop_b[best.j]].loop = 0;
                node_loop_positions[result.loop_b[best.j]].position = best.i;

                result.loop_a[best.i] = result.loop_b[best.j];
                result.loop_b[best.j] = t;
            }
        } else {
            break;
        }
    }
    return result;
}

Solution msls(GraphMatrix graph, s32 iterations) {
    auto best_score = INT32_MAX;
    Solution ls, best_solution;
    for (s32 i = 0; i < iterations; i++) {
        auto init = greedy_loop(graph);
        ls = neighbour_search_edge_cache(graph, init);
        auto new_score = score(graph, ls);
        if (new_score < best_score) {
            best_score = new_score;
            best_solution = ls;
        }
    }
    return best_solution;
}

Solution ils1(GraphMatrix graph, Solution init, s32 time_limit, s32 node_switch, s32 edge_switch) {
    auto current_solution = init;

    std::vector<s32> *loops[2];
    loops[0] = &current_solution.loop_a;
    loops[1] = &current_solution.loop_b;

    Solution best_solution = init;
    auto best_score = score(graph, init);
    
    u64 start = stm_now();
    while (stm_ms(stm_since(start)) < time_limit) {
        for (s32 i = 0; i < edge_switch; i++) {
                // switch two edges from a random loop
                auto& loop = *loops[rand()&1];
                s32 edge_i = rand() % loop.size();
                s32 edge_j = rand() % loop.size();
                while (abs(edge_i - edge_j) < 3)
                    edge_j = rand() % loop.size();
                if (edge_i > edge_j) {
                    s32 t = edge_i;
                    edge_i = edge_j;
                    edge_j = t;
                }
                s32 between = (edge_i+1 + edge_j) / 2;
                for (s32 i = edge_i+1; i <= between; i++) {
                    s32 j = edge_j + edge_i+1 - i;
                    s32 t = loop[i];
                    loop[i] = loop[j];
                    loop[j] = t;
                }
        }
        for (s32 i = 0; i < node_switch; i++) {
                s32 node_i = rand() % current_solution.loop_a.size();
                s32 node_j = rand() % current_solution.loop_b.size();
                s32 t = current_solution.loop_a[node_i];
                current_solution.loop_a[node_i] = current_solution.loop_b[node_j];
                current_solution.loop_b[node_j] = t;
        }

        current_solution = neighbour_search_edge_cache(graph, current_solution);
        
        auto new_score = score(graph, current_solution);
        if (new_score < best_score) {
            best_score = new_score;
            best_solution = current_solution;
        }
    }
    return best_solution;
}

Solution ils2(GraphMatrix graph, Solution init, bool local_search, s32 time_limit, s32 remove_nodes_percentage) {
    auto current_solution = init;

    std::vector<s32> *loops[2];
    loops[0] = &current_solution.loop_a;
    loops[1] = &current_solution.loop_b;

    Solution best_solution = init;
    auto best_score = score(graph, init);
    s32 nodes_to_remove = remove_nodes_percentage*graph.dim/100;
    
    u64 start = stm_now();
    while (stm_ms(stm_since(start)) < time_limit) {
        for (s32 i = 0; i < nodes_to_remove; i++) {
            auto& loop = *loops[rand()&1];
            vector_remove_idx(loop, rand() % loop.size());
        }

        current_solution = greedy_loop(graph, &current_solution);
        if (local_search) {
            current_solution = neighbour_search_edge_cache(graph, current_solution);
        }
        
        auto new_score = score(graph, current_solution);
        if (new_score < best_score) {
            best_score = new_score;
            best_solution = current_solution;
        }
    }
    return best_solution;
}

void greedy_fix(GraphMatrix graph, SolutionList& result, std::vector<s32>& available, std::vector<bool>& avail, s32 loop_count[2]) {
    s32 max_loop_a_size = graph.dim / 2 + (graph.dim & 1);
    s32 max_loop_b_size = graph.dim - max_loop_a_size;

    while (available.size() > 0) {
        s32 best_node = 0;
        s32 best_avail = 0;
        s32 best_cost = INT32_MAX;

        bool can_expand[2] = { loop_count[0] < max_loop_a_size, loop_count[1] < max_loop_b_size };
        for (s32 node = 0; node < graph.dim; node++) {
            if (avail[node]) continue;
            auto l = result.loop[node];
            if (!can_expand[l]) continue;
            s32 next = result[node].next;
            s32 old_cost = graph.get(node, next);

            s32 best_inter = 0;
            s32 best_inter_cost = INT32_MAX;
            for (s32 i = 0; i < available.size(); i++) {
                s32 idx = available[i];
                s32 cost = graph.get(node, idx) + graph.get(next, idx) - old_cost;
                if (cost < best_inter_cost) {
                    best_inter = i;
                    best_inter_cost = cost;
                }
            }

            if (best_inter_cost < best_cost) {
                best_node = node;
                best_avail = best_inter;
                best_cost = best_inter_cost;
            }
        }

        s32 new_node = available[best_avail];
        vector_remove_idx(available, best_avail);
        avail[new_node] = false;

        s32 next = result[best_node].next;
        result[best_node].next = new_node;
        result[new_node].prev = best_node;
        result[new_node].next = next;
        result[next].prev = new_node;
        result.loop[new_node] = result.loop[best_node];
        loop_count[result.loop[best_node]]++;
    }
}

void regret_fix(GraphMatrix graph, SolutionList& result, std::vector<s32>& available, std::vector<bool>& avail, s32 loop_count[2]) {
    s32 max_loop_a_size = graph.dim / 2 + (graph.dim & 1);
    s32 max_loop_b_size = graph.dim - max_loop_a_size;

    while (available.size() > 0) {
        s32 best_node = 0;
        s32 best_avail = 0;
        s32 best_cost = INT32_MAX;

        bool can_expand[2] = { loop_count[0] < max_loop_a_size, loop_count[1] < max_loop_b_size };
        for (s32 node = 0; node < graph.dim; node++) {
            if (avail[node]) continue;
            auto l = result.loop[node];
            if (!can_expand[l]) continue;
            s32 next = result[node].next;
            s32 old_cost = graph.get(node, next);

            s32 best_inter = 0;
            s32 best_inter_cost = INT32_MAX;
            s32 second_inter_cost = INT32_MAX;
            for (s32 i = 0; i < available.size(); i++) {
                s32 idx = available[i];
                s32 cost = graph.get(node, idx) + graph.get(next, idx) - old_cost;
                if (cost < best_inter_cost) {
                    best_inter = i;
                    second_inter_cost = best_inter_cost;
                    best_inter_cost = cost;
                } else if (cost < second_inter_cost) second_inter_cost = cost;
            }
            best_inter_cost = best_inter_cost - second_inter_cost/256;

            if (best_inter_cost < best_cost) {
                best_node = node;
                best_avail = best_inter;
                best_cost = best_inter_cost;
            }
        }

        s32 new_node = available[best_avail];
        vector_remove_idx(available, best_avail);
        avail[new_node] = false;

        s32 next = result[best_node].next;
        result[best_node].next = new_node;
        result[new_node].prev = best_node;
        result[new_node].next = next;
        result[next].prev = new_node;
        result.loop[new_node] = result.loop[best_node];
        loop_count[result.loop[best_node]]++;
    }
}

s32 score(GraphMatrix graph, SolutionList& list) {
    s32 result = 0;
    for (s32 i = 0; i < graph.dim; i++) {
        result += graph.get(i, list[i].next);
    }
    return result;
}

void regret_solve(GraphMatrix graph, SolutionList& list, s32 node_a, s32 node_b) {
    std::vector<s32> available(graph.dim-2);
    std::vector<bool> avail(graph.dim);

    list.nodes.resize(graph.dim);
    list.loop.resize(graph.dim);
    if (node_a > node_b) {
        s32 t = node_a;
        node_a = node_b;
        node_b = t;
    }
    list[node_a].next = list[node_a].prev = node_a;
    list[node_b].next = list[node_b].prev = node_b;
    list.loop[node_a] = 0;
    list.loop[node_b] = 1;
    for (s32 i = 0; i < graph.dim-2; i++) {
        available[i] = i + (i >= node_a) + (i+1 >= node_b);
    }
    avail.flip();
    avail[node_a] = avail[node_b] = false;
    s32 loop_count[2] = { 1, 1 };
    regret_fix(graph, list, available, avail, loop_count);
}

Solution evolve(GraphMatrix graph, s32 time_limit, bool ls) {
    constexpr s32 population_count = 20;
    SolutionList elite[population_count];
    s32 scores[population_count];
    std::unordered_set<s32> score_set(population_count);

    u64 timestart = stm_now();
    s32 added = 0;
    s32 random[population_count];
    for (s32 i = 0; i < population_count; i++) {
        random[i] = rand();
    }
#pragma omp parallel for
    for (s32 i = 0; i < population_count; i++) {
        auto& list = elite[i];
        bool success = false;
        while (!success) {
            s32 node_a = random[i] % (population_count);
            // c std lib rand is not multithreaded by default!
            random[i] = ((random[i] >> 17) ^ (random[i] << 11) ^ (random[i] >> 5)) & 0x7fffffff;
            s32 node_b = random[i] % (population_count-1);
            node_b += (node_b == node_a);
            regret_solve(graph, list, node_a, node_b);
            scores[i] = score(graph, list);
            scores[i] += neighbour_search_edge_cache(graph, list);
#pragma omp critical
            {
                auto [_, inserted] = score_set.insert(scores[i]);
                success = inserted;
            }
        }
    }
    score_set.insert(scores, scores+population_count);
    //printf("init with %zd copies\n", population_count - score_set.size());
    //printf("init took %.3f\n", stm_ms(stm_since(timestart)));

    while (stm_ms(stm_since(timestart)) < time_limit) {
    //s32 iters = 0;
    //while (iters++ < 200) {
        s32 parent_a = rand() % population_count;
        s32 parent_b = rand() % (population_count-1);
        parent_b += (parent_b == parent_a);
        SolutionList child = elite[parent_a];
        auto& pa = elite[parent_a];
        auto& pb = elite[parent_b];

        std::vector<s32> removed_nodes;
        removed_nodes.reserve(graph.dim);
        std::vector<bool> available(graph.dim);
        s32 loop_count[2] = { graph.dim / 2 + (graph.dim & 1), graph.dim / 2 };
        for (s32 i = 0; i < graph.dim; i++) {
            auto a = pa[i];
            auto b = pb[i];
            auto l = child.loop[i];
            //if (loop_count[l] > 1 && (a.next != b.next || a.prev != b.prev)) {
            //if (loop_count[l] > 1 && (a.next != b.next)) {
            if (loop_count[l] > 1 && (a.next != b.next)) {
            //if (rand()%100 < 70) {
                auto n = child[i];
                child[n.prev].next = n.next;
                child[n.next].prev = n.prev;
                removed_nodes.push_back(i);
                available[i] = true;
                loop_count[l]--;
            }
        }

#if 0
        //printf("a: %d\tb: %d\n", loop_count[0], loop_count[1]);
        printf("removing %d nodes\n", removed_nodes.size());
#endif
        //greedy_fix(graph, child, removed_nodes, available, loop_count);
        regret_fix(graph, child, removed_nodes, available, loop_count);
        s32 new_score = score(graph, child);

        if (ls) new_score += neighbour_search_edge_cache(graph, child);
        s32 worst = 0;
        for (s32 i = 0; i < population_count; i++) {
            if (scores[i] > scores[worst]) worst = i;
        }
        if (scores[worst] > new_score) {
            auto [_, inserted] = score_set.insert(new_score);
            if (inserted) {
                score_set.erase(scores[worst]);
                elite[worst] = std::move(child);
                scores[worst] = new_score;
            }
        }
    }

    s32 best = 0;
    for (s32 i = 0; i < population_count; i++) {
        if (scores[i] < scores[best]) best = i;
    }
    return to_visit_list(elite[best]);
}

int cmp_s32(const void *a, const void *b) {
    return (*(s32 *)a > *(s32 *)b) - (*(s32 *)a < *(s32 *)b);
}

Solution evolve_mt(GraphMatrix graph, s32 time_limit, bool ls) {
    constexpr s32 parent_count = 16;
    constexpr s32 population_count = parent_count*(parent_count-1)/2 + parent_count;
    SolutionList elite[population_count];
    SolutionList parents[parent_count];
    struct Score {
        s32 score;
        s32 idx;
    };
    Score scores[population_count];
    std::unordered_set<s32> score_set(population_count);

    u64 timestart = stm_now();
    s32 added = 0;
    s32 random[population_count];
    for (s32 i = 0; i < population_count; i++) {
        random[i] = rand();
    }
#pragma omp parallel for
    for (s32 i = 0; i < population_count; i++) {
        auto& list = elite[i];
        bool success = false;
        scores[i].idx = i;
        s32 node_a = random[i] % (population_count);
        // c std lib rand is not multithreaded by default!
        random[i] = ((random[i] >> 17) ^ (random[i] << 11) ^ (random[i] >> 5)) & 0x7fffffff;
        s32 node_b = random[i] % (population_count-1);
        node_b += (node_b == node_a);
        regret_solve(graph, list, node_a, node_b);
        scores[i].score = score(graph, list);
        scores[i].score += neighbour_search_edge_cache(graph, list);
#pragma omp critical
        {
            auto [_, inserted] = score_set.insert(scores[i].score);
            success = inserted;
        }
        if (!success) {
            scores[i].score = INT32_MAX;
        }
    }
    //score_set.insert(scores, scores+population_count);
    //printf("init with %zd copies\n", population_count - score_set.size());
    //printf("init took %.3f\n", stm_ms(stm_since(timestart)));

    struct Parents { s32 a, b; };
    Parents parent_indices[population_count-parent_count];
    s32 parents_added = 0;
    for (s32 i = 0; i < parent_count; i++) {
        for (s32 j = i+1; j < parent_count; j++) {
            parent_indices[parents_added++] = {i, j};
        }
    }

    while (stm_ms(stm_since(timestart)) < time_limit) {
        //
        // find the best solutions from population
        //
        qsort(scores, population_count, sizeof(Score), cmp_s32);
        for (s32 i = 0; i < parent_count; i++) {
            parents[i] = std::move(elite[scores[i].idx]);
            scores[i].idx = i;
        }
        for (s32 i = 0; i < parent_count; i++) {
            elite[i] = std::move(parents[i]);
        }

        //
        // produce new population from parents
        //
#pragma omp parallel for
        for (s32 child_idx = parent_count; child_idx < population_count; child_idx++) {
            Parents p = parent_indices[child_idx-parent_count];
            s32 parent_a = p.a;
            s32 parent_b = p.b;
            SolutionList child = elite[parent_a];
            //auto& pa = elite[parent_a];
            //auto& pb = elite[parent_b];

            std::vector<s32> removed_nodes;
            removed_nodes.reserve(graph.dim);
            std::vector<bool> available(graph.dim);
            s32 loop_count[2] = { graph.dim / 2 + (graph.dim & 1), graph.dim / 2 };
            for (s32 i = 0; i < graph.dim; i++) {
                //auto a = pa[i];
                //auto b = pb[i];
                auto l = child.loop[i];
                //if (loop_count[l] > 1 && (a.next != b.next || a.prev != b.prev)) {
                //if (loop_count[l] > 1 && ((a.next != b.next) || (rand()%100 < 20))) {
                if (loop_count[l] > 1 && (rand()%100 < 60)) {
                    auto n = child[i];
                    child[n.prev].next = n.next;
                    child[n.next].prev = n.prev;
                    removed_nodes.push_back(i);
                    available[i] = true;
                    loop_count[l]--;
                }
            }

            regret_fix(graph, child, removed_nodes, available, loop_count);
            s32 new_score = score(graph, child);

            if (ls) new_score += neighbour_search_edge_cache(graph, child);

            bool success = false;
#pragma omp critical
            {
                auto [_, inserted] = score_set.insert(new_score);
                success = inserted;
            }

            if (success) {
                elite[child_idx] = std::move(child);
                scores[child_idx].idx = child_idx;
                scores[child_idx].score = new_score;
            } else {
                scores[child_idx].score = INT32_MAX;
            }
        }
    }

    s32 best = 0;
    for (s32 i = 0; i < population_count; i++) {
        if (scores[i].score < scores[best].score) best = i;
    }
    return to_visit_list(elite[best]);
}

ExperimentResult run_experiment(Instance instance, const char *method_name, const char *instance_name, bool random_initial_solution, s32 time_limit, s32 method_switch, s32 node_switch=1, s32 edge_switch=5, s32 percentage=20) {
    Solution min_solution, max_solution;
    ExperimentResult experiment_result = {
        .min = INT32_MAX,
        .max = INT32_MIN
        };
    s32 n = 100;
    s32 total_score = 0;
    u64 min_time = INT64_MAX;
    u64 max_time = 0;
    u64 total_time = 0;
    
    for (int i = 0; i < n; i++) {
        Solution solution, initial_solution;
        u64 start = stm_now();
        switch(method_switch) {
            case 0:
                solution = msls(instance.graph, 100);
            case 1:
                initial_solution = random_initial_solution ? random_solution(instance.graph.dim) : greedy_loop(instance.graph);
                solution = ils1(instance.graph, initial_solution, time_limit, node_switch, edge_switch);
                break;
            case 2:
                initial_solution = random_initial_solution ? random_solution(instance.graph.dim) : greedy_loop(instance.graph);
                solution = ils2(instance.graph, initial_solution, true, time_limit, percentage);
                break;
            case 3:
                initial_solution = random_initial_solution ? random_solution(instance.graph.dim) : greedy_loop(instance.graph);
                solution = ils2(instance.graph, initial_solution, false, time_limit, percentage);
                break;
            case 4:
            // 282.54
                solution = evolve(instance.graph, 282, false);
                break;
            case 5:
                //solution = evolve(instance.graph, 282, true);
                solution = evolve_mt(instance.graph, 282-4, true);
                break;
        }
        u64 elapsed = stm_since(start);
        if (elapsed < min_time) min_time = elapsed;
        if (elapsed > max_time) max_time = elapsed;
        total_time += elapsed;
        auto solution_score = score(instance.graph, solution);
        total_score += solution_score;

        if (solution_score < experiment_result.min) {
            experiment_result.min = solution_score;
            experiment_result.best_solution = solution;
        }
        if (solution_score > experiment_result.max) {
            experiment_result.max = solution_score;
        }
    }

    experiment_result.average = total_score/(f64)n;
    auto average_time = stm_ms(total_time)/(f64)n;
    experiment_result.average_time = average_time;

    std::string result_filepath = "results/"+std::string(instance_name)+"/"+std::string(method_name);
    write_entire_file((result_filepath+"-a.dat").c_str(), experiment_result.best_solution.loop_a);
    write_entire_file((result_filepath+"-b.dat").c_str(), experiment_result.best_solution.loop_b);

    verify_solution(instance, experiment_result.best_solution);
    printf("%s %s: %.2f (%d - %d)\t/\t%.2f (%.2f - %.2f)\n", instance_name, method_name, experiment_result.average, experiment_result.min, experiment_result.max, average_time, stm_ms(min_time), stm_ms(max_time));

    return experiment_result;
}

void run_experiment_for_instance(std::string instance_name) {
    auto instance = parse_file(("data/" + instance_name + ".tsp").c_str());
    write_entire_file(("results/" + instance_name + "/pos.dat").c_str(), instance.positions);
    printf("running experiments for %s\n", instance_name.c_str());
    auto time_limit = run_experiment(instance, "msls", instance_name.c_str(), false, 0, 0).average_time - 20;
    run_experiment(instance, "ils1", instance_name.c_str(), false, time_limit, 1);
    run_experiment(instance, "ils2", instance_name.c_str(), false, time_limit, 2);
    run_experiment(instance, "ils2a", instance_name.c_str(), false, time_limit, 3);
}

#if 0
void globalConvexity(Instance instance, s32 n = 1000) {
    std::vector<Solution> solutions;
    solutions.reserve(n);
    s32 solution_score;
    Solution solution;
    s32 best_score = INT32_MAX;
    s32 best_index;

    printf("solving random+ls\n");
    for(s32 i=0; i<n; i++) {
        if (i%(n/10)==0) {
            printf("%d\n", i);
        }
        solution = random_solution(instance.graph.dim);
        solution = neighbour_search_edge_cache(instance.graph, solution);
        solutions.push_back(solution);
        solution_score = score(instance.graph, solution);
        if (solution_score < best_score) {
            best_score = solution_score;
            best_index = solutions.size()-1;
        }
    }

    s32 dim = instance.graph.dim;
    std::vector<GraphMatrix> solutions_to_matrix;
    solutions_to_matrix.reserve(n);
    s32 solutions_to_list[n][dim] = {0};
    std::vector<s32> scores;
    scores.reserve(n);

    for(s32 k=0; k<n; k++) {
        if (k%(n/10)==0) {
            printf("%d\n", k);
        }
        solution = solutions[k];
        
        GraphMatrix matrix;
        matrix.dim = dim;
        matrix.cells = (s32 *)malloc(sizeof(s32) * dim * dim);
        std::fill_n(matrix.cells, dim*dim, 0);

        for(s32 i=0; i<solution.loop_a.size(); i++) {
            s32 v1 = solution.loop_a[i];
            s32 v2 = solution.loop_a[(i+1) % solution.loop_a.size()];
            matrix.cells[v1*dim+v2]=1;
            matrix.cells[v2*dim+v1]=1;
            solutions_to_list[k][v1]=1;
        }

        for(s32 i=0; i<solution.loop_b.size(); i++) {
            s32 v1 = solution.loop_b[i];
            s32 v2 = solution.loop_b[(i+1) % solution.loop_b.size()];
            matrix.cells[v1*dim+v2]=-1;
            matrix.cells[v2*dim+v1]=-1;
            solutions_to_list[k][v1]=-1;
        }

        scores.push_back(score(instance.graph, solution));
        solutions_to_matrix.push_back(matrix);
    }

    
    std::vector<f64> common_edges_averages;
    std::vector<f64> common_edges_best;
    std::vector<f64> common_vertices_averages;
    std::vector<f64> common_vertices_best;

    common_edges_averages.reserve(n);
    common_edges_best.reserve(n);
    common_vertices_averages.reserve(n);
    common_vertices_best.reserve(n);


    printf("calculating similarity\n");
    for(s32 k=0; k<n; k++){
        if (k%(n/10)==0) {
            printf("%d\n", k);
        }
        solution = solutions[k];

        f64 common_edges_count=0;
        f64 common_vertices_count=0;

        f64 common_edges_count_best=0;
        f64 common_vertices_count_best=0;

        //average
        for (s32 l=0; l<n; l++) {
            for(s32 i=0; i<solution.loop_a.size(); i++) {
                s32 v1 = solution.loop_a[i];
                s32 v2 = solution.loop_a[(i+1) % solution.loop_a.size()];
                s32 edge = solutions_to_matrix[l].get(v1, v2);
                if(edge!=0) {
                // if(edge==1) {
                    if (k!=l) common_edges_count++;
                    if (l==best_index) common_edges_count_best++;
                }
            }

            for(s32 i=0; i<solution.loop_b.size(); i++) {
                s32 v1 = solution.loop_b[i];
                s32 v2 = solution.loop_b[(i+1) % solution.loop_b.size()];
                s32 edge = solutions_to_matrix[l].get(v1, v2);
                if(edge!=0) {
                // if(edge==-1) { 
                    if (k!=l) common_edges_count++;
                    if (l==best_index) common_edges_count_best++;
                }
            }

            s32 aa=0, ab=0, ba=0, bb=0; // aa, ab, ba, bb
            for(s32 i=0; i<dim; i++) {
                s32 vk = solutions_to_list[k][i];
                s32 vl = solutions_to_list[l][i];

                if (vk != 0 && vl != 0) {
                    if (vk == 1 && vl == 1) aa++;
                    if (vk == 1 && vl == -1) ab++;
                    if (vk == -1 && vl == 1) ba++;
                    if (vk == -1 && vl == -1) bb++;
                }
            }
            if (k!=l) common_vertices_count += std::max(aa, ab) + std::max(ba, bb);
            if (l==best_index) common_vertices_count_best = std::max(aa, ab) + std::max(ba, bb);
            
            if (l==best_index) {
                common_edges_best.push_back(common_edges_count_best);
                common_vertices_best.push_back(common_vertices_count_best);
            }
        }

         
        common_edges_averages.push_back(common_edges_count / (float) (n - 1));
        common_vertices_averages.push_back(common_vertices_count / (float) (n-1));
    }   

    write_entire_file("results/scores.dat", scores);
    write_entire_file("results/edge_average.dat", common_edges_averages);
    write_entire_file("results/vertex_average.dat", common_vertices_averages);
    write_entire_file("results/edge_best.dat", common_edges_best);
    write_entire_file("results/vertex_best.dat", common_vertices_best);
}
#endif

int main(int argc, char *argv[]) {
    srand(time(0));
    stm_setup();

    auto instance = parse_file("data/kroA200.tsp");

    //globalConvexity(instance);

    // write_entire_file("results/kroA/pos.dat", instance.positions);

    //run_experiment(instance, "HEA-LS", "kroA", false, 0, 4);
    run_experiment(instance, "HEA+LS", "kroA", false, 0, 5);
    // instance = parse_file("data/kroB200.tsp");
    // write_entire_file("results/kroB/pos.dat", instance.positions);
    // run_experiment(instance, "HEA-LS", "kroB", false, 0, 4);
    // run_experiment(instance, "HEA+LS", "kroB", false, 0, 5);

    return 0;
}
