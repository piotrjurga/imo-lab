#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <iterator>
#include <set>
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
};

ExperimentResult run_experiment(Instance instance, Solver solving_method, const char *method_name, const char *instance_name) {
    Solution min_solution, max_solution;
    ExperimentResult experiment_result = {
        .min = INT32_MAX,
        .max = INT32_MIN
        };
    s32 n = 100;
    s32 total_score = 0;
    for (int i = 0; i < n; i++) {
        auto solution = solving_method(instance.graph);
        auto solution_score = score(instance.graph, solution);
        total_score += solution_score;

        if (solution_score < experiment_result.min) {
            experiment_result.min = solution_score;
            experiment_result.best_solution = solution;
        }
        else if (solution_score > experiment_result.max) {
            experiment_result.max = solution_score;
        }
    }
    experiment_result.average = total_score/(f64)n;

    char result_filepath[128];
    sprintf(result_filepath, "results/%s/%s-a.dat", instance_name, method_name);
    write_entire_file(result_filepath, experiment_result.best_solution.loop_a);
    sprintf(result_filepath, "results/%s/%s-b.dat", instance_name, method_name);
    write_entire_file(result_filepath, experiment_result.best_solution.loop_b);

    verify_solution(instance, experiment_result.best_solution);
    printf("%s %s: %.2f (%d - %d)\n\n", instance_name, method_name, experiment_result.average, experiment_result.min, experiment_result.max);

    return experiment_result;
}

void greedy_experiments() {
    auto kroA100 = parse_file("data/kroA100.tsp");
    write_entire_file("results/kroA100/pos.dat", kroA100.positions);

    auto result1 = run_experiment(kroA100, greedy_simple, "greedy_simple", "kroA100");
    auto result2 = run_experiment(kroA100, greedy_loop, "greedy_loop", "kroA100");
    auto result3 = run_experiment(kroA100, regret_loop, "regret_loop", "kroA100");

    auto kroB100 = parse_file("data/kroB100.tsp");
    write_entire_file("results/kroB100/pos.dat", kroB100.positions);

    auto result4 = run_experiment(kroB100, greedy_simple, "greedy_simple", "kroB100");
    auto result5 = run_experiment(kroB100, greedy_loop, "greedy_loop", "kroB100");
    auto result6 = run_experiment(kroB100, regret_loop, "regret_loop", "kroB100");
}

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

Solution neighbour_search_edge(GraphMatrix graph, Solution init, bool greedy) {
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

ExperimentResult run_experiment_2(Instance instance, bool use_random_solution, LS_Solver local_search_method, const char *method_name, const char *instance_name, bool greedy) {
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
        u64 start = stm_now();
        auto initial_solution = use_random_solution ? random_solution(instance.graph.dim) : greedy_loop(instance.graph);
        auto solution =  local_search_method(instance.graph, initial_solution, greedy);
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
        else if (solution_score > experiment_result.max) {
            experiment_result.max = solution_score;
        }
    }

    experiment_result.average = total_score/(f64)n;
    auto average_time = stm_ms(total_time)/(f64)n;

    char result_filepath[128];
    sprintf(result_filepath, "results/%s/%s-a.dat", instance_name, method_name);
    write_entire_file(result_filepath, experiment_result.best_solution.loop_a);
    sprintf(result_filepath, "results/%s/%s-b.dat", instance_name, method_name);
    write_entire_file(result_filepath, experiment_result.best_solution.loop_b);

    verify_solution(instance, experiment_result.best_solution);
    printf("%s %s: %.2f (%d - %d)\t/\t%.2f (%.2f - %.2f)\n", instance_name, method_name, experiment_result.average, experiment_result.min, experiment_result.max, average_time, stm_ms(min_time), stm_ms(max_time));


    return experiment_result;
}

ExperimentResult run_random_walk_experiment_2(Instance instance, bool use_random_solution, const char *method_name, const char *instance_name, s32 steps) {
    Solution min_solution, max_solution;
    ExperimentResult experiment_result = {
        .min = INT32_MAX,
        .max = INT32_MIN
        };
    s32 n = 100;
    s32 total_score = 0;
    u64 min_time = UINT64_MAX;
    u64 max_time = 0;
    u64 total_time = 0;
    
    for (int i = 0; i < n; i++) {
        u64 start = stm_now();
        auto initial_solution = use_random_solution ? random_solution(instance.graph.dim) : greedy_loop(instance.graph);
        auto solution =  random_walk(instance.graph, initial_solution, steps);
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
        else if (solution_score > experiment_result.max) {
            experiment_result.max = solution_score;
        }
    }
    experiment_result.average = total_score/(f64)n;
    f64 average_time = stm_ms(total_time)/(f64)n;

    char result_filepath[128];
    sprintf(result_filepath, "results/%s/%s-a.dat", instance_name, method_name);
    write_entire_file(result_filepath, experiment_result.best_solution.loop_a);
    sprintf(result_filepath, "results/%s/%s-b.dat", instance_name, method_name);
    write_entire_file(result_filepath, experiment_result.best_solution.loop_b);

    verify_solution(instance, experiment_result.best_solution);
    printf("%s %s: %.2f (%d - %d)\t/\t%.2f (%.2f - %.2f)\n", instance_name, method_name, experiment_result.average, experiment_result.min, experiment_result.max, average_time, stm_ms(min_time), stm_ms(max_time));

    return experiment_result;
}

void run_experiment_2_for_instance(Instance instance, const char *instance_name) {
    auto random_sol = random_solution(instance.graph.dim);
    auto greedy_loop_sol = greedy_loop(instance.graph);

    // losowe rozwiązanie
    run_random_walk_experiment_2(instance, true, "random_loop-random_walk", instance_name, 185000);
    run_experiment_2(instance, true, neighbour_search_node, "random_loop-greedy-neighbour_search_node", instance_name, true);
    run_experiment_2(instance, true, neighbour_search_edge, "random_loop-greedy-neighbour_search_edge", instance_name, true);
    run_experiment_2(instance, true, neighbour_search_node, "random_loop-best-neighbour_search_node", instance_name, false);
    run_experiment_2(instance, true, neighbour_search_edge, "random_loop-best-neighbour_search_edge", instance_name, false);

    // rozbudowa cyklu
    run_random_walk_experiment_2(instance, false, "greedy_loop-random_walk", instance_name, 185000);
    run_experiment_2(instance, false, neighbour_search_node, "greedy_loop-greedy-neighbour_search_node", instance_name, true);
    run_experiment_2(instance, false, neighbour_search_edge, "greedy_loop-greedy-neighbour_search_edge", instance_name, true);
    run_experiment_2(instance, false, neighbour_search_node, "greedy_loop-best-neighbour_search_node", instance_name, false);
    run_experiment_2(instance, false, neighbour_search_edge, "greedy_loop-best-neighbour_search_edge", instance_name, false);
}

void local_search_experiments() {
    auto kroA100 = parse_file("data/kroA100.tsp");
    write_entire_file("results/kroA100/pos.dat", kroA100.positions);
    printf("running experiments for kroA100\n");
    run_experiment_2_for_instance(kroA100, "kroA100");


    auto kroB100 = parse_file("data/kroB100.tsp");
    write_entire_file("results/kroB100/pos.dat", kroB100.positions);
    printf("running experiments for kroB100\n");
    run_experiment_2_for_instance(kroB100, "kroB100");
}

struct Node {
    s32 next;
    s32 prev;
};

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

Solution neighbour_search_edge_cache(GraphMatrix graph, Solution init) {
    // create a sorted list of all moves
    std::multiset<Exchange> moves;

    // add all the edge exchange moves
    std::vector<s32> *loops[2] = { &init.loop_a, &init.loop_b };
    for (s32 loop_i = 0; loop_i < 2; loop_i++) {
        auto& loop = *loops[loop_i];
        s32 loop_count = loop.size();
        Exchange move;
        move.type = (ExchangeType)loop_i;
        for (s32 i = 0; i < loop_count; i++) {
            move.i_from = loop[i];
            move.i_to = loop[(i+1) % loop_count];
            s32 i_cost = graph.get(move.i_from, move.i_to);
            for (s32 j = i+1; j < loop_count; j++) {
                move.j_from = loop[j];
                move.j_to = loop[(j+1) % loop_count];
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
            }
        }
    }

    // add all the node exchange moves
    {
        Exchange move;
        move.type = EXCHANGE_LOOP_AB;

        auto& loop_i = init.loop_a;
        auto& loop_j = init.loop_b;
        s32 loop_i_count = loop_i.size();
        s32 loop_j_count = loop_j.size();
        for (s32 i = 0; i < loop_i_count; i++) {
            move.i_prev = loop_i[(i-1+loop_i_count) % loop_i_count];
            move.i_curr = loop_i[i];
            move.i_next = loop_i[(i+1) % loop_i_count];
            for (s32 j = 0; j < loop_j.size(); j++) {
                move.j_prev = loop_j[(j-1+loop_j_count) % loop_j_count];
                move.j_curr = loop_j[j];
                move.j_next = loop_j[(j+1) % loop_j_count];
                move.delta = node_exchange_delta(graph, move);
                if (move.delta < 0) moves.insert(move);
            }
        }
    }

    auto list = to_linked_list(init);

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

                    assert(list.loop[move.i_curr] != list.loop[move.j_curr]);
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
            s32 iter = 0;
            while (move.j_from != move.i_from) {
                if (iter++ > 1000) exit(0);
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

    return to_visit_list(list);
}

struct NodeLoopPosition {
    s32 loop;
    s32 position;
};

Solution neighbour_search_edge_candidates(GraphMatrix graph, Solution init, s32 max_neighbours) {
    Solution result = init;

    std::vector<s32> *loops[2];
    loops[0] = &result.loop_a;
    loops[1] = &result.loop_b;

    // std::vector<std::vector<std::pair<s32, s32>>> best_neighbours(graph.dim);
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
                        neighbors_visited++;
                        s32 candidate_position = node_loop_positions[candidate_node].position;
                        s32 candidate_cost = candidate.second;

                        s32 candidate_pred = loop[(candidate_position-1 + loop_count) % loop_count];
                        s32 candidate_succ = loop[(candidate_position+1) % loop_count];
                        s32 candidate_pred_cost = graph.get(candidate_pred, candidate_node);
                        s32 candidate_succ_cost = graph.get(candidate_node, candidate_succ);

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

        best_cross_loop = best_node_exchange(graph, result.loop_a, result.loop_b, false); // should this also account for closest neighbours in another loop only?
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

int main() {
    srand(time(0));
    stm_setup();

    auto instance = parse_file("data/kroA200.tsp");
    auto init = greedy_loop(instance.graph);

    printf("greedy loop score = %d\n", score(instance.graph, init));

    s32 start = stm_now();
    auto better_old = neighbour_search_edge(instance.graph, init, false);
    s32 old_time = stm_since(start);
    printf("old implementation score = %d\n", score(instance.graph, better_old));
    printf("old time = %.3f ms\n", stm_ms(old_time));

    start = stm_now();
    auto cache = neighbour_search_edge_cache(instance.graph, init);
    s32 cache_time = stm_since(start);
    printf("cache implementation score = %d\n", score(instance.graph, cache));
    printf("cache time = %.3f ms\n", stm_ms(cache_time));

    start = stm_now();
    auto candidates = neighbour_search_edge_candidates(instance.graph, init, 2);
    s32 candidates_time = stm_since(start);
    printf("candidates implementation score = %d\n", score(instance.graph, candidates));
    printf("candidates time = %.3f ms\n", stm_ms(candidates_time));


    verify_solution(instance, candidates);

    return 0;
}
