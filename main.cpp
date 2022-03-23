#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <set>
#include <string>

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
                    s32 cost = graph.cells[node*graph.dim + node_a] + graph.cells[node*graph.dim + node_b]; 
                    cost -= graph.cells[node_a*graph.dim + node_b];
                    
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
                    s32 cost = graph.cells[node*graph.dim + node_a] + graph.cells[node*graph.dim + node_b]; 
                    cost -= graph.cells[node_a*graph.dim + node_b];

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

void verifySolution(const char* instance_name, const char* method_name, Instance instance, Solution s) {
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
    ExperimentResult experimentResult = {
        .min = INT32_MAX,
        .max = INT32_MIN
        };
    s32 n = 100;
    s32 total_score = 0;
    for (int i = 0; i < n; i++) {
        auto solution = solving_method(instance.graph);
        auto solution_score = score(instance.graph, solution);
        total_score += solution_score;

        if (solution_score < experimentResult.min) {
            experimentResult.min = solution_score;
            experimentResult.best_solution = solution;
        }
        else if (solution_score > experimentResult.max) {
            experimentResult.max = solution_score;
        }
    }
    experimentResult.average = total_score/(f64)n;

    char *result_filepath;
    sprintf(result_filepath, "results/%s/%s-a.dat", instance_name, method_name);
    write_entire_file(result_filepath, experimentResult.best_solution.loop_a);
    sprintf(result_filepath, "results/%s/%s-b.dat", instance_name, method_name);
    write_entire_file(result_filepath, experimentResult.best_solution.loop_b);

    verifySolution(instance_name, method_name, instance, experimentResult.best_solution);
    printf("%s %s: %.2f (%d - %d)\n\n", instance_name, method_name, experimentResult.average, experimentResult.min, experimentResult.max);

    return experimentResult;
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

Solution neighbour_search_node(GraphMatrix graph, Solution init, bool greedy) {
    std::vector<s32> *best_loop = NULL;
    s32 node_i = 0, node_j = 0;
    s32 best_delta = 0;

    Solution sol = init;

    std::vector<s32> *loops[2];
    loops[0] = &sol.loop_a;
    loops[1] = &sol.loop_b;

    while (true) {
        best_delta = 0;

        for (s32 loop_i = 0; loop_i < 2; loop_i++) {
            auto& loop = *loops[loop_i];
            s32 loop_count = loop.size();
            for (s32 i = 0; i < loop_count; i++) {
                s32 i_prev = loop[(i-1+loop_count) % loop_count];
                s32 i_curr = loop[i];
                s32 i_next = loop[(i+1) % loop_count];
                s32 i_cost = graph.cells[i_curr*graph.dim + i_prev] +
                             graph.cells[i_curr*graph.dim + i_next];
                for (s32 j = i+2; j < loop_count; j++) {
                    s32 j_prev = loop[j-1];
                    s32 j_curr = loop[j];
                    s32 j_next = loop[(j+1) % loop_count];
                    s32 j_cost = graph.cells[j_curr*graph.dim + j_prev] +
                                 graph.cells[j_curr*graph.dim + j_next];
                    s32 new_i_cost = graph.cells[j_curr*graph.dim + i_prev] +
                                     graph.cells[j_curr*graph.dim + i_next];
                    s32 new_j_cost = graph.cells[i_curr*graph.dim + j_prev] +
                                     graph.cells[i_curr*graph.dim + j_next];
                    s32 delta = new_i_cost + new_j_cost - i_cost - j_cost;
                    //if (j_curr == i_next || i_curr == j_next) continue;
                    if (j_curr == i_next) {
                        s32 old_cost = graph.cells[i_curr*graph.dim + i_prev] +
                                       graph.cells[j_curr*graph.dim + j_next];
                        s32 new_cost = graph.cells[j_curr*graph.dim + i_prev] +
                                       graph.cells[i_curr*graph.dim + j_next];
                        delta = new_cost - old_cost;
                    } else if (i_curr == j_next) {
                        s32 old_cost = graph.cells[j_curr*graph.dim + j_prev] +
                                       graph.cells[i_curr*graph.dim + i_next];
                        s32 new_cost = graph.cells[i_curr*graph.dim + j_prev] +
                                       graph.cells[j_curr*graph.dim + i_next];
                        delta = new_cost - old_cost;
                    }
                    if (delta < best_delta) {
                        best_delta = delta;
                        node_i = i;
                        node_j = j;
                        best_loop = &loop;
                        if (greedy) {
                            goto apply_found_move;
                        }
                    }
                }
            }
        }

apply_found_move:
        if (best_delta < 0) {
            auto& loop = *best_loop;
#if 0
            printf("delta: %d\n", best_delta);
            printf("switch nodes: (%d, %d)\n", node_i, node_j);
            if (node_i == 29 && node_j == 30) {
                s32 i_prev = loop[(node_i-1+loop.size()) % loop.size()];
                s32 i_curr = loop[node_i];
                s32 i_next = loop[(node_i+1) % loop.size()];

                s32 j_prev = loop[node_j-1];
                s32 j_curr = loop[node_j];
                s32 j_next = loop[(node_j+1) % loop.size()];

                printf("i-1 ->  i  = %d\n", graph.cells[i_prev*graph.dim + i_curr]);
                printf("i-1 ->  j  = %d\n", graph.cells[i_prev*graph.dim + j_curr]);
                printf("j-1 ->  i  = %d\n", graph.cells[j_prev*graph.dim + i_curr]);
                printf("j-1 ->  j  = %d\n", graph.cells[j_prev*graph.dim + j_curr]);
                printf(" i  -> i+1 = %d\n", graph.cells[i_curr*graph.dim + i_next]);
                printf(" i  -> j+1 = %d\n", graph.cells[i_curr*graph.dim + j_next]);
                printf(" j  -> i+1 = %d\n", graph.cells[j_curr*graph.dim + i_next]);
                printf(" j  -> j+1 = %d\n", graph.cells[j_curr*graph.dim + j_next]);
            }
#endif
            s32 t = loop[node_i];
            loop[node_i] = loop[node_j];
            loop[node_j] = t;
        } else {
            break;
        }
    }
    return sol;
}

Solution neighbour_search_edge(GraphMatrix graph, Solution init, bool greedy) {
    std::vector<s32> *best_loop = NULL;
    s32 edge_i = 0, edge_j = 0;
    s32 best_delta = 0;

    Solution sol = init;

    std::vector<s32> *loops[2];
    loops[0] = &sol.loop_a;
    loops[1] = &sol.loop_b;

    while (true) {
        best_delta = 0;

        for (s32 loop_i = 0; loop_i < 2; loop_i++) {
            auto& loop = *loops[loop_i];
            s32 loop_count = loop.size();
            for (s32 i = 0; i < loop_count; i++) {
                s32 i_from = loop[i];
                s32 i_to = loop[(i+1) % loop_count];
                s32 i_cost = graph.cells[i_from*graph.dim + i_to];
                for (s32 j = i+1; j < loop_count; j++) {
                    s32 j_from = loop[j];
                    s32 j_to = loop[(j+1) % loop_count];
                    s32 j_cost = graph.cells[j_from*graph.dim + j_to];

                    s32 new_i_cost = graph.cells[j_from*graph.dim + i_from];
                    s32 new_j_cost = graph.cells[j_to*graph.dim + i_to];
                    s32 delta = new_i_cost + new_j_cost - i_cost - j_cost;
                    if (delta < best_delta) {
                        best_delta = delta;
                        edge_i = i;
                        edge_j = j;
                        best_loop = &loop;
                        if (greedy) {
                            goto apply_found_move;
                        }
                    }
                }
            }
        }

apply_found_move:
        if (best_delta < 0) {
            auto& loop = *best_loop;
            s32 between = (edge_i+1 + edge_j) / 2;
            for (s32 i = edge_i+1; i <= between; i++) {
                s32 j = edge_j + edge_i+1 - i;
                s32 t = loop[i];
                loop[i] = loop[j];
                loop[j] = t;
            }
        } else {
            break;
        }
    }
    return sol;
}


int main() {
    srand(time(0));

    auto kroA100 = parse_file("data/kroA100.tsp");
    write_entire_file("results/kroA100/pos.dat", kroA100.positions);

    auto s1 = greedy_loop(kroA100.graph);
    auto s2 = neighbour_search_edge(kroA100.graph, s1, false);
    auto s3 = neighbour_search_node(kroA100.graph, s1, false);
    printf("greedy score: %d\n", score(kroA100.graph, s1));
    printf("edge search score: %d\n", score(kroA100.graph, s2));
    printf("node search score: %d\n", score(kroA100.graph, s3));

    return 0;
}
