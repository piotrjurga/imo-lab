#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <vector>

typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int8_t   s8;
typedef int16_t  s16;
typedef int32_t  s32;
typedef int64_t  s64;

struct Node {
    std::vector<Node *> next;
    std::vector<s32> cost;
};

struct Graph {
    std::vector<Node> nodes;
};

Graph parse_file(char *filepath) {
    FILE *file = fopen(filepath, "r");
    if (!file) return {};

    // skip first 6 lines
    for (s32 i = 0; i < 6; i++) {
        char line[80];
        getline(line, sizeof(line), file);
    }

    struct Pos2D { s32 x, y; };
    std::vector<Pos2D> positions;
    while (true) {
        s32 node_idx;
        Pos2D pos;
        fscanf(" %d %d %d ", &node_idx, &pos.x, &pos.y);
        positions.push_back(pos);
    }

    Graph result;
    result.nodes.resize(positions.size());
    for (s32 i = 0; i < result.nodes.size(); i++) {
        auto& node = result.nodes[i];
        node.next.resize(result.nodes.size());
        node.cost.resize(result.nodes.size());
        auto pos_a = positions[i];
        for (s32 j = 0; j < result.nodes.size(); j++) {
            node.next[j] = &result.nodes[j];
            auto pos_b = positions[j];
            s32 dx = pos_a.x - pos_b.x;
            s32 dy = pos_a.y - pos_b.y;
            f32 delta = sqrtf(dx*dx + dy*dy);
            node.cost[j] = (s32)(delta + 0.5f);
        }
    }
}

int main() {
    auto graph = parse_file("data/a280.tsp");
}
