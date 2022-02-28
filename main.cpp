#include <stdio.h>
#include <stdint.h>
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

    for (s32 i = 0; i < 6; i++) {
        char line[80];
        getline(line, sizeof(line), file);
    }

    while (true) {
        s32 node_idx, x, y;
        fscanf();
    }
}

int main() {
    auto graph = parse_file("data/a280.tsp");
}
