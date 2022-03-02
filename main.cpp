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

GraphMatrix parse_file(const char *filepath) {
    FILE *file = fopen(filepath, "r");
    if (!file) return {};

    // skip first 6 lines
    skip_lines(file, 6);

    struct Pos2D { s32 x, y; };
    std::vector<Pos2D> positions;
    while (true) {
        s32 node_idx;
        Pos2D pos;
        s32 read_elems = fscanf(file, " %d %d %d ", &node_idx, &pos.x, &pos.y);
        if (read_elems != 3) break;
        positions.push_back(pos);
    }

    GraphMatrix result;
    result.dim = positions.size();
    result.cells = (s32 *)malloc(sizeof(s32) * result.dim * result.dim);

    for (s32 i = 0; i < result.dim; i++) {
        Pos2D a = positions[i];
        for (s32 j = 0; j < result.dim; j++) {
            Pos2D b = positions[j];
            s32 dx = a.x - b.x;
            s32 dy = a.y - b.y;
            f32 delta = sqrtf(dx*dx + dy*dy);
            s32 cost = (s32)(delta + 0.5f);
            result.cells[result.dim * i + j] = cost;
        }
    }

    return result;
}

int main() {
    auto graph = parse_file("data/a280.tsp");

    return 0;
}
