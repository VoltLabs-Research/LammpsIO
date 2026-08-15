#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <lammpsio/reader_registry.hpp>

namespace {

struct Expectation {
    const char* file;
    const char* format;
    int frameCount;
    int atomCount;
};

const Expectation EXPECTED[] = {
    { "dump-ortho-3frames.dump",      "lammps-dump",        3, 4 },
    { "dump-triclinic-scaled.dump",   "lammps-dump",        1, 3 },
    { "dump-mixed-boundaries.dump",   "lammps-dump",        1, 2 },
    { "dump-binary-2frames.bin",      "lammps-dump-binary", 2, 4 },
    { "dump-yaml-2frames.yaml",       "lammps-dump-yaml",   2, 4 },
    { "data-atomic.data",             "lammps-data",        1, 4 },
    { "data-full.data",               "lammps-data",        1, 3 },
    { "data-triclinic.data",          "lammps-data",        1, 2 },
    { "xyz-extended-2frames.xyz",     "extxyz",             2, 3 },
    { "xyz-plain.xyz",                "extxyz",             1, 4 }
};

int failures = 0;

void check(bool condition, const std::string& what) {
    if (condition) return;
    std::printf("  FAIL: %s\n", what.c_str());
    failures++;
}

}

int main(int argc, char** argv) {
    const std::string root = argc > 1 ? argv[1] : "test/fixtures";

    for (const Expectation& expected : EXPECTED) {
        const std::string path = root + "/" + expected.file;
        std::string error;

        const char* format = lammpsio::detectFormat(path.c_str(), error);
        if (!error.empty()) {
            std::printf("%-30s ERROR %s\n", expected.file, error.c_str());
            failures++;
            continue;
        }

        std::vector<lammpsio::FrameIndexEntry> frames;
        if (!lammpsio::scanFrames(path.c_str(), frames, error)) {
            std::printf("%-30s ERROR %s\n", expected.file, error.c_str());
            failures++;
            continue;
        }

        lammpsio::VectorFrameAllocator allocator(lammpsio::PositionPrecision::Float64);
        lammpsio::ParsedFrame frame;
        lammpsio::ReadOptions options;
        options.includeIds = true;
        options.properties = { "*" };
        options.maxThreads = 1;

        if (!lammpsio::readFrame(path.c_str(), options, allocator, frame, error)) {
            std::printf("%-30s ERROR %s\n", expected.file, error.c_str());
            failures++;
            continue;
        }

        std::printf("%-30s %-20s %d frame(s)  %d atoms  pos[0]=(%g, %g, %g)  extras=%zu\n",
                    expected.file, format ? format : "(none)", (int)frames.size(),
                    frame.header.atomCount,
                    allocator.positions()[0], allocator.positions()[1], allocator.positions()[2],
                    frame.extras.size());

        check(format && std::strcmp(format, expected.format) == 0,
              std::string(expected.file) + ": format");
        check((int)frames.size() == expected.frameCount,
              std::string(expected.file) + ": frame count");
        check(frame.header.atomCount == expected.atomCount,
              std::string(expected.file) + ": atom count");
        check(allocator.positions().size() == (size_t)expected.atomCount * 3,
              std::string(expected.file) + ": position buffer size");
        check(allocator.types().size() == (size_t)expected.atomCount,
              std::string(expected.file) + ": type buffer size");
    }

    {
        std::string error;
        const std::string path = root + "/../cpp/smoke.cpp";
        const char* format = lammpsio::detectFormat(path.c_str(), error);
        check(error.empty() && format == nullptr, "a C++ source file is not a trajectory");
    }

    std::printf(failures == 0 ? "\nC++ smoke: OK\n" : "\nC++ smoke: %d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
