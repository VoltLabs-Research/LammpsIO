#pragma once

#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#if defined(_WIN32)
    #ifndef NOMINMAX
        #define NOMINMAX
    #endif
    #ifndef WIN32_LEAN_AND_MEAN
        #define WIN32_LEAN_AND_MEAN
    #endif
    #include <windows.h>
#else
    #include <fcntl.h>
    #include <sys/mman.h>
    #include <sys/stat.h>
    #include <unistd.h>
#endif

#include <lammpsio/external/fast_float.h>

#if defined(_MSC_VER)
    #define UNLIKELY(x) (!!(x))
    #define ALWAYS_INLINE __forceinline
    #define HOT
    #define RESTRICT
#else
    #define UNLIKELY(x) __builtin_expect(!!(x), 0)
    #define ALWAYS_INLINE __attribute__((always_inline)) inline
    #define HOT __attribute__((hot))
    #define RESTRICT __restrict__
#endif

namespace lammpsio {

ALWAYS_INLINE const void* memmemFallback(
    const void* haystack,
    size_t haystackLength,
    const void* needle,
    size_t needleLength
) {
    if (needleLength == 0) return haystack;
    if (haystackLength < needleLength) return nullptr;

    const char* const begin = static_cast<const char*>(haystack);
    const char* const lastStart = begin + (haystackLength - needleLength);
    const char firstByte = *static_cast<const char*>(needle);

    for (const char* candidate = begin; candidate <= lastStart; ++candidate) {
        candidate = static_cast<const char*>(
            memchr(candidate, firstByte, static_cast<size_t>(lastStart - candidate) + 1)
        );
        if (candidate == nullptr) return nullptr;
        if (memcmp(candidate, needle, needleLength) == 0) return candidate;
    }

    return nullptr;
}

#if defined(_MSC_VER)
ALWAYS_INLINE const void* memmem(
    const void* haystack,
    size_t haystackLength,
    const void* needle,
    size_t needleLength
) {
    return memmemFallback(haystack, haystackLength, needle, needleLength);
}
#endif

struct MappedFile {
    const char* data;
    size_t size;
#if defined(_WIN32)
    HANDLE file;
    HANDLE mapping;
#else
    int fd;
#endif
    bool valid;
};

#if defined(_WIN32)

ALWAYS_INLINE MappedFile mapFile(const char* filepath) {
    MappedFile f = {nullptr, 0, INVALID_HANDLE_VALUE, nullptr, false};

    f.file = CreateFileA(
        filepath,
        GENERIC_READ,
        FILE_SHARE_READ,
        nullptr,
        OPEN_EXISTING,
        FILE_ATTRIBUTE_NORMAL | FILE_FLAG_SEQUENTIAL_SCAN,
        nullptr
    );
    if (UNLIKELY(f.file == INVALID_HANDLE_VALUE)) return f;

    LARGE_INTEGER fileSize;
    if (UNLIKELY(!GetFileSizeEx(f.file, &fileSize))) {
        CloseHandle(f.file);
        f.file = INVALID_HANDLE_VALUE;
        return f;
    }

    if (UNLIKELY(fileSize.QuadPart == 0)) {
        CloseHandle(f.file);
        f.file = INVALID_HANDLE_VALUE;
        return f;
    }

    f.mapping = CreateFileMappingA(f.file, nullptr, PAGE_READONLY, 0, 0, nullptr);
    if (UNLIKELY(f.mapping == nullptr)) {
        CloseHandle(f.file);
        f.file = INVALID_HANDLE_VALUE;
        return f;
    }

    const void* view = MapViewOfFile(f.mapping, FILE_MAP_READ, 0, 0, 0);
    if (UNLIKELY(view == nullptr)) {
        CloseHandle(f.mapping);
        CloseHandle(f.file);
        f.mapping = nullptr;
        f.file = INVALID_HANDLE_VALUE;
        return f;
    }

    f.data = static_cast<const char*>(view);
    f.size = static_cast<size_t>(fileSize.QuadPart);
    f.valid = true;
    return f;
}

ALWAYS_INLINE void unmapFile(MappedFile& f) {
    if (f.data) UnmapViewOfFile(f.data);
    if (f.mapping) CloseHandle(f.mapping);
    if (f.file != INVALID_HANDLE_VALUE) CloseHandle(f.file);
    f.data = nullptr;
    f.mapping = nullptr;
    f.file = INVALID_HANDLE_VALUE;
    f.valid = false;
}

#else

ALWAYS_INLINE MappedFile mapFile(const char* filepath) {
    MappedFile f = {nullptr, 0, -1, false};

    f.fd = open(filepath, O_RDONLY);
    if (UNLIKELY(f.fd < 0)) return f;

    struct stat sb;
    if (UNLIKELY(fstat(f.fd, &sb) < 0)) {
        close(f.fd);
        return f;
    }

    f.size = sb.st_size;
    if (UNLIKELY(f.size == 0)) {
        close(f.fd);
        return f;
    }

    f.data = (const char*)mmap(nullptr, f.size, PROT_READ, MAP_PRIVATE | MAP_NORESERVE, f.fd, 0);
    if (UNLIKELY(f.data == MAP_FAILED)) {
        close(f.fd);
        f.data = nullptr;
        return f;
    }

    madvise((void*)f.data, f.size, MADV_SEQUENTIAL | MADV_WILLNEED);

    f.valid = true;
    return f;
}

ALWAYS_INLINE void unmapFile(MappedFile& f) {
    if (f.data) munmap((void*)f.data, f.size);
    if (f.fd >= 0) close(f.fd);
    f.valid = false;
}

#endif


HOT ALWAYS_INLINE double fastAtof(const char* RESTRICT p, const char* RESTRICT end) {
    if (UNLIKELY(p >= end)) return 0.0;
    double result = 0.0;
    auto answer = fast_float::from_chars(p, end, result);
    if (answer.ec != std::errc()) return 0.0;
    return result;
}

HOT ALWAYS_INLINE int fastAtoi(const char* RESTRICT p, const char* RESTRICT end) {
    if (UNLIKELY(p >= end)) return 0;

    int sign = 1;
    if (*p == '-') {
        sign = -1;
        p++;
    } else if (*p == '+') {
        p++;
    }

    int result = 0;
    while (p < end) {
        unsigned int d = (unsigned int)(*p - '0');
        if (d > 9) break;
        result = result * 10 + d;
        p++;
    }

    return sign * result;
}

enum class ColumnDtype { Int32, Float32 };

ALWAYS_INLINE const char* columnDtypeString(ColumnDtype dtype) {
    return dtype == ColumnDtype::Int32 ? "i32" : "f32";
}

HOT ALWAYS_INLINE bool isIntegerToken(const char* RESTRICT p, const char* RESTRICT end) {
    if (UNLIKELY(p >= end)) return false;
    if (*p == '+' || *p == '-') p++;
    if (UNLIKELY(p >= end)) return false;
    while (p < end) {
        unsigned int d = (unsigned int)(*p - '0');
        if (d > 9) return false;
        p++;
    }
    return true;
}

HOT ALWAYS_INLINE const char* skipWhitespace(const char* RESTRICT p, const char* RESTRICT end) {
    while (p < end && *p <= ' ') p++;
    return p;
}

HOT ALWAYS_INLINE const char* findTokenEnd(const char* RESTRICT p, const char* RESTRICT end) {
    while (p < end && *p > ' ') p++;
    return p;
}

ALWAYS_INLINE const char* jumpToNextLine(const char* RESTRICT p, const char* RESTRICT end) {
    if (UNLIKELY(p >= end)) return end;
    const char* nl = (const char*)memchr(p, '\n', (size_t)(end - p));
    return nl ? nl + 1 : end;
}

ALWAYS_INLINE const char* findLineEnd(const char* RESTRICT p, const char* RESTRICT end) {
    if (UNLIKELY(p >= end)) return end;
    const char* nl = (const char*)memchr(p, '\n', (size_t)(end - p));
    return nl ? nl : end;
}

struct alignas(16) BoundingBox {
    float minX, minY, minZ, _pad1;
    float maxX, maxY, maxZ, _pad2;
    
    void init() {
        minX = minY = minZ = 1e30f;
        maxX = maxY = maxZ = -1e30f;
    }
    
    HOT ALWAYS_INLINE void update(float x, float y, float z) {
        if (x < minX) minX = x;
        if (x > maxX) maxX = x;
        if (y < minY) minY = y;
        if (y > maxY) maxY = y;
        if (z < minZ) minZ = z;
        if (z > maxZ) maxZ = z;
    }
    
    void merge(const BoundingBox& other) {
        if (other.minX < minX) minX = other.minX;
        if (other.maxX > maxX) maxX = other.maxX;
        if (other.minY < minY) minY = other.minY;
        if (other.maxY > maxY) maxY = other.maxY;
        if (other.minZ < minZ) minZ = other.minZ;
        if (other.maxZ > maxZ) maxZ = other.maxZ;
    }
};

struct SimulationBox {
    double xlo, xhi;
    double ylo, yhi;
    double zlo, zhi;
    double xy = 0.0, xz = 0.0, yz = 0.0;
};

ALWAYS_INLINE bool isPositionColumn(const char* head, size_t len) {
    if (len == 1) return true;
    if (len == 2) return head[1] == 'u' || head[1] == 's';
    if (len == 3) return head[1] == 's' && head[2] == 'u';
    return false;
}

struct ColumnMapping {
    int idxId = -1;
    int idxType = -1;
    int idxX = -1;
    int idxY = -1;
    int idxZ = -1;
    int maxIdx = 0;
    bool scaledCoords = false;
    
    std::vector<int> extraPropIndices;
    std::vector<std::string> extraPropNames;
    
    void computeMaxIdx() {
        maxIdx = idxType;
        if (idxX > maxIdx) maxIdx = idxX;
        if (idxY > maxIdx) maxIdx = idxY;
        if (idxZ > maxIdx) maxIdx = idxZ;
        if (idxId > maxIdx) maxIdx = idxId;
        for (int idx : extraPropIndices) {
            if (idx > maxIdx) maxIdx = idx;
        }
    }
};

ALWAYS_INLINE void detectDataColumnStyle(int colCount, ColumnMapping& cols) {
    cols.idxId = 0;
    
    if (colCount >= 7) {
        cols.idxType = 2;
        cols.idxX = 4;
        cols.idxY = 5;
        cols.idxZ = 6;
    } else if (colCount == 6) {
        cols.idxType = 1;
        cols.idxX = 3;
        cols.idxY = 4;
        cols.idxZ = 5;
    } else {
        cols.idxType = 1;
        cols.idxX = 2;
        cols.idxY = 3;
        cols.idxZ = 4;
    }
    cols.computeMaxIdx();
}

}
