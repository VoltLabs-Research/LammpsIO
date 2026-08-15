#include <lammpsio/reader_registry.hpp>

#include <string>
#include <vector>

namespace lammpsio {

namespace lammps_dump_text { extern const FormatReader reader; }
namespace lammps_dump_binary { extern const FormatReader reader; }
namespace lammps_dump_yaml { extern const FormatReader reader; }
namespace lammps_data { extern const FormatReader reader; }
namespace extxyz { extern const FormatReader reader; }

namespace {

const std::vector<const FormatReader*>& readerTable() {
    static const std::vector<const FormatReader*> table = {
        &lammps_dump_text::reader,
        &lammps_dump_binary::reader,
        &lammps_dump_yaml::reader,
        &lammps_data::reader,
        &extxyz::reader
    };
    return table;
}

std::string openFailure(const char* path) {
    return std::string("Failed to open file: ") + path;
}

std::string unsupported(const char* path) {
    return std::string("Unsupported trajectory format: ") + path;
}

bool locateFrame(const char* path, int frameIndex, const ScopedMappedFile& mapped,
                 const FormatReader*& reader, FrameIndexEntry& entry, std::string& error) {
    if (!mapped.valid()) {
        error = openFailure(path);
        return false;
    }

    reader = detectReader(mapped.get());
    if (!reader) {
        error = unsupported(path);
        return false;
    }

    std::vector<FrameIndexEntry> frames;
    if (!reader->scan(mapped.get(), frames, error)) return false;

    if (frameIndex < 0 || frameIndex >= (int)frames.size()) {
        error = "Frame " + std::to_string(frameIndex) + " is out of range: " + path +
                " has " + std::to_string(frames.size()) + " frame(s)";
        return false;
    }

    entry = frames[(size_t)frameIndex];
    return true;
}

}

const std::vector<const FormatReader*>& readers() {
    return readerTable();
}

const FormatReader* detectReader(const MappedFile& file) {
    for (const FormatReader* reader : readerTable()) {
        if (reader->sniff(file)) return reader;
    }
    return nullptr;
}

const char* detectFormat(const char* path, std::string& error) {
    ScopedMappedFile mapped(path);
    if (!mapped.valid()) {
        error = openFailure(path);
        return nullptr;
    }

    const FormatReader* reader = detectReader(mapped.get());
    return reader ? reader->id : nullptr;
}

bool scanFrames(const char* path, std::vector<FrameIndexEntry>& frames, std::string& error,
                const char** formatOut) {
    ScopedMappedFile mapped(path);
    if (!mapped.valid()) {
        error = openFailure(path);
        return false;
    }

    const FormatReader* reader = detectReader(mapped.get());
    if (!reader) {
        error = unsupported(path);
        return false;
    }

    if (formatOut) *formatOut = reader->id;
    return reader->scan(mapped.get(), frames, error);
}

bool readHeader(const char* path, int frame, FrameHeader& header, std::string& error) {
    ScopedMappedFile mapped(path);
    const FormatReader* reader = nullptr;
    FrameIndexEntry entry;

    if (!locateFrame(path, frame, mapped, reader, entry, error)) return false;

    if (!reader->readHeader(mapped.get(), entry, header, error)) return false;

    header.format = reader->id;
    return true;
}

bool readFrame(const char* path, const ReadOptions& options, FrameAllocator& allocator,
               ParsedFrame& frame, std::string& error) {
    ScopedMappedFile mapped(path);
    const FormatReader* reader = nullptr;
    FrameIndexEntry entry;

    if (!locateFrame(path, options.frame, mapped, reader, entry, error)) return false;

    if (!reader->readFrame(mapped.get(), entry, options, allocator, frame, error)) return false;

    frame.header.format = reader->id;
    return true;
}

FrameBuffers VectorFrameAllocator::allocate(int atomCount, bool withIds) {
    const size_t count = atomCount > 0 ? (size_t)atomCount : 0;
    FrameBuffers buffers;

    if (precision_ == PositionPrecision::Float64) {
        positions64_.assign(count * 3, 0.0);
        buffers.positions64 = positions64_.data();
    } else {
        positions32_.assign(count * 3, 0.0f);
        buffers.positions32 = positions32_.data();
    }

    types_.assign(count, 0);
    buffers.types = types_.data();

    if (withIds) {
        ids_.assign(count, 0);
        buffers.ids = ids_.data();
    }

    return buffers;
}

}
