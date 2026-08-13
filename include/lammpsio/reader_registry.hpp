#pragma once

// The public C++ surface: everything a consumer needs to read a trajectory, with no
// N-API anywhere. The Node addon is one client of this; a C++ project linking the static
// library is another, and both get the same format coverage from the same code.
//
//   const char* format = lammpsio::detectFormat("run.dump");     // nullptr if unknown
//
//   std::vector<lammpsio::FrameIndexEntry> frames;
//   std::string error;
//   lammpsio::scanFrames("run.dump", frames, error);
//
//   lammpsio::VectorFrameAllocator allocator;                    // float64 positions
//   lammpsio::ParsedFrame frame;
//   lammpsio::ReadOptions options;
//   options.frame = 2;
//   options.includeIds = true;
//   options.properties = { "*" };
//   lammpsio::readFrame("run.dump", options, allocator, frame, error);
//   // allocator.positions() now holds 3 * frame.header.atomCount doubles.

#include <string>
#include <vector>
#include <lammpsio/frame.hpp>

namespace lammpsio {

/** Keeps a memory-mapped file open for as long as it is needed, and closes it once. */
class ScopedMappedFile {
public:
    explicit ScopedMappedFile(const char* path) : file_(mapFile(path)) {}
    ~ScopedMappedFile() { if (file_.valid) unmapFile(file_); }

    ScopedMappedFile(const ScopedMappedFile&) = delete;
    ScopedMappedFile& operator=(const ScopedMappedFile&) = delete;

    bool valid() const { return file_.valid; }
    const MappedFile& get() const { return file_; }

private:
    MappedFile file_;
};

/**
 * The readers, in sniffing order. Exposed so a caller can enumerate what is supported
 * rather than hard-coding a list that drifts.
 */
const std::vector<const FormatReader*>& readers();

/** The reader that claims this file, or nullptr when none does. */
const FormatReader* detectReader(const MappedFile& file);

/**
 * Format id of a file, or nullptr when nothing recognizes it. Sets `error` and returns
 * nullptr when the file itself could not be opened — unreadable and unrecognized are
 * different answers, and a caller usually wants to report them differently.
 */
const char* detectFormat(const char* path, std::string& error);

/** `formatOut`, when given, receives the id of the reader that claimed the file. */
bool scanFrames(const char* path, std::vector<FrameIndexEntry>& frames, std::string& error,
                const char** formatOut = nullptr);

bool readHeader(const char* path, int frame, FrameHeader& header, std::string& error);

bool readFrame(const char* path, const ReadOptions& options, FrameAllocator& allocator,
               ParsedFrame& frame, std::string& error);

/**
 * A FrameAllocator that owns its buffers, for consumers with no special memory
 * requirements. Positions are float64, which is what a geometry pipeline generally wants;
 * pass Float32 to halve the footprint when the values feed a renderer.
 */
class VectorFrameAllocator final : public FrameAllocator {
public:
    explicit VectorFrameAllocator(PositionPrecision precision = PositionPrecision::Float64)
        : precision_(precision) {}

    FrameBuffers allocate(int atomCount, bool withIds) override;
    PositionPrecision positionPrecision() const override { return precision_; }

    /** 3 * atomCount entries, xyz interleaved. Empty when the other precision was used. */
    const std::vector<double>& positions() const { return positions64_; }
    const std::vector<float>& positionsFloat() const { return positions32_; }
    const std::vector<uint16_t>& types() const { return types_; }
    /** Empty unless ids were requested and the format carries them. */
    const std::vector<uint32_t>& ids() const { return ids_; }

private:
    PositionPrecision precision_;
    std::vector<double> positions64_;
    std::vector<float> positions32_;
    std::vector<uint16_t> types_;
    std::vector<uint32_t> ids_;
};

} // namespace lammpsio
