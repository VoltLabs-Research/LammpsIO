#pragma once

#include <string>
#include <vector>
#include <lammpsio/frame.hpp>

namespace lammpsio {

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

const std::vector<const FormatReader*>& readers();

const FormatReader* detectReader(const MappedFile& file);

const char* detectFormat(const char* path, std::string& error);

bool scanFrames(const char* path, std::vector<FrameIndexEntry>& frames, std::string& error,
                const char** formatOut = nullptr);

bool readHeader(const char* path, int frame, FrameHeader& header, std::string& error);

bool readFrame(const char* path, const ReadOptions& options, FrameAllocator& allocator,
               ParsedFrame& frame, std::string& error);

class VectorFrameAllocator final : public FrameAllocator {
public:
    explicit VectorFrameAllocator(PositionPrecision precision = PositionPrecision::Float64)
        : precision_(precision) {}

    FrameBuffers allocate(int atomCount, bool withIds) override;
    PositionPrecision positionPrecision() const override { return precision_; }

    const std::vector<double>& positions() const { return positions64_; }
    const std::vector<float>& positionsFloat() const { return positions32_; }
    const std::vector<uint16_t>& types() const { return types_; }
    const std::vector<uint32_t>& ids() const { return ids_; }

private:
    PositionPrecision precision_;
    std::vector<double> positions64_;
    std::vector<float> positions32_;
    std::vector<uint16_t> types_;
    std::vector<uint32_t> ids_;
};

}
