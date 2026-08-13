// The addon's entry points. Every exported function follows the same shape: map the
// file, pick a reader, hand the work to it, convert the result. Adding a format means
// adding a translation unit and one line to READERS — nothing here changes.

#include <node_api.h>
#include <string>
#include <vector>
#include "frame.hpp"
#include "napi_bridge.hpp"

namespace lammps_dump_text { extern const FormatReader reader; }
namespace lammps_data { extern const FormatReader reader; }
namespace lammps_dump_binary { extern const FormatReader reader; }
namespace lammps_dump_yaml { extern const FormatReader reader; }
namespace extxyz { extern const FormatReader reader; }

namespace {

/**
 * Sniffing order matters: the text dump check is an exact match on the first line, so
 * it is unambiguous and goes first. The data-file check is structural (an atom count
 * plus box bounds) and therefore looser. The XYZ check is the loosest of all — a line
 * holding nothing but an integer — so it goes last. The binary check reads a whole header
 * successfully or not at all, so it is safe anywhere ahead of the text-structural ones.
 */
const FormatReader* const READERS[] = {
    &lammps_dump_text::reader,
    &lammps_dump_binary::reader,
    &lammps_dump_yaml::reader,
    &lammps_data::reader,
    &extxyz::reader
};

/** Keeps a mapped file open for the duration of a call, however the call exits. */
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

struct CallArgs {
    std::string path;
    ReadOptions options;
};

bool readCallArgs(napi_env env, napi_callback_info info, CallArgs& args) {
    size_t argc = 2;
    napi_value values[2] = { nullptr, nullptr };
    napi_get_cb_info(env, info, &argc, values, nullptr, nullptr);

    if (argc < 1) {
        napi_throw_type_error(env, nullptr, "A file path is required");
        return false;
    }

    args.path = napi_bridge::readString(env, values[0]);
    if (argc >= 2) args.options = napi_bridge::readOptions(env, values[1]);
    return true;
}

/**
 * Identifies the format of an already-mapped file, throwing into JS and returning null
 * when the file could not be opened or nothing recognizes it.
 */
const FormatReader* detectOrThrow(napi_env env, const ScopedMappedFile& mapped, const std::string& path) {
    if (!mapped.valid()) {
        napi_throw_error(env, nullptr, ("Failed to open file: " + path).c_str());
        return nullptr;
    }

    for (const FormatReader* reader : READERS) {
        if (reader->sniff(mapped.get())) return reader;
    }

    napi_throw_error(env, nullptr, ("Unsupported trajectory format: " + path).c_str());
    return nullptr;
}

/**
 * Resolves the requested frame index against a scan, so every read shares one story
 * about what "frame 3" means and what happens when it does not exist.
 */
bool selectFrame(napi_env env, const FormatReader& reader, const MappedFile& file,
                 const std::string& path, int index, FrameIndexEntry& entry) {
    std::vector<FrameIndexEntry> frames;
    std::string error;

    if (!reader.scan(file, frames, error)) {
        napi_throw_error(env, nullptr, error.c_str());
        return false;
    }

    if (index < 0 || index >= (int)frames.size()) {
        const std::string message = "Frame " + std::to_string(index) + " is out of range: " +
            path + " has " + std::to_string(frames.size()) + " frame(s)";
        napi_throw_range_error(env, nullptr, message.c_str());
        return false;
    }

    entry = frames[(size_t)index];
    return true;
}

napi_value DetectFormat(napi_env env, napi_callback_info info) {
    CallArgs args;
    if (!readCallArgs(env, info, args)) return nullptr;

    ScopedMappedFile mapped(args.path.c_str());
    if (!mapped.valid()) {
        // Unreadable is not the same as unrecognized, and a caller asking "what is
        // this?" has to be able to tell the two apart.
        napi_throw_error(env, nullptr, ("Failed to open file: " + args.path).c_str());
        return nullptr;
    }

    napi_value result;
    for (const FormatReader* reader : READERS) {
        if (!reader->sniff(mapped.get())) continue;
        napi_create_string_utf8(env, reader->id, NAPI_AUTO_LENGTH, &result);
        return result;
    }

    napi_get_null(env, &result);
    return result;
}

napi_value ScanFrames(napi_env env, napi_callback_info info) {
    CallArgs args;
    if (!readCallArgs(env, info, args)) return nullptr;

    ScopedMappedFile mapped(args.path.c_str());
    const FormatReader* reader = detectOrThrow(env, mapped, args.path);
    if (!reader) return nullptr;

    std::vector<FrameIndexEntry> frames;
    std::string error;
    if (!reader->scan(mapped.get(), frames, error)) {
        napi_throw_error(env, nullptr, error.c_str());
        return nullptr;
    }

    return napi_bridge::buildScanObject(env, reader->id, frames);
}

napi_value ReadHeader(napi_env env, napi_callback_info info) {
    CallArgs args;
    if (!readCallArgs(env, info, args)) return nullptr;

    ScopedMappedFile mapped(args.path.c_str());
    const FormatReader* reader = detectOrThrow(env, mapped, args.path);
    if (!reader) return nullptr;

    FrameIndexEntry entry;
    if (!selectFrame(env, *reader, mapped.get(), args.path, args.options.frame, entry)) return nullptr;

    FrameHeader header;
    std::string error;
    if (!reader->readHeader(mapped.get(), entry, header, error)) {
        napi_throw_error(env, nullptr, error.c_str());
        return nullptr;
    }

    return napi_bridge::buildHeaderObject(env, reader->id, header);
}

napi_value ReadFrame(napi_env env, napi_callback_info info) {
    CallArgs args;
    if (!readCallArgs(env, info, args)) return nullptr;

    ScopedMappedFile mapped(args.path.c_str());
    const FormatReader* reader = detectOrThrow(env, mapped, args.path);
    if (!reader) return nullptr;

    FrameIndexEntry entry;
    if (!selectFrame(env, *reader, mapped.get(), args.path, args.options.frame, entry)) return nullptr;

    napi_bridge::V8FrameAllocator allocator(env);
    ParsedFrame frame;
    std::string error;
    if (!reader->readFrame(mapped.get(), entry, args.options, allocator, frame, error)) {
        napi_throw_error(env, nullptr, error.c_str());
        return nullptr;
    }

    return napi_bridge::buildFrameObject(env, reader->id, frame, allocator);
}

void exportFunction(napi_env env, napi_value exports, const char* name, napi_callback callback) {
    napi_value fn;
    napi_create_function(env, name, NAPI_AUTO_LENGTH, callback, nullptr, &fn);
    napi_set_named_property(env, exports, name, fn);
}

} // namespace

static napi_value Init(napi_env env, napi_value exports) {
    exportFunction(env, exports, "detectFormat", DetectFormat);
    exportFunction(env, exports, "scanFrames", ScanFrames);
    exportFunction(env, exports, "readHeader", ReadHeader);
    exportFunction(env, exports, "readFrame", ReadFrame);
    return exports;
}

NAPI_MODULE(NODE_GYP_MODULE_NAME, Init)
