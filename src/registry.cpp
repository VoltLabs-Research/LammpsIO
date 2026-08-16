#include <node_api.h>
#include <string>
#include <vector>
#include "napi_bridge.hpp"
#include <lammpsio/reader_registry.hpp>

struct CallArgs {
    std::string path;
    lammpsio::ReadOptions options;
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

napi_value throwReaderError(napi_env env, const std::string& error) {
    if (error.rfind("Frame ", 0) == 0 && error.find("out of range") != std::string::npos) {
        napi_throw_range_error(env, nullptr, error.c_str());
    } else {
        napi_throw_error(env, nullptr, error.c_str());
    }
    return nullptr;
}

napi_value DetectFormat(napi_env env, napi_callback_info info) {
    CallArgs args;
    if (!readCallArgs(env, info, args)) return nullptr;

    std::string error;
    const char* format = lammpsio::detectFormat(args.path.c_str(), error);
    if (!error.empty()) return throwReaderError(env, error);

    napi_value result;
    if (format) {
        napi_create_string_utf8(env, format, NAPI_AUTO_LENGTH, &result);
    } else {
        napi_get_null(env, &result);
    }
    return result;
}

napi_value ScanFrames(napi_env env, napi_callback_info info) {
    CallArgs args;
    if (!readCallArgs(env, info, args)) return nullptr;

    std::vector<lammpsio::FrameIndexEntry> frames;
    const char* format = nullptr;
    std::string error;

    if (!lammpsio::scanFrames(args.path.c_str(), frames, error, &format)) {
        return throwReaderError(env, error);
    }

    return napi_bridge::buildScanObject(env, format, frames);
}

napi_value ReadHeader(napi_env env, napi_callback_info info) {
    CallArgs args;
    if (!readCallArgs(env, info, args)) return nullptr;

    lammpsio::FrameHeader header;
    std::string error;

    if (!lammpsio::readHeader(args.path.c_str(), args.options.frame, header, error)) {
        return throwReaderError(env, error);
    }

    return napi_bridge::buildHeaderObject(env, header.format, header);
}

napi_value ReadFrame(napi_env env, napi_callback_info info) {
    CallArgs args;
    if (!readCallArgs(env, info, args)) return nullptr;

    napi_bridge::V8FrameAllocator allocator(env);
    lammpsio::ParsedFrame frame;
    std::string error;

    if (!lammpsio::readFrame(args.path.c_str(), args.options, allocator, frame, error)) {
        return throwReaderError(env, error);
    }

    return napi_bridge::buildFrameObject(env, frame.header.format, frame, allocator);
}

void exportFunction(napi_env env, napi_value exports, const char* name, napi_callback callback) {
    napi_value fn;
    napi_create_function(env, name, NAPI_AUTO_LENGTH, callback, nullptr, &fn);
    napi_set_named_property(env, exports, name, fn);
}

static napi_value Init(napi_env env, napi_value exports) {
    exportFunction(env, exports, "detectFormat", DetectFormat);
    exportFunction(env, exports, "scanFrames", ScanFrames);
    exportFunction(env, exports, "readHeader", ReadHeader);
    exportFunction(env, exports, "readFrame", ReadFrame);
    return exports;
}

NAPI_MODULE(NODE_GYP_MODULE_NAME, Init)
