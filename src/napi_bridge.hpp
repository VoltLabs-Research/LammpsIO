#pragma once

// The only translation unit boundary that touches N-API. Reading options off a JS
// object and turning a ParsedFrame into one, used to be copy-pasted per parser; every
// format now shares this one copy.

#include <node_api.h>
#include <string>
#include <vector>
#include <lammpsio/frame.hpp>

using namespace lammpsio;

namespace napi_bridge {

inline std::string readString(napi_env env, napi_value value) {
    size_t length = 0;
    napi_get_value_string_utf8(env, value, nullptr, 0, &length);
    std::string out(length, '\0');
    napi_get_value_string_utf8(env, value, &out[0], length + 1, &length);
    return out;
}

inline bool isObject(napi_env env, napi_value value) {
    napi_valuetype type;
    napi_typeof(env, value, &type);
    return type == napi_object;
}

inline ReadOptions readOptions(napi_env env, napi_value value) {
    ReadOptions options;
    if (!isObject(env, value)) return options;

    napi_value field;
    if (napi_get_named_property(env, value, "includeIds", &field) == napi_ok) {
        napi_valuetype type;
        napi_typeof(env, field, &type);
        if (type == napi_boolean) napi_get_value_bool(env, field, &options.includeIds);
    }

    if (napi_get_named_property(env, value, "maxThreads", &field) == napi_ok) {
        napi_valuetype type;
        napi_typeof(env, field, &type);
        if (type == napi_number) napi_get_value_int32(env, field, &options.maxThreads);
    }

    if (napi_get_named_property(env, value, "frame", &field) == napi_ok) {
        napi_valuetype type;
        napi_typeof(env, field, &type);
        if (type == napi_number) napi_get_value_int32(env, field, &options.frame);
    }

    if (napi_get_named_property(env, value, "properties", &field) == napi_ok) {
        bool isArray = false;
        napi_is_array(env, field, &isArray);
        if (isArray) {
            uint32_t length = 0;
            napi_get_array_length(env, field, &length);
            for (uint32_t i = 0; i < length; i++) {
                napi_value element;
                napi_get_element(env, field, i, &element);
                options.properties.push_back(readString(env, element));
            }
        }
    }

    return options;
}

inline void setDouble(napi_env env, napi_value target, const char* name, double value) {
    napi_value wrapped;
    napi_create_double(env, value, &wrapped);
    napi_set_named_property(env, target, name, wrapped);
}

inline void setInt(napi_env env, napi_value target, const char* name, int32_t value) {
    napi_value wrapped;
    napi_create_int32(env, value, &wrapped);
    napi_set_named_property(env, target, name, wrapped);
}

inline void setString(napi_env env, napi_value target, const char* name, const char* value) {
    napi_value wrapped;
    napi_create_string_utf8(env, value, NAPI_AUTO_LENGTH, &wrapped);
    napi_set_named_property(env, target, name, wrapped);
}

inline void setBool(napi_env env, napi_value target, const char* name, bool value) {
    napi_value wrapped;
    napi_get_boolean(env, value, &wrapped);
    napi_set_named_property(env, target, name, wrapped);
}

inline napi_value toStringArray(napi_env env, const std::vector<std::string>& values) {
    napi_value array;
    napi_create_array_with_length(env, values.size(), &array);
    for (size_t i = 0; i < values.size(); i++) {
        napi_value element;
        napi_create_string_utf8(env, values[i].c_str(), values[i].size(), &element);
        napi_set_element(env, array, i, element);
    }
    return array;
}

inline napi_value toVec3(napi_env env, double x, double y, double z) {
    napi_value array, element;
    napi_create_array_with_length(env, 3, &array);
    napi_create_double(env, x, &element); napi_set_element(env, array, 0, element);
    napi_create_double(env, y, &element); napi_set_element(env, array, 1, element);
    napi_create_double(env, z, &element); napi_set_element(env, array, 2, element);
    return array;
}

/**
 * Allocates the bulk per-atom buffers directly inside V8-visible ArrayBuffers, so a
 * reader writes atom data once and it reaches JavaScript with no intermediate copy.
 */
class V8FrameAllocator final : public FrameAllocator {
public:
    explicit V8FrameAllocator(napi_env env) : env_(env) {}

    FrameBuffers allocate(int atomCount, bool withIds) override {
        const size_t count = atomCount > 0 ? (size_t)atomCount : 0;
        FrameBuffers buffers;
        void* raw = nullptr;
        napi_value arrayBuffer;

        // float32 for positions: this is what the viewer's GLB and the parquet frames
        // store, so anything wider would be narrowed a moment later anyway.
        napi_create_arraybuffer(env_, count * 3 * sizeof(float), &raw, &arrayBuffer);
        napi_create_typedarray(env_, napi_float32_array, count * 3, arrayBuffer, 0, &positions_);
        buffers.positions32 = (float*)raw;

        napi_create_arraybuffer(env_, count * sizeof(uint16_t), &raw, &arrayBuffer);
        napi_create_typedarray(env_, napi_uint16_array, count, arrayBuffer, 0, &types_);
        buffers.types = (uint16_t*)raw;

        if (withIds) {
            napi_create_arraybuffer(env_, count * sizeof(uint32_t), &raw, &arrayBuffer);
            napi_create_typedarray(env_, napi_uint32_array, count, arrayBuffer, 0, &ids_);
            buffers.ids = (uint32_t*)raw;
            hasIds_ = true;
        }

        return buffers;
    }

    PositionPrecision positionPrecision() const override { return PositionPrecision::Float32; }

    napi_value positions() const { return positions_; }
    napi_value types() const { return types_; }
    napi_value ids() const { return ids_; }
    bool hasIds() const { return hasIds_; }

private:
    napi_env env_;
    napi_value positions_ = nullptr;
    napi_value types_ = nullptr;
    napi_value ids_ = nullptr;
    bool hasIds_ = false;
};

/**
 * Materializes one staged column into the typed array its dtype calls for: Int32Array
 * for i32, Float32Array for f32. The staging doubles are exact for both casts.
 */
inline napi_value toColumnArray(napi_env env, const ExtraColumn& column) {
    const size_t count = column.values.size();
    void* raw = nullptr;
    napi_value arrayBuffer, typedArray;

    if (column.dtype == ColumnDtype::Int32) {
        napi_create_arraybuffer(env, count * sizeof(int32_t), &raw, &arrayBuffer);
        napi_create_typedarray(env, napi_int32_array, count, arrayBuffer, 0, &typedArray);
        int32_t* RESTRICT out = (int32_t*)raw;
        for (size_t i = 0; i < count; i++) out[i] = (int32_t)column.values[i];
    } else {
        napi_create_arraybuffer(env, count * sizeof(float), &raw, &arrayBuffer);
        napi_create_typedarray(env, napi_float32_array, count, arrayBuffer, 0, &typedArray);
        float* RESTRICT out = (float*)raw;
        for (size_t i = 0; i < count; i++) out[i] = (float)column.values[i];
    }

    return typedArray;
}

inline napi_value buildHeaderObject(napi_env env, const char* formatId, const FrameHeader& header) {
    napi_value result, box, pbc, cell, element;
    napi_create_object(env, &result);
    napi_create_object(env, &box);

    setDouble(env, box, "xlo", header.box.xlo);
    setDouble(env, box, "xhi", header.box.xhi);
    setDouble(env, box, "ylo", header.box.ylo);
    setDouble(env, box, "yhi", header.box.yhi);
    setDouble(env, box, "zlo", header.box.zlo);
    setDouble(env, box, "zhi", header.box.zhi);
    setDouble(env, box, "xy", header.box.xy);
    setDouble(env, box, "xz", header.box.xz);
    setDouble(env, box, "yz", header.box.yz);

    napi_create_array_with_length(env, 3, &pbc);
    for (size_t axis = 0; axis < 3; axis++) {
        napi_get_boolean(env, header.periodic[axis], &element);
        napi_set_element(env, pbc, axis, element);
    }

    // The three cell vectors as rows, so cellVectors[0] is a. This is the faithful
    // representation: boxBounds cannot express a lattice that is not upper triangular.
    napi_create_array_with_length(env, 3, &cell);
    for (size_t vector = 0; vector < 3; vector++) {
        napi_set_element(env, cell, vector,
            toVec3(env, header.cell[vector][0], header.cell[vector][1], header.cell[vector][2]));
    }

    // Ordered [name, value] pairs, so a round trip can put the sections back where they were.
    napi_value extras;
    napi_create_array_with_length(env, header.extraSections.size(), &extras);
    for (size_t i = 0; i < header.extraSections.size(); i++) {
        napi_value pair, item;
        napi_create_array_with_length(env, 2, &pair);
        napi_create_string_utf8(env, header.extraSections[i].first.c_str(),
                                header.extraSections[i].first.size(), &item);
        napi_set_element(env, pair, 0, item);
        napi_create_string_utf8(env, header.extraSections[i].second.c_str(),
                                header.extraSections[i].second.size(), &item);
        napi_set_element(env, pair, 1, item);
        napi_set_element(env, extras, i, pair);
    }
    napi_set_named_property(env, result, "extraSections", extras);
    setBool(env, result, "positionsWereScaled", header.positionsWereScaled);

    setString(env, result, "format", formatId);
    setInt(env, result, "timestep", header.timestep);
    setInt(env, result, "natoms", header.atomCount);
    napi_set_named_property(env, result, "boxBounds", box);
    napi_set_named_property(env, result, "cellVectors", cell);
    napi_set_named_property(env, result, "cellOrigin",
        toVec3(env, header.origin[0], header.origin[1], header.origin[2]));
    napi_set_named_property(env, result, "pbc", pbc);
    napi_set_named_property(env, result, "headers", toStringArray(env, header.headers));

    return result;
}

inline napi_value buildFrameObject(napi_env env, const char* formatId,
                                   const ParsedFrame& frame, const V8FrameAllocator& allocator) {
    napi_value result;
    napi_create_object(env, &result);

    napi_set_named_property(env, result, "positions", allocator.positions());
    napi_set_named_property(env, result, "types", allocator.types());
    if (allocator.hasIds() && frame.hasIds) {
        napi_set_named_property(env, result, "ids", allocator.ids());
    }

    if (!frame.extras.empty()) {
        napi_value properties, dtypes;
        napi_create_object(env, &properties);
        napi_create_object(env, &dtypes);

        for (const auto& column : frame.extras) {
            napi_set_named_property(env, properties, column.name.c_str(), toColumnArray(env, column));
            setString(env, dtypes, column.name.c_str(), columnDtypeString(column.dtype));
        }

        napi_set_named_property(env, result, "properties", properties);
        napi_set_named_property(env, result, "propertyDtypes", dtypes);
    }

    napi_set_named_property(env, result, "metadata", buildHeaderObject(env, formatId, frame.header));

    napi_set_named_property(env, result, "min",
        toVec3(env, frame.bbox.minX, frame.bbox.minY, frame.bbox.minZ));
    napi_set_named_property(env, result, "max",
        toVec3(env, frame.bbox.maxX, frame.bbox.maxY, frame.bbox.maxZ));

    // Emitted independently: a data file's Masses section gives both, but an
    // extended-XYZ file gives element symbols and no masses at all.
    napi_value element;

    if (!frame.massesByType.empty()) {
        napi_value masses;
        napi_create_array_with_length(env, frame.massesByType.size(), &masses);
        for (size_t i = 0; i < frame.massesByType.size(); i++) {
            napi_create_double(env, frame.massesByType[i], &element);
            napi_set_element(env, masses, i, element);
        }
        napi_set_named_property(env, result, "massesByType", masses);
    }

    if (!frame.elementHintsByType.empty()) {
        napi_value hints;
        napi_create_array_with_length(env, frame.elementHintsByType.size(), &hints);
        for (size_t i = 0; i < frame.elementHintsByType.size(); i++) {
            if (frame.elementHintsByType[i].empty()) {
                napi_get_null(env, &element);
            } else {
                napi_create_string_utf8(env, frame.elementHintsByType[i].c_str(),
                                        frame.elementHintsByType[i].size(), &element);
            }
            napi_set_element(env, hints, i, element);
        }
        napi_set_named_property(env, result, "elementHintsByType", hints);
    }

    return result;
}

inline napi_value buildScanObject(napi_env env, const char* formatId,
                                  const std::vector<FrameIndexEntry>& frames) {
    napi_value result, list, entry;
    napi_create_object(env, &result);
    napi_create_array_with_length(env, frames.size(), &list);

    for (size_t i = 0; i < frames.size(); i++) {
        napi_create_object(env, &entry);
        setInt(env, entry, "index", frames[i].index);
        // Offsets and lengths can exceed int32 on large trajectories, so they go out as
        // doubles — exact up to 2^53, far past any file size that matters here.
        setDouble(env, entry, "byteOffset", (double)frames[i].byteOffset);
        setDouble(env, entry, "byteLength", (double)frames[i].byteLength);
        setInt(env, entry, "timestep", frames[i].timestep);
        setInt(env, entry, "natoms", frames[i].atomCount);
        napi_set_element(env, list, i, entry);
    }

    setString(env, result, "format", formatId);
    napi_set_named_property(env, result, "frames", list);
    return result;
}

} // namespace napi_bridge
