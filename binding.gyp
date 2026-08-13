{
  "targets": [
    {
      "target_name": "lammps_io",
      "sources": [
        "src/registry.cpp",
        "src/formats/lammps_dump_text.cpp",
        "src/formats/lammps_data.cpp",
        "src/formats/lammps_dump_binary.cpp",
        "src/formats/lammps_dump_yaml.cpp",
        "src/formats/extxyz.cpp"
      ],
      "include_dirs": [
        "src"
      ],
      "cflags!": [
        "-fno-exceptions"
      ],
      "cflags_cc!": [
        "-fno-exceptions"
      ],
      "cflags": [
        "-O3"
      ],
      "cflags_cc": [
        "-O3",
        "-std=c++17"
      ],
      "defines": [
        "NAPI_CPP_EXCEPTIONS"
      ]
    }
  ]
}
