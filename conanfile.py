import os

from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, cmake_layout
from conan.tools.files import copy


class LammpsIOConan(ConanFile):
    """The trajectory reader core, packaged for the C++ side of the ecosystem.

    Same sources as the npm package `@voltstack/lammps-io`: that one adds an N-API layer on
    top and ships a .node module, this one ships a static library. The formats are read by
    the same code either way, which is the point — a format fix lands once.
    """

    name = "lammpsio"
    version = "2.1.3"
    package_type = "static-library"
    license = "MIT"
    description = "Native readers for LAMMPS and LAMMPS-adjacent trajectory formats"
    settings = "os", "arch", "compiler", "build_type"
    exports_sources = "CMakeLists.txt", "include/*", "src/*", "test/cpp/*", "LICENSE"

    def layout(self):
        cmake_layout(self)

    def generate(self):
        toolchain = CMakeToolchain(self)
        toolchain.cache_variables["LAMMPSIO_BUILD_TESTS"] = False
        toolchain.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()
        copy(self, "LICENSE", self.source_folder, os.path.join(self.package_folder, "licenses"))

    def package_info(self):
        self.cpp_info.libs = ["lammpsio"]
        self.cpp_info.includedirs = ["include"]
        if self.settings.os in ("Linux", "FreeBSD"):
            self.cpp_info.system_libs = ["pthread"]
