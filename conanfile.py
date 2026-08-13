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
    # No third-party requirements: fast_float is vendored and the rest is the standard
    # library. Keeping it that way is deliberate — this sits at the bottom of the plugin
    # dependency chain, where every requirement is one every plugin inherits.
    exports_sources = "CMakeLists.txt", "include/*", "src/*", "test/cpp/*", "LICENSE"

    def layout(self):
        cmake_layout(self)

    def generate(self):
        toolchain = CMakeToolchain(self)
        # The smoke test is for developing the library, not for consumers building it.
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
            # The text reader parses in parallel with std::thread.
            self.cpp_info.system_libs = ["pthread"]
