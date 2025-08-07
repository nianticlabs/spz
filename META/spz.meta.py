import metabuild.public as mb

mb.cxx_library(
    name = "spz",
    srcs = [
        mb.cpp_glob("$(project_root)/../src/cc/**"),
    ],
    public_include_directories = [
        "$(project_root)/../src/cc"
    ],
    preferred_linkage="static",
    exported_deps = [
        "zlib//:zlib",
    ],
)

mb.set_root_directory("$(project_root)/../src/pybindings")

mb.cxx_library(
    name = "spz_bindings_pyd",
    product_name = "spz_bindings",
    srcs = [
        mb.cpp_glob("**"),
    ],
    exported_deps = [
        "spz:spz",
        "pybind11//:module",
    ],
    preferred_linkage="shared",
    filter=mb.target.windows,
    product_extension=".pyd",
)

mb.cxx_library(
    name = "spz_bindings_so",
    product_name = "spz_bindings",
    srcs = [
        mb.cpp_glob("**"),
    ],
    exported_deps = [
        "spz:spz",
        "pybind11//:module",
    ],
    preferred_linkage="shared",
    filter=~mb.target.windows,
    product_extension=".so",
    xcode_flags={
        "PRODUCT_BUNDLE_IDENTIFIER": "com.adobe.spz.spz_bindings",
    },
)
