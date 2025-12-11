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
