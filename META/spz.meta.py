import metabuild.public as mb

mb.cxx_library(
    name = "spz",
    srcs = [
        mb.cpp_glob("$(project_root)/../src/cc/**"),
        (mb.target.option("//options:build_extensions"), mb.cpp_glob("$(project_root)/../extensions/cc/**")),
    ],
    public_include_directories = [
        "$(project_root)/../src/cc",
        (mb.target.option("//options:build_extensions"), "$(project_root)/../extensions/cc"),
    ],
    preferred_linkage="static",
    exported_deps = [
        "zlib//:zlib",
    ],
)
