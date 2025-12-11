import metabuild.public as mb

mb.cxx_library(
    name="cli_tools_common",
    project_subpath="cli_tools",
    public_include_directories=[
        "$(project_root)/../src",
    ],
    exported_deps=[
        "//spz:spz",
    ],
)

mb.cxx_binary(
    name="ply_to_spz",
    project_subpath="cli_tools",
    deps=[
        ":cli_tools_common",
    ],
    xcode_product_type="tool",  # don't bundle the app
    srcs=[
        "$(project_root)/../cli_tools/src/ply_to_spz.cpp",
    ],
)

mb.cxx_binary(
    name="spz_to_ply",
    project_subpath="cli_tools",
    deps=[
        ":cli_tools_common",
    ],
    xcode_product_type="tool",  # don't bundle the app
    srcs=[
        "$(project_root)/../cli_tools/src/spz_to_ply.cpp",
    ],
)

mb.cxx_binary(
    name="spz_info",
    project_subpath="cli_tools",
    deps=[
        ":cli_tools_common",
    ],
    xcode_product_type="tool",  # don't bundle the app
    srcs=[
        "$(project_root)/../cli_tools/src/spz_info.cpp",
    ],
)

mb.group(
    name="cli_tools",
    deps=[
        ":ply_to_spz",
        ":spz_to_ply",
        ":spz_info",
    ],
)