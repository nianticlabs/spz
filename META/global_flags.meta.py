import metabuild.public as mb

mb.cxx_library(
    name = "global_flags",
    exported_compiler_flags = [
        (mb.target.linux, [
            "-fPIC", # for dynamic linkage
        ]),
        (~mb.target.windows, [
            "-Wno-shorten-64-to-32"
        ]),
    ],
    xcode_flags={
        "PRODUCT_BUNDLE_IDENTIFIER": "com.nianticlabs.spz",
    },
    exported_preprocessor_macros=[
        "_USE_MATH_DEFINES",
    ],
)
