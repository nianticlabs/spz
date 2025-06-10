import metabuild.public as mb

mb.set_project_name('spz')


def on_http_request_callback(ctx):
    if "artifactory.corp.adobe.com" not in ctx.request.url:
        raise Exception("Unexpected HTTP request URL: " + ctx.request.url)

mb.meta_policy(
    name = "spz_policy",
    on_http_request = on_http_request_callback
)

# build directory will be .build by default, in repo root.
mb.set_output_directory('$(project_root)/../build_mb')

mb.import_rules('//deps')

mb.cxx_library(
    name="universe_flags",
    exported_deps="//global_flags:global_flags",
)

# group(name = 'spz', deps='//components/spz:spz')
mb.group(
    name = 'main',
    deps = [
        (mb.target.windows, '//spz:spz_bindings_pyd'),
        (~mb.target.windows, '//spz:spz_bindings_so'),
    ]
)
