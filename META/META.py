import metabuild.public as mb

mb.set_project_name("spz")


def on_http_request_callback(ctx):
    if "artifactory.corp.adobe.com" not in ctx.request.url:
        raise Exception("Unexpected HTTP request URL: " + ctx.request.url)

mb.meta_policy(
    name = "spz_policy",
    on_http_request = on_http_request_callback
)

mb.import_rules("//deps")

mb.cxx_library(
    name="universe_flags",
    exported_deps="//global_flags:global_flags",
)

mb.group(
    name = "main",
    deps = [
        "//spz:spz",
        (mb.target.option("//options:build_tools"), "//cli_tools:cli_tools"),
    ]
)
