import metabuild.public as mb

mb.meta_policy(
    name = "spz_policy",
    on_http_request = lambda : None,
    on_git_request = lambda : None,
)
