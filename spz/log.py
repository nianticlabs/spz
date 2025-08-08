import boa_toolkit.utils.log as base

root_logger = base.root_logger
logger: base.Logger = root_logger.sub_logger("SPZ")
Channel = base.Channel
ScopedLog = base.ScopedLog
