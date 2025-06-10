# This file is adapted from https://git.corp.adobe.com/3di/python-scaffold

import boa_toolkit.utils.log as base

root_logger = base.root_logger
logger: base.Logger = root_logger.sub_logger("SPZ")
Channel = base.Channel
ScopedLog = base.ScopedLog
