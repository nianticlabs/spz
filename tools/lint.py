"""
Run all the linting tools simultaneously on the code.
"""
# This file is adapted from https://git.corp.adobe.com/3di/python-scaffold

import argparse
import asyncio
import multiprocessing
import sys
from pathlib import Path

from boa_toolkit.parallel import aioworker
from boa_toolkit.utils.log import root_logger

logger = root_logger.sub_logger("spz_tools")

ROOT = Path(__file__).resolve().parents[1]

FORMAT_TARGETS = [
    str(ROOT / "tools"),
    str(ROOT / "spz"),
    str(ROOT / "tests"),
    str(ROOT / "setup.py"),
]

LINT_TARGETS = FORMAT_TARGETS


def parse_args():
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--check", action="store_true")
    common.add_argument("--check-diff", action="store_true")
    common.add_argument("--verbose", action="store_true")

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="action", required=True)

    subparsers.add_parser("black", parents=[common])
    subparsers.add_parser("pylint", parents=[common])
    subparsers.add_parser("isort", parents=[common])
    subparsers.add_parser("all", parents=[common])

    return parser.parse_args()


def start_execute_in_bg(args):
    # Start running a command on the event loop and return
    # the asyncio task.
    return asyncio.create_task(aioworker.execute_in_bg(args, check=False))


async def main(args):
    if args.verbose:
        logger.console_handler.setLevel("DEBUG")

    num_proces = multiprocessing.cpu_count()
    failed = False

    tasks = {}
    if args.action in ["black", "all"]:
        if args.check:
            tasks["black"] = start_execute_in_bg(["black", "--check", *FORMAT_TARGETS])
        else:
            tasks["black"] = start_execute_in_bg(["black", *FORMAT_TARGETS])
    if args.action in ["isort", "all"]:
        if args.check:
            tasks["isort"] = start_execute_in_bg(["isort", "--check-only", *FORMAT_TARGETS])
        else:
            tasks["isort"] = start_execute_in_bg(["isort", *FORMAT_TARGETS])
    if args.action in ["pylint", "all"]:
        tasks["pylint"] = start_execute_in_bg(["pylint", "-j", str(min(16, num_proces)), *LINT_TARGETS])

    if "black" in tasks:
        code, out, err = await tasks["black"]
        if code != 0:
            logger.error(
                "[black] Some of your files are not formatted correctly. "
                "Either run `python tools/lint.py black` or open the faulty files in vscode and run `format document`. "
                "Black output:\n```\n%s```",
                out + err,
            )
            failed = True
        else:
            logger.info("[black] passed.\n%s", out + err)
    if "isort" in tasks:
        code, out, err = await tasks["isort"]
        if code != 0:
            logger.error(
                "[isort] The order of includes in some of your files is not formatted correctly. "
                "Either run `python tools/lint.py isort` or open the faulty files in vscode and run `sort includes`. "
                "isort output:\n```\n%s```",
                out + err,
            )
            failed = True
        else:
            logger.info("[isort] passed.\n%s", out + err)
    if "pylint" in tasks:
        code, out, err = await tasks["pylint"]
        if code != 0:
            logger.error("[pylint] detected some issues with your code```\n%s```", out + err)
            failed = True
        else:
            logger.info("[pylint] passed.\n%s", out + err)

    if failed:
        raise Exception("Linting failed.")


def _main():
    """
    Wrapper for main (sync, exception catching, etc.)
    """
    success = True
    try:
        args = parse_args()
        aioworker.resolve_coroutine_now(main(args))
    except KeyboardInterrupt:
        logger.warning("Keyboard interrupt ...")
        success = False
    except Exception as exc:  # pylint: disable=broad-except
        logger.error("Exception happened %r", exc)
        logger.verbose(exc, exc_info=True)
        success = False

    if not success:
        sys.exit(-1)


if __name__ == "__main__":
    _main()
