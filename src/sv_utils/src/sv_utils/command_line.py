#!/usr/bin/env python

import sys
import os
from sv_utils import common
import contextlib
import io
from typing import Text, List, TextIO, Optional, Iterator, FrozenSet


src_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


def _kebab_to_snake(kebab_str: str) -> str:
    """
    Allow passing commands in kebab case by converting to snake case (to match file names)
    """
    return kebab_str.replace('-', '_')


def _snake_to_kebab(snake_str: str) -> str:
    """
    Convert snake-case files to kebab-case for display
    """
    return snake_str.replace('_', '-')


def _find_sub_modules(package: str) -> Iterator[Text]:
    module = common.dynamic_import(package)
    module_path = module.__path__
    # noinspection PyProtectedMember
    module_folder = module_path[0] if isinstance(module_path, List) else module_path._path[0]
    command_module = os.path.basename(__file__)  # but exclude this module if it happens to be there
    for sub_module_name in os.listdir(module_folder):
        if sub_module_name != command_module:
            yield sub_module_name.rsplit(".py", 1)[0]


def _get_help_summary(package: str, sub_module: str) -> Optional[str]:
    try:
        submodule_arg_parser = common.dynamic_import(f"{package}.{sub_module}.__parse_arguments")
    except ModuleNotFoundError:
        # no __parse_arguments, so not a command-line module. Return None help string
        return None
    string_buffer = io.StringIO()
    with contextlib.redirect_stdout(string_buffer):
        try:
            # noinspection PyUnresolvedReferences
            submodule_arg_parser([sub_module, "--help"])
        except SystemExit:
            pass
        # get the help summary as the part of the parsed help preceded and followed by empty lines
        summary_lines = string_buffer.getvalue().split("\n\n", 2)
        # extra if statements in case the "help" isn't working for a badly-behaved package
        return summary_lines[1] if len(summary_lines) >= 2 else summary_lines[0] if summary_lines \
            else "(No help available)"


def _command_help_func(
        package: str,
        sub_modules: FrozenSet[str],
        file_descriptor: TextIO = sys.stdout,
        num_indent: int = 2,
        max_width: int = 78
):
    print(f"{package} [command] [args...]", file=file_descriptor)
    print("Valid commands are:", file=file_descriptor)
    indent1 = " " * num_indent
    indent2 = " " * (num_indent * 2)
    for command_name in sub_modules:
        help_str = _get_help_summary(package, command_name)
        if help_str is None:
            continue
        print(f"{indent1}{_snake_to_kebab(command_name)}:", file=file_descriptor)
        help_str = help_str.strip()
        if len(indent2 + help_str) <= max_width:
            print(f"{indent2}{help_str}", file=file_descriptor)
        else:
            help_words = help_str.split()
            help_line = indent2 + help_words[0]
            for next_word in help_words[1:]:
                if len(help_line + " " + next_word) <= max_width:
                    help_line += " " + next_word
                else:
                    print(str(help_line), file=file_descriptor)
                    help_line = indent2 + next_word
            print(str(help_line), file=file_descriptor)
    print("", file=file_descriptor)


def _bad_command_func(
        package: str,
        sub_modules: FrozenSet[str],
        argv: List[Text],
        file_descriptor: TextIO = sys.stderr
):
    if len(argv) >= 1:
        print("Bad command: %s" % argv[0], file=file_descriptor)
    else:
        print("No command specified.", file=file_descriptor)

    _command_help_func(package=package, sub_modules=sub_modules, file_descriptor=file_descriptor)
    sys.exit(1)


def main(argv: Optional[List[Text]] = None):
    """
    Dispatch arguments to appropriate module function, with call to
    command_line.py removed, so that the dispatched command will behave the same
    as if it had been called directly.
    Args:
        argv: input arguments to command. If called from command-line (usual use
              case), this will simply be sys.argv
    """
    if argv is None:
        argv = sys.argv
    package = os.path.basename(_kebab_to_snake(argv[0]))
    sub_modules = frozenset(_find_sub_modules(package))
    if len(argv) < 2:
        _bad_command_func(package=package, sub_modules=sub_modules, argv=argv[1:])
    else:
        sub_module_command = argv[1]
        if sub_module_command in {"help", "-help", "--help"}:
            _command_help_func(package=package, sub_modules=sub_modules)
            return
        sub_module_command = _kebab_to_snake(sub_module_command)
        if sub_module_command in sub_modules:
            try:
                call_func = common.dynamic_import(f"{package}.{sub_module_command}.main")
            except ModuleNotFoundError:
                raise ValueError(
                    f"{sub_module_command} does not have a 'main' function, and thus cannot be invoked as a "
                    "command-line program."
                )
            call_func(argv[1:])
        else:
            # bad command
            _bad_command_func(package=package, sub_modules=sub_modules, argv=argv[1:])


if __name__ == "__main__":
    main()
