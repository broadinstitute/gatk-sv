"""
This Script checks the provided WDL files with miniwdl for
syntax, validation, and import errors. See the following
page on the checks. Note that this script does not
check for "warning" or "suggestions"; you may use
`$ miniwdl check WDL_FILENAME` to get warnings and suggestions.

https://miniwdl.readthedocs.io/en/latest/WDL.html#module-WDL.Error
"""

import argparse
import sys
import WDL


def format_error_message(error_type, message, pos=None, tips=None):
    tips = f"Tips: {tips}" if tips else ""
    msg = f"{error_type}: '{message}'"
    if pos:
        msg += f"at {pos.abspath}, Line: {pos.line}, Col: {pos.column}."
    msg += tips
    return msg


def check(wdl_filename, imports_dir, import_max_depth):
    print(f"Validating {wdl_filename} ...", end="")
    error_msg = ""
    try:
        WDL.load(wdl_filename, path=[imports_dir], import_max_depth=import_max_depth)
        print("\033[92m Pass! \033[0m")  # print in Green
        return True
    except WDL.Error.SyntaxError as e:
        error_msg = format_error_message("Syntax Error", e, e.pos)
    except WDL.Error.ValidationError as e:
        error_msg = format_error_message("Validation Error", e, e.node.pos if e.node else None)
    except WDL.Error.MultipleValidationErrors as e:
        for ex in e.exceptions:
            error_msg += format_error_message("Validation Error", ex, ex.node.pos if ex.node else None) + "\n"
    except WDL.Error.ImportError as e:
        error_msg = format_error_message(
            "Import Error", e, e.pos,
            f"This could indicate two issues; either the WDLs imported in {wdl_filename} "
            f"cannot be found in the given imports directory (i.e., `{imports_dir}`), "
            f"or miniwdl identifies errors while parsing the imported WDLs in "
            f"maximum {import_max_depth} import depth. You may separately check the "
            f"WDLs imported in {wdl_filename} and make sure they all pass the "
            f"miniwdl validation.")
    print("\033[91m Fail! Error(s):\033[0m")  # print in red.
    print("\033[93m" + f"\n{error_msg}\n" + "\033[0m")
    return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", help="The WDL files to be checked.")
    parser.add_argument("--imports-dir", required=True, help="The directory to search for WDLs to be imported.")
    parser.add_argument(
        "--max-import-depth",
        default=15,
        help="Set the maximum allowed depth of imports in the WDLs. "
             "For a given WDL, miniwdl first asserts the imports, and "
             "it will error-out if it reaches the given max depth "
             "and the WDL has imports. For instance, with max depth set "
             "to 2, it will error-out validating "
             "Foo.wdl -> (import Bar.wdl -> (import Baz.wdl)).")
    args = parser.parse_args()

    wdl_filenames = args.files
    imports_dir = args.imports_dir
    max_import_depth = int(args.max_import_depth)

    exit_code = 0
    counter = 0
    for wdl_filename in wdl_filenames:
        if not check(wdl_filename, imports_dir, max_import_depth):
            exit_code = 1
            counter += 1

    print("\n\n")
    if counter > 0:
        # print text in bold red.
        print("\033[1m\033[91m" + f"{counter} of {len(wdl_filenames)} checked WDLs failed the miniwdl validation." + "\033[0m")
    else:
        # print text in bold green.
        print("\033[1m\033[92m" + "All the WDLs successfully passed the miniwdl validation." + "\033[0m")
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
