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


def format_error_message(error_type, message, pos):
    return f"{error_type}: '{message}' at {pos.abspath}, Line: {pos.line}, Col: {pos.column}"


def check(wdl_filename):
    try:
        WDL.load(wdl_filename)
        return True
    except WDL.Error.SyntaxError as e:
        print(format_error_message("Syntax Error", e, e.pos))
    except WDL.Error.ValidationError as e:
        print(format_error_message("Validation Error", e, e.node.pos))
    except WDL.Error.MultipleValidationErrors as e:
        for ex in e.exceptions:
            print(format_error_message("Validation Error", ex, ex.node.pos))
    except WDL.Error.ImportError as e:
        print(e)
        print(format_error_message("Import Error", e, e.pos))
    finally:
        return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", help="The WDL files to be checked.")
    args = parser.parse_args()

    wdl_filenames = []
    for filename in args.files:
        if filename.endswith(".wdl"):
            wdl_filenames.append(filename)
    if len(wdl_filenames) != len(args.files):
        print(f"Warning! Of {len(args.files)} provided files, {len(wdl_filenames)} are WDL.")

    exit_code = 0
    counter = 0
    for wdl_filename in wdl_filenames:
        if not check(wdl_filename):
            exit_code = 1
            counter += 1
    print(f"{counter} of {len(wdl_filenames)} checked WDL files failed the miniwdl validation.")
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
