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


def format_error_message(error_type, message, pos, tips=None):
    tips = f"Tips: {tips}" if tips else ""
    return f"{error_type}: '{message}' at {pos.abspath}, Line: {pos.line}, Col: {pos.column}. {tips}"


def check(wdl_filename, imports_dir, import_max_depth):
    try:
        WDL.load(wdl_filename, path=[imports_dir], import_max_depth=import_max_depth)
        return True
    except WDL.Error.SyntaxError as e:
        print(format_error_message("Syntax Error", e, e.pos))
    except WDL.Error.ValidationError as e:
        print(format_error_message("Validation Error", e, e.node.pos))
    except WDL.Error.MultipleValidationErrors as e:
        for ex in e.exceptions:
            print(format_error_message("Validation Error", ex, ex.node.pos))
    except WDL.Error.ImportError as e:
        print(format_error_message(
            "Import Error", e, e.pos,
            f"This could indicate two issues; either the WDLs imported in {wdl_filename} "
            f"cannot be found in the given imports directory (i.e., `{imports_dir}`), "
            f"or miniwdl identifies errors while parsing the imported WDLs in "
            f"maximum {import_max_depth} import depth. You may separately check the "
            f"WDLs imported in {wdl_filename} and make sure they all pass the "
            f"miniwdl validation."))
    finally:
        return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", help="The WDL files to be checked.")
    parser.add_argument("--imports-dir", required=True, help="The directory to search for WDLs to be imported.")
    parser.add_argument("--max-import-depth", default=2, help="Maximum depth of imports to check for validating WDLs.")
    args = parser.parse_args()

    wdl_filenames = args.files
    imports_dir = args.imports_dir

    exit_code = 0
    counter = 0
    for wdl_filename in wdl_filenames:
        # print text in yellow.
        print("\033[93m" + f"\nValidating file {wdl_filename}" + "\033[0m")
        if not check(wdl_filename, imports_dir, args.max_import_depth):
            exit_code = 1
            counter += 1

    print("\n\n")
    if counter > 0:
        # print text in bold red.
        print("\033[1m\033[91m" + f"{counter} of {len(wdl_filenames)} checked WDL files failed the miniwdl validation." + "\033[0m")
    else:
        # print text in bold green.
        print("\033[1m\033[92m" + f"All WDL files are successfully passed the miniwdl validation." + "\033[0m")
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
