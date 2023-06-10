#
# from pathlib import Path

#

#
# from socialgene.utils.printing import (
#     print_test_name,
#     print_message,
#     print_success,
#     print_error,
# )


# def check_hashed_file(neo4j_outdir, subpath, expected_files):
#     print_test_name("Check tigrfam files:")
#     inspect_dir = Path(neo4j_outdir, *subpath)
#     print_message("Looking for files in:")
#     print_message(str(inspect_dir))
#     filenames = [i.name for i in Path.iterdir(inspect_dir)]
#     # separate for-loops so the different cases print out next to each other
#     for i in expected_files:
#         if i in filenames:
#             print_success(f"exists:  {i}")  # extra space aligns print with missing
#     for i in expected_files:
#         if i not in filenames:
#             print_error(f"missing: {i}")
#     for i in filenames:
#         if i not in expected_files:
#             print_error(f"extra:   {i}")


# def check_for_header(neo4j_outdir, expected_headers):
#     print_test_name("Check for expected header files:")
#     import_directory = Path(neo4j_outdir, "neo4j", "import", "headers")
#     print_message("Looking for directories in:")
#     print_message(str(import_directory))
#     for i in Path.iterdir(import_directory):
#         if i.name in expected_headers:
#             print_success(f"exists:  {i.name}")  # extra space aligns print with missing
#         else:
#             print_error(f"missing: {i.name}")
