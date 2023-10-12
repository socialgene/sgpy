# import sqlite3
# from pathlib import Path

# import pandas as pd
# from rich.progress import Progress

# # this assumes that the protein IDs are all unique
# # because this is only used for refseq/genbank data that should be true


# def create_db(tsv_list, outpath, inmem=True):
#     sqliteout = Path(outpath)
#     # files = list(Path("/home/chase/Downloads/a").glob("*.tsv"))
#     if inmem:
#         # build database in-memory then write out to disk at the end
#         conn = sqlite3.connect(":memory:")
#     else:
#         # build database on disk
#         conn = sqlite3.connect(sqliteout)
#     c = conn.cursor()
#     c.execute("""DROP TABLE IF EXISTS protmap""")
#     conn.execute(
#         """
#         CREATE TABLE protmap (
#             hash CHAR(32) NOT NULL,
#             accession VARCHAR(32) NOT NULL,
#             UNIQUE(accession)
#         );
#         """
#     )
#     c.execute(
#         """
#     CREATE INDEX accession_index
#     ON protmap (accession);
#     """
#     )
#     c.execute("""PRAGMA synchronous = OFF""")
#     c.execute("""PRAGMA journal_mode = OFF""")
#     c.execute("""PRAGMA TEMP_STORE = MEMORY""")
#     c.execute("""PRAGMA cache_size = 1000000""")
#     c.execute("""PRAGMA PAGE_SIZE = 1000000""")
#     c.execute("""PRAGMA ignore_check_constraints = true""")
#     with Progress(transient=True) as progress:
#         task = progress.add_task("Progress...", total=len(tsv_list))
#         for i in tsv_list:
#             progress.update(task, advance=1)
#             pd.read_csv(
#                 i, sep="\t", header=None, names=["sha512t24u", "accession"]
#             ).to_sql("protmap", con=conn, if_exists="append", index=False)
#     if inmem:
#         # save in-memory sqlite database to disk
#         db_disk = sqlite3.connect(sqliteout)
#         conn.backup(db_disk)
