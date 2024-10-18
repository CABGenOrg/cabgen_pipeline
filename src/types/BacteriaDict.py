from typing import TypedDict


class BacteriaDict(TypedDict):
    species: str
    assembly_file: str
    sample: str
    others_db_path: str
    poli_db_path: str
    others_outfile_suffix: str
    poli_outfile_suffix: str
