from typing import TypedDict


class SpeciesDict(TypedDict):
    species: str
    assembly: str
    sample: int
    others_db_path: str
    poli_db_path: str
    fastani_db_path: str
