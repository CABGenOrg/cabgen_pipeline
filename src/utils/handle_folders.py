import os
from typing import List
from shutil import rmtree


def delete_folders_and_files(paths: List[str]):
    try:
        for path in paths:
            if os.path.isfile(path):
                os.remove(path)
            elif os.path.isdir(path):
                rmtree(path)
    except Exception as e:
        print(f"Failed to delete {e}.")
