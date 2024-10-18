from typing import List
from src.models.MongoHandler import MongoHandler


def get_fastqc_tasks() -> List[dict]:
    try:
        handler = MongoHandler()

        match_stage = {"ultimaTarefa": "QUA"}
        lookup_stage = {
            "from": "usuarios",
            "localField": "criadoPor",
            "foreignField": "usuario",
            "as": "usuario"
        }
        project_stage = {
            "arquivofastqr1": 1,
            "arquivofastqr2": 1,
            "email": {"$arrayElemAt": ["$usuario.email", 0]},
            "ultimaTarefa": 1
        }

        tasks = handler.search("sequencias", match_stage,
                               lookup_stage, project_stage)
        handler.close()
        return tasks
    except Exception as e:
        print(f"Can't retrive FastQC tasks.\n\n{e}")
        handler.close()
        return []


def get_complete_tasks():
    try:
        handler = MongoHandler()

        match_stage = {"ultimaTarefa": "TODOS"}
        lookup_stage = {
            "from": "usuarios",
            "localField": "criadoPor",
            "foreignField": "usuario",
            "as": "usuario"
        }
        project_stage = {
            "arquivofastqr1": 1,
            "arquivofastqr2": 1,
            "email": {"$arrayElemAt": ["$usuario.email", 0]},
            "ultimaTarefa": 1
        }

        tasks = handler.search("sequencias", match_stage,
                               lookup_stage, project_stage)
        handler.close()
        return tasks
    except Exception as e:
        print(f"Can't retrive complete tasks.\n\n{e}")
        handler.close()
        return []


def get_genomic_tasks() -> List[dict]:
    try:
        handler = MongoHandler()

        match_stage = {"ultimaTarefa": "ENS"}
        lookup_stage = {
            "from": "usuarios",
            "localField": "criadoPor",
            "foreignField": "usuario",
            "as": "usuario"
        }
        project_stage = {
            "arquivofastqr1": 1,
            "arquivofastqr2": 1,
            "email": {"$arrayElemAt": ["$usuario.email", 0]},
            "ultimaTarefa": 1
        }

        tasks = handler.search("sequencias", match_stage,
                               lookup_stage, project_stage)

        return tasks
    except Exception as e:
        print(f"Can't retrive genomic tasks.\n\n{e}")
        handler.close()
        return []
