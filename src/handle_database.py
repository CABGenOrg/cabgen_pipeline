from pymongo import MongoClient


class MongoSaver:
    def __init__(self, sample_id: int):
        self.database_url = "mongodb://localhost:27017/"
        self.database_name = "sbgmi"
        self.query = {"sequenciaId": sample_id}

    def connect(self):
        try:
            client = MongoClient(self.database_url)
            self.db = client[self.database_name]
        except Exception as error:
            raise Exception(f"Could not connect to MongoDB.\n{error}")

    def _get_db(self):
        if self.db is None:
            self.connect()
        return self.db

    def save(self, key: str, value: str):
        try:
            self._get_db()
            bson = {"$set": {key: value}}
            sample = self.query.get("sequenciaId")
            self.db.insert_one(self.query, bson, upsert=True)
        except Exception as error:
            raise Exception(f"Could not update {sample}.\n{error}")
