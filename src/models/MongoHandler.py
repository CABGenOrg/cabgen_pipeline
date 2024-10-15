from pymongo import MongoClient


class MongoHandler:
    def __init__(self, database_url="mongodb://localhost:27017/",
                 db_name="sgbmi"):
        self.database_url = database_url
        self.db_name = db_name
        self.client = MongoClient(self.database_url)
        self.db = self.client[self.db_name]

    def search(self, collection_name: str, match=None, lookup=None,
               project=None):
        collection = self.db[collection_name]
        pipeline = []

        if match:
            pipeline.append({"$match": match})

        if lookup:
            pipeline.append({"$lookup": lookup})

        if project:
            pipeline.append({"$project": project})

        results = [res for res in collection.aggregate(pipeline)]
        return results

    def save(self, collection_name: str, query: dict, bson: dict):
        try:
            collection = self.db[collection_name]

            if "$set" not in bson:
                bson = {"$set": bson}

            collection.update_one(query, bson, upsert=True)
        except Exception as error:
            raise Exception(f"Could not update document.\n\n{error}")

    def close(self):
        self.client.close()
