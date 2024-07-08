import re
from pymongo import MongoClient

def MongoSaver(collection, numero_da_amotra, tipo_de_resultado, imprimir):
    print(f"save_result_mongo::MongoSaver = {collection} - {numero_da_amotra} - {tipo_de_resultado} - {imprimir}")
    if imprimir:
        # Removing unwanted characters
        imprimir = re.sub(r"[\'\"\\R\\t\\n]", "", imprimir)

        # Creating the BSON update operation
        bson = {
            "$set": {
                "sequenciaId": numero_da_amotra,
                tipo_de_resultado: imprimir
            }
        }

        print(f"bson: {bson}")

        # Connecting to MongoDB and performing the update operation
        client = MongoClient('mongodb://localhost:27017/')
        db = client.sgbmi
        collection = db[collection]
        collection.update_one({"sequenciaId": numero_da_amotra}, bson, upsert=True)
