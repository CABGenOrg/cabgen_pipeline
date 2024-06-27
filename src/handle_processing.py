from handle_programs import run_blastx
from typing import TypedDict, List, Tuple
from handle_mutations import find_acineto_mutations, find_ecloacae_mutations, \
    find_kleb_mutations, find_pseudo_mutations

choose_analysis = {"pseudomonasaeruginosa": find_pseudo_mutations,
                   "klebsiellapneumoniae": find_kleb_mutations,
                   "Enterobacter_cloacae_subsp_cloacae":
                   find_ecloacae_mutations,
                   "Acinetobacter_baumannii": find_acineto_mutations}


class BacteriaDict(TypedDict):
    species: str
    assembly_file: str
    sample: str
    others_db_path: str
    poli_db_path: str
    others_outfile_suffix: str
    poli_outfile_suffix: str


def run_blast_and_check_mutations(
        bacteria_dict: BacteriaDict) -> Tuple[List[str], List[str]]:
    species = bacteria_dict["species"]
    assembly_file = bacteria_dict["assembly_file"]
    sample = bacteria_dict["sample"]
    others_db_path = bacteria_dict["others_db_path"]
    poli_db_path = bacteria_dict["poli_db_path"]
    others_outfile_suffix = bacteria_dict["others_outfile_suffix"]
    poli_outfile_suffix = bacteria_dict["poli_outfile_suffix"]

    others_blast_result = run_blastx(
        assembly_file, others_db_path, sample, others_outfile_suffix)
    poli_blast_result = run_blastx(assembly_file, poli_db_path,
                                   sample, poli_outfile_suffix)

    analysis = choose_analysis.get(species, None)
    if not analysis:
        raise Exception(f"Invalid species {species}!")

    others_mutations_result, _ = analysis(others_blast_result)
    _, poli_mutations_result = analysis(poli_blast_result)

    return others_mutations_result, poli_mutations_result
