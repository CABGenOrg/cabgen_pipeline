import re
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
    """
    Runs BLASTx and checks the mutations for a specific bacteria.

    Args:
        bacteria_dict (BacteriaDict): Dictionary that contains the necessary
        information to run the function.

    Returns:
        Tuple[List[str], List[str]]: A tuple of lists that contain the found
        mutations.
    """
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


def get_abricate_result(file_path: str):
    """
    Processes Abricate result file and returns lines with identity > 90 and
    coverage > 90 or containing gene names starting with 'Van'.

    Args:
        file_path (str): The path to the abricate result file, a tab-delimited
        file.

    Returns:
        List[str]: A list of specific lines from the file.
    """
    results = []

    try:
        with open(file_path, "r") as infile:
            lines = [line.strip() for line in infile.readlines()]
    except FileNotFoundError:
        raise FileNotFoundError(f"File {file_path} not open.")

    for line in lines:
        fields = line.split('\t')

        if len(fields) < 11:
            continue

        try:
            # Check the specific data
            coverage = float(fields[9])
            identity = float(fields[10])
            gene = fields[5]
        except ValueError:
            continue

        if (coverage > 90 and identity > 90) or \
                re.search(r'Van.*', gene, re.I):
            results.append(line)

    return results
