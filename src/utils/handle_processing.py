import re
from os import path
from collections import Counter
from typing import List, Tuple, Union
from src.utils.handle_programs import run_blastx
from src.types.SpeciesDict import SpeciesDict
from src.types.BacteriaDict import BacteriaDict
from src.utils.handle_mutations import find_acineto_mutations, \
    find_ecloacae_mutations, find_kleb_mutations, find_pseudo_mutations

choose_analysis = {"pseudomonasaeruginosa": find_pseudo_mutations,
                   "klebsiellapneumoniae": find_kleb_mutations,
                   "Enterobacter_cloacae_subsp_cloacae":
                   find_ecloacae_mutations,
                   "Acinetobacter_baumannii": find_acineto_mutations}


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
    try:
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
    except Exception as e:
        raise Exception(f"Failed to run blast and check mutations.\n\n{e}")


def get_abricate_result(file_path: str) -> List[str]:
    """
    Processes Abricate result file and returns lines with identity > 90 and
    coverage > 90 or containing gene names starting with "Van".

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
        fields = line.split("\t")

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
                re.search(r"Van.*", gene, re.I):
            results.append(line)

    return results


def count_kraken_words(kraken_output: str) -> Tuple[str, str, int, int]:
    """
    Processes Kraken result file and returns the two most common identified
    bacteria species.

    Args:
        kraken_output (str): The path to the Kraken result file.

    Returns:
        Tuple[str, str]: A tuple of the two most common bacteria species.
    """
    try:
        with open(kraken_output) as infile:
            lines = [line.split("\t")[2].split("(")[0].strip()
                     for line in infile.readlines()
                     if len(line.split("\t")) >= 3]

    except FileNotFoundError:
        raise FileNotFoundError(f"File {kraken_output} not found")

    most_common_species = Counter(lines).most_common(2)
    first_most_common = most_common_species[0][0]
    first_count = most_common_species[0][1]
    second_most_common = most_common_species[1][0]
    second_count = most_common_species[1][1]

    return first_most_common, second_most_common, first_count, second_count


def build_species_data(species_info: SpeciesDict) -> dict:
    others_db_path = species_info.get("others_db_path")
    poli_db_path = species_info.get("poli_db_path")
    fastani_path = species_info.get("fastani_db_path")

    species_data = {
        "pseudomonasaeruginosa": {
            "mlst": "paeruginosa",
            "display_name": "Pseudomonas aeruginosa",
            "poli_fasta": f"{poli_db_path}/proteins_pseudo_poli.fasta",
            "others_fasta": (f"{others_db_path}/"
                             "proteins_outrasMut_pseudo.fasta"),
            "run_blast": True
        },
        "escherichiacoli": {
            "mlst": "ecoli",
            "display_name": "Escherichia coli"
        },
        "staphylococcusaureus": {
            "mlst": "saureus",
            "display_name": "Staphylococcus aureus"
        },
        "streptococcuspyogenes": {
            "mlst": "spyogenes",
            "display_name": "Streptococcus pyogenes"
        },
        "pseudomonasputida": {
            "mlst": "pputida",
            "display_name": "Pseudomonas putida"
        },
        "listeriamonocytogenes": {
            "mlst": "lmonocytogenes",
            "display_name": "Listeria monocytogenes"
        },
        "enterococcusfaecalis": {
            "mlst": "efaecalis",
            "display_name": "Enterococcus faecalis"
        },
        "klebsiellaoxytoca": {
            "mlst": "koxytoca",
            "display_name": "Klebsiella oxytoca"
        },
        "enterococcusfaecium": {
            "mlst": "efaecium",
            "display_name": "Enterococcus faecium"
        },
        "klebsiellapneumoniae": {
            "mlst": "kpneumoniae",
            "display_name": "Klebsiella pneumoniae",
            "poli_fasta": f"{poli_db_path}/proteins_kleb_poli.fasta",
            "others_fasta": f"{others_db_path}/proteins_outrasMut_kleb.fasta",
            "run_blast": True,
            "fastani_list": f"{fastani_path}/kleb_database/lista-kleb"
        },
        "enterobacter_species": {
            "mlst": "ecloacae",
            "display_name": "Enterobacter cloacae subsp cloacae",
            "poli_fasta": f"{poli_db_path}/proteins_Ecloacae_poli.fasta",
            "others_fasta": (f"{others_db_path}/"
                             "proteins_outrasMut_Ecloacae.fasta"),
            "run_blast": True,
            "fastani_list": f"{fastani_path}/fastANI/list_entero"
        },
        "acinetobacter_species": {
            "mlst": "abaumannii_2",
            "display_name": "Acinetobacter baumannii",
            "poli_fasta": f"{poli_db_path}/proteins_acineto_poli.fasta",
            "others_fasta": (f"{others_db_path}/"
                             "proteins_outrasMut_acineto.fasta"),
            "run_blast": True,
            "fastani_list": f"{fastani_path}/fastANI_acineto/list-acineto"
        }
    }

    return species_data


def handle_species(species_info: SpeciesDict, species_data: dict) -> \
        Tuple[Union[Tuple[List[str], List[str]], None], Union[str, None],
              Union[str, None]]:
    species = species_info.get("species")
    desired_species_data = species_data.get(species)

    if desired_species_data:
        assembly = species_info.get("assembly")
        sample = species_info.get("sample")

        mlst_species = desired_species_data.get("mlst")
        print_species = desired_species_data.get("display_name")
        poli_fasta = desired_species_data.get("poli_fasta")
        others_fasta = desired_species_data.get("others_fasta")
        output_path = species_info.get("output_path")

        if poli_fasta and others_fasta:
            bacteria_dict: BacteriaDict = {
                "species": species,
                "assembly_file": assembly,
                "sample": str(sample),
                "others_db_path": others_fasta,
                "poli_db_path": poli_fasta,
                "others_outfile_suffix": path.join(output_path, "blastOthers"),
                "poli_outfile_suffix": path.join(output_path, "blastPoli")
            }
            return run_blast_and_check_mutations(bacteria_dict), \
                print_species, mlst_species
    elif "acinetobacter" in species:
        acineto = species_data.get("acinetobacter_species", {})
        #print_species = acineto.get("display_name") 
        mlst_species = acineto.get("mlst")

        return None, None, mlst_species
    elif "enterobacter" in species:
        entero = species_data.get("enterobacter_species", {})
        #print_species = entero.get("display_name")
        mlst_species = entero.get("mlst")

        return None, None, mlst_species

    return None, None, None


def identify_bacteria_species(species_info: SpeciesDict):
    species_data = build_species_data(species_info)

    blast_result, display_name, mlst = handle_species(species_info,
                                                      species_data)
    return blast_result, display_name, mlst


def process_resfinder(abricate_result: List[str]) -> Tuple[List[str],
                                                           List[str]]:
    gene_results = []
    blast_out_results = []
    ref_list = []

    #open reference file add feed the list
    with open("/cabgen/sequences_database/lista_ncbi_ReferenceGeneCatalog160725.txt") as f: #esse caminho precisa ser alterado de acordo com o servidor
        for line in f: 
            line = line.strip() 
            ref_list.append(line.split("\t")) 

    for line in abricate_result:
        blast_lines = line.split("\t")

        id = blast_lines[10]
        gene = blast_lines[5]
        cov_q = blast_lines[9]
        cov_db = blast_lines[6]
    
        fields = line.strip().split("\t") #novo
        name_gene = fields[5].split("_") #novo

        for ref_item in ref_list:
            #print(f"{ref_item[1]}")
            #get correpondence between the two lists
            antibiotic = re.search(re.escape(name_gene[0]), ref_item[0], re.IGNORECASE)
            #print gene and antibiotic reference
            if antibiotic:
                #print(f"{fields[5]} {ref_item[-17]}")
                gene_results.append(f"{gene} (resistance to {ref_item[-17].lower()}) "
                                    f"(allele confidence {id})")
                break
        else:
            gene_results.append(f"{gene} (allele confidence {id})")

        blast_out = (f"{gene} (ID: {id} COV_Q: {cov_q} COV_DB: {cov_db})")
        blast_out_results.append(blast_out)

    return gene_results, blast_out_results


def process_vfdb(abricate_result: List[str]) -> List[str]:
    blast_out_results = []

    for line in abricate_result:
        blast_lines = line.split("\t")

        id = blast_lines[10]
        gene = blast_lines[5]
        cov_q = blast_lines[9]
        cov_db = blast_lines[6]
        sequence = blast_lines[1]
        product = blast_lines[13]

        blast_out = (f"{sequence}: {gene} {product} ID: {id} COV_Q: {cov_q} "
                     f"COV_DB: {cov_db}| ")
        blast_out_results.append(blast_out)

    return blast_out_results


def process_plasmidfinder(abricate_result: List[str]) -> List[str]:
    blast_out_results = []

    for line in abricate_result:
        blast_lines = line.split("\t")

        id = blast_lines[10]
        gene = blast_lines[5]
        cov_q = blast_lines[9]
        cov_db = blast_lines[6]

        blast_out = (f"{gene} (ID: {id} COV_Q: {cov_q} COV_DB: {cov_db})")
        blast_out_results.append(blast_out)

    return blast_out_results


def format_time(seconds: float) -> str:
    hours, rem = divmod(seconds, 3600)
    minutes, seconds = divmod(rem, 60)

    return f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}"


def handle_fastani_species(species_info: SpeciesDict,
                           fastani_species: str) -> Union[Tuple[List[str],
                                                          List[str]], None]:

    species_data = build_species_data(species_info)

    species = fastani_species
    assembly_file = species_info.get("assembly")
    sample = str(species_info.get("sample"))
    others_outfile_suffix = path.join(species_info.get("output_path"),
                                      "blastOthers")
    poli_outfile_suffix = path.join(species_info.get("output_path"),
                                    "blastPoli")

    fastani_species_dict = {"Acinetobacter_baumannii":
                            species_data.get("acinetobacter_species", {}),
                            "Enterobacter_cloacae_subsp_cloacae":
                            species_data.get("enterobacter_species", {})}
    species_dict = fastani_species_dict.get(fastani_species)

    if species_dict:
        others_db_path = species_dict.get("others_fasta", "")
        poli_db_path = species_dict.get("poli_fasta", "")
    else:
        return None

    bacteria_dict: BacteriaDict = {
        "species": species,
        "assembly_file": assembly_file,
        "sample": sample,
        "others_db_path": others_db_path,
        "poli_db_path": poli_db_path,
        "others_outfile_suffix": others_outfile_suffix,
        "poli_outfile_suffix": poli_outfile_suffix
    }
    blast_result = run_blast_and_check_mutations(bacteria_dict)
    return blast_result
