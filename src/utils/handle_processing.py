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

resistance_pattern = re.compile(
    r'(?P<carbapenemase>(blaKPC|blaNDM|blaVIM|blaIMP|blaSPM|blaOXA-23|'
    r'blaOXA-24|blaOXA-25|blaOXA-26|blaOXA-27|blaOXA-48|blaOXA-58|'
    r'blaOXA-72|blaOXA-98|blaOXA-116|blaOXA-117|blaOXA-160|blaOXA-175|'
    r'blaOXA-176|blaOXA-253))|'
    r'(?P<oxa_51_like>(blaOXA-51|blaOXA-64|blaOXA-65|blaOXA-66|blaOXA-69|'
    r'blaOXA-90|blaOXA-125|blaOXA-259|blaOXA-343|blaOXA-407))|'
    r'(?P<esbl>(blaTEM|blaSHV|blaADC|blaCTX-M|blaGES|blaOXA-(?!23|24|25|'
    r'26|27|48|58|72|98|116|117|160|175|176|253|488|486)))|'
    r'(?P<aminoglycosides_fluoroquinolones>(aac\(6\'\)-Ib-cr))|'
    r'(?P<aminoglycosides>(aph|aac|rmt|aad))|'
    r'(?P<chloramphenicol>(cat|cml|cmx|floR))|'
    r'(?P<fluoroquinolones>(qnr|oqx))|'
    r'(?P<sulfonamides>sul)|'
    r'(?P<trimethoprim>dfrA)|'
    r'(?P<tetracycline>tet)|'
    r'(?P<erythromycin>ere)|'
    r'(?P<lincosamides_macrolides_streptogramins>erm)|'
    r'(?P<rifampicin>ARR)|'
    r'(?P<macrolides>(mph|msr))|'
    r'(?P<vancomycin>Van)|'
    r'(?P<clindamycin>lsa)|'
    r'(?P<polymyxin>mcr)',
    re.I
)

resfinder_patterns = {
    'carbapenemase': "(carbapenemase)",
    'oxa_51_like': "(OXA-51-like carbapenemase)",
    'esbl': "(ESBL)",
    'aminoglycosides_fluoroquinolones': ("(resistance to aminoglycosides and "
                                         "fluoroquinolones)"),
    'aminoglycosides': "(resistance to aminoglycosides)",
    'chloramphenicol': "(resistance to chloramphenicol)",
    'fluoroquinolones': "(resistance to fluoroquinolones)",
    'sulfonamides': "(resistance to sulfonamides)",
    'trimethoprim': "(resistance to trimethoprim)",
    'tetracycline': "(resistance to tetracycline)",
    'erythromycin': "(resistance to erythromycin)",
    'lincosamides_macrolides_streptogramins':
        "(resistance to lincosamides, macrolides, and streptogramins)",
    'rifampicin': "(resistance to rifampicin)",
    'macrolides': "(resistance to macrolides)",
    'vancomycin': "(resistance to vancomycin)",
    'clindamycin': "(resistance to clindamycin)",
    'polymyxin': "(resistance to polymyxin)"
}


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


def build_species_data(species_info: SpeciesDict) -> Tuple[dict, dict]:
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
        }
    }

    fastani_species = {
        "klebsiellapneumoniae": {
            "mlst": "kpneumoniae",
            "display_name": "Klebsiella pneumoniae",
            "poli_fasta": f"{poli_db_path}/proteins_kleb_poli.fasta",
            "others_fasta": f"{others_db_path}/proteins_outrasMut_kleb.fasta",
            "fastani_list": f"{fastani_path}/kleb_database/lista-kleb"
        },
        "enterobacter_species": {
            "mlst": "ecloacae",
            "fastani_list": f"{fastani_path}/fastANI/list_entero"
        },
        "acinetobacter_species": {
            "mlst": "abaumannii_2",
            "fastani_list": f"{fastani_path}/fastANI_acineto/list-acineto"
        }
    }

    return species_data, fastani_species


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
        return None, print_species, mlst_species

    return None, None, None


def identify_bacteria_species(species_info: SpeciesDict):
    species_data, _ = build_species_data(species_info)

    blast_result, display_name, mlst = handle_species(species_info,
                                                      species_data)
    return blast_result, display_name, mlst


def process_resfinder(abricate_result: List[str]) -> Tuple[List[str],
                                                           List[str]]:
    gene_results = []
    blast_out_results = []

    for line in abricate_result:
        blast_lines = line.split("\t")

        id = blast_lines[10]
        gene = blast_lines[5]
        cov_q = blast_lines[9]
        cov_db = blast_lines[6]

        blast_out = (f"{gene} (ID: {id} COV_Q: {cov_q} COV_DB: {cov_db})")
        blast_out_results.append(blast_out)

        match = resistance_pattern.search(gene)
        if match:
            for group_name, resistance in resfinder_patterns.items():
                if match.group(group_name):
                    gene_results.append(f"{gene} {resistance} "
                                        f"(allele confidence {id})")
                    break
        else:
            gene_results.append(f"{gene} (allele confidence {id})")
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
