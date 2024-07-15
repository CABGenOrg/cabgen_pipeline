import re
from typing import List, Tuple

sbjct_pattern = re.compile(r">([\w]+)\|", re.I)
sbjct_length_pattern = re.compile(r"Length=(\d+)", re.I)
identities_total_pattern = re.compile(r"Identities = \d+\/(\d+)", re.I)
perc_identity_pattern = re.compile(
    r"Identities = \d+\/\d+\s\((\d+)%\)", re.I)
query_sequence_pattern = re.compile(r"Query\s+\d+\s+(\w+)", re.I)
alignment_sequence_pattern = re.compile(
    r"^\s+[ACDEFGHIKLMNPQRSTVWY\+\s]+$", re.I)
subject_sequence_pattern = re.compile(r"Sbjct\s+(\d+)\s+(\w+)\s+(\d+)", re.I)


def find_mutation(blast_result_path: str, mutations: List[str]) -> List[str]:
    """
    Finds the requested mutations from a BLAST result.

    Args:
        blast_result_path (str): The BLAST result path.
        mutations (List[str]): The list of mutations to find.

    Returns:
        List[str]: A list of the found mutations.
    """
    try:
        with open(blast_result_path, "r") as infile:
            lines = [line.rstrip() for line in infile.readlines()]
    except FileNotFoundError:
        raise FileNotFoundError(f"File {blast_result_path} not open.")

    sbjct = None
    sbjct_length = 0
    identities_total = 0
    perc_identity = 0
    query_sequence = []
    subject_sequence = []
    result = []
    alignment_sequence = []
    position = 0
    test = 0

    for line in lines:
        # Identify the subject
        sbjct_match = sbjct_pattern.search(line)
        if sbjct_match:
            sbjct = sbjct_match.group(1)

        # Identify the length of the subject sequence
        sbjct_length_match = sbjct_length_pattern.search(line)
        if sbjct_length_match:
            sbjct_length = int(sbjct_length_match.group(1))

        # Identify identities
        identities_total_match = identities_total_pattern.search(line)
        perc_identity_match = perc_identity_pattern.search(line)
        if identities_total_match and perc_identity_match:
            identities_total = int(identities_total_match.group(1))
            perc_identity = int(perc_identity_match.group(1))

            # Check for truncation
            if (identities_total < (sbjct_length / 100) * 90) \
                    and (perc_identity > 80):
                mutation = (f"{sbjct} truncation: "
                            f"{identities_total}/{sbjct_length},")
                result.append(mutation)

        # Identify query sequence line
        query_sequence_match = query_sequence_pattern.search(line)
        if query_sequence_match and perc_identity > 90:
            query_sequence = list(query_sequence_match.group(1).strip())

        # Identify alignment line
        alignment_sequence_match = alignment_sequence_pattern.search(line)
        if alignment_sequence_match and perc_identity > 90:
            alignment_sequence = list(alignment_sequence_match.group().strip())

        # Identify subject sequence line
        subject_sequence_match = subject_sequence_pattern.search(line)
        if subject_sequence_match and perc_identity > 90:
            subject_sequence = list(subject_sequence_match.group(2).strip())
            subject_sequence_start = int(subject_sequence_match.group(1))
            subject_sequence_end = int(subject_sequence_match.group(3))

            if subject_sequence_start == 1:
                test = 1

            if test == 1:
                for i in range(60):
                    if i > len(query_sequence) - 1:
                        break
                    query_aa = query_sequence[i].upper()
                    alignment_aa = alignment_sequence[i].upper()
                    subject_aa = subject_sequence[i].upper()

                    if query_aa != alignment_aa:
                        if subject_sequence_start > subject_sequence_end:
                            position = subject_sequence_start - i
                        else:
                            position = subject_sequence_start + i

                        if sbjct in mutations and query_aa != subject_aa:
                            if identities_total > (sbjct_length / 100) * 90:
                                mutation = (f"{sbjct}:{subject_aa}"
                                            f"{position}{query_aa},")
                                result.append(mutation)

    return result


def find_acineto_mutations(blast_result_path: str) -> Tuple[List[str],
                                                            List[str]]:
    """
    Finds the requested mutations for Acineto sp from a BLAST result.

    Args:
        blast_result_path (str): The BLAST result path.

    Returns:
        Tuple[List[str], List[str]]: A tuple with the lists of the found
        mutations.
    """
    other_mutations = ["GyrA", "GyrB", "ParC", "AdeN",
                       "AdeR", "CarO", "OmpA", "AdeL", "AdeS"]
    poli_mutations = ["PmrA", "PmrB", "LpxA", "LpxD", "LpxC"]

    other_result = find_mutation(blast_result_path, other_mutations)
    poli_result = find_mutation(blast_result_path, poli_mutations)
    return other_result, poli_result


def find_ecloacae_mutations(blast_result_path: str) -> Tuple[List[str],
                                                             List[str]]:
    """
    Finds the requested mutations for Enterobacter cloacae from a BLAST result.

    Args:
        blast_result_path (str): The BLAST result path.

    Returns:
        Tuple[List[str], List[str]]: A tuple with the lists of the found
        mutations.
    """
    other_mutations = ["GyrA", "ParC"]
    poli_mutations = ["PmrA", "PmrB", "MgrB", "PhoP", "PhoQ"]

    other_result = find_mutation(blast_result_path, other_mutations)
    poli_result = find_mutation(blast_result_path, poli_mutations)
    return other_result, poli_result


def find_kleb_mutations(blast_result_path: str) -> Tuple[List[str],
                                                         List[str]]:
    """
    Finds the requested mutations for Klebsiella sp from a BLAST result.

    Args:
        blast_result_path (str): The BLAST result path.

    Returns:
        Tuple[List[str], List[str]]: A tuple with the lists of the found
        mutations.
    """
    other_mutations = ["GyrA", "GyrB", "ParC", "AcrR", "RamR"]
    poli_mutations = ["PmrB", "PmrA", "MgrB", "PhoP", "PhoQ"]

    other_result = find_mutation(blast_result_path, other_mutations)
    poli_result = find_mutation(blast_result_path, poli_mutations)
    return other_result, poli_result


def find_pseudo_mutations(blast_result_path: str) -> Tuple[List[str],
                                                           List[str]]:
    """
    Finds the requested mutations for Pseudomonas sp from a BLAST result.

    Args:
        blast_result_path (str): The BLAST result path.

    Returns:
        Tuple[List[str], List[str]]: A tuple with the lists of the found
        mutations.
    """
    other_mutations = ["OprD", "MexT", "AmpC",
                       "AmpR", "GyrA", "GyrB", "ParC", "ParE"]
    poli_mutations = ["PmrA", "PmrB", "PhoQ",
                      "ParR", "ParS", "CrpS", "ColR", "ColS"]

    other_result = find_mutation(blast_result_path, other_mutations)
    poli_result = find_mutation(blast_result_path, poli_mutations)
    return other_result, poli_result
