import re
from typing import List, Tuple


def find_mutation(blast_result_path: str, mutations: List[str]) -> List[str]:
    try:
        with open(blast_result_path, "r") as infile:
            lines = infile.readlines()
    except FileNotFoundError:
        raise FileNotFoundError(f"File {blast_result_path} not open")

    result = []
    sbjct = ""
    tam_ref = 0
    id_porc = 0
    tamanho = []
    space = 0
    aa_1 = []
    aa_alinh = []
    teste = 0

    for line in lines:
        line = line.strip()

        sbjct_match = re.match(r"^>([\w:]+).*", line, re.IGNORECASE)
        if sbjct_match:
            sbjct = sbjct_match.group(1)

        length_match = re.match(r"^Length=(\d*)", line, re.IGNORECASE)
        if length_match:
            tam_ref = int(length_match.group(1))

        identities_match = re.match(
            (r".*Identities\s=\s(\d*/\d*)\s\((\d*)%\),\s."
             r"*Gaps\s=\s(\d*/\d*)\s.*"),
            line, re.IGNORECASE)
        if identities_match:
            identities = identities_match.group(1)
            tamanho = [int(x) for x in identities.split('/')]
            id_porc = int(identities_match.group(2))

            if tamanho[1] < ((tam_ref / 100) * 90) and id_porc > 80:
                mutation = f"{sbjct} truncation: {tamanho[1]}/{tam_ref},"
                result.append(mutation)

        query_match = re.match(r"(Query\s\s(\d*)\s*)(\S*)\s*(\d*)", line,
                               re.IGNORECASE)
        if query_match and id_porc > 90:
            space = len(query_match.group(1))
            aa_1 = list(query_match.group(3))

        align_match = re.match(r"^\s{%d}(.*)" % space, line, re.IGNORECASE)
        if align_match and id_porc > 90:
            aa_alinh = list(align_match.group(1))

        sbjct_match = re.match(r"Sbjct\s\s(\d*)\s*(\S*)\s*(\d*)", line,
                               re.IGNORECASE)
        if sbjct_match and id_porc > 90:
            inicio = int(sbjct_match.group(1))
            final = int(sbjct_match.group(3))
            aa_2 = list(sbjct_match.group(2))

            if inicio == 1:
                teste = 1

            if teste == 1:
                for i in range(min(60, len(aa_1), len(aa_alinh))):
                    if aa_1[i].upper() != aa_alinh[i].upper():
                        position = inicio - i if inicio > final else inicio + i

                        if sbjct in mutations and tamanho[1] > \
                                ((tam_ref / 100) * 90):
                            mutation = f"{sbjct}:{aa_2[i]}{position}{aa_1[i]},"
                            result.append(mutation)

    return result


def find_acineto_mutations(blast_result_path: str) -> Tuple[List[str],
                                                            List[str]]:
    other_mutations = ["GyrA", "GyrB", "ParC", "AdeN",
                       "AdeR", "CarO", "OmpA", "AdeL", "AdeS"]
    poli_mutations = ["PmrA", "PmrB", "LpxA", "LpxD", "LpxC"]

    other_result = find_mutation(blast_result_path, other_mutations)
    poli_result = find_mutation(blast_result_path, poli_mutations)
    return other_result, poli_result


def find_ecloacae_mutations(blast_result_path: str) -> Tuple[List[str],
                                                             List[str]]:
    other_mutations = ["GyrA", "ParC"]
    poli_mutations = ["PmrA", "PmrB", "MgrB", "PhoP", "PhoQ"]

    other_result = find_mutation(blast_result_path, other_mutations)
    poli_result = find_mutation(blast_result_path, poli_mutations)
    return other_result, poli_result


def find_kleb_mutations(blast_result_path: str) -> Tuple[List[str],
                                                         List[str]]:
    other_mutations = ["GyrA", "GyrB", "ParC", "AcrR", "RamR"]
    poli_mutations = ["PmrB", "PmrA", "MgrB", "PhoP", "PhoQ"]

    other_result = find_mutation(blast_result_path, other_mutations)
    poli_result = find_mutation(blast_result_path, poli_mutations)
    return other_result, poli_result


def find_pseudo_mutations(blast_result_path: str) -> Tuple[List[str],
                                                           List[str]]:
    other_mutations = ["OprD", "MexT", "AmpC",
                       "AmpR", "GyrA", "GyrB", "ParC", "ParE"]
    poli_mutations = ["PmrA", "PmrB", "PhoQ",
                      "ParR", "ParS", "CrpS", "ColR", "ColS"]

    other_result = find_mutation(blast_result_path, other_mutations)
    poli_result = find_mutation(blast_result_path, poli_mutations)
    return other_result, poli_result
