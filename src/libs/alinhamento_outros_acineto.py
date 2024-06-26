import subprocess
import re


def outros_acineto(montagem, fasta_outros, sample):
    # Run blastx using the database of genes associated with drug resistance
    blastx_outros = [
        "blastx",
        "-db", fasta_outros,
        "-query", montagem,
        "-evalue", "0.001",
        "-out", f"{sample}_BLASTXoutros"
    ]

    result = subprocess.run(blastx_outros, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"system {blastx_outros} failed: {result.returncode}")

    # Open the blastx output file
    resultado_outros = f"{sample}_BLASTXoutros"
    try:
        with open(resultado_outros, "r") as infile:
            lines = infile.readlines()
    except FileNotFoundError:
        raise FileNotFoundError(f"File {resultado_outros} not open")

    sbjct = ""
    expect = ""
    identities = ""
    gaps = ""
    result = []
    id_porc = 0
    space = 0
    aa_1 = []
    aa_2 = []
    aa_alinh = []
    position = 0
    tamanho = []
    teste = 0
    tam_ref = 0

    for line in lines:
        line = line.strip()
        # Identify the subject, i.e., the sequence that matched the database
        m = re.match(r"^>(\w*)\|.*", line, re.IGNORECASE)
        if m:
            sbjct = m.group(1)

        m = re.match(r"^Length=(\d*)", line, re.IGNORECASE)
        if m:
            tam_ref = int(m.group(1))

        m = re.match(r".*Expect\s=\s(.*),\sMethod.*", line, re.IGNORECASE)
        if m:
            expect = m.group(1)

        m = re.match(
            r".*Identities\s=\s(\d*/\d*)\s\((\d*)%\),\s.*Gaps\s=\s(\d*/\d*)\s.*", line, re.IGNORECASE)
        if m:
            identities = m.group(1)
            tamanho = [int(x) for x in identities.split('/')]
            id_porc = int(m.group(2))
            gaps = m.group(3)

            if tamanho[1] < ((tam_ref / 100) * 90) and id_porc > 80:
                mutation = f"{sbjct} truncation: {tamanho[1]}/{tam_ref},"
                result.append(mutation)

        if re.match(r"(Query\s\s(\d*)\s*)(\S*)\s*(\d*)", line, re.IGNORECASE) and id_porc > 90:
            space = len(re.match(r"(Query\s\s(\d*)\s*)(\S*)\s*(\d*)",
                        line, re.IGNORECASE).group(1))
            aa_1 = list(re.match(r"(Query\s\s(\d*)\s*)(\S*)\s*(\d*)",
                        line, re.IGNORECASE).group(3))

        if re.match(r"^\s{%d}(.*)" % space, line, re.IGNORECASE) and id_porc > 90:
            aa_alinh = list(
                re.match(r"^\s{%d}(.*)" % space, line, re.IGNORECASE).group(1))

        m = re.match(r"Sbjct\s\s(\d*)\s*(\S*)\s*(\d*)", line, re.IGNORECASE)
        if m and id_porc > 90:
            inicio = int(m.group(1))
            final = int(m.group(3))
            aa_2 = list(m.group(2))
            if inicio == 1:
                teste = 1

            if teste == 1:
                for i in range(60):
                    if aa_1[i].upper() != aa_alinh[i].upper():
                        if inicio > final:
                            position = inicio - i
                        else:
                            position = inicio + i

                        if sbjct in ["GyrA", "GyrB", "ParC", "AdeN", "AdeR", "CarO", "OmpA", "AdeL", "AdeS"]:
                            if tamanho[1] > ((tam_ref / 100) * 90):
                                mutation = f"{sbjct}:{aa_2[i]}{position}{aa_1[i]},"
                                result.append(mutation)

    return result
