import subprocess
import re


def poli_pseudo(arquivo_contig, caminho_db_poli, sample):
    # Run blastx using the database of genes associated with drug resistance
    blastx_polim = [
        "blastx",
        "-db", caminho_db_poli,
        "-query", arquivo_contig,
        "-evalue", "0.001",
        "-out", f"{sample}_BLASTXpolimixina"
    ]

    result = subprocess.run(blastx_polim, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"system {blastx_polim} failed: {result.returncode}")

    # Open the blastx output file
    resultado_polim = f"{sample}_BLASTXpolimixina"
    try:
        with open(resultado_polim, "r") as infile:
            lines = infile.readlines()
    except FileNotFoundError:
        raise FileNotFoundError(f"File {resultado_polim} not open")

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
            teste = 2

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

            if id_porc > 80:
                result_entry = f"{sbjct},id:{identities} {id_porc}%,gap:{gaps}|"
                # Uncomment the next line if you want to use `result_entry` in some way
                # result.append(result_entry)

        if re.match(r"(Query\s\s(\d*)\s*)(\S*)\s*(\d*)", line, re.IGNORECASE) and id_porc > 80:
            space = len(re.match(r"(Query\s\s(\d*)\s*)(\S*)\s*(\d*)",
                        line, re.IGNORECASE).group(1))
            aa_1 = list(re.match(r"(Query\s\s(\d*)\s*)(\S*)\s*(\d*)",
                        line, re.IGNORECASE).group(3))

        if re.match(r"^\s{%d}(.*)" % space, line, re.IGNORECASE) and id_porc > 80:
            aa_alinh = list(
                re.match(r"^\s{%d}(.*)" % space, line, re.IGNORECASE).group(1))

        m = re.match(r"Sbjct\s\s(\d*)\s*(\S*)\s*(\d*)", line, re.IGNORECASE)
        if m and id_porc > 80:
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

                        truncation_condition = tamanho[1] < (
                            (tam_ref / 100) * 90)
                        mutation_condition = {
                            "PmrA": "PmrA",
                            "PmrB": "PmrB",
                            "PhoQ": "PhoQ",
                            "ParR": "ParR",
                            "ParS": "ParS",
                            "CrpS": "CrpS",
                            "ColR": "ColR",
                            "ColS": "ColS"
                        }

                        if sbjct in mutation_condition:
                            if truncation_condition:
                                mutation = f"{sbjct} truncation: {tamanho[1]}/{tam_ref},"
                                result.append(mutation)
                                continue
                            else:
                                mutation = f"{sbjct}:{aa_2[i]}{position}{aa_1[i]},"
                                if sbjct == "PmrA" and re.match(r".*T31I.*", mutation, re.IGNORECASE):
                                    provean = re.match(
                                        r"(.*T31I.*),", mutation, re.IGNORECASE).group(1) + '(deleterious),'
                                    print(f"\n\n{provean}\n\n")
                                    result.append(provean)
                                else:
                                    result.append(mutation)

    return result
