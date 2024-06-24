import subprocess
import re


def poli_truncation(arquivo_contig, caminho_db_poli, sample):
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
    linha_sbjct = 0
    posi_subjct = 0
    linha_cabecario = 0
    i = 0

    for line in lines:
        line = line.strip()
        i += 1
        # Identify the subject, i.e., the sequence that matched the database
        m = re.match(r"^>(\w*)\|.*", line, re.IGNORECASE)
        if m:
            sbjct = m.group(1)
            linha_cabecario = i

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

        if re.match(r"(Query\s\s(\d*)\s*)(\S*)\s*(\d*)", line, re.IGNORECASE):
            space = len(re.match(r"(Query\s\s(\d*)\s*)(\S*)\s*(\d*)",
                        line, re.IGNORECASE).group(1))
            aa_1 = list(re.match(r"(Query\s\s(\d*)\s*)(\S*)\s*(\d*)",
                        line, re.IGNORECASE).group(3))

        if re.match(r"^\s{%d}(.*)" % space, line, re.IGNORECASE):
            aa_alinh = list(
                re.match(r"^\s{%d}(.*)" % space, line, re.IGNORECASE).group(1))

        m = re.match(r"Sbjct\s\s(\d*)\s*(\S*)\s*(\d*)", line, re.IGNORECASE)
        if m:
            aa_2 = list(m.group(2))
            linha_sbjct = i
            posi_subjct = int(m.group(1))
            if (linha_sbjct == (linha_cabecario + 9)) and (posi_subjct != 1):
                continue
            else:
                if id_porc < 90:
                    truncation = f"{sbjct} truncada, {id_porc}% presente"
                    result.append(truncation)
                if id_porc > 90:
                    for j in range(60):
                        if aa_1[j].upper() != aa_alinh[j].upper():
                            if posi_subjct > int(m.group(3)):
                                position = posi_subjct - j
                            else:
                                position = posi_subjct + j

                            if sbjct == "PmrB" and tamanho[1] > 360:
                                mutation = f"{sbjct}:{aa_2[j]}{position}{aa_1[j]},"
                                if re.match(r".*(R256G|L16P).*", mutation, re.IGNORECASE):
                                    provean = re.match(
                                        r"(.*(R256G|L16P).*),", mutation, re.IGNORECASE).group(1) + '(deleterious),'
                                    result.append(provean)
                                else:
                                    result.append(mutation)

                            if sbjct == "PmrA" and tamanho[1] > 220:
                                mutation = f"{sbjct}:{aa_2[j]}{position}{aa_1[j]},"
                                result.append(mutation)

                            if sbjct == "MgrB" and tamanho[1] > 40:
                                mutation = f"{sbjct}:{aa_2[j]}{position}{aa_1[j]},"
                                result.append(mutation)

                            if sbjct == "PhoP" and tamanho[1] > 220:
                                mutation = f"{sbjct}:{aa_2[j]}{position}{aa_1[j]},"
                                result.append(mutation)

                            if sbjct == "PhoQ" and tamanho[1] > 480:
                                mutation = f"{sbjct}:{aa_2[j]}{position}{aa_1[j]},"
                                if re.match(r".*(S217R|L203Q|Q424L).*", mutation, re.IGNORECASE):
                                    provean = re.match(
                                        r"(.*(S217R|L203Q|Q424L).*),", mutation, re.IGNORECASE).group(1) + '(deleterious),'
                                    result.append(provean)
                                else:
                                    result.append(mutation)

    return result
