#!/usr/bin/env python3
# -*- coding: utf8 -*-

import re
import sys
import shutil
import pathlib
from os import path, makedirs
from argparse import Namespace
from src.handle_database import MongoSaver
from src.handle_programs import run_command_line
from src.handle_processing import run_blast_and_check_mutations, \
    BacteriaDict, get_abricate_result, count_kraken_words


def pipeline(args: Namespace):
    sys.path[0:0] = ['/opt/pipeline/lib/']

    MLST_result = ''
    # _d = None
    a = ''
    b = ''
    # count2 = 0
    # docker = []
    # especieAb = []
    # especieEc = []
    # especieKp = []
    fasta_outros = ''
    fasta_polimixina = ''
    identificacao = []
    lines = []
    lista_especie = ''
    nearest_sts = ''
    preidentificacao = []
    # resultadoANI = ''
    THREADS = "16"  # new one will have 32

    sys.argv = sys.argv[1:]  # dict(perllib.Array(sys.argv)[1:])
    # use alinhamento_poli_truncation;

    """
    Linha de comando do script perl pipeline_melise_output_gal.pl  <caminho do
    diretorio Output, com resultados de spades/unicycler e prokka
    ex:/home/melise/Output_run18-06.02.21> <nome da amostra na pasta Output
    ex:27563_S12> <caminho onde esta instalado o abricate
    ex:./abricate/bin/abricate> <arquivo da tabela excel formato .xls onde vão
    ser impressos os resultados> <diretorio onde esta instalado o kmer-db
    ex:/home/melise/kmer-db> <caminho do diretorio onde esta instalado o mlst
    ex:/home/melise/identificar_clones> <caminho do DB com as seqs de ptn de
    resistencia mutacional a polimixina ex:/home/melise/resis_poli> <caminho
    do DB com as seqs de ptn de resistencia mutacional a outros antibioticos
    ex:/home/melise/outrasMut> <caminho kraken ex:kraken> <caminho unicycle
    ex:unicycler> <caminho R1> <caminho R2>
    """

    # Guardar o caminho do diretorio "Output", com o resultado da montagem,
    # anotacao e blast para todas as amostra. ex /home/melise/Output
    caminho1 = args.output

    # Guardar o nome da amostra com o _S*
    sample = args.sample
    # para pegar o numero da amostra antes do _S1
    sample1 = str(sample).split('_')
    sample2 = sample1[0]
    print(f"Sample: {sample2}")

    # criar um diretório para a amostra
    makedirs(f"{caminho1}/{sample2}", exist_ok=True)

    # caminho para onde esta instalado o abricate ex: ./abricate/bin/abricate

    caminho_abricate = args.abricate

    # entrar com o caminho da pastar onde esta instalado o mlst.
    # ex: /home/melise/identificar_clones
    mlst_install = args.mlst

    # Caminho para banco de dados com sequências de genes de R a polimixina
    db_polimixina = args.polymyxin_db

    # Caminho para banco de dados com sequências de genes de R mutacionais a
    # outros antibioticos
    db_outrosMut = args.outhers_db
    # entrar com o caminho da pastar onde esta instalado o kraken2 ex: kraken2
    kraken2_install = args.kraken2

    # entrar com o caminho do unicycler
    unicycler = args.unicycler

    # entrar com o caminho do fastANI
    fastANI = args.fastANI

    # Caminho para banco de dados com especies para FastANI
    lista = args.fastani_db

    print('Parametros: ')
    print(f"caminho: {caminho1} ")
    print(f"Sample: {path.basename(sample)} ")
    print(f"SAmple2: {path.basename(sample2)} ")
    print(f"abricate: {caminho_abricate} ")
    print(f"mlst install: {mlst_install}  ")
    print(f"db polimixina: {db_polimixina}  ")
    print(f"db outros mut: {db_outrosMut}  ")
    print(f"kraken2_install: {kraken2_install}  ")
    print(f"unicycler: {unicycler} ")
    print(f"fastANI: {fastANI} ")
    print(f"lista: {lista}  ")

    R1 = args.read1
    R2 = args.read2

    mongo_client = MongoSaver(int(sample2))
    mongo_client.connect()
    ##################################################
    # rodar unicycler
    unicycler_line = (f"{unicycler}  -1 {R1} -2 {R2} "
                      f"-o {caminho1}/{path.basename(sample2)}/"
                      "unicycler --min_fasta_length 500 --mode conservative "
                      f"-t {THREADS} --spades_path "
                      "/opt/SPAdes-3.15.5-Linux/bin/spades.py")
    run_command_line(unicycler_line)

    # arquivo assembly.fasta da amostra

    montagem = (f"{caminho1}/{sample2}"
                "/unicycler/assembly.fasta")

    ###########################################################################
    # rodar prokka
    prokka_line = (f"prokka --outdir {caminho1}/"
                   f"{sample2}/prokka --prefix genome {montagem} "
                   "--force --cpus 0")
    run_command_line(prokka_line)

    # EXCLUSIVO DO PIPELINE OUTPUT

    # variavel para guardar os tipos de resultados que devem ser indicados
    # ex: checkm; especie; contigs; resfinder; VFDB; plasmid; mlst;
    # mutacoes_poli; mutacoes_outra

    # tipo_de_resultado = None
    # o que imprimir desse resultado
    imprimir = None

    ###########################################################################
    # PARA IMPRIMIR O RESULTADO EM UM ARQUIVO TXT PARA GAL

    gal_file = open('resultado_gal.txt', mode="a", encoding='utf-8')
    # MELISE: ESTA USANDO ESSE ARQUIVO PARA O MONGODB? PORQUE O PIPELINE NAO
    # PRECISA DELE

    # printar no arquivo final o nome da amostra na primeira linha
    gal_file.write((f"\nAmostra {path.basename(sample)}\nResultados relevantes"
                    " do sequenciamento do genoma total (WGS):\n"))

    ###########################################################################
    # Rodar o CheckM para saber qualidade da corrida
    # Copiar o arquivo assembly.fasta para a pasta do CheckM checkM_bins
    makedirs("checkM_bins", exist_ok=True)
    shutil.copy(path.join(".", f"{montagem}"), path.join(".", 'checkM_bins'))

    # rodar o CheckM
    checkM_line = ("checkm lineage_wf -x fasta checkM_bins checkM_bins "
                   f"--threads {THREADS} --pplacer_threads {THREADS}")
    checkM_qa_line = (f"checkm qa -o 2 -f checkM_bins/{sample2}"
                      "_resultados --tab_table checkM_bins/lineage.ms "
                      f"checkM_bins --threads {THREADS}")

    run_command_line(checkM_line)
    run_command_line(checkM_qa_line)
    # apagar arquivos gerados, deixando apenas resultados
    shutil.rmtree('checkM_bins/bins', ignore_errors=True)
    shutil.rmtree('checkM_bins/storage', ignore_errors=True)
    # MELISE: NAO ENTENDI PORQUE AQUI NÃO REMOVER OS ARQUIVOS
    pathlib.Path('checkM_bins/assembly.fasta').unlink(missing_ok=True)
    pathlib.Path('checkM_bins/lineage.ms').unlink(missing_ok=True)
    pathlib.Path('checkM_bins/checkm.log').unlink(missing_ok=True)

    # pegar o resultado da contaminacao

    contaminacao = 0

    print('Salvando resultado no mongo relatorios')

    genome_size = 1

    with open(f"checkM_bins/{sample2}_resultados") as IN_check:
        next(IN_check)  # ignore header
        for row in IN_check:
            # remove \n of the line end
            row = row.rstrip("\n")
            # separar as colunas do arquivo em elementos em um array
            lines = row.split("\t")
            # print "$lines[2]\n";
            # printar na tabela OUTPUT as informacoes de qualidade que
            # interessam
            # EXCLUSIVO TABELA OUTPUT
            # print OUT2 "$lines[5]\t$lines[6]\t$lines[8]\t";
            genome_size = lines[8]
            mongo_client.save('checkm_1', lines[5])  # completeness
            mongo_client.save('checkm_2', lines[6])  # contamination
            mongo_client.save('checkm_3', lines[8])  # genome size
            mongo_client.save('checkm_4', lines[11])  # contigs
            contaminacao = lines[6]

    mongo_client.save('sample', sample)
    ###########################################################################
    # Identificar especie usando o kraken

    print('rodar o kraken')
    kraken_line = (f"{kraken2_install}/kraken2 "
                   f"--db {kraken2_install}/minikraken2_v2_8GB"
                   "_201904_UPDATE "
                   f"--use-names --paired {R1} {R2} "
                   f"--output out_kraken --threads {THREADS}")
    run_command_line(kraken_line)

    print("splitting output into %s equal files" % THREADS)
    preffix = "krk"
    splitter_line = (f"split --numeric-suffixes=1 -n l/{THREADS} "
                     f"out_kraken {preffix}")
    run_command_line(splitter_line)

    maior_repeticao, segunda_repeticao, \
        first_count, second_count = count_kraken_words("out_kraken")

    # print "$maior_repeticao\n$segunda_repeticao\n";

    # onde sera guardada o nome da especie
    check_especies = maior_repeticao

    # apagar o arquivo do resultado
    # my @rm2 = ("rm", "out_kraken");

    # system(@rm2) == 0
    #        or die "system @rm2 failes: $?";

    # colocar só o genero e a especie, descartando qualquer outra informação
    # dessa coluna do kraken
    # identificar_especie = ''  # mod 11.05.22
    genero = ''  # mod 11.05.22
    especie = ''  # mod 11.05.22
    # juntar genero e especie para mlst
    # resultado_final_especie = ''  # mod 11.05.22 MELISE: ESSA LINHA NAO ESTA
    # SILENCIADA NO SCRIT PERL
    # resultado que sera impresso
    printar_especies = ''  # mod 11.05.22
    # o que usar para mlst
    especie_mlst = ''  # mod 11.05.22

    # mod 11.05.22
    if (re.findall(re.compile(r'\w+\s\w+', re.I), check_especies)):
        check_especies = check_especies.strip()
        genero, especie = check_especies.split(" ")
        # print "$genero\n$especie\n";
        # Associar o nome da especie ao banco do mlst e gravar o nome que sera
        # dado como resultado final
        # resultado_final_especie = f"{genero}{especie}" #MELISE: ESSA LINHA
        # NAO ESTA SILENCIADA NO SCRIT PERL
        # $printar_especies = $resultado_final_especie;
        # $especie_mlst = "";
    else:
        printar_especies = check_especies
        # MELISE: ESSA LINHA NAO EXISTIA NO SCRIPT EM PERL
        genero = check_especies
    # mod ate aqui 20.05.22

    ###########################################################################
    # Sequencia para verificar mutacoes pontuais usando subrotinas proprias

    # guardar o resultado das mutacoes para polimixina

    result2 = []
    # guardar o resultado das mutacoes para outros antibioticos
    result3 = []
    # guardar resultados dos fatores de virulencia
    # vfdb = []

    resultado_final_especie = f"{genero}{especie}".lower()
    print(f"resultado_final_especie: {resultado_final_especie}")
    if resultado_final_especie == 'pseudomonasaeruginosa':
        especie_mlst = 'paeruginosa'
        printar_especies = 'Pseudomonas aeruginosa'
        fasta_polimixina = (f"{db_polimixina}/"
                            "proteins_pseudo_poli.fasta")
        fasta_outros = (f"{db_outrosMut}"
                            "/proteins_outrasMut_pseudo.fasta")
        bacteria_dict: BacteriaDict = {"species": resultado_final_especie,
                                       "assembly_file": montagem,
                                       "sample": sample,
                                       "others_db_path": fasta_outros,
                                       "poli_db_path": fasta_polimixina,
                                       "others_outfile_suffix": "blastOthers",
                                       "poli_outfile_suffix": "blastPoli"}

        result3, result2 = run_blast_and_check_mutations(bacteria_dict)
    elif resultado_final_especie == 'escherichiacoli':
        especie_mlst = 'ecoli'
        printar_especies = 'Escherichia coli'
    elif resultado_final_especie == 'staphylococcusaureus':
        especie_mlst = 'saureus'
        printar_especies = 'Staphylococcus aureus'
    elif resultado_final_especie == 'streptococcuspyogenes':
        especie_mlst = 'spyogenes'
        printar_especies = 'Streptococcus pyogenes'
    elif resultado_final_especie == 'pseudomonasputida':
        especie_mlst = 'pputida'
        printar_especies = 'pseudomonas putida'
    # modificado 19.11.21
    elif resultado_final_especie == 'Listeriamonocytogenes':
        especie_mlst = 'lmonocytogenes'
        printar_especies = 'Listeria monocytogenes'
    elif resultado_final_especie == 'enterococcusfaecalis':
        especie_mlst = 'efaecalis'
        printar_especies = 'Enterococcus faecalis'
    elif resultado_final_especie == 'klebsiellaoxytoca':
        especie_mlst = 'koxytoca'
        printar_especies = 'Klebsiella oxytoca'
    elif resultado_final_especie == 'enterococcusfaecium':
        especie_mlst = 'efaecium'
        printar_especies = 'Enterococcus faecium'
    elif resultado_final_especie in ('klebsiellapneumoniae',
                                     'acinetobacterbaumannii',
                                     "acinetobacternosocomialis",
                                     "acinetobacterpittii",
                                     "acinetobacterseifertii",
                                     "acinetobacterdijkshoorniae",
                                     "acinetobacterlactucae",
                                     "acinetobactercalcoaceticus",
                                     "enterobactercloacae",
                                     "enterobacterhormaechei",
                                     "enterobacterasburiae",
                                     "enterobacterkobei",
                                     "enterobacterroggenkampii",
                                     "enterobacterludwigii"):
        fastANI_txt = "Para FastANI"
        if resultado_final_especie == 'klebsiellapneumoniae':
            especie_mlst = 'kpneumoniae'
            # $printar_especies = "Klebsiella pneumoniae";
            fasta_polimixina = (f"{db_polimixina}"
                                "/proteins_kleb_poli.fasta")
            fasta_outros = (f"{db_outrosMut}"
                            "/proteins_outrasMut_kleb.fasta")

            bacteria_dict: BacteriaDict = {"species": resultado_final_especie,
                                           "assembly_file": montagem,
                                           "sample": sample,
                                           "others_db_path": fasta_outros,
                                           "poli_db_path": fasta_polimixina,
                                           "others_outfile_suffix":
                                           "blastOthers",
                                           "poli_outfile_suffix": "blastPoli"}

            result3, result2 = run_blast_and_check_mutations(bacteria_dict)
            lista_especie = (f"{lista}"
                             "/kleb_database/lista-kleb")
        elif re.match(r"acinetobacter.*", resultado_final_especie, re.I):
            especie_mlst = 'abaumannii_2' # perform MLST for any species in the genus Acinetobacter
            if resultado_final_especie in ("acinetobacterbaumannii",
                                         "acinetobacternosocomialis",
                                         "acinetobacterpittii",
                                         "acinetobacterseifertii",
                                         "acinetobacterdijkshoorniae",
                                         "acinetobacterlactucae",
                                         "acinetobactercalcoaceticus"):
                # $printar_especies = "Acinetobacter baumannii";
                lista_especie = (f"{lista}"
                                "/fastANI_acineto/list-acineto")
        elif resultado_final_especie in ("enterobactercloacae",
                                         "enterobacterhormaechei",
                                         "enterobacterasburiae",
                                         "enterobacterkobei",
                                         "enterobacterroggenkampii",
                                         "enterobacterludwigii"):
            especie_mlst = "ecloacae"
            lista_especie = (f"{lista}"
                             "/list_entero")
        fastANI_txt = 'Rodar fastANI para subespecie'

        print(fastANI_txt)
        # Abrir o arquivo lista

        fastani_line = (f"{fastANI} -q {caminho1}/"
                        f"{sample2}/unicycler/assembly.fasta "
                        f"--rl {lista_especie} -o {sample2}_out-fastANI "
                        f"--threads {THREADS}")
        run_command_line(fastani_line)

        # abrir output
        # Abrir o arquivo do output de distancia
        # array para guardar especies

        print('resultado do fastANI')
        with open(f"{sample2}_out-fastANI", "r") as IN7:
            # first line only
            especiE = IN7.readline().rstrip("\n").split('\t')
            preidentificacao = especiE[1].split("/")
            identificacao = preidentificacao[-1].split(".")
            printar_especies = identificacao[0]

        if 'Enterobacter_cloacae_subsp_cloacae' == printar_especies:
            fasta_polimixina = (f"{db_polimixina}"
                                "/proteins_Ecloacae_poli.fasta")
            fasta_outros = (f"{db_outrosMut}"
                            "/proteins_outrasMut_Ecloacae.fasta")

            bacteria_dict: BacteriaDict = {"species": printar_especies,
                                           "assembly_file": montagem,
                                           "sample": sample,
                                           "others_db_path": fasta_outros,
                                           "poli_db_path": fasta_polimixina,
                                           "others_outfile_suffix":
                                           "blastOthers",
                                           "poli_outfile_suffix": "blastPoli"}

            result3, result2 = run_blast_and_check_mutations(bacteria_dict)

        if 'Acinetobacter_baumannii' in printar_especies:
            fasta_polimixina = (f"{db_polimixina}"
                                "/proteins_acineto_poli.fasta")
            fasta_outros = (f"{db_outrosMut}"
                            "/proteins_outrasMut_acineto.fasta")

            bacteria_dict: BacteriaDict = {"species": printar_especies,
                                           "assembly_file": montagem,
                                           "sample": sample,
                                           "others_db_path": fasta_outros,
                                           "poli_db_path": fasta_polimixina,
                                           "others_outfile_suffix":
                                           "blastOthers",
                                           "poli_outfile_suffix": "blastPoli"}

            result3, result2 = run_blast_and_check_mutations(bacteria_dict)

    else:
        # mod 20.05.22
        printar_especies = f"{genero} {especie}"  # mod 10.05.22
        especie_mlst = 'Nao disponivel'  # mod 26.08.22
    # mod 20.05.22

    # print "$especie_mlst ";

    print('contaminacao...')
    # printar no arquivo final o nome da especie
    if float(contaminacao) <= 10.:
        # print OUT2 "$printar_especies\t";
        mongo_client.save('especie', printar_especies)
        # para o gal
        gal_file.write("Espécie identificada: "
                       f"{path.basename(printar_especies)}\n")
    else:
        # print OUT2 "$printar_especies\t";
        imprimir = (f"{maior_repeticao} {first_count} "
                    f"{segunda_repeticao} {second_count}")
        mongo_client.save('especie', imprimir)
        # para o gal
        gal_file.write(f"Espécie: CONTAMINAÇÃO {path.basename(imprimir)}\n")

    # else {
    #        	print OUT2 "$maior_repeticao $count_ordenado2{$maior_repeticao}
    # $segunda_repeticao $count_ordenado2{$segunda_repeticao}\t";
    # }

    ###########################################################################
    # Rodar ABRICATE
    # Para resistencia usando o ResFinder (porque so tem resistencia adquirida)
    abricate_out = f"{sample2}_outAbricateRes"
    abricate_line = (f"{caminho_abricate} --db resfinder "
                     f"{caminho1}/{sample2}/prokka"
                     f"/genome.ffn > {abricate_out} --threads {THREADS}")
    run_command_line(abricate_line)

    selected = get_abricate_result(abricate_out)
    select_imprimir = []

    # criar um @ para cada classe de antibioticos
    genes = []

    # print gal
    gal_file.write('Genes de resistência encontrados no WGS:\n')

    # MELISE: ADICIONEI NAS LINHAS DE genes.append "(allele confidence" +
    # lines_blast[10] + ")")
    for n_l in selected:
        # print "$n\n";
        # separar as colunas do arquivo em elementos de um array
        lines_blast = n_l.split("\t")
        # concatenar os resultado
        out_blast = (f"{lines_blast[5]} (ID:{lines_blast[10]} "
                     f"COV_Q:{lines_blast[9]} COV_DB:{lines_blast[6]})")
        select_imprimir.append(out_blast)
        # imprimir no arquivo do gal
        if re.match(r'.*(blaKPC|blaNDM|blaVIM|blaIMP|blaSPM|blaOXA-23|blaOXA-24|blaOXA-25|blaOXA-26|blaOXA-27|blaOXA-48|blaOXA-58|blaOXA-72|blaOXA-98|blaOXA-116|blaOXA-117|blaOXA-160|blaOXA-175|blaOXA-176|blaOXA-253).*', lines_blast[5], re.I):
            # print OUT2 "$lines_blast[5] (carbapenemase)\n";
            genes.append(lines_blast[5] + " (carbapenemase)" +
                         "(allele confidence " + lines_blast[10] + ")")
        elif re.match(r'.*(blaOXA-51|blaOXA-64|blaOXA-65|blaOXA-69|blaOXA-90|blaOXA-259|blaOXA-343).*', lines_blast[5], re.I):
            # print OUT2 "$lines_blast[5] (OXA-51-like carbapenemase)\n";
            genes.append(lines_blast[5] + " (OXA-51-like carbapenemase)" +
                         "(allele confidence " + lines_blast[10] + ")")
        elif re.match(r'.*(blaTEM|blaSHV|blaADC|blaCTX-M|blaGES|blaOXA-(?!23|24|25|26|27|48|58|72|98|116|117|160|175|176|253|488|486)).*', lines_blast[5], re.I):
            # print OUT2 "$lines_blast[5] (ESBL)\n";
            genes.append(lines_blast[5] + " (ESBL)" +
                         "(allele confidence " + lines_blast[10] + ")")
        elif re.match(r"(aac\(6\'\)-Ib-cr).*", lines_blast[5], re.I):
            # print OUT2 "$lines_blast[5] (resistencia a aminoglicosídeos e
            # fluoroquinolonas)\n";
            genes.append(f"{lines_blast[5]} (resistance to aminoglycosides and"
                         " fluoroquinolones) (allele confidence "
                         f"lines_blast[10])")
        elif re.match(r'(aph|aac|rmt|aad).*', lines_blast[5], re.I):
            # print OUT2 "$lines_blast[5] (resistencia a aminoglicosídeos)\n";
            genes.append(f"{lines_blast[5]} (resistance to aminoglycosides)" +
                         "(allele confidence " + lines_blast[10] + ")")
        elif re.match(r'(cat|cml|cmx|floR).*', lines_blast[5], re.I):
            # print OUT2 "$lines_blast[5] (resistencia ao cloranfenicol)\n";
            genes.append(f"{lines_blast[5]} (resistance to chloramphenicol)" +
                         "(allele confidence " + lines_blast[10] + ")")
        elif re.match(r'(qnr|oqx).*', lines_blast[5], re.I):
            # print OUT2 "$lines_blast[5] (resistencia a fluoroquinolonas)\n";
            genes.append(f"{lines_blast[5]} (resistance to fluoroquinolones)" +
                         "(allele confidence " + lines_blast[10] + ")")
        elif re.match(r'sul.*', lines_blast[5], re.I):
            # print OUT2 "$lines_blast[5] (resistencia a sulfonamidas)\n";
            genes.append(f"{lines_blast[5]} (resistance to sulfonamidas)" +
                         "(allele confidence " + lines_blast[10] + ")")
        elif re.match(r'dfrA.*', lines_blast[5], re.I):
            # print OUT2 "$lines_blast[5] (resistencia a trimetoprim)\n";
            genes.append(f"{lines_blast[5]} (resistance to trimetoprim)" +
                         "(allele confidence " + lines_blast[10] + ")")
        elif re.match(r'tet.*', lines_blast[5], re.I):
            # print OUT2 "$lines_blast[5] (resistencia a tetraciclina)\n";
            genes.append(f"{lines_blast[5]} (resistance to tetracycline)" +
                         "(allele confidence " + lines_blast[10] + ")")
        elif re.match(r'ere.*', lines_blast[5], re.I):
            # print OUT2 "$lines_blast[5] (resistencia a eritromicina)\n";
            genes.append(f"{lines_blast[5]} (resistance to eritromicina)" +
                         "(allele confidence " + lines_blast[10] + ")")
        elif re.match(r'erm.*', lines_blast[5], re.I):
            # print OUT2 "$lines_blast[5] (resistencia a lincosamidas,
            # macrolideos e estreptograminas)\n";
            genes.append(f"{lines_blast[5]} (resistance to lincosamides, "
                         "macrolides and streptogramins) (allele confidence "
                         f"{lines_blast[10]})")
        elif re.match(r'ARR.*', lines_blast[5], re.I):
            # print OUT2 "$lines_blast[5] (resistencia a rifampicina)\n";
            genes.append(f"{lines_blast[5]} (resistance to rifampicin)" +
                         "(allele confidence " + lines_blast[10] + ")")
        elif re.match(r'(mph|msr).*', lines_blast[5], re.I):
            # print OUT2 "$lines_blast[5] (resistencia a macrolideos)\n";
            genes.append(f"{lines_blast[5]} (resistance to macrolides)" +
                         "(allele confidence " + lines_blast[10] + ")")
        elif re.match(r'.*Van.*', lines_blast[5], re.I):
            # print OUT2 "$lines_blast[5] (resistencia a vancomicina)\n";
            genes.append(f"{lines_blast[5]} (resistance to vancomycin)" +
                         "(allele confidence " + lines_blast[10] + ")")
        elif re.match(r'.*lsa.*', lines_blast[5], re.I):
            # print OUT2 "$lines_blast[5] (resistencia a vancomicina)\n";
            genes.append(f"{lines_blast[5]} (resistance to clindamycin)" +
                         "(allele confidence " + lines_blast[10] + ")")
        elif re.match(r'.*mcr.*', lines_blast[5], re.I):
            # print OUT2 "$lines_blast[5] (resistencia a vancomicina)\n";
            # COLOCAR O NOME EM INGLES
            genes.append(f"{lines_blast[5]} (resistance to polymyxin)" +
                         "(allele confidence " + lines_blast[10] + ")")
        else:
            # print OUT2 "$lines_blast[5]\n";
            genes.append(f"{lines_blast[5]}" +
                         "(allele confidence " + lines_blast[10] + ")")

    # imprimir resultados com a classe do antimicrobiano

    mongo_client.save("gene", "<br>".join(genes))
    mongo_client.save("resfinder", "<br>".join(select_imprimir))

    ###########################################################################
    # Rodar abricate para VFDB (Virulence factor)
    abricate_out = f"{sample2}_outAbricateVFDB"
    abricate_line = (f"{caminho_abricate} --db vfdb "
                     f"{caminho1}/{sample2}/prokka/"
                     f"genome.ffn > {abricate_out} --threads {THREADS}")

    run_command_line(abricate_line)

    selected = get_abricate_result(abricate_out)
    select_imprimir = []

    # ler o array selected
    for n_l in selected:
        # print "$n\n";
        # separar as colunas do arquivo em elementos de um array
        lines_blast = n_l.split("\t")
        # print OUT2 "$lines_blast[5] ID:$lines_blast[10] COV_Q:$lines_blast[9]
        # COV_DB:$lines_blast[6]\|";
        out_blast = (f"{lines_blast[1]}: {lines_blast[5]} {lines_blast[13]} "
                     f"ID:{lines_blast[10]} COV_Q:{lines_blast[9]} "
                     f"COV_DB:{lines_blast[6]}| ")
        select_imprimir.append(out_blast)

    mongo_client.save('VFDB', "<br>".join(select_imprimir))
    pathlib.Path(abricate_out).unlink(missing_ok=True)

    ###########################################################################
    # Rodar abricate para PlasmidFinder
    abricate_out = f"{sample2}_outAbricatePlasmid"
    abricate_line = (f"{caminho_abricate} --db plasmidfinder "
                     f"{caminho1}/{sample2}/unicycler"
                     f"/assembly.fasta > {abricate_out} --threads {THREADS}")

    run_command_line(abricate_line)

    selected = get_abricate_result(abricate_out)
    select_imprimir = []

    # ler o array selected
    imprimir = 'Not found'
    if not selected:
        # print OUT2 "Nao encontrado\t";
        mongo_client.save('plasmid', imprimir)
        # para o gal
    else:
        imprimir = ""
        for n_l in selected:
            # print "$n\n";
            # separar as colunas do arquivo em elementos de um array
            lines_blast = n_l.split("\t")
            # print OUT2 "$lines_blast[5] ID:$lines_blast[10]
            # COV_Q:$lines_blast[9] COV_DB:$lines_blast[6]\|";
            out_blast = lines_blast[5] + 'ID' + ':' + lines_blast[10] + ' ' + \
                'COV_Q:' + lines_blast[9] + ' ' + \
                'COV_DB:' + lines_blast[6] + '|' + ' '
            select_imprimir.append(out_blast)
            imprimir += f"\n{lines_blast[5]}"
    gal_file.write(f"Plasmídeos encontrados:{path.basename(imprimir)}\n")

    mongo_client.save("plasmid", "<br>".join(select_imprimir))
    pathlib.Path(abricate_out).unlink(missing_ok=True)

    ###########################################################################

    print(f"Rodar o MLST {especie_mlst}")

    MLST_result = (f"{caminho1}/{sample2}"
                   "/unicycler/data.json")
    # se nao tem mlst disponivel, ai tem que avisar
    # mod 26-08-22
    if (especie_mlst == 'Nao disponivel') or (especie_mlst == ''):
        # print OUT2 "Nao disponivel\t";
        imprimir = 'Not available for this species'  # mod 26.08.22
        mongo_client.save('mlst', imprimir)
        # para o gal
        gal_file.write((f"Clone ST {path.basename(imprimir)} "
                        "(determinado por MLST)\n"))
    else:
        # mod 26-08-22
        docker_line = (f"docker run --rm -i -v {mlst_install}/"
                       "mlst_db:/database -v "
                       f"{caminho1}/{sample2}/"
                       "unicycler:/workdir mlst -i assembly.fasta -o . "
                       f"-s {especie_mlst}")
        # rodar o mlst
        run_command_line(docker_line)
    # mod

    ST = None
    print('ler o resultado do mlst')
    mlst_json = pathlib.Path(MLST_result)
    # will create file, if it exists will do nothing
    mlst_json.touch(exist_ok=True)

    with open(mlst_json, "r") as IN3:
        line = IN3.readline().rstrip("\n")  # single line file
        a = re.search(r'.*sequence_type":\s"(\d{1,4})".*', line, re.IGNORECASE)
        b = re.search(r'.*sequence_type":\s"(\d*!,\d*!)".*',
                      line, re.IGNORECASE)
        c = re.search(
            r'.*sequence_type":\s"(\d{1,4}\*)".*', line, re.IGNORECASE)
        for m in (a, b, c):
            if not m:
                continue
            ST = m.group(1)
            # print OUT2 "$ST\t";
            imprimir = ST
            mongo_client.save('mlst', imprimir)
            # para o gal
            gal_file.write((f"Clone ST {path.basename(imprimir)} "
                            "(determinado por MLST)\n"))
        m = re.search(r'nearest_sts":\s"((\d*,)*\d*)".*', line, re.IGNORECASE)
        if m:
            nearest_sts = m.group(1)
            if nearest_sts:
                # print OUT2 "Nearest $nearest_sts\t";
                imprimir = f"Nearest {nearest_sts}"
                mongo_client.save('mlst', imprimir)
                # para o gal
                gal_file.write(f"Clone ST {path.basename(imprimir)}"
                               " (determinado por MLST)\n")
        m = re.search(r'.*sequence_type":\s"(Unknown)".*', line, re.IGNORECASE)
        if m:
            ST = m.group(1)
            # print OUT2 "Unknown\t";
            imprimir = 'Unknown'
            mongo_client.save('mlst', imprimir)
            # para o gal
            gal_file.write((f"Clone ST {path.basename(imprimir)} "
                            "(determinado por MLST)\n"))

    mongo_client.save('mutacoes_poli', "<br>".join(result2))
    gal_file.write("Mutações polimixina: %s" % "<br>".join(result2))
    mongo_client.save('mutacoes_outras', "<br>".join(result3))

    ######################################################################
    print('rodar coverage')

    # figuring out if file is compressed or not
    catcmd = "cat"
    res = run_command_line(f"file {R1}")
    if res and str(res).find("gzip compressed") > -1:
        catcmd = "zcat"

    zcat = f"echo $({catcmd} {R1} | wc -l)/4 | bc"
    res_r1 = run_command_line(zcat)
    n_reads1 = res_r1.rstrip("\n")

    # o mesmo para o arquivo R2
    zcat2 = f"echo $({catcmd} {R2} | wc -l)/4 | bc"
    res_r2 = run_command_line(zcat2)
    n_reads2 = res_r2.rstrip("\n")

    soma_reads = (float(n_reads1) + float(n_reads2))

    # calcular tamanho medio das reads, vou usar só as R1 como base
    zcat3 = (f"{catcmd} {R1} | awk '{{if(NR%4==2) "
             f"{{count++; bases += length}} }} END{{print bases/count}}'")
    res_avg = run_command_line(zcat3)
    average_length2 = res_avg.rstrip("\n")

    gal_file.close()

    coverage = (float(average_length2) * soma_reads) / float(genome_size)

    mongo_client.save('coverage', str(coverage))
