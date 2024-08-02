from argparse import ArgumentParser
from pipeline_website3 import pipeline


def main():
    parser = ArgumentParser(description="CABGen pipeline.")

    parser.add_argument("--polymyxin_db", default="/opt/resis_poli",
                        type=str,
                        help=("Caminho do DB com as seqs de ptn de "
                              "resistência mutacional a polimixina."))
    parser.add_argument("--outhers_db", default="/opt/outrasMut",
                        type=str,
                        help=("Caminho do DB com as seqs de ptn de "
                              "resistência mutacional a outros antibióticos."))
    parser.add_argument("--fastani_db", default="/opt/genomas_enterobacter",
                        type=str,
                        help=("Caminho para os dbs usados "
                              "pelo fastANI."))
    parser.add_argument("--output", default="/sgbmi_node/src/views/uploads/",
                        type=str,
                        help=("Caminho para guardar os resultados de montagem "
                              "e anotacao do genoma."))
    parser.add_argument("--sample", default="",
                        type=str,
                        help=("Numero sequencial da amostra de acordo com o "
                              "CabGen."))
    parser.add_argument("--relatorios_mongo", default="",
                        type=str,
                        help=("Arquivo que armazena os resultados para o mongo"
                              " CONFIRMAR"))

    parser.add_argument("--abricate", default="/opt/abricate/bin/abricate",
                        type=str,
                        help=("Caminho onde esta o bin do programa abricate."))
    parser.add_argument("--mlst",
                        default="/opt/identificar_clones/mlst_run.py",
                        type=str,
                        help=("Docker com o programa, banco de dados e "
                              "dependencias para realizar o MLST."))
    parser.add_argument("--kraken2", default="/opt/kraken2",
                        type=str,
                        help=("Arquivos e bin do Kraken "
                              "usado para identificar especie."))
    parser.add_argument("--unicycler", default="/usr/local/bin",
                        type=str,
                        help=("Arquivos e bin do Unicycler "
                              "usado para montagem do genoma."))
    parser.add_argument("--read1",
                        type=str,
                        help=("Caminho para as reads "
                              "R1 que o usuario carregou."))
    parser.add_argument("--read2",
                        type=str,
                        help=("Caminho para as reads "
                              "R2 que o usuario carregou."))
    parser.add_argument("--fastANI", default="/opt/FastANI/fastANI",
                        type=str,
                        help=("Caminho para o programa fastANI "
                              "usado para identificar subespecie."))

    args = parser.parse_args()
    pipeline(args)


if __name__ == "__main__":
    main()
