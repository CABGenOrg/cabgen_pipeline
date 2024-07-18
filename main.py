from argparse import ArgumentParser


def main():
    parser = ArgumentParser()
    parser.add_argument("--polymyxin_db", default="/opt/resis_poli",
                        type=str,
                        help=("Caminho do DB com as seqs de ptn de "
                              "resistência mutacional a polimixina."))

def main():
    parser = ArgumentParser()
    parser.add_argument("--outros_db", default="/opt/outrasMut",
                        type=str,
                        help=("Caminho do DB com as seqs de ptn de "
                              "resistência mutacional a outros antibióticos."))

def main():
    parser = ArgumentParser()
    parser.add_argument("--output", default="/sgbmi_node/src/views/uploads/",
                        type=str,
                        help=("Caminho para guardar os resultados de montagem "
                              "e anotacao do genoma."))

def main():
    parser = ArgumentParser()
    parser.add_argument("--sample", default="",
                        type=str,
                        help=("Numero sequencial da amostra de acordo com o CabGen."))

def main():
    parser = ArgumentParser()
    parser.add_argument("--abricate", default="/opt/abricate/bin/abricate",
                        type=str,
                        help=("Caminho onde esta o bin do programa abricate."))

def main():
    parser = ArgumentParser()
    parser.add_argument("--relatorios_mongo", default="",
                        type=str,
                        help=("Arquivo que armazena os resultados para o mongo CONFIRMAR"))

def main():
    parser = ArgumentParser()
    parser.add_argument("--mlst", default="/opt/identificar_clones/mlst_run.py",
                        type=str,
                        help=("Docker com o programa, banco de dados e dependencias "
                              "para realizar o MLST."))
    
def main():
    parser = ArgumentParser()
    parser.add_argument("--kraken2", default="/opt/kraken2/kraken2",
                        type=str,
                        help=("Arquivos e bin do Kraken "
                              "usado para identificar especie."))

def main():
    parser = ArgumentParser()
    parser.add_argument("--unicycler", default="/usr/local/bin/unicycler",
                        type=str,
                        help=("Arquivos e bin do Unicycler "
                              "usado para montagem do genoma."))

def main():
    parser = ArgumentParser()
    parser.add_argument("--reads", default="/sgbmi_node/src/views/uploads/",
                        type=str,
                        help=("Caminho para as reads "
                              "R1 e R2 que o usuario carregou."))
def main():
    parser = ArgumentParser()
    parser.add_argument("--fastANI", default="/opt/FastANI/fastANI",
                        type=str,
                        help=("Caminho para o programa fastANI "
                              "usado para identificar subespecie."))




