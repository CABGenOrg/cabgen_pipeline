from argparse import ArgumentParser


def main():
    parser = ArgumentParser()
    parser.add_argument("--polymyxin_db", default="/home/melise/resis_poli",
                        type=str,
                        help=("Caminho do DB com as seqs de ptn de "
                              "resistência mutacional a outros antibióticos."))
