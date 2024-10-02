import re
from os import getenv, path, makedirs, listdir
from src.utils.handle_errors import fatal_error
from src.handle_database import MongoSaver
from src.types.SpeciesDict import SpeciesDict
from src.handle_programs import run_command_line
from src.utils.handle_folders import delete_folders_and_files
from src.handle_processing import count_kraken_words, build_species_data, \
    identify_bacteria_species


class CabgenPipeline:
    def __init__(self, sample: int, read1: str, read2: str, output: str):
        self.sample = int(sample)
        self.read1 = read1
        self.read2 = read2
        self.output = output
        self.threads = 16
        self.mongo_client = MongoSaver(self.sample)
        self.mongo_client.connect()

    def _check_params(self):
        try:
            print("Parameters:")
            params = ["sample", "read1", "read2", "output"]
            for param in params:
                param_value = getattr(self, param, None)
                if param_value:
                    print(f"{param}: {param_value}")
                else:
                    raise ValueError(f"The param {param} is empty.")
        except ValueError as e:
            fatal_error(f"Failed to check parameters.\n\n{e}")

    def _create_dirs(self):
        try:
            sample_directory = path.join(self.output, str(self.sample))
            unicycler_directory = path.join(sample_directory, "/unicycler")
            checkm_directory = path.join(sample_directory, "/checkM_bins")

            self.sample_directory = sample_directory
            self.unicycler_directory = unicycler_directory
            self.checkm_directory = checkm_directory

            dirs_to_create = [sample_directory, unicycler_directory,
                              checkm_directory]

            for dir in dirs_to_create:
                makedirs(dir, exist_ok=True)
        except Exception as e:
            raise Exception(f"Can't create sample directories.\n\n{e}")

    def _load_programs(self):
        self.abricate = getenv("ABRICATE_PATH") or ""
        self.mlst = getenv("MLST_PATH") or ""
        self.polimyxin_db = getenv("POLIMYXIN_DB_PATH") or ""
        self.outhers_db = getenv("OUTHERS_DB_PATH") or ""
        self.kraken2 = getenv("KRAKEN2_PATH") or ""
        self.kraken_db = getenv("KRAKEN_DB_PATH") or ""
        self.unicycler = getenv("UNICYCLER_PATH") or ""
        self.fastani = getenv("FASTANI_PATH") or ""
        self.fastani_db = getenv("FASTANI_DB_PATH") or ""
        self.loaded_programs = ["abricate", "mlst",
                                "polimyxin_db", "outhers_db",
                                "kraken2", "kraken_db", "unicycler",
                                "fastani", "fastani_db"]

    def _check_programs(self):
        try:
            for program in self.loaded_programs:
                program_value = getattr(self, program, None)
                if program_value:
                    print(f"{program}: {program_value}")
                else:
                    raise ValueError(
                        f"The {program} is not defined. Check the env file.")
        except ValueError as e:
            fatal_error(f"Failed to check programs.\n\n{e}")

    def _run_unicycler(self):
        try:
            print("Running Unicycler")
            unicycler_line = (f"{self.unicycler} -1 {self.read1} "
                              f"-2 {self.read2} "
                              f"-o {self.unicycler_directory} "
                              "--min_fasta_length 500 --mode conservative "
                              f"-t {self.threads}")
            program_output = run_command_line(unicycler_line)
            if program_output:
                print(program_output)
        except Exception as e:
            fatal_error(f"Failed to run Unicycler.\n\n{e}")

    def _run_prokka(self):
        try:
            self.assembly_path = path.join(f"{self.output}/{self.sample}",
                                           "/unicycler/assembly.fasta")
            print("Run Prokka")
            prokka_line = (f"prokka --outdir {self.output}/"
                           f"{self.sample}/prokka --prefix genome "
                           f"{self.assembly_path} --force "
                           f"--cpus {self.threads}")
            run_command_line(prokka_line)
        except Exception as e:
            fatal_error(f"Failed to run Prokka.\n\n{e}")

    def _run_checkm(self):
        try:
            print("Run CheckM")
            checkM_line = ("checkm lineage_wf -x fasta "
                           f"{self.unicycler_directory} "
                           f"{self.checkm_directory} --threads {self.threads} "
                           f"--pplacer_threads {self.threads}")
            checkM_qa_line = ("checkm qa -o 2 "
                              f"-f {self.checkm_directory}/{self.sample}"
                              "_resultados "
                              f"--tab_table {self.checkm_directory}/lineage.ms"
                              f" {self.checkm_directory} "
                              f"--threads {self.threads}")

            run_command_line(checkM_line)
            run_command_line(checkM_qa_line)

            files_to_delete = [path.join(self.checkm_directory, file) for file
                               in listdir(self.checkm_directory)
                               if "resultados" not in file]
            delete_folders_and_files(files_to_delete)
        except Exception as e:
            fatal_error(f"Failed to run checkM.\n\n{e}")

    def _process_checkm_result(self):
        try:
            print("Saving CheckM result to MongoDB")
            checkM_results_file = (f"{self.checkm_directory}/"
                                   f"{self.sample}_resultados")
            with open(checkM_results_file) as inp:
                next(inp)
                for row in inp:
                    row = row.rstrip("\n")
                    lines = row.split("\t")
                    self.genome_size = lines[8] or 1
                    self.mongo_client.save('checkm_1', lines[5])
                    self.mongo_client.save('checkm_2', lines[6])
                    self.mongo_client.save('checkm_3', lines[8])
                    self.mongo_client.save('checkm_4', lines[11])
                    self.contamination = lines[6] or 0
                    self.mongo_client.save('sample', str(self.sample))
        except Exception as e:
            fatal_error(f"Failed to process checkM result.\n\n{e}")

    def _run_kraken2(self):
        try:
            print("Run Kraken2")
            kraken_line = (f"{self.kraken2} --db {self.kraken_db} "
                           f"--use-names --paired {self.read1} {self.read2} "
                           f"--output {self.sample_directory}/out_kraken "
                           f"--threads {self.threads}")
            run_command_line(kraken_line)
        except Exception as e:
            fatal_error(f"Failed to run kraken2.\n\n{e}")

    def _process_kraken2_result(self):
        try:
            print(f"Splitting output into {self.threads} equal files.")
            preffix = "krk"
            splitter_line = (f"split --numeric-suffixes=1 -n l/{self.threads} "
                             f"{self.sample_directory}/out_kraken {preffix}")
            run_command_line(splitter_line)

            most_common, second_most_common, \
                first_count, second_count = count_kraken_words(
                    f"{self.sample_directory}/out_kraken")

            self.most_common = most_common
            self.second_most_common = second_most_common
            self.first_count = first_count
            self.second_count = second_count
        except Exception as e:
            fatal_error(f"Failed to process kraken2 result.\n\n{e}")

    def _process_species(self):
        try:
            if (re.findall(re.compile(r'\w+\s\w.*', re.I), self.most_common)):
                check_especies = self.most_common.strip()
                split_especies = check_especies.split(" ")
                genus = split_especies[0]
                species = split_especies[1]
            else:
                genus = self.most_common
                species = ""

            species_final_result = f"{genus}{species}".lower()
            print(f"Final result of the species: {species_final_result}.")
            species_info: SpeciesDict = {"species": species_final_result,
                                         "assembly": self.assembly_path,
                                         "sample": self.sample,
                                         "poli_db_path": self.polimyxin_db,
                                         "others_db_path": self.outhers_db,
                                         "fastani_db_path": self.fastani_db}

            blast_result, display_name, mlst_species = \
                identify_bacteria_species(species_info)

            _, fastani_species = build_species_data(species_info)
            if not blast_result and species in fastani_species.keys() \
                    or "enterobacter" in species \
                    or "acinetobacter" in species:
                display_name, mlst_species = \
                    self._run_fastani(species_final_result)  # type: ignore

            if not display_name:
                display_name = f"{species}"
                mlst_species = "Não disponível"

            self.display_name = display_name
            self.mlst_species = mlst_species
        except Exception as e:
            fatal_error(f"Failed to process species.\n\n{e}")

    def _run_fastani(self, species_info: SpeciesDict):
        try:
            species = species_info.get("species")
            _, fastani_species = build_species_data(species_info)

            fastani_group = "enterobacter_species" if "enterobacter" in \
                species else "acinetobacter_species"
            if species == "klebsiellapneumoniae":
                desired_species_data = fastani_species.get(
                    "klebsiellapneumoniae", {})
            else:
                desired_species_data = fastani_species.get(fastani_group, {})

            mlst_species = desired_species_data.get("mlst")
            fastani_list = desired_species_data.get("fastani_list")
            fastani_output = (f"{self.sample_directory}/{self.sample}_out-"
                              "fastANI")

            fastani_line = (
                f"{self.fastani} -q {self.assembly_path} --rl {fastani_list} "
                f"-o {fastani_output} --threads {self.threads}"
            )
            run_command_line(fastani_line)

            with open(fastani_output, "r") as file:
                species_name = file.readline().strip().split(
                    "\t")[1].split("/")[-1].split(".")[0]

            return species_name, mlst_species
        except Exception as e:
            fatal_error(f"Failed to run FastAni.\n\n{e}")

    def _save_species_result(self):
        try:
            if float(self.contamination) <= 10.:
                self.mongo_client.save('especie', self.display_name)
            else:
                first_repetition = self.most_common
                first_count = self.first_count
                second_repetition = self.second_most_common
                second_count = self.second_count

                species_info = (f"{first_repetition} {first_count} "
                                f"{second_repetition} {second_count}")
                self.mongo_client.save('especie', species_info)
        except Exception as e:
            print(f"Failed to save species result.\n\n{e}")

    def run(self):
        pass
