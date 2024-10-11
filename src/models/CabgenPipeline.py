import re
import sys
from time import time
from dotenv import load_dotenv
from logging import Logger
from os import getenv, path, makedirs, listdir
from src.models.MongoHandler import MongoHandler
from src.types.SpeciesDict import SpeciesDict
from src.utils.handle_programs import run_command_line
from src.utils.handle_folders import delete_folders_and_files
from src.utils.handle_processing import count_kraken_words, \
    build_species_data, identify_bacteria_species, get_abricate_result, \
    process_resfinder, process_vfdb, process_plasmidfinder, format_time

load_dotenv()
uploaded_sequences_path = getenv("UPLOADED_SEQUENCES_PATH") or ""


class CabgenPipeline:
    def __init__(self, sample: int, read1: str, read2: str, output: str,
                 logger: Logger):
        self.sample = int(sample)
        self.read1 = path.join(uploaded_sequences_path, read1)
        self.read2 = path.join(uploaded_sequences_path, read2)
        self.output = output
        self.threads = 3 if not getenv("THREADS") \
            else int(getenv("THREADS"))  # type: ignore
        self.mongo_client = MongoHandler()
        self.logger = logger

    def _check_params(self):
        try:
            self.logger.info("Parameters:")
            params = ["sample", "read1", "read2", "output"]
            for param in params:
                param_value = getattr(self, param, None)
                if param_value:
                    self.logger.info(f"{param}: {param_value}")
                else:
                    self.logger.error(f"The param {param} is empty.")
                    raise ValueError
        except ValueError as e:
            self.logger.error(f"Failed to check parameters.\n\n{e}")
            sys.exit(1)

    def _create_dirs(self):
        try:
            sample_directory = self.output
            unicycler_directory = path.join(sample_directory, "unicycler")
            checkm_directory = path.join(sample_directory, "checkM_bins")

            self.sample_directory = sample_directory
            self.unicycler_directory = unicycler_directory
            self.checkm_directory = checkm_directory

            dirs_to_create = [sample_directory, unicycler_directory,
                              checkm_directory]

            for dir in dirs_to_create:
                makedirs(dir, exist_ok=True)
        except Exception as e:
            self.logger.error(f"Can't create sample directories.\n\n{e}")
            sys.exit(1)

    def _load_programs(self):
        self.fastqc = getenv("FASTQC") or ""
        self.abricate = getenv("ABRICATE_PATH") or ""
        self.mlst = getenv("MLST_PATH") or ""
        self.polimyxin_db = getenv("POLIMYXIN_DB_PATH") or ""
        self.outhers_db = getenv("OUTHERS_DB_PATH") or ""
        self.kraken2 = getenv("KRAKEN2_PATH") or ""
        self.kraken_db = getenv("KRAKEN_DB_PATH") or ""
        self.unicycler = getenv("UNICYCLER_PATH") or ""
        self.fastani = getenv("FASTANI_PATH") or ""
        self.fastani_db = getenv("FASTANI_DB_PATH") or ""
        self.spades = getenv("SPADES_PATH") or ""
        self.loaded_programs = ["abricate", "mlst",
                                "polimyxin_db", "outhers_db",
                                "kraken2", "kraken_db", "unicycler",
                                "fastani", "fastani_db"]

    def _check_programs(self):
        try:
            for program in self.loaded_programs:
                program_value = getattr(self, program, None)
                if program_value:
                    self.logger.info(f"{program}: {program_value}")
                else:
                    self.logger.error(
                        f"The {program} is not defined. Check the env file.")
                    raise ValueError
        except ValueError as e:
            self.logger.error(f"Failed to check programs.\n\n{e}")
            sys.exit(1)

    def _run_fastqc(self):
        try:
            self.logger.info("Running FastQC")
            fastqc_output_path = getenv("FASTQC_OUTPUT_PATH") or ""

            if not fastqc_output_path:
                self.logger.error("FastQC output path is not defined in .env.")
                raise ValueError

            fastqc_line = (f"{self.fastqc} --quiet {self.read1} {self.read2} "
                           f"--outdir {fastqc_output_path}")
            run_command_line(fastqc_line)
        except Exception as e:
            self.logger.error(f"Failed to run FASTQC.\n\n{e}")
            sys.exit(1)

    def _run_unicycler(self):
        try:
            self.logger.info("Running Unicycler")
            if self.spades:
                unicycler_line = (f"{self.unicycler} -1 {self.read1} "
                                  f"-2 {self.read2} "
                                  f"-o {self.unicycler_directory} "
                                  "--min_fasta_length 500 --mode conservative "
                                  f"-t {self.threads} "
                                  f"--spades_path {self.spades}")
            else:
                unicycler_line = (f"{self.unicycler} -1 {self.read1} "
                                  f"-2 {self.read2} "
                                  f"-o {self.unicycler_directory} "
                                  "--min_fasta_length 500 --mode conservative "
                                  f"-t {self.threads}")
            program_output = run_command_line(unicycler_line)
            if program_output:
                self.logger.info(program_output)
        except Exception as e:
            self.logger.error(f"Failed to run Unicycler.\n\n{e}")
            sys.exit(1)

    def _run_prokka(self):
        try:
            self.assembly_path = path.join(f"{self.output}/{self.sample}",
                                           "/unicycler/assembly.fasta")
            self.logger.info("Run Prokka")
            prokka_line = (f"prokka --outdir {self.output}/"
                           f"{self.sample}/prokka --prefix genome "
                           f"{self.assembly_path} --force "
                           f"--cpus {self.threads}")
            run_command_line(prokka_line)
        except Exception as e:
            self.logger.error(f"Failed to run Prokka.\n\n{e}")
            sys.exit(1)

    def _run_checkm(self):
        try:
            self.logger.info("Run CheckM")
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
            self.logger.error(f"Failed to run checkM.\n\n{e}")
            sys.exit(1)

    def _process_checkm_result(self):
        try:
            self.logger.info("Saving CheckM result to MongoDB")
            checkM_results_file = (f"{self.checkm_directory}/"
                                   f"{self.sample}_resultados")
            with open(checkM_results_file) as inp:
                next(inp)
                for row in inp:
                    row = row.rstrip("\n")
                    lines = row.split("\t")
                    self.genome_size = lines[8] or 1

                    query = {"sequenciaId": self.sample}
                    self.mongo_client.save(
                        "relatorios", query, {"checkm_1": lines[5]})
                    self.mongo_client.save(
                        "relatorios", query, {"checkm_2": lines[6]})
                    self.mongo_client.save(
                        "relatorios", query, {"checkm_3": lines[8]})
                    self.mongo_client.save(
                        "relatorios", query, {"checkm_4": lines[11]})
                    self.contamination = lines[6] or 0
                    self.mongo_client.save(
                        "relatorios", query, {"sample": str(self.sample)})
        except Exception as e:
            self.logger.error(f"Failed to process checkM result.\n\n{e}")
            sys.exit(1)

    def _run_kraken2(self):
        try:
            self.logger.info("Run Kraken2")
            kraken_line = (f"{self.kraken2} --db {self.kraken_db} "
                           f"--use-names --paired {self.read1} {self.read2} "
                           f"--output {self.sample_directory}/out_kraken "
                           f"--threads {self.threads}")
            run_command_line(kraken_line)
        except Exception as e:
            self.logger.error(f"Failed to run kraken2.\n\n{e}")
            sys.exit(1)

    def _process_kraken2_result(self):
        try:
            self.logger.info(f"Splitting output into {self.threads} equal "
                             "files.")
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
            self.logger.error(f"Failed to process kraken2 result.\n\n{e}")
            sys.exit(1)

    def _process_species(self):
        try:
            if (re.findall(re.compile(r"\w+\s\w.*", re.I), self.most_common)):
                check_especies = self.most_common.strip()
                split_especies = check_especies.split(" ")
                genus = split_especies[0]
                species = split_especies[1]
            else:
                genus = self.most_common
                species = ""

            species_final_result = f"{genus}{species}".lower()
            self.logger.info(f"Final result of the species: "
                             f"{species_final_result}.")
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

            if blast_result:
                self.others_mutations_result = blast_result[0]
                self.poli_mutations_result = blast_result[1]
            else:
                self.others_mutations_result = []
                self.poli_mutations_result = []

            if not display_name:
                display_name = f"{species}"
                mlst_species = "Não disponível"

            self.display_name = display_name
            self.mlst_species = mlst_species
        except Exception as e:
            self.logger.error(f"Failed to process species.\n\n{e}")
            sys.exit(1)

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
            self.logger.error(f"Failed to run FastAni.\n\n{e}")
            sys.exit(1)

    def _save_species_result(self):
        try:
            query = {"sequenciaId": self.sample}
            if float(self.contamination) <= 10.:
                self.mongo_client.save(
                    "relatorios", query, {"especie": self.display_name})
            else:
                first_repetition = self.most_common
                first_count = self.first_count
                second_repetition = self.second_most_common
                second_count = self.second_count

                species_info = (f"{first_repetition} {first_count} "
                                f"{second_repetition} {second_count}")
                self.mongo_client.save(
                    "relatorios", query, {"especie": species_info})
        except Exception as e:
            self.logger.error(f"Failed to save species result.\n\n{e}")

    def _run_abricate(self, db: str):
        try:
            if db.lower() == "resfinder":
                self.logger.info("Run Abricate - ResFinder ")
                self.abricate_res_out = path.join(
                    self.sample_directory, f"{self.sample}_outAbricateRes")
            elif db.lower() == "vfdb":
                self.logger.info("Run Abricate - VFDB ")
                self.abricate_vfdb_out = path.join(
                    self.sample_directory, f"{self.sample}_outAbricateVFDB")
            elif db.lower() == "plasmidfinder":
                self.logger.info("Run Abricate - PlasmidFinder ")
                self.abricate_plasmid_out = path.join(
                    self.sample_directory, f"{self.sample}_outAbricatePlasmid")
            else:
                self.logger.error("Abricate database invalid.")

            abricate_line = (f"{self.abricate} --db {db} "
                             f"{self.sample_directory}/prokka/genome.ffn "
                             f"> {self.abricate_res_out} "
                             f"--threads {self.threads}")
            self.logger.info(f"{abricate_line}")
            run_command_line(abricate_line)
        except Exception as e:
            self.logger.error(
                f"Failed to run Abricate with resfinder DB.\n\n{e}")

    def _process_resfinder_result(self):
        try:
            abricate_result = get_abricate_result(
                self.abricate_res_out)
            gene_results, blast_out_results = process_resfinder(
                abricate_result)

            query = {"sequenciaId": self.sample}
            if not gene_results:
                self.mongo_client.save(
                    "relatorios", query, {"gene": "Not found"})
            else:
                self.mongo_client.save("relatorios", query,
                                       {"gene": "<br>".join(gene_results)})
                self.mongo_client.save("relatorios", query,
                                       {"resfinder":
                                        "<br>".join(blast_out_results)})
        except Exception as e:
            self.logger.error(
                f"Failed to process Abricate resfinder result.\n\n{e}")

    def _process_vfdb_result(self):
        try:
            abricate_result = get_abricate_result(
                self.abricate_vfdb_out)

            query = {"sequenciaId": self.sample}
            if abricate_result:
                blast_out_results = process_vfdb(abricate_result)
                self.mongo_client.save(
                    "relatorios", query,
                    {"VFDB": "<br>".join(blast_out_results)})
            else:
                self.mongo_client.save(
                    "relatorios", query, {"VFDB": "Not Found"})
        except Exception as e:
            self.logger.error(
                f"Failed to process Abricate VFDB result.\n\n{e}")

    def _process_plasmid_result(self):
        try:
            abricate_result = get_abricate_result(
                self.abricate_plasmid_out)

            query = {"sequenciaId": self.sample}
            if abricate_result:
                blast_out_results = process_plasmidfinder(abricate_result)
                self.mongo_client.save("relatorios", query, {"plasmid":
                                       "<br>".join(blast_out_results)})
            else:
                self.mongo_client.save(
                    "relatorios", query, {"plasmid": "Not Found"})
        except Exception as e:
            self.logger.error(
                f"Failed to process Abricate PlasmidFinder result.\n\n{e}")

    def _process_abricate_result(self, db: str):
        try:
            if db.lower() == "resfinder":
                self._process_resfinder_result()
            elif db.lower() == "vfdb":
                self._process_vfdb_result()
            elif db.lower() == "plasmidfinder":
                self._process_plasmid_result()
            else:
                self.logger.error("Abricate database invalid.")
                raise ValueError
        except Exception as e:
            self.logger.error(f"Failed to process Abricate result.\n\n{e}")

    def _run_mlst(self):
        try:
            self.logger.info(f"Run MLST for {self.mlst_species}")
            self.mlst_result_path = path.join(self.sample_directory,
                                              "mlst.csv")
            mlst_line = (f"{self.mlst} --threads {self.threads} "
                         "--exclude abaumannii --csv "
                         f"{self.assembly_path} > {self.mlst_result_path}")

            run_command_line(mlst_line)
        except Exception as e:
            self.logger.error(f"Failed to run MLST.\n\n{e}")

    def _process_mlst(self):
        try:
            with open(self.mlst_result_path, "r") as inp:
                line = inp.readline().rstrip("\n")
                out_mlst = line.split(",")
                scheme_mlst = out_mlst[1]
                st = out_mlst[2]

                query = {"sequenciaId": self.sample}
                if st != "-":
                    result = st
                    self.mongo_client.save(
                        "relatorios", query, {"mlst": result})
                    self.logger.info(f"Scheme used {scheme_mlst}")
                elif st == "-":
                    result = "New ST"
                    self.mongo_client.save(
                        "relatorios", query, {"mlst": result})
                    self.logger.info(f"Scheme used {scheme_mlst}")
                elif scheme_mlst == "-":
                    result = "Not available for this specie"
                    self.mongo_client.save(
                        "relatorios", query, {"mlst": result})
                    self.logger.info(f"Scheme used {scheme_mlst}")

            self.mongo_client.save("relatorios", query,
                                   {"mutacoes_poli": "<br>".join(
                                       self.poli_mutations_result)})
            self.mongo_client.save("relatorios", query,
                                   {"mutacoes_outras": "<br>".join(
                                       self.others_mutations_result)})
        except Exception as e:
            self.logger.error(f"Failed to process MLST result.\n\n{e}")

    def _run_coverage(self):
        try:
            self.logger.info("Run coverage")
            catcmd = "cat"
            res = run_command_line(f"file {self.read1}")
            if res and (str(res).find("gzip compressed") > -1 or
                        str(res).find("gzip compatible") > -1):
                catcmd = "zcat"

            zcat = f"echo $({catcmd} {self.read1} | wc -l)/4 | bc"
            res_r1 = run_command_line(zcat)
            n_reads1 = res_r1.rstrip("\n")

            zcat2 = f"echo $({catcmd} {self.read2} | wc -l)/4 | bc"
            res_r2 = run_command_line(zcat2)
            n_reads2 = res_r2.rstrip("\n")

            reads_sum = (float(n_reads1) + float(n_reads2))

            zcat3 = (f"{catcmd} {self.read1} | awk '{{if(NR%4==2) "
                     "{{count++; bases += length}} }} "
                     "END{{print bases/count}}'")
            res_avg = run_command_line(zcat3)
            average_length = res_avg.rstrip("\n")

            pre_coverage = (float(average_length) *
                            reads_sum) / float(self.genome_size)

            coverage = round(pre_coverage, 2)
            query = {"sequenciaId": self.sample}
            self.mongo_client.save(
                "relatorios", query, {"coverage": str(coverage)})
        except Exception as e:
            self.logger.error(f"Failed to run coverage.\n\n{e}")
            sys.exit(1)

    def _run_only_fastqc(self):
        try:
            self._run_fastqc()
            query = {"_id": self.sample}
            bson = {"$currentDate": {"ultimaActualizacao": True},
                    "$set": {"estado": "QUAL", "ultimaTarefa": ""}}
            self.mongo_client.save("sequencias", query, bson)
        except Exception as e:
            self.logger.error(
                f"Failed to run CABGen only FastQC pipeline.\n\n{e}")
            sys.exit(1)

    def _run_only_genomic(self):
        try:
            self._run_unicycler()
            self._run_prokka()
            self._run_checkm()
            self._process_checkm_result()
            self._run_kraken2()
            self._process_kraken2_result()
            self._process_species()
            self._save_species_result()
            abricate_dbs = ["resfinder", "vfdb", "plasmidfinder"]
            for db in abricate_dbs:
                self._run_abricate(db)
                self._process_abricate_result(db)
            self._run_mlst()
            self._process_mlst()
            self._run_coverage()

            query = {"_id": self.sample}
            bson = {"$currentDate": {"ultimaActualizacao": True},
                    "$set": {"estado": "ENSA", "ultimaTarefa": "",
                             "arquivofasta": f"{self.sample}.fasta"}}
            self.mongo_client.save("sequencias", query, bson)
        except Exception as e:
            self.logger.error(
                f"Failed to run CABGen only genomic pipeline.\n\n{e}")
            sys.exit(1)

    def _run_complete(self):
        try:
            self._run_fastqc()
            self._run_only_genomic()
        except Exception as e:
            self.logger.error(
                f"Failed to run CABGen complete pipeline.\n\n{e}")
            sys.exit(1)

    def run(self, only_fastqc=False, only_genomic=False, complete=False):
        try:
            start_time = time()
            # Starting the pipeline dependencies
            self._check_params()
            self._create_dirs()
            self._load_programs()
            self._check_programs()

            # Starting pipeline run
            if only_fastqc:
                self._run_only_fastqc()
            elif only_genomic:
                self._run_only_genomic()
            elif complete:
                self._run_complete()
            else:
                self.logger.error("Invalid pipeline choice!")
                raise ValueError

            self.mongo_client.close()
            runtime = format_time(time() - start_time)
            self.logger.info(f"Total runtime: {runtime}")
        except Exception as e:
            self.mongo_client.close()
            self.logger.error(f"Failed to run CABGen pipeline.\n\n{e}")
            sys.exit(1)
