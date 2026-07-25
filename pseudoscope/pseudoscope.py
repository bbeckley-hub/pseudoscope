#!/usr/bin/env python3
"""
PseudoScope - Unified Orchestrator for P. aeruginosa Analysis
Author: Brown Beckley <brownbeckley94@gmail.com>
GitHub: bbeckley-hub
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-07-25
Version: 1.2.0 (HPC‑friendly with temporary directory isolation + signal handlers)
"""

import os
import sys
import glob
import argparse
import subprocess
import shutil
import random
import tempfile
import logging
import traceback
import signal
import atexit
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

__version__ = "1.2.0"


class Color:
    RESET = '\033[0m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    BRIGHT_BLACK = '\033[90m'
    BRIGHT_RED = '\033[91m'
    BRIGHT_GREEN = '\033[92m'
    BRIGHT_YELLOW = '\033[93m'
    BRIGHT_BLUE = '\033[94m'
    BRIGHT_MAGENTA = '\033[95m'
    BRIGHT_CYAN = '\033[96m'
    BRIGHT_WHITE = '\033[97m'


class PseudoScopeOrchestrator:
    """Main orchestrator: runs all analyses inside temporary directories with signal handling."""

    def __init__(self, quiet: bool = False, keep_temp: bool = False):
        self.base_dir = Path(__file__).parent
        self.quiet = quiet
        self.keep_temp = keep_temp
        self.logger = None
        self.user_output_dir = None
        self.temp_dirs = []          # track temp dirs for cleanup
        self.running = True

        self.setup_colors()
        self.quotes = self._get_scientific_quotes()
        self.quote_colors = [
            Color.BRIGHT_CYAN, Color.BRIGHT_GREEN, Color.BRIGHT_YELLOW,
            Color.BRIGHT_MAGENTA, Color.BRIGHT_BLUE, Color.BRIGHT_RED,
            Color.CYAN, Color.GREEN, Color.YELLOW, Color.MAGENTA
        ]

        self.output_dirs = {
            'qc': 'fasta_qc_results',
            'mlst': 'mlst_results',
            'past': 'past_results',
            'abricate': 'pseudo_abricate_results',
            'amr': 'pseudo_amrfinder_results',
            'summary': 'GENIUS_PSEUDOMONAS_ULTIMATE_GENE_CENTRIC_REPORTS',
            'viz': 'GENIUS_PSEUDOMONAS_VISUAL_DASHBOARD',
            'sample_centric': 'GENIUS_PSEUDOMONAS_SAMPLE_CENTRIC_REPORTS'
        }

        # Files required for gene‑centric summary
        self.summary_html_files = {
            'mlst_summary.html': 'mlst_results/mlst_summary.html',
            'past_summary.html': 'past_results/past_summary.html',
            'PseudoScope_FASTA_QC_summary.html': 'fasta_qc_results/PseudoScope_FASTA_QC_summary.html',
            'pseudo_amrfinder_summary_report.html': 'pseudo_amrfinder_results/pseudo_amrfinder_summary_report.html',
            'pseudo_ncbi_summary_report.html': 'pseudo_abricate_results/pseudo_ncbi_summary_report.html',
            'pseudo_card_summary_report.html': 'pseudo_abricate_results/pseudo_card_summary_report.html',
            'pseudo_resfinder_summary_report.html': 'pseudo_abricate_results/pseudo_resfinder_summary_report.html',
            'pseudo_vfdb_summary_report.html': 'pseudo_abricate_results/pseudo_vfdb_summary_report.html',
            'pseudo_argannot_summary_report.html': 'pseudo_abricate_results/pseudo_argannot_summary_report.html',
            'pseudo_plasmidfinder_summary_report.html': 'pseudo_abricate_results/pseudo_plasmidfinder_summary_report.html',
            'pseudo_megares_summary_report.html': 'pseudo_abricate_results/pseudo_megares_summary_report.html',
            'pseudo_ecoh_summary_report.html': 'pseudo_abricate_results/pseudo_ecoh_summary_report.html',
            'pseudo_bacmet2_summary_report.html': 'pseudo_abricate_results/pseudo_bacmet2_summary_report.html',
            'pseudo_ecoli_vf_summary_report.html': 'pseudo_abricate_results/pseudo_ecoli_vf_summary_report.html',
            'mutation_summary.html': 'pseudo_amrfinder_results/mutation_summary.html',
        }

        self.warnings = []
        self.errors = []

        # Register signal handlers
        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)
        atexit.register(self._cleanup_temp_dirs)

    def _signal_handler(self, signum, frame):
        self.print_warning(f"Received signal {signum}. Cleaning up and exiting...")
        self.running = False
        self._cleanup_temp_dirs()
        sys.exit(1)

    def _cleanup_temp_dirs(self):
        if self.keep_temp:
            self.print_info("Keeping temporary directories (--keep-temp set).")
            return
        for d in self.temp_dirs:
            if d and Path(d).exists():
                shutil.rmtree(d, ignore_errors=True)
                self.logger.info(f"Removed temporary directory: {d}")
        self.temp_dirs = []

    def setup_colors(self) -> None:
        self.color_info = Color.CYAN
        self.color_success = Color.BRIGHT_GREEN
        self.color_warning = Color.BRIGHT_YELLOW
        self.color_error = Color.BRIGHT_RED
        self.color_reset = Color.RESET

    def setup_logging(self, output_dir: Path) -> None:
        log_file = output_dir / "pseudoscope_run.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[logging.FileHandler(log_file, mode='w')]
        )
        self.logger = logging.getLogger("PseudoScope")
        self.logger.info(f"Logging to {log_file}")
        self.user_output_dir = output_dir

    def _print(self, message: str, color: str = Color.RESET, bold: bool = False, force: bool = False) -> None:
        if self.quiet and not force:
            return
        style = Color.BOLD if bold else ''
        print(f"{style}{color}{message}{Color.RESET}")

    def print_header(self, title: str, subtitle: str = "") -> None:
        if not self.quiet:
            print()
            print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'='*80}{Color.RESET}")
            print(f"{Color.BOLD}{Color.BRIGHT_CYAN}{' ' * 20}{title}{Color.RESET}")
            if subtitle:
                print(f"{Color.DIM}{Color.WHITE}{' ' * 22}{subtitle}{Color.RESET}")
            print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'='*80}{Color.RESET}")
            print()
        if self.logger:
            self.logger.info(f"=== {title} ===")
            if subtitle:
                self.logger.info(f"   {subtitle}")

    def print_info(self, message: str) -> None:
        if not self.quiet:
            print(f"{self.color_info}[INFO]{Color.RESET} {message}")
        if self.logger:
            self.logger.info(message)

    def print_success(self, message: str) -> None:
        if not self.quiet:
            print(f"{self.color_success}✓{Color.RESET} {message}")
        if self.logger:
            self.logger.info(f"SUCCESS: {message}")

    def print_warning(self, message: str) -> None:
        if not self.quiet:
            print(f"{self.color_warning}⚠️{Color.RESET} {message}")
        if self.logger:
            self.logger.warning(message)
        self.warnings.append(message)

    def print_error(self, message: str) -> None:
        print(f"{self.color_error}✗{Color.RESET} {message}")
        if self.logger:
            self.logger.error(message)
        self.errors.append(message)

    def _get_scientific_quotes(self) -> List[Dict[str, str]]:
        return [
            {"quote": "Pseudomonas aeruginosa: the master of adaptation.", "author": "Unknown"},
            {"quote": "The important thing is not to stop questioning.", "author": "Albert Einstein"},
            {"quote": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"},
            {"quote": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"},
            {"quote": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie"},
            {"quote": "Research is what I'm doing when I don't know what I'm doing.", "author": "Wernher von Braun"},
            {"quote": "The science of today is the technology of tomorrow.", "author": "Edward Teller"},
            {"quote": "Science knows no country, because knowledge belongs to humanity.", "author": "Louis Pasteur"},
            {"quote": "In science, there are no shortcuts to truth.", "author": "Karl Popper"},
            {"quote": "Every great advance in science has issued from a new audacity of imagination.", "author": "John Dewey"},
            {"quote": "The most beautiful thing we can experience is the mysterious.", "author": "Albert Einstein"},
            {"quote": "Microbes are the dark matter of the biological world.", "author": "Jack Gilbert"},
            {"quote": "Antibiotic resistance is a ticking time bomb.", "author": "Sally Davies"},
            {"quote": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson"},
            {"quote": "We are all connected; to each other, biologically; to the earth, chemically; to the rest of the universe, atomically.", "author": "Neil deGrasse Tyson"},
            {"quote": "The history of life is written in the genes.", "author": "Sydney Brenner"},
            {"quote": "Knowledge is power. Information is liberating.", "author": "Kofi Annan"},
            {"quote": "One microgram of truth can neutralise tonnes of lies.", "author": "Carl Sagan"},
            {"quote": "The universe is under no obligation to make sense to you.", "author": "Neil deGrasse Tyson"},
            {"quote": "In the middle of difficulty lies opportunity.", "author": "Albert Einstein"},
            {"quote": "Curiosity is the engine of achievement.", "author": "Ken Robinson"},
            {"quote": "The secret of getting ahead is getting started.", "author": "Mark Twain"},
            {"quote": "Do not let what you cannot do interfere with what you can do.", "author": "John Wooden"},
            {"quote": "The only source of knowledge is experience.", "author": "Albert Einstein"},
            {"quote": "The greatest enemy of knowledge is not ignorance, it is the illusion of knowledge.", "author": "Stephen Hawking"},
            {"quote": "It is not the strongest of the species that survives, but the one most responsive to change.", "author": "Charles Darwin"},
            {"quote": "Science is a way of thinking much more than it is a body of knowledge.", "author": "Carl Sagan"},
            {"quote": "If I have seen further it is by standing on the shoulders of Giants.", "author": "Isaac Newton"},
            {"quote": "The beauty of a living thing is not the atoms that go into it, but the way those atoms are put together.", "author": "Carl Sagan"},
            {"quote": "Genes are the stories, and the genome is the library.", "author": "Sam Kean"},
            {"quote": "The microbiome is the last great frontier of human biology.", "author": "Rob Knight"},
        ]

    def display_random_quote(self) -> None:
        if self.quiet:
            return
        if not self.quotes:
            return
        quote_data = random.choice(self.quotes)
        quote = quote_data["quote"]
        author = quote_data["author"]
        quote_color = random.choice(self.quote_colors)
        print()
        print(f"{Color.DIM}{Color.WHITE}{'─' * 80}{Color.RESET}")
        print(f"{quote_color}   \"{quote}\"{Color.RESET}")
        print(f"{Color.BOLD}{Color.WHITE}   — {author}{Color.RESET}")
        print(f"{Color.DIM}{Color.WHITE}{'─' * 80}{Color.RESET}")
        print()

    # ----------------------------------------------------------------------
    # Setup checks
    # ----------------------------------------------------------------------
    def check_abricate_databases(self) -> bool:
        try:
            result = subprocess.run(["abricate", "--list"], capture_output=True, text=True)
            return result.returncode == 0 and "ncbi" in result.stdout.lower()
        except FileNotFoundError:
            return False

    def check_amr_database(self) -> bool:
        amr_module = self.base_dir / "modules" / "amr_module"
        db_root = amr_module / "data" / "amrfinder_db"
        if not db_root.exists():
            return False
        for item in db_root.iterdir():
            if item.is_dir() and item.name.startswith("20"):
                if (item / "version.txt").exists() or any(item.glob("*.hmm")):
                    return True
        return False

    def update_amr_database(self) -> bool:
        amr_module = self.base_dir / "modules" / "amr_module"
        script = amr_module / "p_amrfinder.py"
        if not script.exists():
            self.print_error(f"AMR script not found: {script}")
            return False
        self.print_info("Updating AMRfinderPlus database...")
        cmd = [sys.executable, str(script), "--update-db"]
        result = subprocess.run(cmd, cwd=amr_module, capture_output=True, text=True)
        if result.stdout:
            self.logger.info(f"AMR update STDOUT:\n{result.stdout}")
        if result.stderr:
            self.logger.warning(f"AMR update STDERR:\n{result.stderr}")
        if result.returncode == 0:
            self.print_success("AMR database updated successfully.")
            ver_cmd = [sys.executable, str(script), "--db-version"]
            ver_result = subprocess.run(ver_cmd, cwd=amr_module, capture_output=True, text=True)
            if ver_result.returncode == 0 and ver_result.stdout:
                self.print_info(f"Database version: {ver_result.stdout.strip()}")
            return True
        else:
            self.print_error("AMR database update failed.")
            return False

    # ----------------------------------------------------------------------
    # File handling
    # ----------------------------------------------------------------------
    def find_fasta_files(self, input_path: str) -> List[Path]:
        if '*' in input_path or '?' in input_path:
            matched = glob.glob(input_path)
            files = [Path(f) for f in matched if Path(f).is_file() and
                     f.lower().endswith(('.fna', '.fasta', '.fa', '.fn')) and
                     not Path(f).name.startswith('.')]
            self.print_info(f"Found {len(files)} FASTA files from pattern")
            return sorted(files)

        path = Path(input_path)
        if path.is_file() and path.suffix.lower() in ('.fna', '.fasta', '.fa', '.fn'):
            return [path]

        if path.is_dir():
            files = []
            for ext in ['*.fna', '*.fasta', '*.fa', '*.fn']:
                files.extend(path.glob(ext))
            files = [f for f in files if f.is_file() and not f.name.startswith('.')]
            self.print_info(f"Found {len(files)} FASTA files in directory {path}")
            return sorted(files)

        self.print_error(f"Input path not found: {input_path}")
        return []

    def get_file_pattern(self, fasta_files: List[Path]) -> str:
        if not fasta_files:
            return "*.fna"
        exts = set(f.suffix.lower() for f in fasta_files)
        if len(exts) == 1:
            return f"*{list(exts)[0]}"
        return "*"

    # ----------------------------------------------------------------------
    # Core temporary directory runner
    # ----------------------------------------------------------------------
    def run_module_in_temp(
        self,
        module_name: str,
        fasta_files: List[Path],
        cmd_parts: List[str],
        result_subdir: Optional[str] = None,
        extra_files: Optional[List[Tuple[str, str]]] = None,
    ) -> Tuple[bool, str]:
        module_orig = self.base_dir / "modules" / module_name
        if not module_orig.exists():
            return False, f"Module directory not found: {module_orig}"

        temp_dir = tempfile.mkdtemp(prefix=f"pseudoscope_{module_name}_")
        self.temp_dirs.append(temp_dir)
        self.logger.info(f"Temporary directory for {module_name}: {temp_dir}")

        try:
            shutil.copytree(module_orig, Path(temp_dir) / module_name, dirs_exist_ok=True)

            for f in fasta_files:
                shutil.copy2(f, Path(temp_dir) / f.name)

            script_name = Path(cmd_parts[1]).name if len(cmd_parts) > 1 else None
            if not script_name:
                return False, f"Could not determine script name from command: {cmd_parts}"

            script_path = Path(temp_dir) / module_name / script_name
            if not script_path.exists():
                py_files = list((Path(temp_dir) / module_name).glob("*.py"))
                if py_files:
                    script_path = py_files[0]
                else:
                    return False, f"Could not locate script in {module_name}"

            abs_cmd = [cmd_parts[0], str(script_path)] + cmd_parts[2:]

            self.logger.info(f"Running {module_name}: {' '.join(abs_cmd)}")
            result = subprocess.run(abs_cmd, cwd=temp_dir, capture_output=True, text=True)

            if result.stdout:
                self.logger.info(f"{module_name} STDOUT:\n{result.stdout}")
            if result.stderr:
                self.logger.warning(f"{module_name} STDERR:\n{result.stderr}")

            if result.returncode != 0:
                self.logger.error(f"{module_name} failed with return code {result.returncode}")
                return False, f"Module {module_name} failed (rc={result.returncode})"

            if result_subdir:
                src = Path(temp_dir) / result_subdir
                if src.exists():
                    dst = self.user_output_dir / result_subdir
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)
                    self.logger.info(f"Results copied to {dst}")
                else:
                    self.logger.warning(f"Expected result subdir not found: {src}")

            if extra_files:
                for src_rel, dst_name in extra_files:
                    src = Path(temp_dir) / src_rel
                    if src.exists():
                        dst = self.user_output_dir / dst_name
                        shutil.copy2(src, dst)
                        self.logger.info(f"Copied {dst_name} to output directory")
                    else:
                        self.logger.warning(f"Extra file not found: {src}")

            return True, f"{module_name} completed successfully"

        except Exception as e:
            self.logger.error(f"Exception in {module_name}: {e}\n{traceback.format_exc()}")
            return False, f"Exception: {e}"

        finally:
            if not self.keep_temp and self.running:
                shutil.rmtree(temp_dir, ignore_errors=True)
                self.logger.info(f"Removed temporary directory: {temp_dir}")

    # ----------------------------------------------------------------------
    # Module runners
    # ----------------------------------------------------------------------
    def run_qc(self, fasta_files: List[Path], final_out: Path, threads: int) -> Tuple[bool, str]:
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "p_qc.py", pattern]
        return self.run_module_in_temp("pa_qc_module", fasta_files, cmd, "fasta_qc_results")

    def run_mlst(self, fasta_files: List[Path], final_out: Path, threads: int) -> Tuple[bool, str]:
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "p_mlst.py", "-i", pattern, "-o", "mlst_results", "--batch"]
        return self.run_module_in_temp("mlst_module", fasta_files, cmd, "mlst_results")

    def run_past(self, fasta_files: List[Path], final_out: Path, threads: int,
                 min_pident: Optional[int] = None, min_coverage: Optional[int] = None) -> Tuple[bool, str]:
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "p_past.py", "-i", pattern, "-o", "past_results", "--cpus", str(threads)]
        if min_pident is not None:
            cmd.extend(["--min-pident", str(min_pident)])
        if min_coverage is not None:
            cmd.extend(["--min-coverage", str(min_coverage)])
        return self.run_module_in_temp("past_module", fasta_files, cmd, "past_results")

    def run_abricate(self, fasta_files: List[Path], final_out: Path, threads: int,
                     minid: Optional[int] = None, mincov: Optional[int] = None) -> Tuple[bool, str]:
        if not self.check_abricate_databases():
            msg = "ABRicate databases not found. Please run 'abricate --setupdb' first."
            return False, msg
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "p_abricate.py", pattern]
        if minid is not None:
            cmd.extend(["--minid", str(minid)])
        if mincov is not None:
            cmd.extend(["--mincov", str(mincov)])
        return self.run_module_in_temp("abricate_module", fasta_files, cmd, "pseudo_abricate_results")

    def run_amr(self, fasta_files: List[Path], final_out: Path, threads: int,
                min_identity: Optional[float] = None, min_coverage: Optional[float] = None,
                skip_mutations: bool = False) -> Tuple[bool, str]:
        if not self.check_amr_database():
            msg = "AMR database not found. Please run 'pseudoscope --update-amr-db' first."
            return False, msg
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "p_amrfinder.py", pattern]
        if min_identity is not None:
            cmd.extend(["--min-identity", str(min_identity)])
        if min_coverage is not None:
            cmd.extend(["--min-coverage", str(min_coverage)])
        if skip_mutations:
            cmd.append("--skip-mutations")
        return self.run_module_in_temp(
            "amr_module",
            fasta_files,
            cmd,
            result_subdir="pseudo_amrfinder_results",
            extra_files=[("pseudo_amrfinder_results/mutation_summary.html", "mutation_summary.html")]
        )

    def run_summary(self, final_out: Path) -> Tuple[bool, str]:
        module_orig = self.base_dir / "modules" / "gene_centric_module"
        if not module_orig.exists():
            return False, f"Summary module not found: {module_orig}"

        temp_dir = tempfile.mkdtemp(prefix="pseudoscope_summary_")
        self.temp_dirs.append(temp_dir)
        self.logger.info(f"Temporary directory for summary: {temp_dir}")

        try:
            shutil.copytree(module_orig, Path(temp_dir) / "gene_centric_module", dirs_exist_ok=True)

            for target_name, source_rel in self.summary_html_files.items():
                src = final_out / source_rel
                if src.exists():
                    shutil.copy2(src, Path(temp_dir) / target_name)
                    self.logger.info(f"Copied {target_name} to temporary summary directory")
                else:
                    self.logger.warning(f"Required summary file not found: {src}")

            script_path = Path(temp_dir) / "gene_centric_module" / "p_gene_centric.py"
            if not script_path.exists():
                return False, f"Summary script not found: {script_path}"

            cmd = [sys.executable, str(script_path), "-i", "."]
            self.logger.info(f"Running summary: {' '.join(cmd)}")
            result = subprocess.run(cmd, cwd=temp_dir, capture_output=True, text=True)

            if result.stdout:
                self.logger.info(f"Summary STDOUT:\n{result.stdout}")
            if result.stderr:
                self.logger.warning(f"Summary STDERR:\n{result.stderr}")

            if result.returncode != 0:
                self.logger.error(f"Summary failed with return code {result.returncode}")
                return False, f"Summary failed (rc={result.returncode})"

            src = Path(temp_dir) / self.output_dirs['summary']
            dst = final_out / self.output_dirs['summary']
            if src.exists():
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(src, dst)
                self.logger.info(f"Reports copied to {dst}")
            else:
                self.logger.warning("Summary reports not found")

            return True, "Summary completed successfully"

        except Exception as e:
            self.logger.error(f"Summary exception: {e}\n{traceback.format_exc()}")
            return False, f"Exception: {e}"

        finally:
            if not self.keep_temp and self.running:
                shutil.rmtree(temp_dir, ignore_errors=True)
                self.logger.info(f"Removed temporary directory: {temp_dir}")

    def run_sample_centric(self, final_out: Path) -> Tuple[bool, str]:
        module_orig = self.base_dir / "modules" / "sample_centric_module"
        if not module_orig.exists():
            return False, f"Sample‑centric module not found: {module_orig}"

        temp_dir = tempfile.mkdtemp(prefix="pseudoscope_sample_centric_")
        self.temp_dirs.append(temp_dir)
        self.logger.info(f"Temporary directory for sample‑centric: {temp_dir}")

        try:
            shutil.copytree(module_orig, Path(temp_dir) / "sample_centric_module", dirs_exist_ok=True)

            # Copy required files from final_out to temp_dir
            # 1. QC TSV/HTML
            qc_src = final_out / self.output_dirs['qc']
            if qc_src.exists():
                for f in qc_src.glob("PseudoScope_FASTA_QC_summary.*"):
                    shutil.copy2(f, Path(temp_dir) / f.name)
                    self.logger.info(f"Copied {f.name} to sample‑centric temp")
            # 2. MLST files
            mlst_src = final_out / self.output_dirs['mlst']
            if mlst_src.exists():
                for f in mlst_src.glob("mlst_summary.*"):
                    shutil.copy2(f, Path(temp_dir) / f.name)
                    self.logger.info(f"Copied {f.name} to sample‑centric temp")
            # 3. PAST HTML if exists
            past_src = final_out / self.output_dirs['past']
            if past_src.exists():
                for f in past_src.glob("past_summary.html"):
                    shutil.copy2(f, Path(temp_dir) / f.name)
                    self.logger.info(f"Copied {f.name} to sample‑centric temp")
            # 4. AMRfinder TSV and HTML + mutation TSV
            amr_src = final_out / self.output_dirs['amr']
            if amr_src.exists():
                for f in amr_src.glob("pseudo_amrfinder_summary.*"):
                    shutil.copy2(f, Path(temp_dir) / f.name)
                    self.logger.info(f"Copied {f.name} to sample‑centric temp")
                # mutation summary TSV (without prefix)
                mut_tsv = amr_src / "mutation_summary.tsv"
                if mut_tsv.exists():
                    shutil.copy2(mut_tsv, Path(temp_dir) / mut_tsv.name)
                    self.logger.info(f"Copied {mut_tsv.name} to sample‑centric temp")
            # 5. ABRicate TSVs and HTMLs
            abr_src = final_out / self.output_dirs['abricate']
            if abr_src.exists():
                for f in abr_src.glob("pseudo_*_abricate_summary.tsv"):
                    shutil.copy2(f, Path(temp_dir) / f.name)
                    self.logger.info(f"Copied {f.name} to sample‑centric temp")
                for f in abr_src.glob("pseudo_*_summary_report.html"):
                    shutil.copy2(f, Path(temp_dir) / f.name)
                    self.logger.info(f"Copied {f.name} to sample‑centric temp")

            script_path = Path(temp_dir) / "sample_centric_module" / "p_sample_centric.py"
            if not script_path.exists():
                return False, f"Sample‑centric script not found: {script_path}"

            cmd = [sys.executable, str(script_path), "-i", "."]
            self.logger.info(f"Running sample‑centric: {' '.join(cmd)}")
            result = subprocess.run(cmd, cwd=temp_dir, capture_output=True, text=True)

            if result.stdout:
                self.logger.info(f"Sample‑centric STDOUT:\n{result.stdout}")
            if result.stderr:
                self.logger.warning(f"Sample‑centric STDERR:\n{result.stderr}")

            if result.returncode != 0:
                self.logger.error(f"Sample‑centric failed with return code {result.returncode}")
                return False, f"Sample‑centric failed (rc={result.returncode})"

            src = Path(temp_dir) / self.output_dirs['sample_centric']
            dst = final_out / self.output_dirs['sample_centric']
            if src.exists():
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(src, dst)
                self.logger.info(f"Sample‑centric reports copied to {dst}")
            else:
                self.logger.warning("Sample‑centric reports not found")

            return True, "Sample‑centric completed successfully"

        except Exception as e:
            self.logger.error(f"Sample‑centric exception: {e}\n{traceback.format_exc()}")
            return False, f"Exception: {e}"

        finally:
            if not self.keep_temp and self.running:
                shutil.rmtree(temp_dir, ignore_errors=True)
                self.logger.info(f"Removed temporary directory: {temp_dir}")

    def run_viz(self, final_out: Path) -> Tuple[bool, str]:
        module_orig = self.base_dir / "modules" / "viz_module"
        if not module_orig.exists():
            return False, f"Visualisation module not found: {module_orig}"

        temp_dir = tempfile.mkdtemp(prefix="pseudoscope_viz_")
        self.temp_dirs.append(temp_dir)
        self.logger.info(f"Temporary directory for visualisation: {temp_dir}")

        try:
            shutil.copytree(module_orig, Path(temp_dir) / "viz_module", dirs_exist_ok=True)

            ultimate_dir = final_out / self.output_dirs['summary']
            if ultimate_dir.exists():
                for csv_file in ultimate_dir.glob("*.csv"):
                    shutil.copy2(csv_file, Path(temp_dir) / csv_file.name)
                    self.logger.info(f"Copied {csv_file.name} to temporary viz directory")
            else:
                self.logger.warning("Ultimate reports directory not found")

            qc_dir = final_out / self.output_dirs['qc']
            if qc_dir.exists():
                qc_tsv = qc_dir / "PseudoScope_FASTA_QC_summary.tsv"
                qc_html = qc_dir / "PseudoScope_FASTA_QC_summary.html"
                if qc_tsv.exists():
                    shutil.copy2(qc_tsv, Path(temp_dir) / qc_tsv.name)
                if qc_html.exists():
                    shutil.copy2(qc_html, Path(temp_dir) / qc_html.name)

            script_path = Path(temp_dir) / "viz_module" / "p_visualizer.py"
            if not script_path.exists():
                return False, f"Visualizer script not found: {script_path}"

            cmd = [sys.executable, str(script_path), "-i", "."]
            self.logger.info(f"Running visualisation: {' '.join(cmd)}")
            result = subprocess.run(cmd, cwd=temp_dir, capture_output=True, text=True)

            if result.stdout:
                self.logger.info(f"Visualisation STDOUT:\n{result.stdout}")
            if result.stderr:
                self.logger.warning(f"Visualisation STDERR:\n{result.stderr}")

            if result.returncode != 0:
                self.logger.error(f"Visualisation failed with return code {result.returncode}")
                return False, f"Visualisation failed (rc={result.returncode})"

            dashboard = Path(temp_dir) / "genius_pseudomonas_visual_dashboard.html"
            dst_dir = final_out / self.output_dirs['viz']
            dst_dir.mkdir(parents=True, exist_ok=True)
            if dashboard.exists():
                shutil.copy2(dashboard, dst_dir / dashboard.name)
                self.logger.info(f"Dashboard copied to {dst_dir / dashboard.name}")
            else:
                self.logger.warning("Dashboard not found")

            return True, "Visualisation completed successfully"

        except Exception as e:
            self.logger.error(f"Visualisation exception: {e}\n{traceback.format_exc()}")
            return False, f"Exception: {e}"

        finally:
            if not self.keep_temp and self.running:
                shutil.rmtree(temp_dir, ignore_errors=True)
                self.logger.info(f"Removed temporary directory: {temp_dir}")

    # ----------------------------------------------------------------------
    # Cleanup top‑level HTML files
    # ----------------------------------------------------------------------
    def _cleanup_top_level_html(self, final_out: Path) -> None:
        """Remove any .html files directly in the output directory (not in subdirectories)."""
        for item in final_out.iterdir():
            if item.is_file() and item.suffix.lower() == '.html':
                try:
                    item.unlink()
                    self.logger.info(f"Removed top-level HTML file: {item.name}")
                except Exception as e:
                    self.logger.warning(f"Could not remove {item.name}: {e}")

    # ----------------------------------------------------------------------
    # Full pipeline
    # ----------------------------------------------------------------------
    def run_complete_analysis(
        self,
        input_path: str,
        output_dir: str,
        threads: int = 1,
        skip: Dict[str, bool] = None,
        skip_summary: bool = False,
        skip_sample_centric: bool = False,
        skip_viz: bool = False,
        update_amr_only: bool = False,
        abricate_minid: Optional[int] = None,
        abricate_mincov: Optional[int] = None,
        past_min_pident: Optional[int] = None,
        past_min_coverage: Optional[int] = None,
        amr_min_identity: Optional[float] = None,
        amr_min_coverage: Optional[float] = None,
        amr_skip_mutations: bool = False,
    ) -> None:
        if update_amr_only:
            self.update_amr_database()
            return

        if skip is None:
            skip = {}

        start_time = datetime.now()
        self.print_header("PSEUDOSCOPE", "P. aeruginosa Genomic Analysis Pipeline")
        self.print_info(f"Version {__version__}")

        final_out = Path(output_dir)
        final_out.mkdir(parents=True, exist_ok=True)
        self.setup_logging(final_out)

        fasta_files = self.find_fasta_files(input_path)
        if not fasta_files:
            self.print_error("No FASTA files found. Exiting.")
            return

        self.print_success(f"Found {len(fasta_files)} FASTA files")
        for sub in self.output_dirs.values():
            (final_out / sub).mkdir(exist_ok=True)

        plan = [
            ("QC", not skip.get('qc', False)),
            ("MLST", not skip.get('mlst', False)),
            ("PAST", not skip.get('past', False)),
            ("ABRicate", not skip.get('abricate', False)),
            ("AMR", not skip.get('amr', False)),
            ("Ultimate Reporter (Gene‑centric)", not skip_summary),
            ("Sample‑centric Reporter", not skip_sample_centric),
            ("Visualisation", not skip_viz),
        ]
        self.print_header("ANALYSIS PLAN", "Modules to be executed")
        for name, enabled in plan:
            status = f"{Color.BRIGHT_GREEN}✅ ENABLED{Color.RESET}" if enabled else f"{Color.YELLOW}⏸️  SKIPPED{Color.RESET}"
            self._print(f"   {status} - {name}", force=True)
            if self.logger:
                self.logger.info(f"Plan: {name} {'ENABLED' if enabled else 'SKIPPED'}")

        # First batch: QC, MLST, PAST (parallel)
        batch1_tasks = []
        if not skip.get('qc', False):
            batch1_tasks.append(("QC", self.run_qc))
        if not skip.get('mlst', False):
            batch1_tasks.append(("MLST", self.run_mlst))
        if not skip.get('past', False):
            batch1_tasks.append(("PAST", lambda f, o, t: self.run_past(f, o, t, past_min_pident, past_min_coverage)))

        if batch1_tasks:
            self.print_info(f"Running {len(batch1_tasks)} analyses in parallel...")
            results = {}
            with ThreadPoolExecutor(max_workers=len(batch1_tasks)) as executor:
                future_to_name = {executor.submit(run, fasta_files, final_out, threads): name
                                  for name, run in batch1_tasks}
                for future in as_completed(future_to_name):
                    name = future_to_name[future]
                    try:
                        success, log = future.result()
                        results[name] = (success, log)
                    except Exception as e:
                        self.logger.error(f"Exception in {name}: {e}")
                        results[name] = (False, f"Exception: {e}")
            for name, _ in batch1_tasks:
                success, log = results.get(name, (False, "No result"))
                self.print_header(f"{name} Analysis")
                lines = log.strip().split('\n')
                for line in lines[:5]:
                    self._print(line, force=True)
                if len(lines) > 5:
                    for line in lines[5:8]:
                        self._print(line, force=True)
                if success:
                    self.print_success(f"{name} completed")
                else:
                    self.print_error(f"{name} failed")
                self.display_random_quote()
        else:
            self.print_info("No first‑batch modules selected.")

        # Second batch: ABRicate
        if not skip.get('abricate', False):
            self.print_header("ABRICATE ANALYSIS", "Comprehensive Resistance & Virulence")
            success, log = self.run_abricate(fasta_files, final_out, threads, abricate_minid, abricate_mincov)
            lines = log.strip().split('\n')
            for line in lines[:8]:
                self._print(line, force=True)
            if success:
                self.print_success("ABRicate completed")
            else:
                self.print_error("ABRicate failed")
            self.display_random_quote()
        else:
            self.print_info("Skipping ABRicate analysis.")

        # Third batch: AMR
        if not skip.get('amr', False):
            self.print_header("AMR ANALYSIS", "Antimicrobial Resistance Gene Detection")
            success, log = self.run_amr(fasta_files, final_out, threads, amr_min_identity, amr_min_coverage, amr_skip_mutations)
            lines = log.strip().split('\n')
            for line in lines[:8]:
                self._print(line, force=True)
            if success:
                self.print_success("AMR completed")
            else:
                self.print_error("AMR failed")
            self.display_random_quote()
        else:
            self.print_info("Skipping AMR analysis.")

        # Gene‑centric summary
        if not skip_summary:
            self.print_header("ULTIMATE REPORTER (Gene‑centric)", "Gene‑centric Integration")
            success, log = self.run_summary(final_out)
            lines = log.strip().split('\n')
            for line in lines[:8]:
                self._print(line, force=True)
            if success:
                self.print_success("Ultimate reporter completed")
            else:
                self.print_warning("Ultimate reporter had issues")
            self.display_random_quote()
        else:
            self.print_info("Skipping gene‑centric reporter.")

        # Sample‑centric summary
        if not skip_sample_centric:
            self.print_header("SAMPLE‑CENTRIC REPORTER", "Interactive Isolate Boxes")
            success, log = self.run_sample_centric(final_out)
            lines = log.strip().split('\n')
            for line in lines[:8]:
                self._print(line, force=True)
            if success:
                self.print_success("Sample‑centric reporter completed")
            else:
                self.print_warning("Sample‑centric reporter had issues")
            self.display_random_quote()
        else:
            self.print_info("Skipping sample‑centric reporter.")

        # Visualisation
        if not skip_viz and not skip_summary:
            self.print_header("VISUALISATION", "Interactive Dashboards")
            success, log = self.run_viz(final_out)
            lines = log.strip().split('\n')
            for line in lines[:8]:
                self._print(line, force=True)
            if success:
                self.print_success("Visualisation completed")
            else:
                self.print_warning("Visualisation had issues")
            self.display_random_quote()
        elif skip_summary:
            self.print_info("Skipping visualisation because summary was skipped.")

        # Clean up top‑level HTML files
        self._cleanup_top_level_html(final_out)

        elapsed = datetime.now() - start_time
        self.print_header("ANALYSIS COMPLETE", f"Time: {str(elapsed).split('.')[0]}")
        self.print_success(f"All results in: {final_out}")
        for subdir in sorted(final_out.iterdir()):
            if subdir.is_dir():
                cnt = len(list(subdir.glob("*")))
                self.print_info(f"  📁 {subdir.name} ({cnt} files)")
        self.display_random_quote()

        # Print citation
        self.print_info("\nPlease cite:")
        self.print_info("  Beckley B, Amarh V. PseudoScope: a species‑specific bioinformatics suite for rapid and accessible Pseudomonas aeruginosa genomic analysis. Github 2026 https://github.com/brown-beckley/pseudoscope")

        # Print any warnings/errors
        if self.warnings:
            self.print_warning("\nWarnings encountered:")
            for w in self.warnings:
                self.print_info(f"  ⚠️ {w}")
        if self.errors:
            self.print_error("\nErrors encountered:")
            for e in self.errors:
                self.print_info(f"  ✗ {e}")


# -----------------------------------------------------------------------------
# Coloured help function
# -----------------------------------------------------------------------------
def print_colored_help() -> None:
    usage = f"{Color.BRIGHT_CYAN}usage: pseudoscope [OPTIONS]{Color.RESET}"
    description = f"{Color.BRIGHT_WHITE}PseudoScope: Complete P. aeruginosa genomic analysis pipeline (v{__version__}){Color.RESET}"

    options = [
        (f"{Color.BRIGHT_GREEN}-i INPUT, --input INPUT{Color.RESET}", f"Input FASTA file, directory, or glob pattern (e.g., \"{Color.BRIGHT_YELLOW}*.fna{Color.RESET}\")"),
        (f"{Color.BRIGHT_GREEN}-o OUTPUT, --output OUTPUT{Color.RESET}", "Output directory for results"),
        (f"{Color.BRIGHT_GREEN}-t THREADS, --threads THREADS{Color.RESET}", f"Number of threads (default: {Color.BRIGHT_YELLOW}2{Color.RESET})"),
        (f"{Color.BRIGHT_GREEN}--quiet{Color.RESET}", "Suppress all non‑error console output"),
        (f"{Color.BRIGHT_GREEN}--keep-temp{Color.RESET}", "Do not delete temporary directories (for debugging)"),
        (f"{Color.BRIGHT_GREEN}--version{Color.RESET}", "Show version and exit"),
        (f"{Color.BRIGHT_GREEN}--update-amr-db{Color.RESET}", "Update AMRfinderPlus database and exit"),
        (f"{Color.BRIGHT_GREEN}--skip-qc{Color.RESET}", "Skip FASTA QC"),
        (f"{Color.BRIGHT_GREEN}--skip-mlst{Color.RESET}", "Skip MLST"),
        (f"{Color.BRIGHT_GREEN}--skip-past{Color.RESET}", "Skip PAST serotyping"),
        (f"{Color.BRIGHT_GREEN}--skip-abricate{Color.RESET}", "Skip ABRicate"),
        (f"{Color.BRIGHT_GREEN}--skip-amr{Color.RESET}", "Skip AMRfinderPlus"),
        (f"{Color.BRIGHT_GREEN}--skip-summary{Color.RESET}", "Skip gene‑centric ultimate summary report"),
        (f"{Color.BRIGHT_GREEN}--skip-sample-centric{Color.RESET}", "Skip sample‑centric interactive report"),
        (f"{Color.BRIGHT_GREEN}--skip-viz{Color.RESET}", "Skip visualisation dashboard"),
        (f"{Color.BRIGHT_GREEN}--abricate-minid INT{Color.RESET}", "ABRicate minimum %identity (default: 80)"),
        (f"{Color.BRIGHT_GREEN}--abricate-mincov INT{Color.RESET}", "ABRicate minimum %coverage (default: 80)"),
        (f"{Color.BRIGHT_GREEN}--past-min-pident INT{Color.RESET}", "PAST minimum percent identity (default: 95)"),
        (f"{Color.BRIGHT_GREEN}--past-min-coverage INT{Color.RESET}", "PAST minimum percent coverage (default: 95)"),
        (f"{Color.BRIGHT_GREEN}--amr-min-identity FLOAT{Color.RESET}", "AMR minimum identity (0-1, default: module default)"),
        (f"{Color.BRIGHT_GREEN}--amr-min-coverage FLOAT{Color.RESET}", "AMR minimum coverage (0-1, default: module default)"),
        (f"{Color.BRIGHT_GREEN}--amr-skip-mutations{Color.RESET}", "Disable point mutation reporting in AMR (enabled by default)"),
        (f"{Color.BRIGHT_GREEN}-h, --help{Color.RESET}", "Show this help message and exit"),
    ]

    epilog = f"""
{Color.BRIGHT_YELLOW}EXAMPLES:{Color.RESET}
  # Basic analysis on all .fna files in current directory
  pseudoscope -i "{Color.BRIGHT_CYAN}*.fna{Color.RESET}" -o results

  # Use a directory containing FASTA files
  pseudoscope -i {Color.BRIGHT_CYAN}genomes/{Color.RESET} -o results --threads 4

  # Skip QC and visualisation, only run MLST and PAST
  pseudoscope -i "{Color.BRIGHT_CYAN}*.fasta{Color.RESET}" -o results --skip-qc --skip-summary --skip-viz

  # Custom thresholds for PAST and ABRicate
  pseudoscope -i "{Color.BRIGHT_CYAN}*.fna{Color.RESET}" -o results --past-min-pident 90 --past-min-coverage 85 --abricate-minid 85 --abricate-mincov 80

  # AMR with custom identity/coverage and disable mutation reporting
  pseudoscope -i "{Color.BRIGHT_CYAN}*.fna{Color.RESET}" -o results --amr-min-identity 0.95 --amr-min-coverage 0.9 --amr-skip-mutations

  # Update AMR database (mandatory before first run)
  pseudoscope --update-amr-db

  # Show version
  pseudoscope --version

  # Quiet mode (minimal console output)
  pseudoscope -i "{Color.BRIGHT_CYAN}*.fna{Color.RESET}" -o results --quiet

  # Keep temporary directories for debugging
  pseudoscope -i "{Color.BRIGHT_CYAN}*.fna{Color.RESET}" -o results --keep-temp

{Color.BRIGHT_YELLOW}REQUIRED BEFORE FIRST ANALYSIS:{Color.RESET}
  1. {Color.BRIGHT_GREEN}abricate --setupdb{Color.RESET}   (setup ABRicate databases)
  2. {Color.BRIGHT_GREEN}pseudoscope --update-amr-db{Color.RESET}   (download AMRfinderPlus database)

{Color.BRIGHT_YELLOW}SUPPORTED FASTA FORMATS:{Color.RESET} {Color.BRIGHT_CYAN}.fna, .fasta, .fa, .fn{Color.RESET}

{Color.BRIGHT_YELLOW}CITATION:{Color.RESET}
  Beckley B, Amarh V. PseudoScope: a species‑specific bioinformatics suite for rapid and accessible
  Pseudomonas aeruginosa genomic analysis. Github 2026 https://github.com/brown-beckley/pseudoscope
    """

    print()
    print(usage)
    print()
    print(description)
    print()
    print(f"{Color.BRIGHT_YELLOW}OPTIONS:{Color.RESET}")
    for opt, desc in options:
        print(f"  {opt:<38} {desc}")
    print(epilog)


def main() -> None:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-i', '--input', help=argparse.SUPPRESS)
    parser.add_argument('-o', '--output', help=argparse.SUPPRESS)
    parser.add_argument('-t', '--threads', type=int, default=2, help=argparse.SUPPRESS)
    parser.add_argument('--quiet', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--keep-temp', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--version', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--update-amr-db', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-qc', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-mlst', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-past', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-abricate', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-amr', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-summary', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-sample-centric', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-viz', action='store_true', help=argparse.SUPPRESS)
    # Thresholds
    parser.add_argument('--abricate-minid', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--abricate-mincov', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--past-min-pident', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--past-min-coverage', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--amr-min-identity', type=float, help=argparse.SUPPRESS)
    parser.add_argument('--amr-min-coverage', type=float, help=argparse.SUPPRESS)
    parser.add_argument('--amr-skip-mutations', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('-h', '--help', action='store_true', help=argparse.SUPPRESS)

    args = parser.parse_args()

    if args.help:
        print_colored_help()
        sys.exit(0)

    if args.version:
        print(f"pseudoscope {__version__}")
        sys.exit(0)

    if args.update_amr_db:
        orch = PseudoScopeOrchestrator(quiet=args.quiet, keep_temp=args.keep_temp)
        orch.run_complete_analysis("", "", update_amr_only=True)
        sys.exit(0)

    if not args.input or not args.output:
        parser.error("When not using --update-amr-db, both -i/--input and -o/--output are required.")

    skip = {
        'qc': args.skip_qc,
        'mlst': args.skip_mlst,
        'past': args.skip_past,
        'abricate': args.skip_abricate,
        'amr': args.skip_amr,
    }

    orch = PseudoScopeOrchestrator(quiet=args.quiet, keep_temp=args.keep_temp)
    try:
        orch.run_complete_analysis(
            input_path=args.input,
            output_dir=args.output,
            threads=args.threads,
            skip=skip,
            skip_summary=args.skip_summary,
            skip_sample_centric=args.skip_sample_centric,
            skip_viz=args.skip_viz,
            abricate_minid=args.abricate_minid,
            abricate_mincov=args.abricate_mincov,
            past_min_pident=args.past_min_pident,
            past_min_coverage=args.past_min_coverage,
            amr_min_identity=args.amr_min_identity,
            amr_min_coverage=args.amr_min_coverage,
            amr_skip_mutations=args.amr_skip_mutations,
        )
    except KeyboardInterrupt:
        orch.print_error("Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        orch.print_error(f"Critical error: {e}")
        if not args.quiet:
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()