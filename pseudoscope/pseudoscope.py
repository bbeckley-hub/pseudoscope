#!/usr/bin/env python3
"""
PseudoScope - Unified Orchestrator for P. aeruginosa Analysis
Author: Brown Beckley <brownbeckley94@gmail.com>
GitHub: bbeckley-hub
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-04-26
Version: 1.0.0
"""

import os
import sys
import glob
import argparse
import subprocess
import shutil
import random
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

__version__ = "1.0.0"


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
    def __init__(self, quiet: bool = False, verbose: bool = False):
        self.base_dir = Path(__file__).parent
        self.quiet = quiet
        self.verbose = verbose
        self.setup_colors()
        self.quotes = self._get_scientific_quotes()
        self.quote_colors = [
            Color.BRIGHT_CYAN, Color.BRIGHT_GREEN, Color.BRIGHT_YELLOW,
            Color.BRIGHT_MAGENTA, Color.BRIGHT_BLUE, Color.BRIGHT_RED,
            Color.CYAN, Color.GREEN, Color.YELLOW, Color.MAGENTA
        ]

        # Final output subdirectories (user's output folder)
        self.output_dirs = {
            'qc': 'fasta_qc_results',
            'mlst': 'mlst_results',
            'past': 'past_results',
            'abricate': 'pseudo_abricate_results',
            'amr': 'pseudo_amrfinder_results',
            'summary': 'GENIUS_PSEUDOMONAS_ULTIMATE_REPORTS',
            'viz': 'GENIUS_PSEUDOMONAS_VISUAL_DASHBOARD'
        }

        # HTML/TSV files required by the summary module
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
            'pseudo_ecoli_vf_summary_report.html': 'pseudo_abricate_results/pseudo_ecoli_vf_summary_report.html'
        }

    def setup_colors(self):
        self.color_info = Color.CYAN
        self.color_success = Color.BRIGHT_GREEN
        self.color_warning = Color.BRIGHT_YELLOW
        self.color_error = Color.BRIGHT_RED
        self.color_highlight = Color.BRIGHT_CYAN
        self.color_banner = Color.BRIGHT_MAGENTA
        self.color_module = Color.BRIGHT_BLUE
        self.color_sample = Color.GREEN
        self.color_file = Color.YELLOW
        self.color_reset = Color.RESET

    def _print(self, message: str, color: str = Color.RESET, bold: bool = False, force: bool = False):
        if self.quiet and not force:
            return
        style = Color.BOLD if bold else ''
        print(f"{style}{color}{message}{Color.RESET}")

    def print_header(self, title: str, subtitle: str = ""):
        if self.quiet:
            return
        print()
        print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'='*80}{Color.RESET}")
        print(f"{Color.BOLD}{Color.BRIGHT_CYAN}{' ' * 20}{title}{Color.RESET}")
        if subtitle:
            print(f"{Color.DIM}{Color.WHITE}{' ' * 22}{subtitle}{Color.RESET}")
        print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'='*80}{Color.RESET}")
        print()

    def print_info(self, message: str):
        if not self.quiet:
            print(f"{self.color_info}[INFO]{Color.RESET} {message}")

    def print_success(self, message: str):
        if not self.quiet:
            print(f"{self.color_success}✓{Color.RESET} {message}")

    def print_warning(self, message: str):
        if not self.quiet:
            print(f"{self.color_warning}⚠️{Color.RESET} {message}")

    def print_error(self, message: str):
        print(f"{self.color_error}✗{Color.RESET} {message}")

    def print_command(self, cmd_parts: List[str], cwd: Path):
        """Show a short, readable version of the command when verbose."""
        if not self.verbose:
            return
        # Try to make path relative to base_dir
        rel_parts = []
        for p in cmd_parts:
            try:
                rel = str(Path(p).relative_to(self.base_dir))
                rel_parts.append(rel)
            except ValueError:
                rel_parts.append(p)
        cmd_str = " ".join(rel_parts)
        self._print(f"[RUN] {cmd_str}", color=Color.DIM, force=True)

    def _get_scientific_quotes(self):
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

    def display_random_quote(self):
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
    # Mandatory setup checks
    # ----------------------------------------------------------------------
    def check_abricate_databases(self) -> bool:
        """Check if abricate databases are installed."""
        try:
            result = subprocess.run(["abricate", "--list"], capture_output=True, text=True)
            if result.returncode == 0 and "ncbi" in result.stdout.lower():
                return True
            return False
        except FileNotFoundError:
            return False

    def check_amr_database(self) -> bool:
        """Check if AMR database is present (look for latest dated folder)."""
        amr_module = self.base_dir / "modules" / "amr_module"
        db_root = amr_module / "data" / "amrfinder_db"
        if not db_root.exists():
            return False
        # Look for a folder starting with 20
        for item in db_root.iterdir():
            if item.is_dir() and item.name.startswith("20"):
                # Also check that version.txt or some file exists
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
        if result.returncode == 0:
            self.print_success("AMR database updated successfully.")
            # Show version
            ver_cmd = [sys.executable, str(script), "--db-version"]
            ver_result = subprocess.run(ver_cmd, cwd=amr_module, capture_output=True, text=True)
            if ver_result.returncode == 0:
                self.print_info(f"Database version: {ver_result.stdout.strip()}")
            return True
        else:
            self.print_error("AMR database update failed.")
            if result.stderr:
                print(result.stderr)
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

    def copy_fasta_to_module(self, fasta_files: List[Path], module_path: Path):
        module_path.mkdir(parents=True, exist_ok=True)
        for f in fasta_files:
            target = module_path / f.name
            if not target.exists():
                shutil.copy2(f, target)

    def cleanup_module(self, module_path: Path, fasta_files: List[Path]):
        try:
            for f in fasta_files:
                temp = module_path / f.name
                if temp.exists():
                    temp.unlink()
            for d in self.output_dirs.values():
                dir_path = module_path / d
                if dir_path.exists():
                    shutil.rmtree(dir_path)
            for html in module_path.glob("*.html"):
                html.unlink()
        except Exception as e:
            self.print_warning(f"Cleanup issue in {module_path.name}: {e}")

    # ----------------------------------------------------------------------
    # Module runners
    # ----------------------------------------------------------------------
    def run_qc(self, fasta_files: List[Path], final_out: Path, threads: int) -> Tuple[bool, str]:
        module_dir = self.base_dir / "modules" / "pa_qc_module"
        script = module_dir / "p_qc.py"
        if not script.exists():
            return False, f"QC script missing: {script}"
        self.copy_fasta_to_module(fasta_files, module_dir)
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, str(script), pattern]
        log = f"Running QC on {pattern}\n"
        self.print_command(cmd, module_dir)
        result = subprocess.run(cmd, cwd=module_dir, capture_output=True, text=True)
        if result.returncode != 0:
            log += "⚠️ QC had warnings/errors\n"
            if result.stderr:
                log += f"stderr: {result.stderr[:500]}\n"
        else:
            log += "✓ QC completed\n"
        src = module_dir / self.output_dirs['qc']
        dst = final_out / self.output_dirs['qc']
        if src.exists():
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(src, dst)
            log += f"✓ Results copied to {dst}\n"
        else:
            log += f"⚠️ QC results not found\n"
        return result.returncode == 0, log

    def run_mlst(self, fasta_files: List[Path], final_out: Path, threads: int) -> Tuple[bool, str]:
        module_dir = self.base_dir / "modules" / "mlst_module"
        script = module_dir / "p_mlst.py"
        if not script.exists():
            return False, f"MLST script missing: {script}"
        self.copy_fasta_to_module(fasta_files, module_dir)
        pattern = self.get_file_pattern(fasta_files)
        out_sub = self.output_dirs['mlst']
        cmd = [sys.executable, str(script), "-i", pattern, "-o", out_sub, "--batch"]
        log = f"Running MLST on {pattern}\n"
        self.print_command(cmd, module_dir)
        result = subprocess.run(cmd, cwd=module_dir, capture_output=True, text=True)
        if result.returncode != 0:
            log += "⚠️ MLST had warnings/errors\n"
            if result.stderr:
                log += f"stderr: {result.stderr[:500]}\n"
        else:
            log += "✓ MLST completed\n"
        src = module_dir / out_sub
        dst = final_out / out_sub
        if src.exists():
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(src, dst)
            log += f"✓ Results copied to {dst}\n"
        else:
            log += f"⚠️ MLST results not found\n"
        return result.returncode == 0, log

    def run_past(self, fasta_files: List[Path], final_out: Path, threads: int) -> Tuple[bool, str]:
        module_dir = self.base_dir / "modules" / "past_module"
        script = module_dir / "p_past.py"
        if not script.exists():
            return False, f"PAST script missing: {script}"
        self.copy_fasta_to_module(fasta_files, module_dir)
        pattern = self.get_file_pattern(fasta_files)
        out_sub = self.output_dirs['past']
        # Use shell to handle quoted glob pattern
        shell_cmd = f'{sys.executable} {script} -i "{pattern}" -o {out_sub} --cpus {threads}'
        log = f"Running PAST on {pattern}\n"
        self.print_command([sys.executable, str(script), "-i", pattern, "-o", out_sub, "--cpus", str(threads)], module_dir)
        result = subprocess.run(shell_cmd, cwd=module_dir, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            log += "⚠️ PAST had warnings/errors\n"
            if result.stderr:
                log += f"stderr: {result.stderr[:500]}\n"
        else:
            log += "✓ PAST completed\n"
        src = module_dir / out_sub
        dst = final_out / out_sub
        if src.exists():
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(src, dst)
            log += f"✓ Results copied to {dst}\n"
        else:
            log += f"⚠️ PAST results not found\n"
        return result.returncode == 0, log

    def run_abricate(self, fasta_files: List[Path], final_out: Path, threads: int) -> Tuple[bool, str]:
        # Mandatory abricate database check
        if not self.check_abricate_databases():
            msg = ("ABRicate databases not found. Please run 'abricate --setupdb' first.\n"
                   "Then re-run your analysis.")
            return False, msg
        module_dir = self.base_dir / "modules" / "abricate_module"
        script = module_dir / "p_abricate.py"
        if not script.exists():
            return False, f"ABRicate script missing: {script}"
        self.copy_fasta_to_module(fasta_files, module_dir)
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, str(script), pattern]
        log = f"Running ABRicate on {pattern}\n"
        self.print_command(cmd, module_dir)
        result = subprocess.run(cmd, cwd=module_dir, capture_output=True, text=True)
        if result.returncode != 0:
            log += "⚠️ ABRicate had warnings/errors\n"
            if result.stderr:
                log += f"stderr: {result.stderr[:500]}\n"
        else:
            log += "✓ ABRicate completed\n"
        src = module_dir / "pseudo_abricate_results"
        dst = final_out / self.output_dirs['abricate']
        if src.exists():
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(src, dst)
            log += f"✓ Results copied to {dst}\n"
        else:
            log += f"⚠️ ABRicate results not found\n"
        return result.returncode == 0, log

    def run_amr(self, fasta_files: List[Path], final_out: Path, threads: int) -> Tuple[bool, str]:
        # Mandatory AMR database check
        if not self.check_amr_database():
            msg = ("AMR database not found. Please run 'pseudoscope --update-amr-db' first.\n"
                   "Then re-run your analysis.")
            return False, msg
        module_dir = self.base_dir / "modules" / "amr_module"
        script = module_dir / "p_amrfinder.py"
        if not script.exists():
            return False, f"AMR script missing: {script}"
        self.copy_fasta_to_module(fasta_files, module_dir)
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, str(script), pattern]
        log = f"Running AMRfinderPlus on {pattern}\n"
        self.print_command(cmd, module_dir)
        result = subprocess.run(cmd, cwd=module_dir, capture_output=True, text=True)
        if result.returncode != 0:
            log += "⚠️ AMR had warnings/errors\n"
            if result.stderr:
                log += f"stderr: {result.stderr[:500]}\n"
        else:
            log += "✓ AMR completed\n"
        src = module_dir / "pseudo_amrfinder_results"
        dst = final_out / self.output_dirs['amr']
        if src.exists():
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(src, dst)
            log += f"✓ Results copied to {dst}\n"
        else:
            log += f"⚠️ AMR results not found\n"
        return result.returncode == 0, log

    def run_summary(self, final_out: Path) -> Tuple[bool, str]:
        module_dir = self.base_dir / "modules" / "summary_module"
        script = module_dir / "p_ultimate.py"
        if not script.exists():
            return False, f"Summary script missing: {script}"
        # Copy required summary files
        log = "Copying summary files:\n"
        for target_name, source_rel in self.summary_html_files.items():
            src = final_out / source_rel
            if src.exists():
                shutil.copy2(src, module_dir / target_name)
                log += f"  ✓ {target_name}\n"
            else:
                log += f"  ✗ {target_name} (missing)\n"
        # Run ultimate reporter
        cmd = [sys.executable, str(script), "-i", "."]
        self.print_command(cmd, module_dir)
        result = subprocess.run(cmd, cwd=module_dir, capture_output=True, text=True)
        log += result.stdout
        if result.stderr:
            log += f"stderr: {result.stderr[:500]}\n"
        if result.returncode == 0:
            log += "✓ Ultimate reporter completed\n"
        else:
            log += "⚠️ Ultimate reporter had issues\n"
        src = module_dir / self.output_dirs['summary']
        dst = final_out / self.output_dirs['summary']
        if src.exists():
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(src, dst)
            log += f"✓ Reports copied to {dst}\n"
        else:
            log += f"⚠️ Reports not found\n"
        return result.returncode == 0, log

    def run_viz(self, final_out: Path) -> Tuple[bool, str]:
        module_dir = self.base_dir / "modules" / "viz_module"
        script = module_dir / "p_visualizer.py"
        if not script.exists():
            return False, f"Visualizer script missing: {script}"
        ultimate_dir = final_out / self.output_dirs['summary']
        log = "Copying CSV files:\n"
        if ultimate_dir.exists():
            for csv_file in ultimate_dir.glob("*.csv"):
                shutil.copy2(csv_file, module_dir / csv_file.name)
                log += f"  ✓ {csv_file.name}\n"
        else:
            log += f"⚠️ Ultimate reports directory not found\n"
        cmd = [sys.executable, str(script), "-i", "."]
        self.print_command(cmd, module_dir)
        result = subprocess.run(cmd, cwd=module_dir, capture_output=True, text=True)
        log += result.stdout
        if result.stderr:
            log += f"stderr: {result.stderr[:500]}\n"
        if result.returncode == 0:
            log += "✓ Visualizer completed\n"
        else:
            log += "⚠️ Visualizer had issues\n"
        dashboard = module_dir / "genius_pseudomonas_visual_dashboard.html"
        dst = final_out / self.output_dirs['viz'] / "genius_pseudomonas_visual_dashboard.html"
        if dashboard.exists():
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(dashboard, dst)
            log += f"✓ Dashboard copied to {dst}\n"
        else:
            log += f"⚠️ Dashboard not found\n"
        return result.returncode == 0, log

    # ----------------------------------------------------------------------
    # Main orchestration
    # ----------------------------------------------------------------------
    def run_complete_analysis(self, input_path: str, output_dir: str, threads: int = 1,
                              skip: Dict[str, bool] = None, skip_summary: bool = False,
                              skip_viz: bool = False, update_amr_only: bool = False):
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

        fasta_files = self.find_fasta_files(input_path)
        if not fasta_files:
            self.print_error("No FASTA files found. Exiting.")
            return

        self.print_success(f"Found {len(fasta_files)} FASTA files")
        for sub in self.output_dirs.values():
            (final_out / sub).mkdir(exist_ok=True)

        # Show plan
        plan = [
            ("QC", not skip.get('qc', False)),
            ("MLST", not skip.get('mlst', False)),
            ("PAST", not skip.get('past', False)),
            ("ABRicate", not skip.get('abricate', False)),
            ("AMR", not skip.get('amr', False)),
            ("Ultimate Reporter", not skip_summary),
            ("Visualisation", not skip_viz),
        ]
        self.print_header("ANALYSIS PLAN", "Modules to be executed")
        for name, enabled in plan:
            status = f"{Color.BRIGHT_GREEN}✅ ENABLED{Color.RESET}" if enabled else f"{Color.YELLOW}⏸️  SKIPPED{Color.RESET}"
            self._print(f"   {status} - {name}", force=True)

        # First batch: QC, MLST, PAST (parallel)
        batch1_tasks = []
        if not skip.get('qc', False):
            batch1_tasks.append(("QC", self.run_qc))
        if not skip.get('mlst', False):
            batch1_tasks.append(("MLST", self.run_mlst))
        if not skip.get('past', False):
            batch1_tasks.append(("PAST", self.run_past))

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
                        results[name] = (False, f"Exception: {e}")
            for name, _ in batch1_tasks:
                success, log = results.get(name, (False, "No result"))
                self.print_header(f"{name} Analysis")
                # Show first few lines of log (without the "use --verbose" message)
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
            success, log = self.run_abricate(fasta_files, final_out, threads)
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
            success, log = self.run_amr(fasta_files, final_out, threads)
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

        # Ultimate summary
        if not skip_summary:
            self.print_header("ULTIMATE REPORTER", "Gene‑centric Integration")
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
            self.print_info("Skipping ultimate reporter.")

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

        # Cleanup
        for module_name in ['pa_qc_module', 'mlst_module', 'past_module', 'abricate_module', 'amr_module', 'summary_module', 'viz_module']:
            mod_path = self.base_dir / "modules" / module_name
            if mod_path.exists():
                self.cleanup_module(mod_path, fasta_files)

        elapsed = datetime.now() - start_time
        self.print_header("ANALYSIS COMPLETE", f"Time: {str(elapsed).split('.')[0]}")
        self.print_success(f"All results in: {final_out}")
        for subdir in sorted(final_out.iterdir()):
            if subdir.is_dir():
                cnt = len(list(subdir.glob("*")))
                self.print_info(f"  📁 {subdir.name} ({cnt} files)")
        self.display_random_quote()


# -----------------------------------------------------------------------------
# Coloured help function
# -----------------------------------------------------------------------------
def print_colored_help():
    """Print a coloured help message."""
    usage = f"{Color.BRIGHT_CYAN}usage: pseudoscope [OPTIONS]{Color.RESET}"
    description = f"{Color.BRIGHT_WHITE}PseudoScope: Complete P. aeruginosa genomic analysis pipeline{Color.RESET}"

    options = [
        (f"{Color.BRIGHT_GREEN}-i INPUT, --input INPUT{Color.RESET}", f"Input FASTA file, directory, or glob pattern (e.g., \"{Color.BRIGHT_YELLOW}*.fna{Color.RESET}\")"),
        (f"{Color.BRIGHT_GREEN}-o OUTPUT, --output OUTPUT{Color.RESET}", "Output directory for results"),
        (f"{Color.BRIGHT_GREEN}-t THREADS, --threads THREADS{Color.RESET}", f"Number of threads (default: {Color.BRIGHT_YELLOW}2{Color.RESET})"),
        (f"{Color.BRIGHT_GREEN}--quiet{Color.RESET}", "Suppress all non‑error output"),
        (f"{Color.BRIGHT_GREEN}--verbose{Color.RESET}", "Show full command output from modules"),
        (f"{Color.BRIGHT_GREEN}--version{Color.RESET}", "Show version and exit"),
        (f"{Color.BRIGHT_GREEN}--update-amr-db{Color.RESET}", "Update AMRfinderPlus database and exit"),
        (f"{Color.BRIGHT_GREEN}--skip-qc{Color.RESET}", "Skip FASTA QC"),
        (f"{Color.BRIGHT_GREEN}--skip-mlst{Color.RESET}", "Skip MLST"),
        (f"{Color.BRIGHT_GREEN}--skip-past{Color.RESET}", "Skip PAST serotyping"),
        (f"{Color.BRIGHT_GREEN}--skip-abricate{Color.RESET}", "Skip ABRicate"),
        (f"{Color.BRIGHT_GREEN}--skip-amr{Color.RESET}", "Skip AMRfinderPlus"),
        (f"{Color.BRIGHT_GREEN}--skip-summary{Color.RESET}", "Skip ultimate summary report"),
        (f"{Color.BRIGHT_GREEN}--skip-viz{Color.RESET}", "Skip visualisation dashboard"),
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

  # Update AMR database (mandatory before first run)
  pseudoscope --update-amr-db

  # Show version
  pseudoscope --version

  # Quiet mode (minimal output)
  pseudoscope -i "{Color.BRIGHT_CYAN}*.fna{Color.RESET}" -o results --quiet

  # Verbose mode (show all module commands)
  pseudoscope -i "{Color.BRIGHT_CYAN}*.fna{Color.RESET}" -o results --verbose

{Color.BRIGHT_YELLOW}REQUIRED BEFORE FIRST ANALYSIS:{Color.RESET}
  1. {Color.BRIGHT_GREEN}abricate --setupdb{Color.RESET}   (setup ABRicate databases)
  2. {Color.BRIGHT_GREEN}pseudoscope --update-amr-db{Color.RESET}   (download AMRfinderPlus database)

{Color.BRIGHT_YELLOW}SUPPORTED FASTA FORMATS:{Color.RESET} {Color.BRIGHT_CYAN}.fna, .fasta, .fa, .fn{Color.RESET}
    """

    print()
    print(usage)
    print()
    print(description)
    print()
    print(f"{Color.BRIGHT_YELLOW}OPTIONS:{Color.RESET}")
    for opt, desc in options:
        print(f"  {opt:<35} {desc}")
    print(epilog)


def main():
    # Create parser with add_help=False to handle help manually
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-i', '--input', help=argparse.SUPPRESS)
    parser.add_argument('-o', '--output', help=argparse.SUPPRESS)
    parser.add_argument('-t', '--threads', type=int, default=2, help=argparse.SUPPRESS)
    parser.add_argument('--quiet', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--verbose', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--version', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--update-amr-db', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-qc', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-mlst', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-past', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-abricate', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-amr', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-summary', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-viz', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('-h', '--help', action='store_true', help=argparse.SUPPRESS)

    args = parser.parse_args()

    # Show coloured help and exit
    if args.help:
        print_colored_help()
        sys.exit(0)

    if args.version:
        print(f"pseudoscope {__version__}")
        sys.exit(0)

    if args.update_amr_db:
        orch = PseudoScopeOrchestrator(quiet=args.quiet, verbose=args.verbose)
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

    orch = PseudoScopeOrchestrator(quiet=args.quiet, verbose=args.verbose)
    try:
        orch.run_complete_analysis(
            input_path=args.input,
            output_dir=args.output,
            threads=args.threads,
            skip=skip,
            skip_summary=args.skip_summary,
            skip_viz=args.skip_viz
        )
    except KeyboardInterrupt:
        orch.print_error("Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        orch.print_error(f"Critical error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()