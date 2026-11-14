#!/usr/bin/env python3
"""
PseudoScope - Unified Orchestrator for P. aeruginosa Analysis
Author: Brown Beckley <brownbeckley94@gmail.com>
GitHub: bbeckley-hub
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-11-07
Version: 1.3.0

"""

import os
import sys
import glob
import argparse
import subprocess
import shutil
import tempfile
import logging
import traceback
import signal
import atexit
import random
import urllib.request
import zipfile
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

__version__ = "1.3.0"

AUTHOR_NAME = "Brown Beckley"
AUTHOR_EMAIL = "brownbeckley94@gmail.com"
AUTHOR_AFFILIATION = (
    "University of Ghana Medical School - Department of Medical Biochemistry"
)


class Color:
    RESET = '\033[0m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    BRIGHT_RED = '\033[91m'
    BRIGHT_GREEN = '\033[92m'
    BRIGHT_YELLOW = '\033[93m'
    BRIGHT_BLUE = '\033[94m'
    BRIGHT_MAGENTA = '\033[95m'
    BRIGHT_CYAN = '\033[96m'
    BRIGHT_WHITE = '\033[97m'


class ColoredHelpFormatter(argparse.HelpFormatter):
    """Colored help output for readability."""

    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    CYAN = '\033[96m'
    BOLD = '\033[1m'
    RESET = '\033[0m'

    def _format_usage(self, usage, actions, groups, prefix):
        s = super()._format_usage(usage, actions, groups, prefix)
        return f"{self.GREEN}{self.BOLD}usage: {self.RESET}{s}" if s else s

    def _format_action(self, action):
        action_str = super()._format_action(action)
        if not action_str:
            return action_str
        out = []
        for line in action_str.split('\n'):
            if line.strip():
                if line.lstrip().startswith('-'):
                    parts = line.split('  ', 1)
                    if len(parts) == 2:
                        out.append(
                            f"  {self.CYAN}{parts[0].strip()}{self.RESET}  {parts[1]}"
                        )
                    else:
                        out.append(f"  {self.CYAN}{line.strip()}{self.RESET}")
                else:
                    out.append(f"  {self.YELLOW}{line}{self.RESET}")
            else:
                out.append(line)
        return '\n'.join(out)

    def _format_text(self, text):
        return f"{self.BLUE}{text}{self.RESET}" if text else text

    def start_section(self, heading):
        super().start_section(f"{self.BOLD}{self.GREEN}{heading}{self.RESET}")


# ----------------------------------------------------------------------
# Canonical output directory names
# ----------------------------------------------------------------------

DIR_FASTA_QC = "fasta_qc_results"
DIR_MLST = "mlst_results"
DIR_PAST = "past_results"
DIR_ABRICATE = "abricate_results"
DIR_AMR = "amrfinder_results"
DIR_MGE = "mge_results"
DIR_FINAL = "PseudoScope_final_report"
DIR_GENE_CENTRIC = "PSEUDOSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS"
DIR_SAMPLE_CENTRIC = "PSEUDOSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS"
DIR_VIZ = "PSEUDOSCOPE_VISUAL_DASHBOARD"

ALL_MODULE_DIRS = [
    DIR_FASTA_QC, DIR_MLST, DIR_PAST, DIR_ABRICATE, DIR_AMR, DIR_MGE,
]

CITATION_LINES = [
    "Beckley B, Amarh V. PseudoScope: a species-specific bioinformatics suite for rapid and",
    "accessible Pseudomonas aeruginosa genomic analysis. GitHub 2026.",
    "https://github.com/bbeckley-hub/pseudoscope",
]


class PseudoScopeOrchestrator:
    """Main orchestrator for the PseudoScope pipeline.

    Runs analysis modules in isolated temporary directories, collects
    batch summaries into a consolidated final report, and prunes
    intermediate files according to the selected report level.

    Full subprocess stdout/stderr is captured into the run log. On
    non-zero exit codes the complete output is also echoed to the
    console so failures are actionable without opening the log.
    """

    def __init__(self, quiet: bool = False, keep_temp: bool = False,
                 tmpdir: Optional[str] = None):
        self.base_dir = Path(__file__).parent
        self.quiet = quiet
        self.keep_temp = keep_temp
        self.logger = None
        self.user_output_dir: Optional[Path] = None
        self.temp_dirs: List[str] = []
        self.module_times: Dict[str, float] = {}
        self.running = True
        self.report_level = 'summary'

        if tmpdir:
            self.tmp_root = Path(tmpdir)
            self.tmp_root.mkdir(parents=True, exist_ok=True)
        else:
            self.tmp_root = None

        self.setup_colors()
        self.quotes = self._get_scientific_quotes()
        self.quote_colors = [
            Color.BRIGHT_CYAN, Color.BRIGHT_GREEN, Color.BRIGHT_YELLOW,
            Color.BRIGHT_MAGENTA, Color.BRIGHT_BLUE, Color.BRIGHT_RED,
            Color.CYAN, Color.GREEN, Color.YELLOW, Color.MAGENTA,
        ]

        self.output_dirs = {
            'qc': DIR_FASTA_QC,
            'mlst': DIR_MLST,
            'past': DIR_PAST,
            'abricate': DIR_ABRICATE,
            'amr': DIR_AMR,
            'mge': DIR_MGE,
            'summary': DIR_GENE_CENTRIC,
            'sample_centric': DIR_SAMPLE_CENTRIC,
            'viz': DIR_VIZ,
        }

        self.module_subtitles = {
            'QC': 'Sequence Quality Control & Statistics',
            'MLST': 'Multi-Locus Sequence Typing (Pasteur scheme)',
            'PAST': 'O-Antigen Serotyping via pasty',
            'ABRicate': 'Comprehensive Resistance & Virulence Screening',
            'AMR': 'Antimicrobial Resistance Gene Detection',
            'MGE': 'Mobile Genetic Element Profiling',
            'Ultimate Reporter': 'Gene-centric Cross-Genome Integration',
            'Sample-centric Reporter': 'Interactive Isolate-centric Report',
            'Visualisation': 'Interactive Dashboards & Plots',
        }

        self.summary_html_files = {
            'mlst_summary.html':
                f'{DIR_MLST}/mlst_summary.html',
            'past_summary.html':
                f'{DIR_PAST}/past_summary.html',
            'PseudoScope_FASTA_QC_summary.html':
                f'{DIR_FASTA_QC}/PseudoScope_FASTA_QC_summary.html',
            'pseudo_amrfinder_summary_report.html':
                f'{DIR_AMR}/pseudo_amrfinder_summary_report.html',
            'mutation_summary.html':
                f'{DIR_AMR}/mutation_summary.html',
            'pseudo_ncbi_summary_report.html':
                f'{DIR_ABRICATE}/pseudo_ncbi_summary_report.html',
            'pseudo_card_summary_report.html':
                f'{DIR_ABRICATE}/pseudo_card_summary_report.html',
            'pseudo_resfinder_summary_report.html':
                f'{DIR_ABRICATE}/pseudo_resfinder_summary_report.html',
            'pseudo_vfdb_summary_report.html':
                f'{DIR_ABRICATE}/pseudo_vfdb_summary_report.html',
            'pseudo_argannot_summary_report.html':
                f'{DIR_ABRICATE}/pseudo_argannot_summary_report.html',
            'pseudo_plasmidfinder_summary_report.html':
                f'{DIR_ABRICATE}/pseudo_plasmidfinder_summary_report.html',
            'pseudo_megares_summary_report.html':
                f'{DIR_ABRICATE}/pseudo_megares_summary_report.html',
            'pseudo_ecoh_summary_report.html':
                f'{DIR_ABRICATE}/pseudo_ecoh_summary_report.html',
            'pseudo_bacmet2_summary_report.html':
                f'{DIR_ABRICATE}/pseudo_bacmet2_summary_report.html',
            'pseudo_ecoli_vf_summary_report.html':
                f'{DIR_ABRICATE}/pseudo_ecoli_vf_summary_report.html',
            'pseudoscope_mge_summary.html':
                f'{DIR_MGE}/pseudoscope_mge_summary.html',
        }

        self.warnings: List[str] = []
        self.errors: List[str] = []

        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)
        atexit.register(self._cleanup_temp_dirs)

    # ------------------------------------------------------------------
    # Signal handling and temp management
    # ------------------------------------------------------------------
    def _signal_handler(self, signum, frame):
        self.print_warning(f"Received signal {signum}. Cleaning up and exiting...")
        self.running = False
        self._cleanup_temp_dirs()
        sys.exit(1)

    def _cleanup_temp_dirs(self):
        if self.keep_temp:
            self.print_info("Keeping temporary directories (--keep-temp set).")
            return
        for d in list(self.temp_dirs):
            if d and Path(d).exists():
                shutil.rmtree(d, ignore_errors=True)
                if self.logger:
                    self.logger.info(f"Removed temporary directory: {d}")
        self.temp_dirs = []

    def _mkdtemp(self, prefix: str) -> str:
        """Create a temp dir honouring --tmpdir, PSEUDOSCOPE_TMPDIR, TMPDIR."""
        if self.tmp_root:
            return tempfile.mkdtemp(prefix=prefix, dir=str(self.tmp_root))
        env_tmp = os.environ.get("PSEUDOSCOPE_TMPDIR")
        if env_tmp:
            Path(env_tmp).mkdir(parents=True, exist_ok=True)
            return tempfile.mkdtemp(prefix=prefix, dir=env_tmp)
        return tempfile.mkdtemp(prefix=prefix)

    # ------------------------------------------------------------------
    # Console helpers
    # ------------------------------------------------------------------
    def setup_colors(self):
        self.color_info = Color.CYAN
        self.color_success = Color.BRIGHT_GREEN
        self.color_warning = Color.BRIGHT_YELLOW
        self.color_error = Color.BRIGHT_RED

    def setup_logging(self, output_dir: Path):
        log_file = output_dir / "pseudoscope_run.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[logging.FileHandler(log_file, mode='w')],
        )
        self.logger = logging.getLogger("PseudoScope")
        self.logger.info(f"Logging to {log_file}")
        self.user_output_dir = output_dir

    def _print(self, message: str, color: str = Color.RESET, bold: bool = False,
               force: bool = False):
        if self.quiet and not force:
            return
        style = Color.BOLD if bold else ''
        print(f"{style}{color}{message}{Color.RESET}")

    def print_header(self, title: str, subtitle: str = "",
                     show_author: bool = False):
        """Banner-style section header."""
        if not self.quiet:
            width = 80
            top = f"{Color.BOLD}{Color.BRIGHT_BLUE}{'═' * width}{Color.RESET}"
            bottom = f"{Color.BOLD}{Color.BRIGHT_BLUE}{'═' * width}{Color.RESET}"
            print()
            print(top)
            title_pad = max(0, (width - len(title)) // 2)
            print(f"{Color.BOLD}{Color.BRIGHT_CYAN}"
                  f"{' ' * title_pad}{title}"
                  f"{Color.RESET}")
            if subtitle:
                sub_pad = max(0, (width - len(subtitle)) // 2)
                print(f"{Color.DIM}{Color.BRIGHT_WHITE}"
                      f"{' ' * sub_pad}{subtitle}"
                      f"{Color.RESET}")
            if show_author:
                author_line = f"Author: {AUTHOR_NAME} <{AUTHOR_EMAIL}>"
                affil_line = f"{AUTHOR_AFFILIATION}"
                a_pad = max(0, (width - len(author_line)) // 2)
                f_pad = max(0, (width - len(affil_line)) // 2)
                print()
                print(f"{Color.BRIGHT_YELLOW}{' ' * a_pad}{author_line}{Color.RESET}")
                print(f"{Color.DIM}{Color.BRIGHT_WHITE}"
                      f"{' ' * f_pad}{affil_line}"
                      f"{Color.RESET}")
            print(bottom)
            print()
        if self.logger:
            self.logger.info(f"=== {title} ===")
            if subtitle:
                self.logger.info(f"   {subtitle}")
            if show_author:
                self.logger.info(f"Author: {AUTHOR_NAME} <{AUTHOR_EMAIL}>")
                self.logger.info(f"Affiliation: {AUTHOR_AFFILIATION}")

    def display_module_header(self, title: str, subtitle: str = ""):
        """Module header used before each analysis block."""
        if self.quiet:
            return
        width = 80
        line = f"{Color.BOLD}{Color.BRIGHT_MAGENTA}{'─' * width}{Color.RESET}"
        print(line)
        print(f"{Color.BOLD}{Color.BRIGHT_CYAN}  ▸ {title}{Color.RESET}")
        if subtitle:
            print(f"{Color.DIM}{Color.WHITE}    {subtitle}{Color.RESET}")
        print(line)

    def print_info(self, m):
        if not self.quiet:
            print(f"{self.color_info}[i]{Color.RESET} {m}")
        if self.logger:
            self.logger.info(m)

    def print_success(self, m):
        if not self.quiet:
            print(f"{self.color_success}✓{Color.RESET} {m}")
        if self.logger:
            self.logger.info(f"SUCCESS: {m}")

    def print_warning(self, m):
        if not self.quiet:
            print(f"{self.color_warning}⚠{Color.RESET} {m}")
        if self.logger:
            self.logger.warning(m)
        self.warnings.append(m)

    def print_error(self, m):
        print(f"{self.color_error}✗{Color.RESET} {m}")
        if self.logger:
            self.logger.error(m)
        self.errors.append(m)

    def _get_scientific_quotes(self) -> List[Dict[str, str]]:
        return [
            {"quote": "Pseudomonas aeruginosa: the master of adaptation.", "author": "Unknown"},
            {"quote": "The important thing is not to stop questioning.", "author": "Albert Einstein"},
            {"quote": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"},
            {"quote": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"},
            {"quote": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie"},
            {"quote": "The science of today is the technology of tomorrow.", "author": "Edward Teller"},
            {"quote": "In science, there are no shortcuts to truth.", "author": "Karl Popper"},
            {"quote": "Science knows no country, because knowledge belongs to humanity.", "author": "Louis Pasteur"},
            {"quote": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson"},
            {"quote": "Antibiotic resistance is a ticking time bomb.", "author": "Sally Davies"},
            {"quote": "Microbes are the dark matter of the biological world.", "author": "Jack Gilbert"},
            {"quote": "It is not the strongest of the species that survives, but the one most responsive to change.", "author": "Charles Darwin"},
            {"quote": "The greatest enemy of knowledge is not ignorance, it is the illusion of knowledge.", "author": "Stephen Hawking"},
            {"quote": "One microgram of truth can neutralise tonnes of lies.", "author": "Carl Sagan"},
        ]

    def display_random_quote(self):
        if self.quiet or not self.quotes:
            return
        q = random.choice(self.quotes)
        col = random.choice(self.quote_colors)
        print()
        print(f"{Color.DIM}{'─' * 80}{Color.RESET}")
        print(f"{col}{Color.BOLD}   \"{q['quote']}\"{Color.RESET}")
        print(f"{Color.DIM}{Color.BRIGHT_WHITE}   — {q['author']}{Color.RESET}")
        print(f"{Color.DIM}{'─' * 80}{Color.RESET}")
        print()

    def _run_timed(self, label: str, fn, *args, **kwargs):
        """Execute fn, log per-module wall-clock time, return its result."""
        start = datetime.now()
        result = fn(*args, **kwargs)
        elapsed_secs = (datetime.now() - start).total_seconds()
        elapsed_str = str(datetime.now() - start).split('.')[0]
        self.module_times[label] = elapsed_secs
        ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        if not self.quiet:
            print(f"{Color.DIM}[t]{Color.RESET} "
                  f"{Color.DIM}[{ts}]{Color.RESET} "
                  f"{label} finished in {Color.BRIGHT_WHITE}{elapsed_str}{Color.RESET}")
        if self.logger:
            self.logger.info(f"TIMING: {label} = {elapsed_str}")
        return result

    def _log_subprocess_result(self, name: str, result) -> bool:
        """Log full subprocess output; echo to console on failure."""
        if result.stdout:
            self.logger.info(f"=== {name} STDOUT ===\n{result.stdout}")
        if result.stderr:
            self.logger.warning(f"=== {name} STDERR ===\n{result.stderr}")
        if result.returncode != 0:
            print(f"\n{Color.BRIGHT_RED}--- {name} FAILED (rc={result.returncode}) ---"
                  f"{Color.RESET}")
            if result.stdout:
                print(f"{Color.DIM}STDOUT:{Color.RESET}\n{result.stdout}")
            if result.stderr:
                print(f"{Color.DIM}STDERR:{Color.RESET}\n{result.stderr}")
            print(f"{Color.BRIGHT_RED}--- END {name} ---{Color.RESET}\n")
            return False
        return True

    # ------------------------------------------------------------------
    # Module path resolution
    # ------------------------------------------------------------------
    def get_module_path(self, module_name: str) -> Path:
        """Resolve module directory: env override, sys.prefix, then base.

        Checks both PSEUDOSCOPE_<NAME>_MODULE_PATH and PSEUDOSCOPE_<NAME>_PATH
        so pull_mlst_database()'s PSEUDOSCOPE_MLST_MODULE_PATH override is
        honoured by callers of get_module_path("mlst_module").
        """
        for suffix in ("MODULE_PATH", "PATH"):
            env_override = os.environ.get(
                f'PSEUDOSCOPE_{module_name.upper()}_{suffix}'
            )
            if env_override:
                return Path(env_override)
        if hasattr(sys, 'prefix'):
            share_path = (Path(sys.prefix) / "share" / "pseudoscope"
                          / "modules" / module_name)
            if share_path.exists():
                return share_path
        return self.base_dir / "modules" / module_name

    def _get_system_module_path(self, module_name: str) -> Path:
        if hasattr(sys, 'prefix'):
            share_path = (Path(sys.prefix) / "share" / "pseudoscope"
                          / "modules" / module_name)
            if share_path.exists():
                return share_path
        return self.base_dir / "modules" / module_name

    def get_shared_data_root(self) -> Path:
        """Writable user-local data root for read-only module installs."""
        root = Path(os.environ.get(
            "PSEUDOSCOPE_SHARED_DATA",
            Path.home() / ".local" / "share" / "pseudoscope",
        ))
        root.mkdir(parents=True, exist_ok=True)
        return root

    def get_mlst_db_dir(self, for_update: bool = False) -> Path:
        """Resolve the MLST database directory.

        Resolution order:
          1. PSEUDOSCOPE_MLST_DB env override
          2. User-local ~/.local/share/pseudoscope/mlst_module/db
          3. Packaged module db/

        When for_update=True, always return a writable path.
        """
        env_db = os.environ.get("PSEUDOSCOPE_MLST_DB")
        if env_db:
            p = Path(env_db)
            if for_update:
                p.mkdir(parents=True, exist_ok=True)
                return p
            if (p / "pubmlst" / "paeruginosa").exists():
                return p

        user_module = Path(os.environ.get(
            "PSEUDOSCOPE_MLST_MODULE_PATH",
            Path.home() / ".local" / "share" / "pseudoscope" / "mlst_module",
        ))
        user_db = user_module / "db"

        if for_update:
            user_db.mkdir(parents=True, exist_ok=True)
            return user_db

        if (user_db / "pubmlst" / "paeruginosa").exists():
            self.print_info(f"Using user-local MLST database: {user_db}")
            return user_db

        module_db = self.get_module_path("mlst_module") / "db"
        if module_db.exists() and (module_db / "pubmlst" / "paeruginosa").exists():
            return module_db

        self.print_warning(
            "No populated MLST database found. Run `pseudoscope --pull-mlst-db`."
        )
        user_db.mkdir(parents=True, exist_ok=True)
        return user_db

    # ------------------------------------------------------------------
    # MGE database resolution + management
    # ------------------------------------------------------------------
    def get_mge_data_dir(self, for_update: bool = False) -> Optional[Path]:
        """Return the mobileOG-db data directory.

        Resolution order:
          1. PSEUDOSCOPE_MGE_DATA_DIR env override
          2. User-local ~/.local/share/pseudoscope/mge_module/mobileOG-db/
             (populated by `pseudoscope --pull-mge-db`)
          3. Packaged <module>/mobileOG-db/ (Docker / Singularity)
        """
        env_dir = os.environ.get("PSEUDOSCOPE_MGE_DATA_DIR")
        if env_dir:
            env_path = Path(env_dir)
            if for_update:
                env_path.mkdir(parents=True, exist_ok=True)
                return env_path
            if self._mge_db_has_data(env_path):
                return env_path

        user_module = Path(os.environ.get(
            "PSEUDOSCOPE_MGE_MODULE_PATH",
            self.get_shared_data_root() / "mge_module",
        ))
        user_dir = user_module / "mobileOG-db"

        if for_update:
            user_dir.mkdir(parents=True, exist_ok=True)
            return user_dir
        if self._mge_db_has_data(user_dir):
            if self.logger:
                self.logger.info(f"Using user-local mobileOG-db: {user_dir}")
            return user_dir

        module_db = self.get_module_path("mge_module") / "mobileOG-db"
        if self._mge_db_has_data(module_db):
            return module_db

        return None

    @staticmethod
    def _mge_db_has_data(db_dir: Optional[Path]) -> bool:
        """Return True if the directory contains the beatrix DIAMOND index."""
        if db_dir is None or not db_dir.exists():
            return False
        candidates = [
            db_dir / "beatrix-1-6_v1_all" / "mobileOG-db-beatrix-1.6.dmnd",
            db_dir / "mobileOG-db-beatrix-1.6.dmnd",
            db_dir / "beatrix-1-6_v1_all.dmnd",
        ]
        return any(c.exists() for c in candidates)

    def pull_mge_database(self) -> bool:
        """Download the mobileOG-db payload from the GitHub release.

        Docker and Singularity users do not need this — the data is baked
        into the image and get_mge_data_dir() finds it there automatically.
        Conda users run `pseudoscope --pull-mge-db` once after install.
        """
        module_path = self.get_module_path("mge_module")

        if module_path.exists() and os.access(module_path, os.W_OK):
            target_root = module_path
            self.print_info(f"Using writable module directory: {target_root}")
        else:
            target_root = self.get_shared_data_root() / "mge_module"
            target_root.mkdir(parents=True, exist_ok=True)
            self.print_info(
                f"Module dir is read-only; using user-local path: {target_root}"
            )

        repo_url = os.environ.get(
            "PSEUDOSCOPE_MGE_PULL_URL",
            "https://github.com/bbeckley-hub/pseudoscope/releases/download/"
            "v1.3.0/pseudoscope-mge-data-v1.3.0.zip",
        )
        temp_zip = Path(tempfile.mktemp(suffix=".zip"))
        extract_dir = None

        try:
            self.print_info(f"Downloading mobileOG-db from {repo_url}")
            self.print_info(
                "Payload is approximately 300 MB, this may take a few minutes..."
            )
            urllib.request.urlretrieve(repo_url, temp_zip)
            self.print_success("Download complete")

            extract_dir = Path(self._mkdtemp(prefix="pseudoscope_mge_extract_"))
            with zipfile.ZipFile(temp_zip, "r") as zip_ref:
                zip_ref.extractall(extract_dir)
            if self.logger:
                self.logger.info(f"Extracted to {extract_dir}")

            beatrix_dir = None
            for root, dirs, _ in os.walk(extract_dir):
                if "beatrix-1-6_v1_all" in dirs:
                    beatrix_dir = Path(root) / "beatrix-1-6_v1_all"
                    break

            if beatrix_dir is None:
                self.print_error(
                    "Could not find 'beatrix-1-6_v1_all' directory "
                    "in downloaded zip"
                )
                return False

            dmnd = beatrix_dir / "mobileOG-db-beatrix-1.6.dmnd"
            if not dmnd.exists():
                self.print_error(
                    f"DIAMOND index missing from archive: {dmnd.name}"
                )
                return False

            target_beatrix = target_root / "mobileOG-db" / "beatrix-1-6_v1_all"
            if target_beatrix.exists():
                shutil.rmtree(target_beatrix)
            target_beatrix.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(beatrix_dir, target_beatrix)

            dmnd_target = target_beatrix / "mobileOG-db-beatrix-1.6.dmnd"
            if not dmnd_target.exists():
                self.print_error(
                    f"Install verification failed: {dmnd_target} not present"
                )
                return False

            size_mb = dmnd_target.stat().st_size / (1024 * 1024)
            self.print_success(
                f"mobileOG-db installed to {target_root / 'mobileOG-db'} "
                f"(DIAMOND index: {size_mb:.0f} MB)"
            )

            os.environ["PSEUDOSCOPE_MGE_DATA_DIR"] = str(
                target_root / "mobileOG-db"
            )
            os.environ["PSEUDOSCOPE_MGE_MODULE_PATH"] = str(target_root)
            return True

        except Exception as e:
            self.print_error(f"Failed to download/extract mobileOG-db: {e}")
            return False
        finally:
            if temp_zip.exists():
                temp_zip.unlink()
            if extract_dir and extract_dir.exists():
                shutil.rmtree(extract_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Exec-bit restoration
    # ------------------------------------------------------------------
    @staticmethod
    def _restore_exec_bits(root: Path, logger=None) -> int:
        """Restore +x on files under root that should be runnable.

        Zip extraction, some conda package builds, and a few NFS
        configurations strip the Unix executable bit. Without this,
        calling bin/mlst or scripts/*.pl from inside a module fails
        with PermissionError (Errno 13).

        Adds exec bits without touching read/write bits, so an
        already-755 file is unchanged. Returns count of files touched.
        """
        if not root.exists():
            return 0
        executable_exts = {".pl", ".sh", ".py", ".rb"}
        count = 0
        for f in root.rglob("*"):
            if not f.is_file():
                continue
            try:
                parent_is_bin = f.parent.name in ("bin", "scripts")
                is_script = f.suffix.lower() in executable_exts
                if not (parent_is_bin or is_script):
                    continue
                current = f.stat().st_mode
                new_mode = current | 0o111
                if new_mode != current:
                    f.chmod(new_mode)
                    count += 1
            except Exception:
                pass
        if logger and count > 0:
            logger.info(
                f"Restored execute bit on {count} file(s) under {root}"
            )
        return count

    # ------------------------------------------------------------------
    # MLST database: pull from GitHub
    # ------------------------------------------------------------------
    def pull_mlst_database(self) -> bool:
        """Download a pre-built P. aeruginosa MLST scheme from the GitHub repo.

        Writes into the packaged module directory if writable, otherwise
        ~/.local/share/pseudoscope/mlst_module/. After copy, prunes every
        scheme folder except paeruginosa so the end state matches the
        PubMLST update path.
        """
        module_path = self.get_module_path("mlst_module")

        if module_path.exists() and os.access(module_path, os.W_OK):
            target_path = module_path
            self.print_info(f"Using writable module directory: {target_path}")
        else:
            target_path = Path(os.environ.get(
                "PSEUDOSCOPE_MLST_MODULE_PATH",
                Path.home() / ".local" / "share" / "pseudoscope" / "mlst_module",
            ))
            target_path.mkdir(parents=True, exist_ok=True)
            self.print_info(
                f"System MLST module is read-only; using user-local path: {target_path}"
            )

        system_path = self._get_system_module_path("mlst_module")

        repo_url = "https://github.com/bbeckley-hub/mlst/archive/refs/heads/master.zip"
        temp_zip = Path(self._mkdtemp(prefix="pseudoscope_mlst_zip_")) / "mlst.zip"
        extract_dir = None

        try:
            self.print_info("Downloading MLST data from GitHub...")
            urllib.request.urlretrieve(repo_url, temp_zip)
            self.print_info(f"Downloaded {temp_zip.stat().st_size / 1024:.1f} KB")

            extract_dir = Path(self._mkdtemp(prefix="pseudoscope_mlst_extract_"))
            with zipfile.ZipFile(temp_zip, 'r') as zf:
                zf.extractall(extract_dir)

            candidates = list(extract_dir.glob("mlst-*"))
            if not candidates:
                self.print_error(f"Extracted repository root not found in {extract_dir}")
                return False
            src_root = candidates[0]
            if self.logger:
                self.logger.info(f"Extracted repository root: {src_root}")

            for dirname in ('bin', 'db', 'perl5', 'scripts'):
                src = src_root / dirname
                if src.exists():
                    dst = target_path / dirname
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)
                    if self.logger:
                        self.logger.info(f"Copied {dirname}/ to {dst}")
                else:
                    self.print_warning(f"'{dirname}' not found in repository, skipping.")

            self._restore_exec_bits(target_path / "bin", logger=self.logger)
            self._restore_exec_bits(target_path / "scripts", logger=self.logger)

            if target_path != system_path:
                for script in ('p_mlst.py',):
                    src_script = system_path / script
                    if src_script.exists():
                        shutil.copy2(src_script, target_path / script)
                        if self.logger:
                            self.logger.info(f"Copied {script} to {target_path / script}")

            target_pubmlst = target_path / "db" / "pubmlst"
            db_check = target_pubmlst / "paeruginosa"
            if not db_check.exists():
                self.print_error(f"MLST scheme missing after pull: {db_check}")
                self.print_error(
                    "Check that the GitHub repo contains a 'paeruginosa' scheme "
                    "in db/pubmlst/."
                )
                return False

            allele_files = sorted(db_check.glob("*.tfa"))
            if not allele_files:
                self.print_error(f"No allele files (*.tfa) found in {db_check}")
                return False
            self.print_info(f"Verified {len(allele_files)} allele loci in {db_check}")

            removed = 0
            if target_pubmlst.exists():
                for item in target_pubmlst.iterdir():
                    if item.is_dir() and item.name != "paeruginosa":
                        shutil.rmtree(item, ignore_errors=True)
                        removed += 1
                        if self.logger:
                            self.logger.info(f"Pruned unused MLST scheme: {item.name}")
            if removed:
                self.print_info(f"Pruned {removed} unused MLST scheme(s)")

            os.environ['PSEUDOSCOPE_MLST_MODULE_PATH'] = str(target_path)
            os.environ['PSEUDOSCOPE_MLST_DB'] = str(target_path / "db")

            self.print_success(f"MLST module pulled to {target_path}")
            self.print_info(f"Scheme: {db_check}")
            self.print_info(
                "To refresh from PubMLST later (requires API key): "
                "pseudoscope --update-mlst-db"
            )
            return True

        except Exception as e:
            self.print_error(f"Failed to download/extract MLST data: {e}")
            if self.logger:
                self.logger.error(traceback.format_exc())
            return False
        finally:
            if temp_zip.exists():
                temp_zip.unlink()
            if extract_dir and extract_dir.exists():
                shutil.rmtree(extract_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # MLST database: update from PubMLST
    # ------------------------------------------------------------------
    def update_mlst_database(self) -> bool:
        """Refresh the P. aeruginosa scheme from PubMLST via mlstdb."""
        module_path = self.get_module_path("mlst_module")
        if module_path.exists() and os.access(module_path, os.W_OK):
            target_module = module_path
        else:
            target_module = Path(os.environ.get(
                "PSEUDOSCOPE_MLST_MODULE_PATH",
                Path.home() / ".local" / "share" / "pseudoscope" / "mlst_module",
            ))
            target_module.mkdir(parents=True, exist_ok=True)
            self.print_info(
                f"System MLST module read-only; using user-local path: {target_module}"
            )

        db_dir = target_module / "db"
        script_dir = target_module / "bin"
        db_dir.mkdir(parents=True, exist_ok=True)
        script_dir.mkdir(parents=True, exist_ok=True)

        os.environ['PSEUDOSCOPE_MLST_MODULE_PATH'] = str(target_module)
        os.environ['PSEUDOSCOPE_MLST_DB'] = str(db_dir)

        try:
            sys.path.insert(0, str(self.base_dir / "modules" / "mlst_module"))
            from p_mlst import update_mlst_database as _update
        except ImportError as e:
            self.print_error(f"Could not import update_mlst_database from p_mlst.py: {e}")
            return False

        self.print_info(f"Updating MLST database at: {db_dir}")
        return _update(db_dir, script_dir, scheme="paeruginosa")

    # ------------------------------------------------------------------
    # AMR database management
    # ------------------------------------------------------------------
    def check_amr_database(self) -> bool:
        amr_module = self.get_module_path("amr_module")
        db_root = amr_module / "data" / "amrfinder_db"
        if not db_root.exists():
            return False
        for item in db_root.iterdir():
            if item.is_dir() and item.name.startswith("20"):
                if (item / "version.txt").exists() or any(item.glob("*.hmm")):
                    return True
        return False

    def update_amr_database(self, force: bool = False) -> bool:
        """Update AMRfinderPlus database. If force=True, wipe existing first."""
        amr_module = self.get_module_path("amr_module")
        script = amr_module / "p_amrfinder.py"
        if not script.exists():
            self.print_error(f"AMR script not found: {script}")
            return False

        if force:
            db_root = amr_module / "data" / "amrfinder_db"
            if db_root.exists():
                self.print_info(f"Force mode: removing existing database at {db_root}")
                shutil.rmtree(db_root, ignore_errors=True)

        self.print_info("Updating AMRfinderPlus database...")
        cmd = [sys.executable, str(script), "--update-db"]
        result = subprocess.run(cmd, cwd=amr_module, capture_output=True, text=True)
        if not self._log_subprocess_result("AMR database update", result):
            self.print_error("AMR database update failed.")
            return False
        self.print_success("AMR database updated successfully.")
        ver_cmd = [sys.executable, str(script), "--db-version"]
        ver = subprocess.run(ver_cmd, cwd=amr_module, capture_output=True, text=True)
        if ver.returncode == 0 and ver.stdout:
            self.print_info(f"Database version: {ver.stdout.strip()}")
        return True

    def ensure_amr_database(self) -> bool:
        if self.check_amr_database():
            self.print_info("AMR database already present.")
            return True
        self.print_warning("AMR database not found. Attempting automatic update...")
        return self.update_amr_database(force=False)

    # ------------------------------------------------------------------
    # FASTA discovery
    # ------------------------------------------------------------------
    def find_fasta_files(self, input_path: str) -> List[Path]:
        if '*' in input_path or '?' in input_path:
            matched = glob.glob(input_path)
            files = [Path(f) for f in matched if Path(f).is_file()
                     and f.lower().endswith(('.fna', '.fasta', '.fa', '.fsa'))
                     and not Path(f).name.startswith('.')]
            self.print_info(f"Found {len(files)} FASTA files from pattern")
            return sorted(files)

        p = Path(input_path)
        if p.is_file() and p.suffix.lower() in ('.fna', '.fasta', '.fa', '.fsa'):
            return [p]
        if p.is_dir():
            files = []
            for ext in ('*.fna', '*.fasta', '*.fa', '*.fsa',
                        '*.FNA', '*.FASTA', '*.FA', '*.FSA'):
                files.extend(p.glob(ext))
            files = [f for f in files if f.is_file() and not f.name.startswith('.')]
            self.print_info(f"Found {len(files)} FASTA files in {p}")
            return sorted(set(files))
        self.print_error(f"Input path not found: {input_path}")
        return []

    def get_file_pattern(self, fasta_files: List[Path]) -> str:
        if not fasta_files:
            return "*.fna"
        exts = set(f.suffix.lower() for f in fasta_files)
        if len(exts) == 1:
            return f"*{list(exts)[0]}"
        return "*"

    # ------------------------------------------------------------------
    # Copy helpers
    # ------------------------------------------------------------------
    def _copy_module_results(self, src: Path, dst: Path,
                             skip_dirs: Optional[set] = None) -> None:
        """Copy module results to the user output dir, respecting report_level.

        Top-level batch files are always copied. Per-sample folders are
        only copied when report_level == 'full', so large batches never
        incur the copy-then-prune overhead in summary / final mode.
        """
        if not src.exists():
            return
        if dst.exists():
            shutil.rmtree(dst)
        dst.mkdir(parents=True, exist_ok=True)

        skip_dirs = skip_dirs or set()
        for item in src.iterdir():
            if item.name in skip_dirs:
                continue
            if item.is_file():
                shutil.copy2(item, dst / item.name)
            elif item.is_dir() and self.report_level == 'full':
                shutil.copytree(item, dst / item.name)

    # ------------------------------------------------------------------
    # Generic temp-dir module runner
    # ------------------------------------------------------------------
    def run_module_in_temp(
        self,
        module_name: str,
        fasta_files: List[Path],
        cmd_parts: List[str],
        result_subdir: Optional[str] = None,
        src_subdir: Optional[str] = None,
        extra_files: Optional[List[Tuple[str, str]]] = None,
        ignore_patterns: Optional[List[str]] = None,
        env_overrides: Optional[Dict[str, str]] = None,
    ) -> Tuple[bool, str]:
        """Copy module to a temp dir, run it there, copy results back.

        result_subdir : destination folder name inside the user output
                        directory (also used as the source folder name
                        when src_subdir is not provided).
        src_subdir    : name of the folder the module actually writes to
                        inside its temp work dir. Used when the module's
                        internal folder name differs from the friendlier
                        output name — e.g. ABRicate writes to
                        pseudo_abricate_results/ but we publish it as
                        abricate_results/.

        cmd_parts[2:] are passed verbatim to the module script. Most
        p_*.py modules expect a single glob pattern argument (e.g.
        "*.fna") and expand it internally. MGEFinder.py is an exception
        and receives actual filenames — see run_mge().
        """
        module_orig = self.get_module_path(module_name)
        if not module_orig.exists():
            return False, f"Module directory not found: {module_orig}"

        temp_dir = self._mkdtemp(prefix=f"pseudoscope_{module_name}_")
        self.temp_dirs.append(temp_dir)
        if self.logger:
            self.logger.info(f"Temporary directory for {module_name}: {temp_dir}")

        try:
            copytree_kwargs = {'dirs_exist_ok': True}
            if ignore_patterns:
                copytree_kwargs['ignore'] = shutil.ignore_patterns(*ignore_patterns)
            shutil.copytree(
                module_orig, Path(temp_dir) / module_name, **copytree_kwargs
            )
            self._restore_exec_bits(Path(temp_dir) / module_name, logger=self.logger)

            for f in fasta_files:
                shutil.copy2(f, Path(temp_dir) / f.name)

            script_name = Path(cmd_parts[1]).name if len(cmd_parts) > 1 else None
            if not script_name:
                return False, f"Could not determine script name from: {cmd_parts}"
            script_path = Path(temp_dir) / module_name / script_name
            if not script_path.exists():
                py_files = list((Path(temp_dir) / module_name).glob("*.py"))
                if py_files:
                    script_path = py_files[0]
                else:
                    return False, f"Could not locate script in {module_name}"

            abs_cmd = [cmd_parts[0], str(script_path)] + cmd_parts[2:]

            env = os.environ.copy()
            if env_overrides:
                env.update(env_overrides)

            if self.logger:
                self.logger.info(f"Running {module_name}: {' '.join(abs_cmd)}")
            result = subprocess.run(abs_cmd, cwd=temp_dir, capture_output=True,
                                    text=True, env=env)

            if not self._log_subprocess_result(module_name, result):
                return False, f"Module {module_name} failed (rc={result.returncode})"

            if result_subdir:
                actual_src = Path(temp_dir) / (src_subdir or result_subdir)
                if actual_src.exists():
                    dst = self.user_output_dir / result_subdir
                    self._copy_module_results(actual_src, dst)
                    if self.logger:
                        self.logger.info(f"Results copied to {dst}")
                else:
                    if self.logger:
                        self.logger.warning(
                            f"Expected result subdir not found: {actual_src}"
                        )

            if extra_files:
                for src_rel, dst_name in extra_files:
                    src = Path(temp_dir) / src_rel
                    if src.exists():
                        dst = self.user_output_dir / dst_name
                        shutil.copy2(src, dst)
                        if self.logger:
                            self.logger.info(f"Copied {dst_name} to output directory")
                    else:
                        if self.logger:
                            self.logger.warning(f"Extra file not found: {src}")

            return True, f"{module_name} completed successfully"

        except Exception as e:
            if self.logger:
                self.logger.error(f"Exception in {module_name}: {e}\n{traceback.format_exc()}")
            return False, f"Exception: {e}"
        finally:
            if not self.keep_temp and self.running:
                shutil.rmtree(temp_dir, ignore_errors=True)
                if self.logger:
                    self.logger.info(f"Removed temporary directory: {temp_dir}")

    # ------------------------------------------------------------------
    # Module runners
    # ------------------------------------------------------------------
    def run_qc(self, fasta_files, final_out, threads):
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "p_qc.py", pattern, "--cpus", str(threads)]
        return self.run_module_in_temp("pa_qc_module", fasta_files, cmd, DIR_FASTA_QC)

    def run_mlst(self, fasta_files, final_out, threads):
        pattern = self.get_file_pattern(fasta_files)
        db_dir = self.get_mlst_db_dir()
        cmd = [sys.executable, "p_mlst.py", "-i", pattern,
               "-o", "mlst_results", "--batch"]
        return self.run_module_in_temp(
            "mlst_module", fasta_files, cmd, DIR_MLST,
            env_overrides={"MLST_DATADIR": str(db_dir)},
        )

    def run_past(self, fasta_files, final_out, threads,
                 min_pident=None, min_coverage=None):
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "p_past.py", "-i", pattern,
               "-o", "past_results", "--cpus", str(threads)]
        if min_pident is not None:
            cmd.extend(["--min-pident", str(min_pident)])
        if min_coverage is not None:
            cmd.extend(["--min-coverage", str(min_coverage)])
        return self.run_module_in_temp("past_module", fasta_files, cmd, DIR_PAST)

    def run_abricate(self, fasta_files, final_out, threads,
                     minid=None, mincov=None):
        if not self.check_abricate_databases():
            return False, "ABRicate databases not found. Run 'abricate --setupdb' first."
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "p_abricate.py", pattern, "--cpus", str(threads)]
        if minid is not None:
            cmd.extend(["--minid", str(minid)])
        if mincov is not None:
            cmd.extend(["--mincov", str(mincov)])
        return self.run_module_in_temp(
            "abricate_module", fasta_files, cmd,
            result_subdir=DIR_ABRICATE,
            src_subdir="pseudo_abricate_results",
        )

    def run_amr(self, fasta_files, final_out, threads,
                min_identity=None, min_coverage=None,
                skip_mutations=False, force_update=False):
        if not self.check_amr_database():
            if not self.ensure_amr_database():
                return False, "AMR database not found and could not be auto-updated."
        if force_update:
            self.print_info("Forcing AMR database update before analysis...")
            self.update_amr_database(force=True)

        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "p_amrfinder.py", pattern]
        if min_identity is not None:
            cmd.extend(["--min-identity", str(min_identity)])
        if min_coverage is not None:
            cmd.extend(["--min-coverage", str(min_coverage)])
        if skip_mutations:
            cmd.append("--skip-mutations")
        return self.run_module_in_temp(
            "amr_module", fasta_files, cmd,
            result_subdir=DIR_AMR,
            src_subdir="pseudo_amrfinder_results",
        )

    def run_mge(self, fasta_files, final_out, threads):
        """Run the mobileOG-based MGE profiler.

        The resolved mobileOG-db payload is copied into the temp working
        directory (no symlinks) so the pipeline works on filesystems that
        do not support symlinks. Actual filenames are passed to
        MGEFinder.py — it does not expand glob patterns.
        """
        module_orig = self.get_module_path("mge_module")
        if not module_orig.exists():
            return False, f"mge_module not found: {module_orig}"

        mge_data = self.get_mge_data_dir()
        if mge_data is None:
            return False, (
                "mobileOG-db not found. Run `pseudoscope --pull-mge-db` once, "
                "or use the Docker/Singularity image (which has it bundled)."
            )
        if self.logger:
            self.logger.info(f"Resolved mobileOG-db: {mge_data}")

        temp_dir = Path(self._mkdtemp(prefix="pseudoscope_mge_"))
        self.temp_dirs.append(str(temp_dir))
        if self.logger:
            self.logger.info(f"Temporary directory for mge: {temp_dir}")

        try:
            shutil.copytree(
                module_orig, temp_dir, dirs_exist_ok=True,
                ignore=shutil.ignore_patterns(
                    "mge_results", "annotations",
                    "*.fna", "*.fasta", "*.fa", "*.fsa",
                    "mobileOG-db", "__pycache__", "*.pyc",
                ),
            )
            self._restore_exec_bits(temp_dir, logger=self.logger)

            beatrix_src = mge_data / "beatrix-1-6_v1_all"
            if not beatrix_src.exists() and mge_data.name == "beatrix-1-6_v1_all":
                beatrix_src = mge_data
            if not beatrix_src.exists():
                return False, f"beatrix-1-6_v1_all not found under {mge_data}"

            beatrix_dst = temp_dir / "mobileOG-db" / "beatrix-1-6_v1_all"
            beatrix_dst.parent.mkdir(parents=True, exist_ok=True)
            if beatrix_dst.exists():
                shutil.rmtree(beatrix_dst)
            shutil.copytree(beatrix_src, beatrix_dst)
            if self.logger:
                self.logger.info(f"Copied {beatrix_src} -> {beatrix_dst}")

            for f in fasta_files:
                shutil.copy2(f, temp_dir / f.name)

            self.print_info(
                f"Running MGE analysis on {len(fasta_files)} genome(s)..."
            )
            script = temp_dir / "MGEFinder.py"
            if not script.exists():
                return False, "MGEFinder.py not found in mge_module"

            mge_inputs = [f.name for f in fasta_files]
            cmd = ([sys.executable, str(script), "-i"] + mge_inputs +
                   ["-d", "mge_results",
                    "-db_dir", "mobileOG-db/beatrix-1-6_v1_all"])
            if self.logger:
                self.logger.info(
                    f"Running MGE: MGEFinder.py -i <{len(mge_inputs)} filenames> "
                    f"-d mge_results -db_dir mobileOG-db/beatrix-1-6_v1_all"
                )

            result = subprocess.run(cmd, cwd=temp_dir, capture_output=True,
                                    text=True)
            if not self._log_subprocess_result("MGE", result):
                return False, f"MGE analysis failed (rc={result.returncode})"

            src = temp_dir / "mge_results"
            if not src.exists():
                return False, f"MGE results directory not found: {src}"

            dst = self.user_output_dir / DIR_MGE
            self._copy_module_results(src, dst, skip_dirs={"annotations"})

            expected = [
                "pseudoscope_mge_summary.tsv",
                "pseudoscope_mge_hits.tsv",
                "pseudoscope_mge_summary.json",
                "pseudoscope_mge_summary.html",
            ]
            missing = [f for f in expected if not (dst / f).exists()]
            if missing:
                self.print_warning(
                    f"MGE ran but these summary files are missing: "
                    f"{', '.join(missing)}. Reporters may not include MGE data."
                )

            self.print_success("MGE profiling completed successfully.")
            return True, "MGE completed successfully"

        except Exception as e:
            if self.logger:
                self.logger.error(traceback.format_exc())
            return False, f"MGE exception: {e}"
        finally:
            if not self.keep_temp and self.running:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def run_summary(self, final_out):
        """Run the gene-centric ultimate reporter."""
        module_orig = self.get_module_path("gene_centric_module")
        if not module_orig.exists():
            return False, f"Summary module not found: {module_orig}"

        temp_dir = self._mkdtemp(prefix="pseudoscope_summary_")
        self.temp_dirs.append(temp_dir)
        if self.logger:
            self.logger.info(f"Temporary directory for summary: {temp_dir}")

        try:
            shutil.copytree(module_orig, Path(temp_dir) / "gene_centric_module",
                            dirs_exist_ok=True)
            self._restore_exec_bits(Path(temp_dir) / "gene_centric_module",
                                    logger=self.logger)

            for target_name, source_rel in self.summary_html_files.items():
                src = final_out / source_rel
                if src.exists():
                    shutil.copy2(src, Path(temp_dir) / target_name)
                    if self.logger:
                        self.logger.info(f"Copied {target_name}")
                else:
                    if self.logger:
                        self.logger.warning(
                            f"Required summary file not found: {src}"
                        )

            script = Path(temp_dir) / "gene_centric_module" / "p_gene_centric.py"
            if not script.exists():
                return False, f"Summary script not found: {script}"

            cmd = [sys.executable, str(script), "-i", "."]
            if self.logger:
                self.logger.info(f"Running summary: {' '.join(cmd)}")
            result = subprocess.run(cmd, cwd=temp_dir,
                                    capture_output=True, text=True)
            if not self._log_subprocess_result("Summary", result):
                return False, f"Summary failed (rc={result.returncode})"

            src = Path(temp_dir) / DIR_GENE_CENTRIC
            dst = final_out / DIR_GENE_CENTRIC
            if src.exists():
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(src, dst)
                if self.logger:
                    self.logger.info(f"Reports copied to {dst}")
            else:
                return False, "Ultimate gene-centric reports directory not found"
            return True, "Summary completed successfully"

        except Exception as e:
            if self.logger:
                self.logger.error(f"Summary exception: {e}\n{traceback.format_exc()}")
            return False, f"Exception: {e}"
        finally:
            if not self.keep_temp and self.running:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def run_sample_centric(self, final_out):
        """Run the sample-centric interactive reporter."""
        module_orig = self.get_module_path("sample_centric_module")
        if not module_orig.exists():
            return False, f"Sample-centric module not found: {module_orig}"

        temp_dir = self._mkdtemp(prefix="pseudoscope_sample_centric_")
        self.temp_dirs.append(temp_dir)
        if self.logger:
            self.logger.info(
                f"Temporary directory for sample-centric: {temp_dir}"
            )

        try:
            shutil.copytree(module_orig, Path(temp_dir) / "sample_centric_module",
                            dirs_exist_ok=True)
            self._restore_exec_bits(Path(temp_dir) / "sample_centric_module",
                                    logger=self.logger)

            copy_specs = [
                (DIR_FASTA_QC, "PseudoScope_FASTA_QC_summary.*"),
                (DIR_MLST, "mlst_summary.*"),
                (DIR_PAST, "past_summary.html"),
                (DIR_AMR, "pseudo_amrfinder_summary.*"),
                (DIR_AMR, "mutation_summary.tsv"),
                (DIR_ABRICATE, "pseudo_*_abricate_summary.tsv"),
                (DIR_ABRICATE, "pseudo_*_summary_report.html"),
                (DIR_MGE, "pseudoscope_mge_hits.tsv"),
                (DIR_MGE, "pseudoscope_mge_summary.html"),
                (DIR_MGE, "pseudoscope_mge_summary.json"),
            ]
            for sub, pattern in copy_specs:
                src_dir = final_out / sub
                if not src_dir.exists():
                    continue
                for f in src_dir.glob(pattern):
                    if f.is_file():
                        shutil.copy2(f, Path(temp_dir) / f.name)
                        if self.logger:
                            self.logger.info(f"Copied {f.name}")

            script = (Path(temp_dir) / "sample_centric_module"
                      / "p_sample_centric.py")
            if not script.exists():
                return False, f"Sample-centric script not found: {script}"

            cmd = [sys.executable, str(script), "-i", "."]
            if self.logger:
                self.logger.info(f"Running sample-centric: {' '.join(cmd)}")
            result = subprocess.run(cmd, cwd=temp_dir,
                                    capture_output=True, text=True)
            if not self._log_subprocess_result("Sample-centric", result):
                return False, f"Sample-centric failed (rc={result.returncode})"

            src = Path(temp_dir) / DIR_SAMPLE_CENTRIC
            dst = final_out / DIR_SAMPLE_CENTRIC
            if src.exists():
                if dst.exists():
                    shutil.rmtree(dst)
                shutil.copytree(src, dst)
                if self.logger:
                    self.logger.info(f"Sample-centric reports copied to {dst}")
            else:
                return False, "Sample-centric reports directory not found"
            return True, "Sample-centric completed successfully"

        except Exception as e:
            if self.logger:
                self.logger.error(
                    f"Sample-centric exception: {e}\n{traceback.format_exc()}"
                )
            return False, f"Exception: {e}"
        finally:
            if not self.keep_temp and self.running:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def run_viz(self, final_out):
        """Run the visualisation module (dashboard + static plots)."""
        module_orig = self.get_module_path("viz_module")
        if not module_orig.exists():
            return False, f"Visualisation module not found: {module_orig}"

        temp_dir = self._mkdtemp(prefix="pseudoscope_viz_")
        self.temp_dirs.append(temp_dir)
        if self.logger:
            self.logger.info(f"Temporary directory for visualisation: {temp_dir}")

        try:
            shutil.copytree(module_orig, Path(temp_dir) / "viz_module",
                            dirs_exist_ok=True)
            self._restore_exec_bits(Path(temp_dir) / "viz_module",
                                    logger=self.logger)

            ultimate_dir = final_out / DIR_GENE_CENTRIC
            if ultimate_dir.exists():
                for csv_file in ultimate_dir.glob("*.csv"):
                    shutil.copy2(csv_file, Path(temp_dir) / csv_file.name)
                    if self.logger:
                        self.logger.info(
                            f"Copied {csv_file.name} to viz temp dir"
                        )

            qc_dir = final_out / DIR_FASTA_QC
            if qc_dir.exists():
                for name in ("PseudoScope_FASTA_QC_summary.tsv",
                             "PseudoScope_FASTA_QC_summary.html"):
                    p = qc_dir / name
                    if p.exists():
                        shutil.copy2(p, Path(temp_dir) / p.name)

            script = Path(temp_dir) / "viz_module" / "p_visualizer.py"
            if not script.exists():
                return False, f"Visualizer script not found: {script}"

            cmd = [sys.executable, str(script), "-i", "."]
            if self.logger:
                self.logger.info(f"Running visualisation: {' '.join(cmd)}")
            result = subprocess.run(cmd, cwd=temp_dir,
                                    capture_output=True, text=True)
            if not self._log_subprocess_result("Visualisation", result):
                return False, f"Visualisation failed (rc={result.returncode})"

            dashboard = Path(temp_dir) / "genius_pseudomonas_visual_dashboard.html"
            dst_dir = final_out / DIR_VIZ
            dst_dir.mkdir(parents=True, exist_ok=True)
            if dashboard.exists():
                shutil.copy2(dashboard, dst_dir / dashboard.name)
                if self.logger:
                    self.logger.info(
                        f"Dashboard copied to {dst_dir / dashboard.name}"
                    )
            else:
                return False, "Dashboard HTML not found after visualisation run"
            return True, "Visualisation completed successfully"

        except Exception as e:
            if self.logger:
                self.logger.error(
                    f"Visualisation exception: {e}\n{traceback.format_exc()}"
                )
            return False, f"Exception: {e}"
        finally:
            if not self.keep_temp and self.running:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def check_abricate_databases(self) -> bool:
        try:
            result = subprocess.run(["abricate", "--list"],
                                    capture_output=True, text=True)
            return result.returncode == 0 and "ncbi" in result.stdout.lower()
        except FileNotFoundError:
            return False

    # ------------------------------------------------------------------
    # Final report consolidation + pruning
    # ------------------------------------------------------------------
    def copy_summary_results_to_final_directory(self, output_dir: Path) -> None:
        """Copy the top-level deliverables into PseudoScope_final_report/.

        Contents:
          * Gene-centric ultimate reports
          * Sample-centric interactive reports
          * Visualisation dashboard
          * MGE summary files for quick reference
        """
        try:
            self.print_info("Consolidating results into final report directory...")
            final_dir = output_dir / DIR_FINAL
            final_dir.mkdir(parents=True, exist_ok=True)

            for sub in [DIR_GENE_CENTRIC, DIR_SAMPLE_CENTRIC, DIR_VIZ]:
                src = output_dir / sub
                if src.exists() and src.is_dir():
                    dst = final_dir / sub
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)
                    self.print_info(f"  ▸ Copied {sub}/")

            mge_dir = output_dir / DIR_MGE
            if mge_dir.exists():
                for name in ("pseudoscope_mge_summary.html",
                             "pseudoscope_mge_hits.tsv",
                             "pseudoscope_mge_summary.json"):
                    src = mge_dir / name
                    if src.exists():
                        shutil.copy2(src, final_dir / name)
                        self.print_info(f"  ▸ Copied {name}")

            self.print_success(f"Final report ready: {final_dir}")

        except Exception as e:
            self.print_error(f"Error consolidating final report: {e}")

    def _cleanup_top_level_duplicates(self, output_path: Path) -> None:
        """Remove copies that were consolidated into PseudoScope_final_report/."""
        removed = 0
        for dirname in [DIR_GENE_CENTRIC, DIR_SAMPLE_CENTRIC, DIR_VIZ]:
            d = output_path / dirname
            if d.exists() and d.is_dir():
                shutil.rmtree(d, ignore_errors=True)
                removed += 1
        for item in output_path.iterdir():
            if item.is_file() and item.suffix.lower() == '.html':
                try:
                    item.unlink()
                    removed += 1
                except Exception:
                    pass
        if removed:
            self._print(
                f"{Color.DIM}Removed {removed} duplicate top-level item(s){Color.RESET}",
                force=True,
            )

    def _finalize_output_by_report_level(self, output_path: Path) -> None:
        """Prune intermediate files based on the selected report level.

        final   : Only PseudoScope_final_report/ is kept.
        summary : PseudoScope_final_report/ + batch files in module dirs.
                  Per-sample folders are never copied in the first place.
        full    : Everything is kept, including per-sample HTMLs.
        """
        if self.report_level == 'full':
            self.print_info("Report level: full (all intermediate files kept)")
            return

        final_dir = output_path / DIR_FINAL
        final_exists = final_dir.exists() and any(final_dir.iterdir())
        if not final_exists:
            self.print_warning(
                "Final report directory is missing or empty - skipping prune "
                "to avoid losing all output."
            )
            return

        if self.report_level == 'final':
            removed = 0
            for module_dir_name in ALL_MODULE_DIRS:
                module_dir = output_path / module_dir_name
                if not module_dir.exists():
                    continue
                shutil.rmtree(module_dir, ignore_errors=True)
                removed += 1
            self.print_success(
                f"Report level: final - removed {removed} module directories"
            )
            return

        if self.report_level == 'summary':
            removed = 0
            for module_dir_name in ALL_MODULE_DIRS:
                module_dir = output_path / module_dir_name
                if not module_dir.exists():
                    continue
                for item in list(module_dir.iterdir()):
                    if item.is_dir():
                        shutil.rmtree(item, ignore_errors=True)
                        removed += 1
            if removed:
                self.print_success(
                    f"Report level: summary - removed {removed} per-sample folders"
                )
            else:
                self.print_success(
                    "Report level: summary - only batch summaries retained"
                )

    # ------------------------------------------------------------------
    # Full pipeline
    # ------------------------------------------------------------------
    def run_complete_analysis(
        self,
        input_path: str,
        output_dir: str,
        threads: int = 1,
        skip: Optional[Dict[str, bool]] = None,
        skip_summary: bool = False,
        skip_sample_centric: bool = False,
        skip_viz: bool = False,
        update_amr_only: bool = False,
        update_mlst_only: bool = False,
        pull_mlst_only: bool = False,
        pull_mge_only: bool = False,
        amr_force_update: bool = False,
        clean_output: bool = False,
        abricate_minid: Optional[int] = None,
        abricate_mincov: Optional[int] = None,
        past_min_pident: Optional[int] = None,
        past_min_coverage: Optional[int] = None,
        amr_min_identity: Optional[float] = None,
        amr_min_coverage: Optional[float] = None,
        amr_skip_mutations: bool = False,
        report_level: str = 'summary',
    ) -> bool:
        """Execute the full PseudoScope pipeline."""
        self.report_level = report_level

        if update_amr_only:
            temp_log = Path(self._mkdtemp(prefix="pseudoscope_db_"))
            self.temp_dirs.append(str(temp_log))
            self.setup_logging(temp_log)
            ok = self.update_amr_database(force=amr_force_update)
            if not self.keep_temp:
                shutil.rmtree(temp_log, ignore_errors=True)
            return ok
        if update_mlst_only:
            temp_log = Path(self._mkdtemp(prefix="pseudoscope_mlst_db_"))
            self.temp_dirs.append(str(temp_log))
            self.setup_logging(temp_log)
            ok = self.update_mlst_database()
            if not self.keep_temp:
                shutil.rmtree(temp_log, ignore_errors=True)
            return ok
        if pull_mlst_only:
            temp_log = Path(self._mkdtemp(prefix="pseudoscope_mlst_pull_"))
            self.temp_dirs.append(str(temp_log))
            self.setup_logging(temp_log)
            ok = self.pull_mlst_database()
            if not self.keep_temp:
                shutil.rmtree(temp_log, ignore_errors=True)
            return ok
        if pull_mge_only:
            temp_log = Path(self._mkdtemp(prefix="pseudoscope_mge_pull_"))
            self.temp_dirs.append(str(temp_log))
            self.setup_logging(temp_log)
            ok = self.pull_mge_database()
            if not self.keep_temp:
                shutil.rmtree(temp_log, ignore_errors=True)
            return ok

        if skip is None:
            skip = {}

        start = datetime.now()
        self.print_header(
            "PSEUDOSCOPE",
            "P. aeruginosa Genomic Analysis Pipeline",
            show_author=True,
        )
        self.print_info(f"Version {__version__}")

        final_out = Path(output_dir)
        if clean_output and final_out.exists():
            self.print_info(f"Cleaning output directory: {final_out}")
            shutil.rmtree(final_out)
        final_out.mkdir(parents=True, exist_ok=True)
        self.setup_logging(final_out)

        fasta_files = self.find_fasta_files(input_path)
        if not fasta_files:
            self.print_error("No FASTA files found. Exiting.")
            return False

        self.print_success(f"Found {len(fasta_files)} FASTA files")
        for sub in self.output_dirs.values():
            (final_out / sub).mkdir(exist_ok=True)

        if not skip.get('mge', False) and self.get_mge_data_dir() is None:
            self.print_warning(
                "MGE profiling is enabled but mobileOG-db is not installed. "
                "Run `pseudoscope --pull-mge-db` first, or add --skip-mge."
            )

        plan = [
            ("QC",                                  not skip.get('qc', False)),
            ("MLST",                                not skip.get('mlst', False)),
            ("PAST",                                not skip.get('past', False)),
            ("ABRicate",                            not skip.get('abricate', False)),
            ("AMR",                                 not skip.get('amr', False)),
            ("MGE",                                 not skip.get('mge', False)),
            ("Ultimate Reporter (Gene-centric)",    not skip_summary),
            ("Sample-centric Reporter",             not skip_sample_centric and not skip_summary),
            ("Visualisation",                       not skip_viz and not skip_summary),
        ]
        self.print_header("ANALYSIS PLAN", "Modules to be executed")
        for name, enabled in plan:
            status = (f"{Color.BRIGHT_GREEN}✅ ENABLED{Color.RESET}" if enabled
                      else f"{Color.YELLOW}⏸️  SKIPPED{Color.RESET}")
            self._print(f"   {status} - {name}", force=True)
            if self.logger:
                self.logger.info(f"Plan: {name} {'ENABLED' if enabled else 'SKIPPED'}")
        self.print_info(f"Report level: {self.report_level}")

        # ---------- Batch 1 (parallel): QC, MLST, PAST ----------
        batch1 = []
        if not skip.get('qc', False):
            batch1.append(("QC", self.run_qc))
        if not skip.get('mlst', False):
            batch1.append(("MLST", self.run_mlst))
        if not skip.get('past', False):
            batch1.append(
                ("PAST",
                 lambda f, o, t: self.run_past(f, o, t, past_min_pident, past_min_coverage))
            )

        if batch1:
            self.display_module_header(
                "PARALLEL BATCH 1",
                f"QC, MLST, PAST — {len(batch1)} module(s) running concurrently",
            )
            self.print_info(f"Running {len(batch1)} analyses in parallel...")
            results = {}
            with ThreadPoolExecutor(max_workers=len(batch1)) as executor:
                futures = {
                    executor.submit(self._run_timed, name, run,
                                    fasta_files, final_out, threads): name
                    for name, run in batch1
                }
                for fut in as_completed(futures):
                    name = futures[fut]
                    try:
                        results[name] = fut.result()
                    except Exception as e:
                        if self.logger:
                            self.logger.error(f"Exception in {name}: {e}")
                        results[name] = (False, f"Exception: {e}")
            print()
            for name, _ in batch1:
                success, log = results.get(name, (False, "No result"))
                if success:
                    self.print_success(f"{name} completed")
                else:
                    self.print_error(f"{name} failed")
            self.display_random_quote()

        # ---------- Batch 2: ABRicate ----------
        if not skip.get('abricate', False):
            self.display_module_header("ABRICATE", self.module_subtitles['ABRicate'])
            success, log = self._run_timed(
                "ABRicate", self.run_abricate,
                fasta_files, final_out, threads,
                abricate_minid, abricate_mincov,
            )
            if success:
                self.print_success("ABRicate completed")
            else:
                self.print_error("ABRicate failed")
            self.display_random_quote()

        # ---------- Batch 3: AMR ----------
        if not skip.get('amr', False):
            self.display_module_header("AMRFinderPlus", self.module_subtitles['AMR'])
            success, log = self._run_timed(
                "AMR", self.run_amr,
                fasta_files, final_out, threads,
                amr_min_identity, amr_min_coverage,
                amr_skip_mutations, amr_force_update,
            )
            if success:
                self.print_success("AMR completed")
            else:
                self.print_error("AMR failed")
            self.display_random_quote()

        # ---------- Batch 4: MGE ----------
        if not skip.get('mge', False):
            self.display_module_header("MGE PROFILING", self.module_subtitles['MGE'])
            success, log = self._run_timed(
                "MGE", self.run_mge, fasta_files, final_out, threads,
            )
            if success:
                self.print_success("MGE completed")
            else:
                self.print_warning("MGE had issues")
            self.display_random_quote()

        # ---------- Gene-centric summary ----------
        if not skip_summary:
            self.display_module_header(
                "ULTIMATE REPORTER (Gene-centric)",
                self.module_subtitles['Ultimate Reporter'],
            )
            success, log = self._run_timed(
                "Ultimate Reporter", self.run_summary, final_out,
            )
            if success:
                self.print_success("Ultimate reporter completed")
            else:
                self.print_warning("Ultimate reporter had issues")
            self.display_random_quote()
        else:
            self.print_info("Skipping gene-centric reporter.")

        # ---------- Sample-centric ----------
        if not skip_sample_centric and not skip_summary:
            self.display_module_header(
                "SAMPLE-CENTRIC REPORTER",
                self.module_subtitles['Sample-centric Reporter'],
            )
            success, log = self._run_timed(
                "Sample-centric Reporter", self.run_sample_centric, final_out,
            )
            if success:
                self.print_success("Sample-centric reporter completed")
            else:
                self.print_warning("Sample-centric reporter had issues")
            self.display_random_quote()
        elif skip_summary and not skip_sample_centric:
            self.print_warning(
                "Sample-centric reporter requires the gene-centric summary CSVs/JSON."
            )
            self.print_warning(
                "Since --skip-summary was used, sample-centric will be skipped."
            )
        elif skip_sample_centric:
            self.print_info("Skipping sample-centric reporter.")

        # ---------- Visualisation ----------
        if not skip_viz and not skip_summary:
            self.display_module_header(
                "VISUALISATION", self.module_subtitles['Visualisation'],
            )
            success, log = self._run_timed(
                "Visualisation", self.run_viz, final_out,
            )
            if success:
                self.print_success("Visualisation completed")
            else:
                self.print_warning("Visualisation had issues")
            self.display_random_quote()
        elif skip_summary:
            self.print_info("Skipping visualisation because summary was skipped.")
        elif skip_viz:
            self.print_info("Skipping visualisation.")

        # ---------- Consolidation + pruning ----------
        self.copy_summary_results_to_final_directory(final_out)
        self._cleanup_top_level_duplicates(final_out)
        self._finalize_output_by_report_level(final_out)

        elapsed = datetime.now() - start
        elapsed_str = str(elapsed).split('.')[0]

        self.print_header(
            "ANALYSIS COMPLETE",
            f"Total time: {elapsed_str}  |  "
            f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            show_author=True,
        )

        if self.module_times:
            self.display_module_header(
                "TIMING SUMMARY", "Per-module wall-clock time"
            )
            for label, secs in sorted(self.module_times.items(), key=lambda x: -x[1]):
                mins = int(secs // 60)
                rem = secs % 60
                time_str = f"{mins}m {rem:.1f}s" if mins else f"{rem:.1f}s"
                self._print(
                    f"   {Color.BRIGHT_CYAN}{label:<28}{Color.RESET} "
                    f"{Color.BRIGHT_WHITE}{time_str}{Color.RESET}",
                    force=True,
                )
            self._print(f"   {Color.DIM}{'─' * 40}{Color.RESET}", force=True)
            self._print(
                f"   {Color.BOLD}Total (pipeline){Color.RESET}{' ' * 14} "
                f"{Color.BRIGHT_GREEN}{elapsed_str}{Color.RESET}",
                force=True,
            )

        self.print_success(f"All results in: {final_out}")
        for subdir in sorted(final_out.iterdir()):
            if subdir.is_dir():
                cnt = len(list(subdir.rglob("*")))
                self.print_info(
                    f"  📁 {Color.BRIGHT_CYAN}{subdir.name}{Color.RESET} "
                    f"{Color.DIM}({cnt} items){Color.RESET}"
                )
        self.display_random_quote()

        self.print_info("Please cite:")
        for line in CITATION_LINES:
            self.print_info(f"  {line}")

        if self.warnings:
            print()
            if not self.quiet:
                print(f"{self.color_warning}Warnings encountered:{Color.RESET}")
            for w in self.warnings:
                self.print_info(f"  - {w}")
        if self.errors:
            print()
            if not self.quiet:
                print(f"{self.color_error}Errors encountered:{Color.RESET}")
            for e in self.errors:
                self.print_info(f"  - {e}")

        return True


# -----------------------------------------------------------------------------
# Colored help
# -----------------------------------------------------------------------------
def print_colored_help():
    usage = f"{Color.BRIGHT_CYAN}usage: pseudoscope [OPTIONS]{Color.RESET}"
    description = (f"{Color.BRIGHT_WHITE}PseudoScope: Complete P. aeruginosa genomic "
                   f"analysis pipeline (v{__version__}){Color.RESET}")

    options = [
        (f"{Color.BRIGHT_GREEN}-i INPUT, --input INPUT{Color.RESET}", 'Input FASTA file, directory, or glob pattern (e.g. "*.fna")'),
        (f"{Color.BRIGHT_GREEN}-o OUTPUT, --output OUTPUT{Color.RESET}", "Output directory for results"),
        (f"{Color.BRIGHT_GREEN}-t THREADS, --threads THREADS{Color.RESET}", f"Number of threads (default: {Color.BRIGHT_YELLOW}2{Color.RESET})"),
        (f"{Color.BRIGHT_GREEN}--quiet{Color.RESET}", "Suppress all non-error console output"),
        (f"{Color.BRIGHT_GREEN}--keep-temp{Color.RESET}", "Do not delete temporary directories (for debugging)"),
        (f"{Color.BRIGHT_GREEN}--tmpdir PATH{Color.RESET}", "Override temporary directory root (HPC/Docker)"),
        (f"{Color.BRIGHT_GREEN}--clean-output{Color.RESET}", "Delete output directory before analysis (fresh run)"),
        (f"{Color.BRIGHT_GREEN}--report-level {{final,summary,full}}{Color.RESET}", "Controls which files are kept after the run (default: summary)"),
        (f"{Color.BRIGHT_GREEN}-v, --version{Color.RESET}", "Show version and exit"),
        (f"{Color.BRIGHT_GREEN}--update-amr-db{Color.RESET}", "Update AMRfinderPlus database (incremental) and exit"),
        (f"{Color.BRIGHT_GREEN}--force-update-amr-db{Color.RESET}", "Force complete AMR database overwrite and exit"),
        (f"{Color.BRIGHT_GREEN}--pull-mlst-db{Color.RESET}", "Download pre-built P. aeruginosa MLST scheme from GitHub (no API key) and exit"),
        (f"{Color.BRIGHT_GREEN}--update-mlst-db{Color.RESET}", "Update P. aeruginosa MLST database from PubMLST (API key required) and exit"),
        (f"{Color.BRIGHT_GREEN}--pull-mge-db{Color.RESET}", "Download the mobileOG-db payload (~300 MB) and exit (conda users only)"),
        (f"{Color.BRIGHT_GREEN}--skip-qc{Color.RESET}", "Skip FASTA QC"),
        (f"{Color.BRIGHT_GREEN}--skip-mlst{Color.RESET}", "Skip MLST"),
        (f"{Color.BRIGHT_GREEN}--skip-past{Color.RESET}", "Skip PAST serotyping"),
        (f"{Color.BRIGHT_GREEN}--skip-abricate{Color.RESET}", "Skip ABRicate"),
        (f"{Color.BRIGHT_GREEN}--skip-amr{Color.RESET}", "Skip AMRfinderPlus"),
        (f"{Color.BRIGHT_GREEN}--skip-mge{Color.RESET}", "Skip MGE profiling"),
        (f"{Color.BRIGHT_GREEN}--skip-summary{Color.RESET}", "Skip gene-centric ultimate summary report"),
        (f"{Color.BRIGHT_GREEN}--skip-sample-centric{Color.RESET}", "Skip sample-centric interactive report"),
        (f"{Color.BRIGHT_GREEN}--skip-visualization{Color.RESET}", "Skip visualisation dashboard"),
        (f"{Color.BRIGHT_GREEN}--abricate-minid INT{Color.RESET}", "ABRicate minimum %identity (default: 80)"),
        (f"{Color.BRIGHT_GREEN}--abricate-mincov INT{Color.RESET}", "ABRicate minimum %coverage (default: 80)"),
        (f"{Color.BRIGHT_GREEN}--past-min-pident INT{Color.RESET}", "PAST minimum percent identity (default: 95)"),
        (f"{Color.BRIGHT_GREEN}--past-min-coverage INT{Color.RESET}", "PAST minimum percent coverage (default: 95)"),
        (f"{Color.BRIGHT_GREEN}--amr-min-identity FLOAT{Color.RESET}", "AMR minimum identity (0-1)"),
        (f"{Color.BRIGHT_GREEN}--amr-min-coverage FLOAT{Color.RESET}", "AMR minimum coverage (0-1)"),
        (f"{Color.BRIGHT_GREEN}--amr-skip-mutations{Color.RESET}", "Disable point mutation reporting in AMR"),
        (f"{Color.BRIGHT_GREEN}--amr-force-update{Color.RESET}", "Force AMR database update before analysis"),
        (f"{Color.BRIGHT_GREEN}-h, --help{Color.RESET}", "Show this help and exit"),
    ]

    epilog = f"""
{Color.BRIGHT_YELLOW}EXAMPLES:{Color.RESET}
  # Basic analysis on all .fna files
  pseudoscope -i "{Color.BRIGHT_CYAN}*.fna{Color.RESET}" -o results

  # Directory input, 4 threads
  pseudoscope -i {Color.BRIGHT_CYAN}genomes/{Color.RESET} -o results --threads 4

  # Report level control (manage disk usage on large batches)
  pseudoscope -i "{Color.BRIGHT_CYAN}*.fna{Color.RESET}" -o results --report-level summary   # default
  pseudoscope -i "{Color.BRIGHT_CYAN}*.fna{Color.RESET}" -o results --report-level final     # final reports only
  pseudoscope -i "{Color.BRIGHT_CYAN}*.fna{Color.RESET}" -o results --report-level full      # keep everything

  # HPC/Docker — override temp location
  pseudoscope -i "{Color.BRIGHT_CYAN}*.fna{Color.RESET}" -o results --tmpdir /scratch/user/pseudoscope_tmp

  # Custom thresholds
  pseudoscope -i "{Color.BRIGHT_CYAN}*.fna{Color.RESET}" -o results \\
      --past-min-pident 90 --past-min-coverage 85 \\
      --abricate-minid 85 --abricate-mincov 80 \\
      --amr-min-identity 0.95 --amr-min-coverage 0.9

  # Database management
  pseudoscope --update-amr-db
  pseudoscope --force-update-amr-db
  pseudoscope --pull-mlst-db
  pseudoscope --update-mlst-db
  pseudoscope --pull-mge-db

  # Keep temp dirs for debugging
  pseudoscope -i "{Color.BRIGHT_CYAN}*.fna{Color.RESET}" -o results --keep-temp

{Color.BRIGHT_YELLOW}REQUIRED BEFORE FIRST ANALYSIS:{Color.RESET}
  1. {Color.BRIGHT_GREEN}abricate --setupdb{Color.RESET}                (ABRicate databases)
  2. {Color.BRIGHT_GREEN}pseudoscope --update-amr-db{Color.RESET}      (AMRfinderPlus database)
  3. {Color.BRIGHT_GREEN}pseudoscope --pull-mlst-db{Color.RESET}       (MLST Pasteur scheme, from GitHub)
  4. {Color.BRIGHT_GREEN}pseudoscope --pull-mge-db{Color.RESET}        (mobileOG-db payload, ~300 MB — conda users only)

{Color.BRIGHT_YELLOW}REPORT LEVELS:{Color.RESET}
  {Color.YELLOW}final{Color.RESET}   - Only {DIR_FINAL}/ is kept. All module dirs are deleted.
            Best for archival and very large batches.
  {Color.YELLOW}summary{Color.RESET} - (default) {DIR_FINAL}/ + batch TSVs/HTMLs/JSONs in each
            module directory. Per-sample folders are removed.
  {Color.YELLOW}full{Color.RESET}    - Everything is kept, including all per-sample HTMLs.
            Best for small batches or when debugging.

{Color.BRIGHT_YELLOW}MLST DATABASE MANAGEMENT:{Color.RESET}
  --pull-mlst-db     : Download pre-built P. aeruginosa MLST scheme from the
                       GitHub repo. No API key required.
  --update-mlst-db   : Refresh the MLST scheme directly from PubMLST using
                       mlstdb. Requires a one-time PubMLST API key setup:
                         1. Register at https://pubmlst.org/site-accounts
                         2. Generate an API key in your profile
                         3. Run: mlstdb connect --db pubmlst --api-key
                         4. Paste your key when prompted

{Color.BRIGHT_YELLOW}MGE DATABASE MANAGEMENT:{Color.RESET}
  --pull-mge-db      : Download the mobileOG-db payload (~300 MB) from the
                       GitHub release. Required once after conda install
                       if you want MGE profiling. Docker and Singularity
                       users do NOT need this — the data is already baked
                       into the image.

                       Location: ~/.local/share/pseudoscope/mge_module/mobileOG-db/
                       Env override: PSEUDOSCOPE_MGE_DATA_DIR
                       Override URL: PSEUDOSCOPE_MGE_PULL_URL

{Color.BRIGHT_YELLOW}HPC / DOCKER NOTES:{Color.RESET}
  • All module writes go to a temporary directory (TMPDIR, --tmpdir, or PSEUDOSCOPE_TMPDIR).
  • MLST database resolved from PSEUDOSCOPE_MLST_DB, then
    ~/.local/share/pseudoscope/mlst_module/db, then packaged db/.
  • MGE database resolved from PSEUDOSCOPE_MGE_DATA_DIR, then
    ~/.local/share/pseudoscope/mge_module/mobileOG-db, then packaged
    mge_module/mobileOG-db/.
  • Module paths can be overridden with PSEUDOSCOPE_<MODULE>_MODULE_PATH
    or PSEUDOSCOPE_<MODULE>_PATH env vars.

{Color.BRIGHT_YELLOW}SUPPORTED FASTA FORMATS:{Color.RESET} {Color.BRIGHT_CYAN}.fna, .fasta, .fa, .fsa{Color.RESET}

{Color.BRIGHT_YELLOW}AUTHOR:{Color.RESET}
  {AUTHOR_NAME} <{AUTHOR_EMAIL}>
  {AUTHOR_AFFILIATION}

{Color.BRIGHT_YELLOW}CITATION:{Color.RESET}
  Beckley B, Amarh V. PseudoScope: a species-specific bioinformatics suite for rapid and
  accessible Pseudomonas aeruginosa genomic analysis. GitHub 2026.
  https://github.com/bbeckley-hub/pseudoscope
    """

    print()
    print(usage)
    print()
    print(description)
    print()
    print(f"{Color.BRIGHT_YELLOW}OPTIONS:{Color.RESET}")
    for opt, desc in options:
        print(f"  {opt:<52} {desc}")
    print(epilog)


def main():
    parser = argparse.ArgumentParser(add_help=False,
                                     formatter_class=ColoredHelpFormatter)
    parser.add_argument('-i', '--input', help=argparse.SUPPRESS)
    parser.add_argument('-o', '--output', help=argparse.SUPPRESS)
    parser.add_argument('-t', '--threads', type=int, default=2,
                        help=argparse.SUPPRESS)
    parser.add_argument('--quiet', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--keep-temp', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--tmpdir', help=argparse.SUPPRESS)
    parser.add_argument('--clean-output', action='store_true',
                        help=argparse.SUPPRESS)
    parser.add_argument('--report-level',
                        choices=['final', 'summary', 'full'],
                        default='summary', help=argparse.SUPPRESS)
    parser.add_argument('-v', '--version', action='store_true',
                        help=argparse.SUPPRESS)

    parser.add_argument('--update-amr-db', action='store_true',
                        help=argparse.SUPPRESS)
    parser.add_argument('--force-update-amr-db', action='store_true',
                        help=argparse.SUPPRESS)
    parser.add_argument('--pull-mlst-db', action='store_true',
                        help=argparse.SUPPRESS)
    parser.add_argument('--update-mlst-db', action='store_true',
                        help=argparse.SUPPRESS)
    parser.add_argument('--pull-mge-db', action='store_true',
                        help=argparse.SUPPRESS)

    parser.add_argument('--skip-qc', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-mlst', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-past', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-abricate', action='store_true',
                        help=argparse.SUPPRESS)
    parser.add_argument('--skip-amr', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-mge', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--skip-summary', action='store_true',
                        help=argparse.SUPPRESS)
    parser.add_argument('--skip-sample-centric', action='store_true',
                        help=argparse.SUPPRESS)
    parser.add_argument('--skip-visualization', action='store_true',
                        help=argparse.SUPPRESS)

    parser.add_argument('--abricate-minid', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--abricate-mincov', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--past-min-pident', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--past-min-coverage', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--amr-min-identity', type=float, help=argparse.SUPPRESS)
    parser.add_argument('--amr-min-coverage', type=float, help=argparse.SUPPRESS)
    parser.add_argument('--amr-skip-mutations', action='store_true',
                        help=argparse.SUPPRESS)
    parser.add_argument('--amr-force-update', action='store_true',
                        help=argparse.SUPPRESS)

    parser.add_argument('-h', '--help', action='store_true',
                        help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.help:
        print_colored_help()
        sys.exit(0)

    if args.version:
        print(f"pseudoscope {__version__}")
        sys.exit(0)

    orch = PseudoScopeOrchestrator(
        quiet=args.quiet,
        keep_temp=args.keep_temp,
        tmpdir=args.tmpdir,
    )

    if args.update_amr_db or args.force_update_amr_db:
        ok = orch.run_complete_analysis(
            "", "", update_amr_only=True,
            amr_force_update=args.force_update_amr_db,
        )
        sys.exit(0 if ok else 1)

    if args.pull_mlst_db:
        ok = orch.run_complete_analysis("", "", pull_mlst_only=True)
        sys.exit(0 if ok else 1)

    if args.update_mlst_db:
        ok = orch.run_complete_analysis("", "", update_mlst_only=True)
        sys.exit(0 if ok else 1)

    if args.pull_mge_db:
        ok = orch.run_complete_analysis("", "", pull_mge_only=True)
        sys.exit(0 if ok else 1)

    if not args.input or not args.output:
        parser.error(
            "When not using a database-management flag, both -i/--input and "
            "-o/--output are required."
        )

    if (args.report_level == 'final'
            and args.skip_summary
            and args.skip_sample_centric):
        parser.error(
            "--report-level final requires at least one reporter to run "
            "(--skip-summary and --skip-sample-centric both set would leave "
            "the final report empty and pruning would delete all output)."
        )

    skip = {
        'qc':       args.skip_qc,
        'mlst':     args.skip_mlst,
        'past':     args.skip_past,
        'abricate': args.skip_abricate,
        'amr':      args.skip_amr,
        'mge':      args.skip_mge,
    }

    try:
        ok = orch.run_complete_analysis(
            input_path=args.input,
            output_dir=args.output,
            threads=args.threads,
            skip=skip,
            skip_summary=args.skip_summary,
            skip_sample_centric=args.skip_sample_centric,
            skip_viz=args.skip_visualization,
            clean_output=args.clean_output,
            abricate_minid=args.abricate_minid,
            abricate_mincov=args.abricate_mincov,
            past_min_pident=args.past_min_pident,
            past_min_coverage=args.past_min_coverage,
            amr_min_identity=args.amr_min_identity,
            amr_min_coverage=args.amr_min_coverage,
            amr_skip_mutations=args.amr_skip_mutations,
            amr_force_update=args.amr_force_update,
            report_level=args.report_level,
        )
        sys.exit(0 if ok else 1)
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