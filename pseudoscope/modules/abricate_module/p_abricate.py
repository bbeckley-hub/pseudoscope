#!/usr/bin/env python3
"""
PseudoScope ABRicate Module - Comprehensive Gene Tracking for Pseudomonas aeruginosa
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: (Updated: 2026-10-02)

REFACTORED (2026-11-12) to match  consolidation pattern:
  - Per-sample output reduced to 2 files:
      {genome}_comprehensive_abricate_report.html
      abricate_raw.tsv (all databases merged)
  - Batch-level JSON reduced to ONE file:
      pseudo_abricate_master_summary.json
  - All raw per-database .txt files now live inside a tempfile.TemporaryDirectory
    and are automatically deleted once the merged TSV has been written.
"""

import subprocess
import sys
import os
import glob
import logging
import tempfile
import random
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any, Optional
import argparse
import re
from datetime import datetime
import psutil
import json

from collections import defaultdict, Counter

# pandas is only needed for the optional Excel export
try:
    import pandas as pd
except ImportError:
    pd = None


class PseudoAbricateExecutor:
    """ABRicate executor for Pseudomonas aeruginosa with comprehensive HTML reporting - MAXIMUM SPEED"""

    def __init__(self, cpus: Optional[int] = None, minid: int = 80, mincov: int = 80):
        self.minid = minid
        self.mincov = mincov
        self.logger = self._setup_logging()
        self.available_ram = self._get_available_ram()
        self.cpus = self._calculate_optimal_cpus(cpus)

        self.required_databases = [
            'ncbi', 'card', 'resfinder', 'vfdb', 'argannot',
            'plasmidfinder', 'megares', 'ecoh', 'bacmet2', 'ecoli_vf'
        ]

        # ─── Gene categories  ─────────────────────────────────
        self.critical_carbapenemases = {
            'blaOXA-2', 'blaOXA-10', 'blaOXA-14', 'blaOXA-17', 'blaOXA-19', 'blaOXA-28',
            'blaOXA-35', 'blaOXA-45', 'blaOXA-50', 'blaOXA-198',
            'OXA-2', 'OXA-10', 'OXA-14', 'OXA-17', 'OXA-19', 'OXA-28', 'OXA-35',
            'OXA-45', 'OXA-50', 'OXA-198',
            'blaIMP', 'blaVIM', 'blaNDM', 'blaKPC', 'blaGES', 'blaSPM', 'blaAIM', 'blaDIM',
            'IMP', 'VIM', 'NDM', 'KPC', 'GES', 'SPM', 'AIM', 'DIM'
        }

        self.critical_esbls = {
            'blaPER', 'blaVEB', 'blaBEL', 'blaGES', 'blaTEM', 'blaSHV', 'blaCTX-M',
            'PER', 'VEB', 'BEL', 'GES', 'TEM', 'SHV', 'CTX-M'
        }

        self.critical_colistin = {
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8',
            'mcr-9', 'mcr-10', 'pmrA', 'pmrB', 'phoP', 'phoQ', 'mgrB', 'lpxA', 'lpxC',
            'lpxD', 'arnT', 'eptA', 'eptB', 'basS', 'basR'
        }

        self.critical_aminoglycoside = {
            'armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD', 'rmtE', 'rmtF', 'rmtG', 'rmtH', 'npmA',
            'aac(3)', 'aac(6\')', 'ant(2\")', 'ant(4\')', 'aph(3\')', 'aph(6)',
            'aacC1', 'aacC2', 'aacC4', 'aacA4', 'aacA7', 'aadA1', 'aadA2', 'aadA5', 'aadA7',
            'strA', 'strB', 'aphA1', 'aphA2', 'aphA3', 'aphA6', 'aac3', 'aac6', 'aadA', 'aadB'
        }

        self.high_risk_resistance = {
            'mexA', 'mexB', 'mexC', 'mexD', 'mexE', 'mexF', 'mexX', 'mexY',
            'oprM', 'oprN', 'oprJ', 'adeA', 'adeB', 'adeC',
            'tetA', 'tetB', 'tetC', 'tetD', 'tetE', 'tetG', 'tetH', 'tetK', 'tetL',
            'tetM', 'tetO', 'tetQ', 'tetS', 'tetX',
            'sul1', 'sul2', 'sul3', 'sul4',
            'dfrA1', 'dfrA5', 'dfrA7', 'dfrA8', 'dfrA12', 'dfrA14', 'dfrA17', 'dfrA19',
            'dfrA20', 'dfrA21', 'dfrB1', 'dfrB2', 'dfrB3', 'dfrB4', 'dfrB5', 'dfrB6', 'dfrB7',
            'catA1', 'catA2', 'catB2', 'catB3', 'catB8', 'catI', 'catII', 'catIII',
            'cmlA', 'cmlA1', 'cmlA5', 'cmlA6', 'cmlA7', 'floR',
            'ermA', 'ermB', 'ermC', 'ermF', 'ermG', 'ermX', 'ermY',
            'mphA', 'mphB', 'mphC', 'mphD', 'mphE', 'msrA', 'msrB', 'msrC', 'msrD',
            'qnrA1', 'qnrB1', 'qnrB2', 'qnrB4', 'qnrB6', 'qnrB10', 'qnrB19', 'qnrS1',
            'qnrS2', 'qnrVC1', 'qnrVC4', 'aac(6\')-Ib-cr',
            'qepA', 'qepA1', 'qepA2', 'qepA3', 'qepA4',
            'fosA', 'fosB', 'fosC', 'fosX',
            'arr-2', 'arr-3', 'arr-4', 'arr-5', 'arr-6', 'arr-7',
        }

        self.virulence_genes = {
            'exoU', 'exoS', 'exoT', 'exoY',
            'pscC', 'pscD', 'pscE', 'pscF', 'pscG', 'pscH', 'pscI', 'pscJ', 'pscK', 'pscL',
            'pcrV', 'pcrG', 'pcrH', 'pcrD', 'pcrC', 'pcrB', 'pcrR',
            'exsA', 'exsB', 'exsC', 'exsD',
            'algD', 'algU', 'alg8', 'alg44', 'algK', 'algE', 'mucA', 'mucB', 'mucC', 'mucD',
            'lasA', 'lasB', 'lasI', 'lasR', 'rhlA', 'rhlB', 'rhlI', 'rhlR',
            'aprA', 'aprB', 'aprC', 'aprD', 'aprE', 'aprF',
            'plcH', 'plcN', 'plcB',
            'phzA', 'phzB', 'phzC', 'phzD', 'phzE', 'phzF', 'phzG', 'phzM', 'phzS',
            'pvdA', 'pvdD', 'pvdE', 'pvdF', 'pvdG', 'pvdH', 'pvdI', 'pvdJ', 'pvdL',
            'pvdN', 'pvdO', 'pvdP', 'pvdQ',
            'fpvA', 'fpvB', 'fpvC', 'fpvD', 'fpvE', 'fpvF', 'fpvG', 'fpvH', 'fpvI',
            'fpvJ', 'fpvK',
            'fliC', 'fliD', 'fliE', 'fliF', 'fliG', 'fliH', 'fliI', 'fliJ', 'fliK',
            'fliL', 'fliM', 'fliN', 'fliO', 'fliP', 'fliQ', 'fliR',
            'fleN', 'fleQ', 'fleR',
            'pilA', 'pilB', 'pilC', 'pilD', 'pilE', 'pilF', 'pilG', 'pilH', 'pilI',
            'pilJ', 'pilK', 'pilL', 'pilM', 'pilN', 'pilO', 'pilP', 'pilQ', 'pilR',
            'pilS', 'pilT', 'pilU',
            'pslA', 'pslB', 'pslC', 'pslD', 'pslE', 'pslF', 'pslG', 'pslH', 'pslI',
            'pslJ', 'pslK', 'pslL', 'pslM', 'pslN', 'pslO', 'pslP',
            'pelA', 'pelB', 'pelC', 'pelD', 'pelE', 'pelF', 'pelG',
            'lpxA', 'lpxC', 'lpxD', 'lpxL', 'lpxO', 'lpxP',
            'lasI', 'lasR', 'rhlI', 'rhlR', 'pqsA', 'pqsB', 'pqsC', 'pqsD', 'pqsE',
            'pqsH', 'mvfR'
        }

        self.critical_resistance_genes = self.critical_carbapenemases.union(
            self.critical_esbls
        ).union(self.critical_colistin).union(self.critical_aminoglycoside)

        self.critical_virulence_genes = {'exoU', 'exoS', 'exoT', 'exoY'}

        self.science_quotes = [
            "Pseudomonas aeruginosa: the master of adaptation. - Unknown",
            "The important thing is not to stop questioning. Curiosity has its own reason for existing. - Albert Einstein",
            "Science is not only a disciple of reason but also one of romance and passion. - Stephen Hawking",
            "Somewhere, something incredible is waiting to be known. - Carl Sagan",
            "In science, there are no shortcuts to truth. - Karl Popper",
            "Science knows no country, because knowledge belongs to humanity. - Louis Pasteur",
            "The science of today is the technology of tomorrow. - Edward Teller",
            "Nothing in life is to be feared, it is only to be understood. - Marie Curie",
            "Research is what I'm doing when I don't know what I'm doing. - Wernher von Braun"
        ]

        self.ascii_art = r"""
██████╗ ███████╗███████╗██╗   ██╗██████╗  ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██╔════╝██║   ██║██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
██████╔╝███████╗█████╗  ██║   ██║██║  ██║██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗
██╔═══╝ ╚════██║██╔══╝  ██║   ██║██║  ██║██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝
██║     ███████║███████╗╚██████╔╝██████╔╝╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝     ╚══════╝╚══════╝ ╚═════╝ ╚═════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
"""

        self.metadata = {
            "tool_name": "PseudoScope ABRicate for Pseudomonas aeruginosa",
            "version": "1.3.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School - Department of Medical Biochemistry",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "min_identity": self.minid,
            "min_coverage": self.mincov
        }

    # ─────────────────────────── setup helpers ───────────────────────────
    def _setup_logging(self):
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)

    def _get_available_ram(self) -> float:
        try:
            return psutil.virtual_memory().available / (1024 ** 3)
        except Exception as e:
            self.logger.warning(f"Could not detect RAM: {e}")
            return 8

    def _calculate_optimal_cpus(self, user_cpus: int = None) -> int:
        if user_cpus is not None:
            self._log_resource_info(user_cpus)
            return user_cpus
        try:
            total_physical_cores = psutil.cpu_count(logical=False) or os.cpu_count() or 2
            if total_physical_cores <= 4:
                optimal_cpus = total_physical_cores
            elif total_physical_cores <= 8:
                optimal_cpus = total_physical_cores - 1
            elif total_physical_cores <= 16:
                optimal_cpus = max(8, total_physical_cores - 1)
            elif total_physical_cores <= 32:
                optimal_cpus = max(16, total_physical_cores - 2)
            else:
                optimal_cpus = min(32, int(total_physical_cores * 0.95))
            optimal_cpus = max(1, min(optimal_cpus, total_physical_cores))
            self._log_resource_info(optimal_cpus, total_physical_cores)
            return optimal_cpus
        except Exception as e:
            self.logger.warning(f"Could not detect CPU cores, using maximum available: {e}")
            return os.cpu_count() or 4

    def _log_resource_info(self, cpus: int, total_cores: int = None):
        self.logger.info(f"Available RAM: {self.available_ram:.1f} GB")
        if total_cores:
            self.logger.info(f"System CPU cores: {total_cores}")
            utilization = (cpus / total_cores) * 100
            self.logger.info(f"Using CPU cores: {cpus} ({utilization:.1f}% of available cores)")
        else:
            self.logger.info(f"Using user-specified CPU cores: {cpus}")
        if cpus == 1:
            self.logger.info("Performance: single-core")
        elif cpus <= 4:
            self.logger.info("Performance: multi-core (small system)")
        elif cpus <= 8:
            self.logger.info("Performance: high-speed mode")
        else:
            self.logger.info("Performance: MAXIMUM SPEED MODE")

    def check_abricate_installed(self) -> bool:
        try:
            result = subprocess.run(['abricate', '--version'],
                                    capture_output=True, text=True, check=True)
            version_line = result.stdout.strip()
            self.logger.info("ABRicate version: %s", version_line)
            version_match = re.search(r'(\d+\.\d+\.\d+)', version_line)
            if version_match:
                version_str = version_match.group(1)
                if version_str >= "1.2.0":
                    self.logger.info("ABRicate version meets requirement (>=1.2.0)")
                    return True
                self.logger.error("ABRicate version too old: %s. Required >=1.2.0", version_str)
                return False
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.logger.error(
                "ABRicate not found. Install with: conda install -c bioconda abricate"
            )
            return False

    def setup_abricate_databases(self):
        self.logger.info("Setting up ABRicate databases for P. aeruginosa analysis...")
        available_dbs = []
        missing_dbs = []
        try:
            check_result = subprocess.run(['abricate', '--list'],
                                          capture_output=True, text=True, check=True)
            for db in self.required_databases:
                if db in check_result.stdout:
                    self.logger.info("Database available: %s", db)
                    available_dbs.append(db)
                else:
                    self.logger.warning("Database not available: %s", db)
                    missing_dbs.append(db)
            for db in missing_dbs:
                self.logger.info("Attempting to setup database: %s", db)
                try:
                    subprocess.run(['abricate', '--setupdb', '--db', db],
                                   capture_output=True, text=True, check=True)
                    self.logger.info("Database setup completed: %s", db)
                    available_dbs.append(db)
                except subprocess.CalledProcessError as e:
                    self.logger.error("Failed to setup database %s: %s", db, e.stderr)
            self.required_databases = available_dbs
            self.logger.info("Using databases: %s", ", ".join(self.required_databases))
        except subprocess.CalledProcessError as e:
            self.logger.error("Error checking ABRicate databases: %s", e.stderr)
        except Exception as e:
            self.logger.error("Unexpected error setting up databases: %s", e)

    # ─────────────────────── abricate per-db runner ──────────────────────
    def run_abricate_single_db(self, genome_file: str, database: str,
                               temp_dir: str) -> Dict[str, Any]:
        """Run ABRicate for one genome against one database. Writes into temp_dir only."""
        genome_name = Path(genome_file).stem
        raw_file = os.path.join(temp_dir, f"abricate_{database}.txt")
        cmd = [
            'abricate', genome_file, '--db', database,
            '--minid', str(self.minid), '--mincov', str(self.mincov),
        ]
        self.logger.info("Running ABRicate: %s --db %s (minid=%d, mincov=%d)",
                         genome_name, database, self.minid, self.mincov)
        try:
            with open(raw_file, 'w') as outfile:
                subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE,
                               text=True, check=True)
            hits = self._parse_abricate_output(raw_file)
            with open(raw_file, 'r') as f:
                raw_content = f.read()
            return {
                'database': database,
                'genome': genome_name,
                'hits': hits,
                'hit_count': len(hits),
                'raw_content': raw_content,
                'status': 'success',
            }
        except subprocess.CalledProcessError as e:
            self.logger.error("ABRicate failed for %s on %s: %s",
                              database, genome_name, e.stderr)
            return {
                'database': database,
                'genome': genome_name,
                'hits': [],
                'hit_count': 0,
                'raw_content': '',
                'status': 'failed',
            }

    def _parse_abricate_output(self, abricate_file: str) -> List[Dict]:
        hits = []
        try:
            with open(abricate_file, 'r') as f:
                lines = f.readlines()
            if not lines:
                return hits
            headers = []
            data_lines = []
            for line in lines:
                if line.startswith('#FILE') and not headers:
                    headers = line.strip().replace('#', '').split('\t')
                elif line.strip() and not line.startswith('#'):
                    data_lines.append(line.strip())
            if not headers:
                self.logger.warning("No headers found in %s", abricate_file)
                return hits
            expected_columns = len(headers)
            for line_num, line in enumerate(data_lines, 1):
                parts = line.split('\t')
                if len(parts) > expected_columns:
                    combined_parts = parts[:expected_columns - 1]
                    combined_parts.append('\t'.join(parts[expected_columns - 1:]))
                    parts = combined_parts
                elif len(parts) < expected_columns:
                    parts.extend([''] * (expected_columns - len(parts)))
                if len(parts) == expected_columns:
                    hit = {}
                    for i, header in enumerate(headers):
                        hit[header] = parts[i] if i < len(parts) else ''
                    processed_hit = {
                        'file': hit.get('FILE', ''),
                        'sequence': hit.get('SEQUENCE', ''),
                        'start': hit.get('START', ''),
                        'end': hit.get('END', ''),
                        'strand': hit.get('STRAND', ''),
                        'gene': hit.get('GENE', ''),
                        'coverage': hit.get('COVERAGE', ''),
                        'coverage_map': hit.get('COVERAGE_MAP', ''),
                        'gaps': hit.get('GAPS', ''),
                        'coverage_percent': hit.get('%COVERAGE', ''),
                        'identity_percent': hit.get('%IDENTITY', ''),
                        'database': hit.get('DATABASE', ''),
                        'accession': hit.get('ACCESSION', ''),
                        'product': hit.get('PRODUCT', ''),
                        'resistance': hit.get('RESISTANCE', ''),
                    }
                    hits.append(processed_hit)
                else:
                    self.logger.warning(
                        "Line %d has %d parts, expected %d: %s",
                        line_num, len(parts), expected_columns, line[:100] + "..."
                    )
        except Exception as e:
            self.logger.error("Error parsing %s: %s", abricate_file, e)
        self.logger.info("Parsed %d hits from %s", len(hits), abricate_file)
        return hits

    def _write_merged_raw_tsv(self, results: Dict[str, Dict], output_dir: str):
        """Merge every per-database ABRicate output into ONE abricate_raw.tsv."""
        output_file = os.path.join(output_dir, "abricate_raw.tsv")
        header_written = False
        data_lines = 0
        with open(output_file, 'w') as f:
            for db, result in results.items():
                raw = result.get('raw_content', '') or ''
                if not raw:
                    continue
                for line in raw.splitlines():
                    if not line.strip():
                        continue
                    if line.startswith('#'):
                        if not header_written:
                            f.write(line + '\n')
                            header_written = True
                    else:
                        f.write(line + '\n')
                        data_lines += 1
        self.logger.info("Merged raw TSV: %s (%d hits across %d databases)",
                         output_file, data_lines, len(results))

    # ───────────────────── resistance/virulence analysis ─────────────────
    def analyze_pseudo_resistance(self, all_hits: List[Dict]) -> Dict[str, Any]:
        analysis = {
            'carbapenemase_status': 'negative',
            'esbl_status': 'negative',
            'colistin_resistance': 'negative',
            'critical_carbapenemase_genes': [],
            'critical_esbl_genes': [],
            'critical_colistin_genes': [],
            'critical_aminoglycoside_genes': [],
            'high_risk_resistance_genes': [],
            'critical_virulence_genes': [],
            'high_risk_virulence_genes': [],
            'moderate_risk_genes': [],
            'resistance_classes': {},
            'total_critical_resistance': 0,
            'total_high_risk_resistance': 0,
            'total_critical_virulence': 0,
            'total_high_risk_virulence': 0,
            'total_hits': len(all_hits),
        }
        for hit in all_hits:
            gene = hit['gene']
            gene_base = gene.split('-')[0] if '-' in gene else gene
            gene_lower = gene.lower()

            if any(carba in gene_lower for carba in ['oxa', 'imp', 'vim', 'ndm', 'kpc', 'ges', 'spm', 'aim', 'dim']):
                if any(crit in gene_lower for crit in [g.lower() for g in self.critical_carbapenemases]):
                    analysis['carbapenemase_status'] = 'positive'
                    analysis['critical_carbapenemase_genes'].append(self._hit_record(hit, 'CARBAPENEMASE'))
            elif any(esbl in gene_lower for esbl in ['per', 'veb', 'bel', 'ges', 'tem', 'shv', 'ctx']):
                if any(crit in gene_lower for crit in [g.lower() for g in self.critical_esbls]):
                    analysis['esbl_status'] = 'positive'
                    analysis['critical_esbl_genes'].append(self._hit_record(hit, 'ESBL'))
            elif any(mcr in gene_lower for mcr in ['mcr']) or any(pm in gene_lower for pm in ['pmr', 'pho', 'mgr', 'lpx', 'arn', 'ept', 'bas']):
                if any(crit in gene_lower for crit in [g.lower() for g in self.critical_colistin]):
                    analysis['colistin_resistance'] = 'positive'
                    analysis['critical_colistin_genes'].append(self._hit_record(hit, 'COLISTIN-RES'))
            elif any(ag in gene_lower for ag in ['arm', 'rmt', 'npm', 'aac', 'ant', 'aph', 'str']):
                if any(crit in gene_lower for crit in [g.lower() for g in self.critical_aminoglycoside]):
                    analysis['critical_aminoglycoside_genes'].append(self._hit_record(hit, 'CRITICAL_AMINOGLYCOSIDE'))
            elif any(hr in gene_lower for hr in [g.lower() for g in self.high_risk_resistance]):
                analysis['high_risk_resistance_genes'].append(self._hit_record(hit, 'HIGH-RISK'))

            if any(crit_vf in gene_base for crit_vf in self.critical_virulence_genes):
                analysis['critical_virulence_genes'].append(self._hit_record(hit, 'CRITICAL-VIRULENCE'))
            elif any(vf in gene_lower for vf in [g.lower() for g in self.virulence_genes]):
                is_high = any(hr in gene_lower for hr in ['exo', 'psc', 'pcr', 'exs', 'las', 'rhl', 'apr', 'plc', 'tox'])
                if is_high:
                    analysis['high_risk_virulence_genes'].append(self._hit_record(hit, 'HIGH-VIRULENCE'))
                else:
                    analysis['moderate_risk_genes'].append(self._hit_record(hit, 'MODERATE'))

            resistance_class = self._classify_resistance(hit['product'], gene)
            if resistance_class:
                if resistance_class not in analysis['resistance_classes']:
                    analysis['resistance_classes'][resistance_class] = []
                if gene not in [g['gene'] for g in analysis['resistance_classes'][resistance_class]]:
                    analysis['resistance_classes'][resistance_class].append({
                        'gene': gene, 'product': hit['product']
                    })

        analysis['total_critical_resistance'] = (
            len(analysis['critical_carbapenemase_genes'])
            + len(analysis['critical_esbl_genes'])
            + len(analysis['critical_colistin_genes'])
            + len(analysis['critical_aminoglycoside_genes'])
        )
        analysis['total_high_risk_resistance'] = len(analysis['high_risk_resistance_genes'])
        analysis['total_critical_virulence'] = len(analysis['critical_virulence_genes'])
        analysis['total_high_risk_virulence'] = len(analysis['high_risk_virulence_genes'])
        return analysis

    @staticmethod
    def _hit_record(hit: Dict, risk_level: str) -> Dict:
        return {
            'gene': hit['gene'],
            'product': hit['product'],
            'database': hit['database'],
            'coverage': hit['coverage_percent'],
            'identity': hit['identity_percent'],
            'risk_level': risk_level,
        }

    def _classify_resistance(self, product: str, gene: str) -> str:
        p = product.lower()
        g = gene.lower()
        if any(t in p or t in g for t in ['carbapenem', 'oxa', 'imp', 'vim', 'ndm', 'kpc', 'spm', 'aim']):
            return 'Carbapenem resistance'
        if any(t in p or t in g for t in ['beta-lactam', 'esbl', 'per', 'veb', 'bel', 'ges', 'tem', 'shv', 'ctx']):
            return 'Beta-lactam resistance'
        if any(t in p or t in g for t in ['aminoglycoside', 'aac', 'ant', 'aph', 'str', 'arm', 'rmt']):
            return 'Aminoglycoside resistance'
        if any(t in p or t in g for t in ['tetracycline', 'tet']):
            return 'Tetracycline resistance'
        if any(t in p or t in g for t in ['sulfonamide', 'sul']):
            return 'Sulfonamide resistance'
        if any(t in p or t in g for t in ['trimethoprim', 'dfr']):
            return 'Trimethoprim resistance'
        if any(t in p or t in g for t in ['chloramphenicol', 'cat', 'cml', 'flor']):
            return 'Chloramphenicol resistance'
        if any(t in p or t in g for t in ['macrolide', 'erm', 'mph', 'msr']):
            return 'Macrolide resistance'
        if any(t in p or t in g for t in ['quinolone', 'qnr', 'qep']):
            return 'Quinolone resistance'
        if any(t in p or t in g for t in ['colistin', 'polymyxin', 'mcr', 'pmr', 'pho', 'mgr', 'lpx']):
            return 'Polymyxin resistance'
        if any(t in p or t in g for t in ['fosfomycin', 'fos']):
            return 'Fosfomycin resistance'
        if any(t in p or t in g for t in ['rifampicin', 'arr']):
            return 'Rifampicin resistance'
        if any(t in p or t in g for t in ['efflux', 'mex', 'opr']):
            return 'Efflux pumps'
        if any(t in p or t in g for t in ['virulence', 'toxin', 'exo', 'exs', 'psc', 'las', 'rhl', 'alg', 'pil', 'fli']):
            return 'Virulence factors'
        return 'Other resistance'

    # ─────────────────────── shared HTML helpers ────────────────────────
    def _get_shared_css(self) -> str:
        return """
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            background: linear-gradient(135deg, #2A1408 0%, #5C2E0B 40%, #8B4513 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff; padding: 20px; min-height: 100vh;
        }
        .container { max-width: 1400px; margin: 0 auto; }
        .header { text-align: center; margin-bottom: 30px; }
        .ascii-container {
            background: rgba(0, 0, 0, 0.75); padding: 20px; border-radius: 15px;
            margin-bottom: 20px; box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(210, 180, 140, 0.4);
        
            text-align: center;
        }
        .ascii-art {
            font-family: 'Courier New', monospace; font-size: 10px; line-height: 1.1;
            white-space: pre; color: #E8C39E;
            text-shadow: 0 0 10px rgba(232, 195, 158, 0.4); overflow-x: auto;
        
            display: inline-block;
            text-align: left;
        }
        .card {
            background: rgba(255, 255, 255, 0.96); color: #2c2c2c;
            padding: 25px; margin: 20px 0; border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.2);
        }
        .gene-table {
            width: 100%; border-collapse: collapse; margin: 20px 0;
            background: white; border-radius: 8px; overflow: hidden;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
        }
        .gene-table th, .gene-table td {
            padding: 15px; text-align: left; border-bottom: 1px solid #e0e0e0;
        }
        .gene-table th {
            background: linear-gradient(135deg, #A0522D 0%, #5C2E0B 100%);
            color: white; font-weight: 600;
        }
        tr:hover { background-color: #f8f9fa; }
        .summary-stats {
            display: flex; justify-content: space-around; margin: 20px 0; flex-wrap: wrap;
        }
        .stat-card {
            background: linear-gradient(135deg, #A0522D 0%, #5C2E0B 100%);
            color: white; padding: 20px; border-radius: 12px; text-align: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            margin: 10px; flex: 1; min-width: 200px;
        }
        .quote-container {
            background: rgba(255, 255, 255, 0.08); color: white;
            padding: 20px; border-radius: 12px; margin: 20px 0; text-align: center;
            font-style: italic; border-left: 4px solid #8B4513;
        }
.authorship {
        margin-top: 15px;
        padding: 15px;
        background: rgba(255,255,255,0.08);
        border-radius: 8px;
        font-size: 12px;
    }
                .footer {
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0,0,0,0.35);
            border-radius: 10px;
            font-size: 14px;
        }
        .footer a { color: #E8C39E; text-decoration: none; }
        .footer a:hover { text-decoration: underline; }
        .present { background-color: #d4edda; }
        .high-risk { background-color: #fff3cd; }
        .critical { background-color: #f5e6d3; font-weight: bold; }
        .risk-badge {
            display: inline-block; background: #8B4513; color: white;
            padding: 5px 10px; border-radius: 15px; margin: 2px; font-size: 0.9em;
        }
        .warning-badge {
            display: inline-block; background: #ffc107; color: black;
            padding: 5px 10px; border-radius: 15px; margin: 2px; font-size: 0.9em;
        }
        .safe-badge {
            display: inline-block; background: #28a745; color: white;
            padding: 5px 10px; border-radius: 15px; margin: 2px; font-size: 0.9em;
        }
        .percentage-badge {
            display: inline-block; background: #e3f2fd; color: #1976d2;
            padding: 2px 8px; border-radius: 12px; font-size: 0.85em;
            margin-left: 5px; border: 1px solid #bbdefb;
        }
        .frequency-cell { font-weight: 600; color: #2c3e50; }
        """

    def _get_quote_js(self) -> str:
        return f"""
        <script>
            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            function rotateQuote() {{
                document.getElementById('science-quote').innerHTML = quotes[currentQuote];
                currentQuote = (currentQuote + 1) % quotes.length;
            }}
            setInterval(rotateQuote, 10000);
            document.addEventListener('DOMContentLoaded', rotateQuote);
        </script>
        """

    def _get_footer_html(self, title: str) -> str:
        """Return the shared pasty-style footer block for HTML reports."""
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        clean_title = title
        for prefix in ("PseudoScope - ", "PseudoScope ", "PSEUDOSCOPE - "):
            if clean_title.startswith(prefix):
                clean_title = clean_title[len(prefix):]
                break
        return f"""
    <div class="footer">
        <p><strong>PSEUDOSCOPE</strong> - {clean_title}</p>
        <p class="timestamp">Generated: {current_time}</p>
        <div class="authorship">
            <p><strong>Technical Support &amp; Inquiries:</strong></p>
            <p>Maintainer: Brown Beckley</p>
            <p>GitHub: <a href="https://github.com/bbeckley-hub" style="color:#E8C39E;">bbeckley-hub</a></p>
            <p>Email: <a href="mailto:brownbeckley94@gmail.com" style="color:#E8C39E;">brownbeckley94@gmail.com</a></p>
            <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
            <p style="margin-top:10px;font-size:11px;opacity:0.8;">
                Powered by ABRicate (Seemann) | minID={self.minid}% | minCov={self.mincov}%
            </p>
        </div>
    </div>
    """

    def _render_hit_table(self, title: str, genes: List[Dict], row_class: str, badge_class: str) -> str:
        if not genes:
            return ""
        html = f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #8B4513; padding-bottom: 10px;">{title}</h2>
            <table class="gene-table">
                <thead>
                    <tr>
                        <th>Gene</th><th>Product</th><th>Database</th>
                        <th>Coverage</th><th>Identity</th><th>Risk Level</th>
                    </tr>
                </thead>
                <tbody>
        """
        for g in genes:
            product = g['product']
            if len(product) > 1000:
                product = product[:1000] + "..."
            html += f"""
                    <tr class="{row_class}">
                        <td><strong>{g['gene']}</strong></td>
                        <td title="{g['product']}">{product}</td>
                        <td>{g['database']}</td>
                        <td>{g['coverage']}%</td>
                        <td>{g['identity']}%</td>
                        <td><span class="{badge_class}">{g['risk_level']}</span></td>
                    </tr>
            """
        html += "</tbody></table></div>"
        return html

    # ─────────────────────── comprehensive HTML report ─────────────────
    def create_comprehensive_html_report(self, genome_name: str, results: Dict, output_dir: str):
        all_hits = []
        for db_result in results.values():
            all_hits.extend(db_result['hits'])
        analysis = self.analyze_pseudo_resistance(all_hits)
        random_quote = random.choice(self.science_quotes)

        # Critical alerts
        alerts_html = ""
        critical_alerts = []
        if analysis['carbapenemase_status'] == 'positive':
            critical_alerts.append("🔴 CARBAPENEMASE DETECTED")
        if analysis['esbl_status'] == 'positive':
            critical_alerts.append("🟡 ESBL DETECTED")
        if analysis['colistin_resistance'] == 'positive':
            critical_alerts.append("🔴 COLISTIN RESISTANCE DETECTED")
        if critical_alerts:
            badges = "".join(f'<span class="risk-badge">{a}</span>' for a in critical_alerts)
            alerts_html = f"""
        <div class="card" style="border-left: 4px solid #8B4513;">
            <h2 style="color: #5C2E0B;">⚠️ CRITICAL RESISTANCE ALERTS</h2>
            <div style="margin: 10px 0;">{badges}</div>
        </div>
        """

        # Resistance classes block
        resistance_classes_html = ""
        if analysis['resistance_classes']:
            blocks = ""
            for class_name, genes in analysis['resistance_classes'].items():
                gene_list = ", ".join(sorted(g['gene'] for g in genes))
                blocks += f"""
                <div style="margin: 10px 0; padding: 10px; background: #f8f9fa; border-radius: 8px;">
                    <strong style="color: #8B4513;">{class_name}</strong> ({len(genes)} genes)
                    <br><span style="color: #666; font-size: 0.9em;">{gene_list}</span>
                </div>
                """
            resistance_classes_html = f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #8B4513; padding-bottom: 10px;">🧪 Resistance Classes Detected</h2>
            {blocks}
        </div>
        """

        tables_html = ""
        tables_html += self._render_hit_table(
            "🔴 CRITICAL Carbapenemase Genes", analysis['critical_carbapenemase_genes'],
            "critical", "risk-badge")
        tables_html += self._render_hit_table(
            "🟡 Critical ESBL Genes", analysis['critical_esbl_genes'],
            "critical", "warning-badge")
        tables_html += self._render_hit_table(
            "🔴 Colistin Resistance Genes", analysis['critical_colistin_genes'],
            "critical", "risk-badge")
        tables_html += self._render_hit_table(
            "🔴 Critical Aminoglycoside Genes", analysis['critical_aminoglycoside_genes'],
            "critical", "risk-badge")
        tables_html += self._render_hit_table(
            "🟡 High-Risk Resistance Genes", analysis['high_risk_resistance_genes'],
            "high-risk", "warning-badge")
        tables_html += self._render_hit_table(
            "🟢 Critical Virulence Genes", analysis['critical_virulence_genes'],
            "high-risk", "warning-badge")
        tables_html += self._render_hit_table(
            "🟢 High-Risk Virulence Genes", analysis['high_risk_virulence_genes'],
            "high-risk", "warning-badge")
        tables_html += self._render_hit_table(
            "Moderate Risk Genes", analysis['moderate_risk_genes'],
            "present", "safe-badge")

        # DB summary rows
        db_summary_rows = ""
        for db, result in results.items():
            status_icon = "OK" if result['status'] == 'success' else "FAIL"
            db_summary_rows += f"""
                    <tr>
                        <td>{db}</td>
                        <td>{result['hit_count']}</td>
                        <td>{status_icon} {result['status']}</td>
                    </tr>
            """

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PseudoScope ABRicate - Comprehensive Report - {genome_name}</title>
    <style>{self._get_shared_css()}</style>
    {self._get_quote_js()}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.ascii_art}</div>
            </div>
            <div class="quote-container">
                <div id="science-quote" style="font-size: 1.1em;"></div>
            </div>
        </div>

        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #8B4513; padding-bottom: 10px;">📊 P. aeruginosa AMR Summary</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total Genes</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['total_hits']}</p>
                </div>
                <div class="stat-card">
                    <h3>Critical Resistance</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['total_critical_resistance']}</p>
                </div>
                <div class="stat-card">
                    <h3>Virulence Factors</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['total_high_risk_virulence'] + analysis['total_critical_virulence']}</p>
                </div>
            </div>
            <p><strong>Genome:</strong> {genome_name}</p>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>Thresholds:</strong> minID={self.minid}%, minCov={self.mincov}%</p>
        </div>

        {alerts_html}
        {resistance_classes_html}
        {tables_html}

        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #8B4513; padding-bottom: 10px;">🗃️ Database Results Summary</h2>
            <table class="gene-table">
                <thead>
                    <tr><th>Database</th><th>Hits</th><th>Status</th></tr>
                </thead>
                <tbody>
                    {db_summary_rows}
                </tbody>
            </table>
        </div>

        {self._get_footer_html("PseudoScope - Comprehensive ABRicate Report")}
    </div>
</body>
</html>"""

        html_file = os.path.join(output_dir, f"{genome_name}_comprehensive_abricate_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info("Comprehensive HTML report generated: %s", html_file)

    # ───────────────────── batch-level TSV + HTML summaries ─────────────
    def create_database_summaries(self, all_results: Dict[str, Any], output_base: str):
        """Write per-database TSV and HTML batch summaries (NO per-database JSON)."""
        self.logger.info("Creating P. aeruginosa database summary files and HTML reports...")
        db_results: Dict[str, List[Dict]] = {}
        for genome_name, genome_result in all_results.items():
            for db, db_result in genome_result['results'].items():
                db_results.setdefault(db, [])
                for hit in db_result['hits']:
                    h = hit.copy()
                    h['genome'] = genome_name
                    db_results[db].append(h)

        for db, hits in db_results.items():
            if not hits:
                self.logger.info("No hits for database %s, skipping summary", db)
                continue
            # TSV summary
            summary_file = os.path.join(output_base, f"pseudo_{db}_abricate_summary.tsv")
            headers = list(hits[0].keys())
            with open(summary_file, 'w') as f:
                f.write('\t'.join(headers) + '\n')
                for hit in hits:
                    row = [str(hit.get(h, '')) for h in headers]
                    f.write('\t'.join(row) + '\n')
            self.logger.info("Created %s summary: %s (%d hits)", db, summary_file, len(hits))
            # HTML summary
            self._create_database_summary_html(db, hits, output_base)

    def _create_database_summary_html(self, database: str, hits: List[Dict], output_base: str):
        unique_genomes = sorted(set(h['genome'] for h in hits))
        total_genomes = len(unique_genomes)

        genes_per_genome: Dict[str, set] = {}
        for hit in hits:
            genes_per_genome.setdefault(hit['genome'], set()).add(hit['gene'])

        genes_by_genome_rows = ""
        for genome in unique_genomes:
            genes = genes_per_genome.get(genome, set())
            gene_list = ", ".join(sorted(genes))
            genes_by_genome_rows += f"""
                    <tr class="present">
                        <td><strong>{genome}</strong></td>
                        <td>{len(genes)}</td>
                        <td>{gene_list}</td>
                    </tr>
            """

        gene_frequency: Dict[str, set] = {}
        for hit in hits:
            gene_frequency.setdefault(hit['gene'], set()).add(hit['genome'])

        frequency_rows = ""
        for gene, genomes in sorted(gene_frequency.items(), key=lambda x: len(x[1]), reverse=True):
            genome_list = ", ".join(sorted(genomes))
            count = len(genomes)
            pct = (count / total_genomes * 100) if total_genomes else 0
            frequency_rows += f"""
                    <tr>
                        <td><strong>{gene}</strong></td>
                        <td class="frequency-cell">{count}
                            <span class="percentage-badge">({pct:.1f}%)</span>
                        </td>
                        <td>{genome_list}</td>
                    </tr>
            """

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PseudoScope ABRicate - Batch Summary - {database.upper()}</title>
    <style>{self._get_shared_css()}</style>
    {self._get_quote_js()}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.ascii_art}</div>
            </div>
            <div class="quote-container">
                <div id="science-quote" style="font-size: 1.1em;"></div>
            </div>
        </div>

        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #8B4513; padding-bottom: 10px;">📊 {database.upper()} Database Batch Summary</h2>
            <div class="summary-stats">
                <div class="stat-card"><h3>Total Hits</h3><p style="font-size: 2em; margin: 0;">{len(hits)}</p></div>
                <div class="stat-card"><h3>Genomes</h3><p style="font-size: 2em; margin: 0;">{total_genomes}</p></div>
                <div class="stat-card"><h3>Unique Genes</h3><p style="font-size: 2em; margin: 0;">{len(gene_frequency)}</p></div>
            </div>
            <p><strong>Analysis date:</strong> {self.metadata['analysis_date']}</p>
        </div>

        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #8B4513; padding-bottom: 10px;">🔍 Genes by Genome</h2>
            <table class="gene-table">
                <thead><tr><th>Genome</th><th>Gene Count</th><th>Genes Detected</th></tr></thead>
                <tbody>{genes_by_genome_rows}</tbody>
            </table>
        </div>

        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #8B4513; padding-bottom: 10px;">📈 Gene Frequency</h2>
            <table class="gene-table">
                <thead><tr><th>Gene</th><th>Frequency</th><th>Genomes</th></tr></thead>
                <tbody>{frequency_rows}</tbody>
            </table>
        </div>

        {self._get_footer_html(f"PseudoScope - {database.upper()} Batch Summary")}
    </div>
</body>
</html>"""
        html_file = os.path.join(output_base, f"pseudo_{database}_summary_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info("Database summary HTML report: %s", html_file)

    # ───────────────────── database detailed helper (for master) ────────
    def _compute_database_detailed_summary(self, database: str, hits: List[Dict]) -> Dict[str, Any]:
        """Compute the detailed per-database summary that used to live in its own JSON.

        This is now folded into the master JSON so that we only emit ONE JSON
        at the summary level.
        """
        unique_genomes = list(set(h['genome'] for h in hits))
        unique_genes = list(set(h['gene'] for h in hits))

        gene_frequency: Dict[str, set] = {}
        gene_details: Dict[str, Dict[str, Any]] = {}
        for hit in hits:
            gene = hit['gene']
            if gene not in gene_frequency:
                gene_frequency[gene] = set()
                gene_details[gene] = {
                    'product': hit['product'],
                    'databases': set(),
                    'coverage_sum': 0.0,
                    'identity_sum': 0.0,
                    'count': 0,
                }
            gene_frequency[gene].add(hit['genome'])
            gene_details[gene]['databases'].add(hit['database'])
            try:
                cov = float(str(hit['coverage_percent']).replace('%', ''))
                idd = float(str(hit['identity_percent']).replace('%', ''))
                gene_details[gene]['coverage_sum'] += cov
                gene_details[gene]['identity_sum'] += idd
                gene_details[gene]['count'] += 1
            except (ValueError, TypeError):
                pass

        for gene, d in gene_details.items():
            n = d['count'] or 1
            d['avg_coverage'] = round(d['coverage_sum'] / n, 2)
            d['avg_identity'] = round(d['identity_sum'] / n, 2)
            del d['coverage_sum'], d['identity_sum'], d['count']
            d['databases'] = sorted(d['databases'])

        genomes_summary: Dict[str, Dict[str, Any]] = {}
        for genome in unique_genomes:
            genome_hits = [h for h in hits if h['genome'] == genome]
            genomes_summary[genome] = {
                'total_hits': len(genome_hits),
                'unique_genes': len(set(h['gene'] for h in genome_hits)),
                'genes': sorted(set(h['gene'] for h in genome_hits)),
            }

        risk_genes: Dict[str, List[str]] = {
            'critical_carbapenemase': [],
            'critical_esbl': [],
            'critical_colistin': [],
            'critical_aminoglycoside': [],
            'high_risk_resistance': [],
            'virulence_factors': [],
        }
        for gene in unique_genes:
            gene_lower = gene.lower()
            if any(c in gene_lower for c in ['oxa', 'imp', 'vim', 'ndm', 'kpc', 'ges', 'spm', 'aim', 'dim']):
                if any(c.lower() in gene_lower for c in self.critical_carbapenemases):
                    risk_genes['critical_carbapenemase'].append(gene)
                    continue
            if any(e in gene_lower for e in ['per', 'veb', 'bel', 'ges', 'tem', 'shv', 'ctx']):
                if any(c.lower() in gene_lower for c in self.critical_esbls):
                    risk_genes['critical_esbl'].append(gene)
                    continue
            if 'mcr' in gene_lower or any(p in gene_lower for p in ['pmr', 'pho', 'mgr', 'lpx', 'arn', 'ept', 'bas']):
                if any(c.lower() in gene_lower for c in self.critical_colistin):
                    risk_genes['critical_colistin'].append(gene)
                    continue
            if any(a in gene_lower for a in ['arm', 'rmt', 'npm', 'aac', 'ant', 'aph', 'str']):
                if any(c.lower() in gene_lower for c in self.critical_aminoglycoside):
                    risk_genes['critical_aminoglycoside'].append(gene)
                    continue
            if any(h.lower() in gene_lower for h in self.high_risk_resistance):
                risk_genes['high_risk_resistance'].append(gene)
                continue
            if any(v.lower() in gene_lower for v in self.virulence_genes):
                risk_genes['virulence_factors'].append(gene)

        return {
            'unique_genes': len(unique_genes),
            'gene_analysis': {
                'most_frequent_genes': {
                    g: len(gs) for g, gs in sorted(
                        gene_frequency.items(), key=lambda x: len(x[1]), reverse=True
                    )
                },
                'gene_details': gene_details,
            },
            'genome_analysis': genomes_summary,
            'risk_genes': risk_genes,
        }

    # ───────────────────────── master JSON summary ─────────────────────
    def create_master_json_summary(self, all_results: Dict[str, Any], output_base: str):
        """Write the ONE master JSON summarising every database & genome.

        Per-database JSONs are intentionally NOT emitted — all of their useful
        data is embedded under ``database_summaries[<db>]`` here.
        """
        self.logger.info("Creating master JSON summary (single consolidated JSON)...")
        master_report = {
            "metadata": {
                **self.metadata,
                "analysis_type": "P. aeruginosa Comprehensive Antimicrobial Resistance & Virulence Analysis",
                "databases_analyzed": self.required_databases,
                "total_genomes": len(all_results),
                "analysis_parameters": {
                    "minimum_identity": self.minid,
                    "minimum_coverage": self.mincov,
                },
            },
            "overall_summary": {
                "total_genomes": len(all_results),
                "total_hits": 0,
                "critical_carbapenemase_genes": 0,
                "critical_esbl_genes": 0,
                "critical_colistin_genes": 0,
                "critical_aminoglycoside_genes": 0,
                "carbapenemase_positive_genomes": 0,
                "esbl_positive_genomes": 0,
                "colistin_resistant_genomes": 0,
            },
            "database_summaries": {},
            "genome_summaries": {},
            "risk_assessment": {
                "critical_carbapenemase_genes_found": [],
                "critical_esbl_genes_found": [],
                "critical_colistin_genes_found": [],
                "critical_aminoglycoside_genes_found": [],
                "high_risk_genes_found": [],
                "virulence_genes_found": [],
                "resistance_classes_found": defaultdict(list),
            },
            "detailed_results": {},
        }

        for genome_name, genome_result in all_results.items():
            all_genome_hits = []
            for db_result in genome_result['results'].values():
                all_genome_hits.extend(db_result['hits'])
            analysis = self.analyze_pseudo_resistance(all_genome_hits)

            master_report["overall_summary"]["total_hits"] += analysis['total_hits']
            master_report["overall_summary"]["critical_carbapenemase_genes"] += len(analysis['critical_carbapenemase_genes'])
            master_report["overall_summary"]["critical_esbl_genes"] += len(analysis['critical_esbl_genes'])
            master_report["overall_summary"]["critical_colistin_genes"] += len(analysis['critical_colistin_genes'])
            master_report["overall_summary"]["critical_aminoglycoside_genes"] += len(analysis['critical_aminoglycoside_genes'])
            if analysis['carbapenemase_status'] == 'positive':
                master_report["overall_summary"]["carbapenemase_positive_genomes"] += 1
            if analysis['esbl_status'] == 'positive':
                master_report["overall_summary"]["esbl_positive_genomes"] += 1
            if analysis['colistin_resistance'] == 'positive':
                master_report["overall_summary"]["colistin_resistant_genomes"] += 1

            master_report["genome_summaries"][genome_name] = {
                "total_hits": analysis['total_hits'],
                "critical_carbapenemase": len(analysis['critical_carbapenemase_genes']),
                "critical_esbl": len(analysis['critical_esbl_genes']),
                "critical_colistin": len(analysis['critical_colistin_genes']),
                "critical_aminoglycoside": len(analysis['critical_aminoglycoside_genes']),
                "high_risk_genes": analysis['total_high_risk_resistance'],
                "virulence_factors": analysis['total_high_risk_virulence'] + analysis['total_critical_virulence'],
                "carbapenemase_status": analysis['carbapenemase_status'],
                "esbl_status": analysis['esbl_status'],
                "colistin_resistance": analysis['colistin_resistance'],
                "database_hits": {db: r['hit_count'] for db, r in genome_result['results'].items()},
            }

            for key_analysis, key_master in [
                ('critical_carbapenemase_genes', 'critical_carbapenemase_genes_found'),
                ('critical_esbl_genes', 'critical_esbl_genes_found'),
                ('critical_colistin_genes', 'critical_colistin_genes_found'),
                ('critical_aminoglycoside_genes', 'critical_aminoglycoside_genes_found'),
                ('high_risk_resistance_genes', 'high_risk_genes_found'),
            ]:
                for g in analysis[key_analysis]:
                    if g['gene'] not in master_report["risk_assessment"][key_master]:
                        master_report["risk_assessment"][key_master].append(g['gene'])

            for g in analysis['high_risk_virulence_genes'] + analysis['critical_virulence_genes']:
                if g['gene'] not in master_report["risk_assessment"]["virulence_genes_found"]:
                    master_report["risk_assessment"]["virulence_genes_found"].append(g['gene'])

            for class_name, genes in analysis['resistance_classes'].items():
                for g in genes:
                    if g['gene'] not in [x['gene'] for x in master_report["risk_assessment"]["resistance_classes_found"][class_name]]:
                        master_report["risk_assessment"]["resistance_classes_found"][class_name].append(g)

            master_report["detailed_results"][genome_name] = {
                "total_hits": genome_result['total_hits'],
                "database_results": genome_result['results'],
            }

        # ─── Per-database summary block (enriched; replaces separate JSONs) ───
        for db in self.required_databases:
            db_hits = []
            for genome_name, genome_result in all_results.items():
                if db in genome_result['results']:
                    for hit in genome_result['results'][db]['hits']:
                        h = hit.copy()
                        h['genome'] = genome_name
                        db_hits.append(h)

            if db_hits:
                unique_genes = len(set(h['gene'] for h in db_hits))
                unique_genomes = len(set(h['genome'] for h in db_hits))
                detail = self._compute_database_detailed_summary(db, db_hits)
                master_report["database_summaries"][db] = {
                    "total_hits": len(db_hits),
                    "unique_genes": unique_genes,
                    "unique_genomes": unique_genomes,
                    "hits_per_genome": round(len(db_hits) / unique_genomes, 2) if unique_genomes else 0,
                    "most_frequent_genes": dict(Counter(h['gene'] for h in db_hits).most_common(10)),
                    **detail,
                }
            else:
                master_report["database_summaries"][db] = {
                    "total_hits": 0,
                    "unique_genes": 0,
                    "unique_genomes": 0,
                    "hits_per_genome": 0,
                    "most_frequent_genes": {},
                    "gene_analysis": {"most_frequent_genes": {}, "gene_details": {}},
                    "genome_analysis": {},
                    "risk_genes": {},
                }

        for key in ("critical_carbapenemase_genes_found", "critical_esbl_genes_found",
                    "critical_colistin_genes_found", "critical_aminoglycoside_genes_found",
                    "high_risk_genes_found", "virulence_genes_found"):
            master_report["risk_assessment"][key].sort()
        master_report["risk_assessment"]["resistance_classes_found"] = dict(
            master_report["risk_assessment"]["resistance_classes_found"]
        )

        master_report["performance_metrics"] = {
            "cpu_cores_used": self.cpus,
            "available_ram_gb": round(self.available_ram, 2),
            "processing_mode": "MAXIMUM SPEED" if self.cpus > 8 else "STANDARD",
        }

        # Optional Excel export (kept because it's a different format; the JSON is still the ONE json)
        if pd is not None:
            try:
                excel_data = []
                for genome_name, genome_result in all_results.items():
                    all_genome_hits = []
                    for db_result in genome_result['results'].values():
                        all_genome_hits.extend(db_result['hits'])
                    a = self.analyze_pseudo_resistance(all_genome_hits)
                    excel_data.append({
                        'Genome': genome_name,
                        'Total_Hits': genome_result['total_hits'],
                        'Critical_Carbapenemase': len(a['critical_carbapenemase_genes']),
                        'Critical_ESBL': len(a['critical_esbl_genes']),
                        'Critical_Colistin': len(a['critical_colistin_genes']),
                        'Critical_Aminoglycoside': len(a['critical_aminoglycoside_genes']),
                        'High_Risk_Resistance': len(a['high_risk_resistance_genes']),
                        'Virulence_Factors': len(a['high_risk_virulence_genes']) + len(a['critical_virulence_genes']),
                        'Carbapenemase_Status': a['carbapenemase_status'],
                        'ESBL_Status': a['esbl_status'],
                        'Colistin_Resistance': a['colistin_resistance'],
                    })
                df = pd.DataFrame(excel_data)
                excel_file = os.path.join(output_base, "pseudo_abricate_master_summary.xlsx")
                df.to_excel(excel_file, index=False)
                master_report["excel_summary_file"] = excel_file
                self.logger.info("Created Excel summary: %s", excel_file)
            except Exception as e:
                self.logger.warning("Could not create Excel summary: %s", e)

        master_file = os.path.join(output_base, "pseudo_abricate_master_summary.json")
        with open(master_file, 'w') as f:
            json.dump(master_report, f, indent=4, default=str)
        self.logger.info("Created master JSON summary: %s", master_file)
        return master_file

    # ─────────────────────────── per-genome pipeline ────────────────────
    def process_single_genome(self, genome_file: str,
                              output_base: str = "pseudo_abricate_results") -> Dict[str, Any]:
        """Process one P. aeruginosa genome. Emits exactly 2 files:
            1. {genome}_comprehensive_abricate_report.html
            2. abricate_raw.tsv
        """
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        self.logger.info("=== PROCESSING P. AERUGINOSA GENOME: %s ===", genome_name)
        os.makedirs(results_dir, exist_ok=True)

        results: Dict[str, Dict] = {}

        # All raw per-database .txt files live here and vanish on exit
        with tempfile.TemporaryDirectory() as temp_dir:
            for db in self.required_databases:
                result = self.run_abricate_single_db(genome_file, db, temp_dir)
                results[db] = result
                icon = "OK" if result['status'] == 'success' else "FAIL"
                self.logger.info("%s %s: %d hits", icon, db, result['hit_count'])

            # Single merged raw TSV escapes the temp dir
            self._write_merged_raw_tsv(results, results_dir)

        # Free raw text from memory before building HTML
        for db_result in results.values():
            db_result.pop('raw_content', None)

        # Single comprehensive HTML report
        self.create_comprehensive_html_report(genome_name, results, results_dir)

        return {
            'genome': genome_name,
            'results': results,
            'total_hits': sum(r['hit_count'] for r in results.values()),
        }

    def process_multiple_genomes(self, genome_pattern: str,
                                 output_base: str = "pseudo_abricate_results") -> Dict[str, Any]:
        if not self.check_abricate_installed():
            raise RuntimeError("ABRicate not properly installed")
        self.setup_abricate_databases()

        fasta_patterns = [genome_pattern, f"{genome_pattern}.fasta",
                          f"{genome_pattern}.fa", f"{genome_pattern}.fna",
                          f"{genome_pattern}.faa"]
        genome_files = []
        for pattern in fasta_patterns:
            genome_files.extend(glob.glob(pattern))
        genome_files = list(set(genome_files))

        if not genome_files:
            raise FileNotFoundError(f"No FASTA files found matching pattern: {genome_pattern}")

        self.logger.info("Found %d P. aeruginosa genomes: %s",
                         len(genome_files), [Path(f).name for f in genome_files])
        os.makedirs(output_base, exist_ok=True)

        all_results: Dict[str, Any] = {}
        if len(genome_files) > 1 and self.cpus > 1:
            self.logger.info("Using parallel processing with %d CPU cores", self.cpus)
            with ThreadPoolExecutor(max_workers=self.cpus) as executor:
                future_to_genome = {
                    executor.submit(self.process_single_genome, g, output_base): g
                    for g in genome_files
                }
                for future in as_completed(future_to_genome):
                    genome = future_to_genome[future]
                    try:
                        result = future.result()
                        all_results[Path(genome).stem] = result
                        self.logger.info("Completed: %s (%d total hits)",
                                         result['genome'], result['total_hits'])
                    except Exception as e:
                        self.logger.error("Failed: %s - %s", genome, e)
        else:
            for genome in genome_files:
                try:
                    result = self.process_single_genome(genome, output_base)
                    all_results[Path(genome).stem] = result
                    self.logger.info("Completed: %s (%d total hits)",
                                     result['genome'], result['total_hits'])
                except Exception as e:
                    self.logger.error("Failed: %s - %s", genome, e)

        self.create_database_summaries(all_results, output_base)
        master_json = self.create_master_json_summary(all_results, output_base)

        self.logger.info("=== P. AERUGINOSA ANALYSIS COMPLETE ===")
        self.logger.info("Processed %d genomes", len(all_results))
        self.logger.info("Results saved to: %s", output_base)
        self.logger.info("Master JSON summary: %s", master_json)
        return all_results


def main():
    """Command line interface for P. aeruginosa ABRicate analysis"""
    parser = argparse.ArgumentParser(
        description='PseudoScope ABRicate Analysis - Consolidated Output for Pseudomonas aeruginosa',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python p_abricate.py "*.fna"
  python p_abricate.py "PA_*.fasta"
  python p_abricate.py "*.fa" --output my_pseudo_results
  python p_abricate.py "*.fna" --cpus 4
  python p_abricate.py "*.fna" --minid 85 --mincov 80

OUTPUTS (consolidated):
  Per genome (2 files each):
    - {genome}_comprehensive_abricate_report.html
    - abricate_raw.tsv  (all databases merged)

  Batch-level:
    - pseudo_{db}_abricate_summary.tsv     (per database)
    - pseudo_{db}_summary_report.html      (per database)
    - pseudo_abricate_master_summary.json  (the ONE summary JSON)
    - pseudo_abricate_master_summary.xlsx  (optional convenience export)
        """
    )
    parser.add_argument('pattern', help='File pattern for P. aeruginosa genomes (e.g., "*.fasta")')
    parser.add_argument('--cpus', '-c', type=int, default=None,
                        help='Number of CPU cores to use (default: auto-detect)')
    parser.add_argument('--output', '-o', default='pseudo_abricate_results',
                        help='Output directory (default: pseudo_abricate_results)')
    parser.add_argument('--minid', type=int, default=80,
                        help='Minimum percent identity (default: 80)')
    parser.add_argument('--mincov', type=int, default=80,
                        help='Minimum percent coverage (default: 80)')
    args = parser.parse_args()

    if not 0 <= args.minid <= 100:
        parser.error("--minid must be between 0 and 100")
    if not 0 <= args.mincov <= 100:
        parser.error("--mincov must be between 0 and 100")

    print("\n" + "=" * 80)
    print(r"""
██████╗ ███████╗███████╗██╗   ██╗██████╗  ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██╔════╝██║   ██║██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
██████╔╝███████╗█████╗  ██║   ██║██║  ██║██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗
██╔═══╝ ╚════██║██╔══╝  ██║   ██║██║  ██║██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝
██║     ███████║███████╗╚██████╔╝██████╔╝╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝     ╚══════╝╚══════╝ ╚═════╝ ╚═════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
""")
    print("PseudoScope ABRicate - Consolidated Output Version")
    print("=" * 80)

    executor = PseudoAbricateExecutor(cpus=args.cpus, minid=args.minid, mincov=args.mincov)
    try:
        results = executor.process_multiple_genomes(args.pattern, args.output)

        executor.logger.info("\n" + "=" * 50)
        executor.logger.info("PseudoScope ABRicate FINAL SUMMARY")
        executor.logger.info("=" * 50)
        total_carb = total_esbl = total_col = 0
        for genome_name, result in results.items():
            all_hits = []
            for db_result in result['results'].values():
                all_hits.extend(db_result['hits'])
            a = executor.analyze_pseudo_resistance(all_hits)
            executor.logger.info(
                "%s: %d hits | Carb=%d ESBL=%d Colistin=%d",
                genome_name, result['total_hits'],
                len(a['critical_carbapenemase_genes']),
                len(a['critical_esbl_genes']),
                len(a['critical_colistin_genes']),
            )
            total_carb += len(a['critical_carbapenemase_genes'])
            total_esbl += len(a['critical_esbl_genes'])
            total_col += len(a['critical_colistin_genes'])

        executor.logger.info("\n" + "=" * 50)
        executor.logger.info("CONSOLIDATED OUTPUT STRUCTURE")
        executor.logger.info("=" * 50)
        executor.logger.info("%s/", args.output)
        executor.logger.info("├── [SAMPLE_DIRECTORIES]/")
        executor.logger.info("│   ├── abricate_raw.tsv                       (all DBs merged)")
        executor.logger.info("│   └── [SAMPLE]_comprehensive_abricate_report.html")
        executor.logger.info("├── pseudo_*_abricate_summary.tsv              (per DB)")
        executor.logger.info("├── pseudo_*_summary_report.html               (per DB)")
        executor.logger.info("└── pseudo_abricate_master_summary.json        (ONE summary JSON)")
        executor.logger.info("    pseudo_abricate_master_summary.xlsx        (optional)")

        executor.logger.info("\n" + "=" * 50)
        executor.logger.info("PERFORMANCE SUMMARY")
        executor.logger.info("=" * 50)
        executor.logger.info("CPU cores utilized: %d", executor.cpus)
        executor.logger.info("Available RAM: %.1f GB", executor.available_ram)
        executor.logger.info("Genomes processed: %d", len(results))
        executor.logger.info("Total critical carbapenemase genes: %d", total_carb)
        executor.logger.info("Total critical ESBL genes: %d", total_esbl)
        executor.logger.info("Total critical colistin genes: %d", total_col)

        executor.logger.info("\n%s", random.choice(executor.science_quotes))

    except Exception as e:
        executor.logger.error("P. aeruginosa analysis failed: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()