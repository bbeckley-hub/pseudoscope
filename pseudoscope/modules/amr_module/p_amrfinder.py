#!/usr/bin/env python3
"""
PseudoScope AMRfinderPlus - P. aeruginosa AMR Analysis with Mutation Reporting
Comprehensive AMR analysis for Pseudomonas aeruginosa with beautiful HTML reporting - INTERACTIVE ENHANCED VERSION
DYNAMIC DATABASE with auto-update capability
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2025-12-28 / Updated 2026-06-24
Version: 1.1.0
"""

import subprocess
import sys
import os
import glob
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any, Optional
import argparse
import re
from datetime import datetime
import psutil
import math
import json
from collections import defaultdict
import csv
import base64
import random

class PseudoAMRfinderPlus:
    """AMRfinderPlus executor for P. aeruginosa with DYNAMIC database detection and update capability"""

    def __init__(self, cpus: int = None):
        self.logger = self._setup_logging()
        self.module_dir = os.path.dirname(os.path.abspath(__file__))
        self.available_ram = self._get_available_ram()
        self.cpus = self._calculate_optimal_cpus(cpus)

        self.bundled_amrfinder = os.path.join(self.module_dir, "bin", "amrfinder")
        self.bundled_update = os.path.join(self.module_dir, "bin", "amrfinder_update")
        self.bundled_database = self._get_latest_database()

        if self.bundled_database is None:
            self.logger.warning("No AMRfinderPlus database found. Please run with --update-db to download.")

        db_version = self._get_database_version() if self.bundled_database else "Unknown"

        self.metadata = {
            "tool_name": "PseudoScope AMRfinderPlus",
            "version": "1.1.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School - Department of Medical Biochemistry",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "amrfinder_version": "4.2.7",
            "database_version": db_version
        }

        # P. aeruginosa specific gene sets - COMPREHENSIVE
        self.critical_carbapenemases = {
            'blaOXA-2', 'blaOXA-10', 'blaOXA-14', 'blaOXA-17', 'blaOXA-19', 'blaOXA-28', 'blaOXA-35', 'blaOXA-45', 'blaOXA-50', 'blaOXA-198',
            'OXA-2', 'OXA-10', 'OXA-14', 'OXA-17', 'OXA-19', 'OXA-28', 'OXA-35', 'OXA-45', 'OXA-50', 'OXA-198',
            'blaIMP', 'blaVIM', 'blaNDM', 'blaKPC', 'blaGES', 'blaSPM', 'blaAIM', 'blaDIM',
            'IMP', 'VIM', 'NDM', 'KPC', 'GES', 'SPM', 'AIM', 'DIM'
        }

        self.critical_esbls = {
            'blaPER', 'blaVEB', 'blaBEL', 'blaGES', 'blaTEM', 'blaSHV', 'blaCTX-M',
            'PER', 'VEB', 'BEL', 'GES', 'TEM', 'SHV', 'CTX-M'
        }

        self.critical_colistin = {
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10',
            'pmrA', 'pmrB', 'phoP', 'phoQ', 'mgrB', 'lpxA', 'lpxC', 'lpxD', 'arnT', 'eptA', 'eptB', 'basS', 'basR'
        }

        self.critical_aminoglycoside = {
            'armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD', 'rmtE', 'rmtF', 'rmtG', 'rmtH', 'npmA',
            'aac(3)', 'aac(6\')', 'ant(2")', 'ant(4\')', 'aph(3\')', 'aph(6)',
            'aacC1', 'aacC2', 'aacC4', 'aacA4', 'aacA7', 'aadA1', 'aadA2', 'aadA5', 'aadA7',
            'strA', 'strB', 'aphA1', 'aphA2', 'aphA3', 'aphA6', 'aac3', 'aac6', 'aadA', 'aadB'
        }

        self.high_risk_resistance = {
            'mexA', 'mexB', 'mexC', 'mexD', 'mexE', 'mexF', 'mexX', 'mexY',
            'oprM', 'oprN', 'oprJ',
            'adeA', 'adeB', 'adeC',
            'tetA', 'tetB', 'tetC', 'tetD', 'tetE', 'tetG', 'tetH', 'tetK', 'tetL', 'tetM', 'tetO', 'tetQ', 'tetS', 'tetX',
            'sul1', 'sul2', 'sul3', 'sul4',
            'dfrA1', 'dfrA5', 'dfrA7', 'dfrA8', 'dfrA12', 'dfrA14', 'dfrA17', 'dfrA19', 'dfrA20', 'dfrA21',
            'dfrB1', 'dfrB2', 'dfrB3', 'dfrB4', 'dfrB5', 'dfrB6', 'dfrB7',
            'catA1', 'catA2', 'catB2', 'catB3', 'catB8', 'catI', 'catII', 'catIII',
            'cmlA', 'cmlA1', 'cmlA5', 'cmlA6', 'cmlA7', 'floR',
            'ermA', 'ermB', 'ermC', 'ermF', 'ermG', 'ermX', 'ermY',
            'mphA', 'mphB', 'mphC', 'mphD', 'mphE', 'msrA', 'msrB', 'msrC', 'msrD',
            'qnrA1', 'qnrB1', 'qnrB2', 'qnrB4', 'qnrB6', 'qnrB10', 'qnrB19', 'qnrS1', 'qnrS2', 'qnrVC1', 'qnrVC4',
            'aac(6\')-Ib-cr', 'qepA', 'qepA1', 'qepA2', 'qepA3', 'qepA4',
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
            'pvdA', 'pvdD', 'pvdE', 'pvdF', 'pvdG', 'pvdH', 'pvdI', 'pvdJ', 'pvdL', 'pvdN', 'pvdO', 'pvdP', 'pvdQ',
            'fpvA', 'fpvB', 'fpvC', 'fpvD', 'fpvE', 'fpvF', 'fpvG', 'fpvH', 'fpvI', 'fpvJ', 'fpvK',
            'fliC', 'fliD', 'fliE', 'fliF', 'fliG', 'fliH', 'fliI', 'fliJ', 'fliK', 'fliL', 'fliM', 'fliN', 'fliO', 'fliP', 'fliQ', 'fliR',
            'fleN', 'fleQ', 'fleR',
            'pilA', 'pilB', 'pilC', 'pilD', 'pilE', 'pilF', 'pilG', 'pilH', 'pilI', 'pilJ', 'pilK', 'pilL', 'pilM', 'pilN', 'pilO', 'pilP', 'pilQ', 'pilR', 'pilS', 'pilT', 'pilU',
            'pslA', 'pslB', 'pslC', 'pslD', 'pslE', 'pslF', 'pslG', 'pslH', 'pslI', 'pslJ', 'pslK', 'pslL', 'pslM', 'pslN', 'pslO', 'pslP',
            'pelA', 'pelB', 'pelC', 'pelD', 'pelE', 'pelF', 'pelG',
            'lpxA', 'lpxC', 'lpxD', 'lpxL', 'lpxO', 'lpxP',
            'lasI', 'lasR', 'rhlI', 'rhlR', 'pqsA', 'pqsB', 'pqsC', 'pqsD', 'pqsE', 'pqsH', 'mvfR'
        }

        self.high_risk_genes = self.critical_carbapenemases.union(
            self.critical_esbls
        ).union(self.critical_colistin).union(self.critical_aminoglycoside).union(self.high_risk_resistance)

        self.critical_risk_genes = self.critical_carbapenemases.union(self.critical_colistin)

        self.ascii_art = r"""
██████╗ ███████╗███████╗██╗   ██╗██████╗  ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██╔════╝██║   ██║██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
██████╔╝███████╗█████╗  ██║   ██║██║  ██║██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔═══╝ ╚════██║██╔══╝  ██║   ██║██║  ██║██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║     ███████║███████╗╚██████╔╝██████╔╝╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝     ╚══════╝╚══════╝ ╚═════╝ ╚═════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
"""

        self.science_quotes = [
            {"text": "Pseudomonas aeruginosa: the master of adaptation.", "author": "Unknown"},
            {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"},
            {"text": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"},
            {"text": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"},
            {"text": "In science, there are no shortcuts to truth.", "author": "Karl Popper"},
            {"text": "Science knows no country, because knowledge belongs to humanity.", "author": "Louis Pasteur"},
            {"text": "The science of today is the technology of tomorrow.", "author": "Edward Teller"},
            {"text": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie"},
            {"text": "Research is what I'm doing when I don't know what I'm doing.", "author": "Wernher von Braun"}
        ]

    def _setup_logging(self):
        """Setup logging - must be called first in __init__"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)

    def _get_available_ram(self) -> int:
        """Get available RAM in GB"""
        try:
            ram_gb = psutil.virtual_memory().available / (1024 ** 3)
            return ram_gb
        except Exception as e:
            self.logger.warning(f"Could not detect RAM: {e}")
            return 8

    def _calculate_optimal_cpus(self, user_cpus: int = None) -> int:
        """Calculate optimal number of CPU cores for MAXIMUM SPEED"""
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
                optimal_cpus = max(16, total_physical_cores - 1)
            else:
                optimal_cpus = min(32, int(total_physical_cores * 0.95))
            optimal_cpus = max(1, min(optimal_cpus, total_physical_cores))
            self._log_resource_info(optimal_cpus, total_physical_cores)
            return optimal_cpus
        except Exception as e:
            self.logger.warning(f"Could not detect CPU cores, using maximum available: {e}")
            return os.cpu_count() or 4

    def _log_resource_info(self, cpus: int, total_cores: int = None):
        """Log resource allocation information"""
        self.logger.info(f"Available RAM: {self.available_ram:.1f} GB")
        if total_cores:
            self.logger.info(f"System CPU cores: {total_cores}")
            utilization = (cpus / total_cores) * 100
            self.logger.info(f"Using CPU cores: {cpus} ({utilization:.1f}% of available cores)")
        else:
            self.logger.info(f"Using user-specified CPU cores: {cpus}")
        if cpus == 1:
            self.logger.info("💡 Performance: Single-core (max speed for 1-core systems)")
        elif cpus <= 4:
            self.logger.info("💡 Performance: Multi-core (max speed for small systems)")
        elif cpus <= 8:
            self.logger.info("💡 Performance: High-speed multi-core mode")
        elif cpus <= 16:
            self.logger.info("💡 Performance: Ultra-speed multi-core mode 🚀")
        elif cpus <= 32:
            self.logger.info("💡 Performance: MAXIMUM SPEED MULTI-CORE MODE 🚀🔥")
        else:
            self.logger.info("💡 Performance: EXTREME SPEED MULTI-CORE MODE 🚀🔥💨")
        self.logger.info("📝 STRATEGY: Processing MULTIPLE P. aeruginosa samples concurrently with optimal core allocation for maximum throughput")

    def _get_latest_database(self) -> Optional[str]:
        """Find the latest dated database folder in data/amrfinder_db/ (starts with 20)"""
        db_root = os.path.join(self.module_dir, "data", "amrfinder_db")
        if not os.path.exists(db_root):
            self.logger.warning(f"Database root directory not found: {db_root}")
            return None
        candidates = []
        for item in os.listdir(db_root):
            full_path = os.path.join(db_root, item)
            if os.path.isdir(full_path) and item.startswith('20'):
                candidates.append(item)
        if not candidates:
            self.logger.warning("No database folder starting with '20' found.")
            return None
        latest = sorted(candidates)[-1]
        latest_path = os.path.join(db_root, latest)
        self.logger.info(f"Using latest database: {latest_path}")
        return latest_path

    def _get_database_version(self) -> str:
        """Read version.txt from the database folder or fallback to folder name"""
        if not self.bundled_database:
            return "Unknown"
        version_file = os.path.join(self.bundled_database, "version.txt")
        if os.path.exists(version_file):
            with open(version_file, 'r') as f:
                return f.read().strip()
        return os.path.basename(self.bundled_database)

    def update_database(self) -> bool:
        """Download the latest AMRfinderPlus database using bundled amrfinder_update"""
        if not os.path.exists(self.bundled_update):
            self.logger.error(f"amrfinder_update not found at {self.bundled_update}")
            return False
        if not os.access(self.bundled_update, os.X_OK):
            self.logger.warning("amrfinder_update not executable, fixing permissions...")
            os.chmod(self.bundled_update, 0o755)
        db_dir = os.path.join(self.module_dir, "data", "amrfinder_db")
        os.makedirs(db_dir, exist_ok=True)
        self.logger.info("Updating AMRfinderPlus database...")
        try:
            cmd = [self.bundled_update, "--database", db_dir]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            self.logger.info("Database update completed successfully.")
            self.bundled_database = self._get_latest_database()
            if self.bundled_database:
                self.metadata['database_version'] = self._get_database_version()
                self.logger.info(f"New database version: {self.metadata['database_version']}")
                return True
            else:
                self.logger.error("Database update succeeded but no database folder found.")
                return False
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Database update failed: {e}")
            self.logger.error(f"STDERR: {e.stderr}")
            return False

    def check_amrfinder_installed(self) -> bool:
        """Check if bundled AMRfinderPlus is available and database exists"""
        try:
            if not os.path.exists(self.bundled_amrfinder):
                self.logger.error(f"Bundled AMRfinderPlus not found at: {self.bundled_amrfinder}")
                return False
            if not os.access(self.bundled_amrfinder, os.X_OK):
                self.logger.warning(f"Bundled AMRfinderPlus not executable, fixing permissions...")
                os.chmod(self.bundled_amrfinder, 0o755)
            result = subprocess.run(
                [self.bundled_amrfinder, '--version'],
                capture_output=True,
                text=True,
                check=True
            )
            version_line = result.stdout.strip()
            self.logger.info(f"Bundled AMRfinderPlus version: {version_line}")
            if self.bundled_database and os.path.exists(self.bundled_database):
                self.logger.info(f"✅ Bundled database found: {self.bundled_database}")
                db_version = self._get_database_version()
                self.logger.info(f"✅ Database version: {db_version}")
                return True
            else:
                self.logger.warning(f"⚠️ Bundled database not found at expected location.")
                self.logger.info("Please run with --update-db to download the latest database.")
                return False
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            self.logger.error(f"Bundled AMRfinderPlus check failed: {e}")
            return False

    def run_amrfinder_single_genome(self, genome_file: str, output_dir: str,
                                    min_identity: float = None, min_coverage: float = None,
                                    report_mutations: bool = True) -> Dict[str, Any]:
        """Run AMRfinderPlus on a single P. aeruginosa genome - USING BUNDLED BINARY with DYNAMIC DATABASE"""
        genome_name = Path(genome_file).stem
        output_file = os.path.join(output_dir, f"{genome_name}_amrfinder.txt")

        if not os.path.exists(self.bundled_amrfinder):
            self.logger.error(f"Bundled AMRfinderPlus not found at: {self.bundled_amrfinder}")
            return {
                'genome': genome_name,
                'output_file': output_file,
                'hits': [],
                'hit_count': 0,
                'status': 'failed',
                'error': 'AMRfinder binary not found'
            }

        cmd = [
            self.bundled_amrfinder,
            '-n', genome_file,
            '-O', 'Pseudomonas_aeruginosa',
            '--output', output_file,
            '--plus'
        ]

        if self.bundled_database and os.path.exists(self.bundled_database):
            cmd.extend(['--database', self.bundled_database])
            self.logger.info(f"Using bundled database: {self.bundled_database}")
        else:
            self.logger.warning("Using default AMRfinderPlus database location")

        if min_identity is not None:
            cmd.extend(['--ident_min', str(min_identity)])
            self.logger.info(f"Using minimum identity: {min_identity}")
        if min_coverage is not None:
            cmd.extend(['--coverage_min', str(min_coverage)])
            self.logger.info(f"Using minimum coverage: {min_coverage}")

        mut_file = None
        if report_mutations:
            mut_file = os.path.join(output_dir, f"{genome_name}_mutations.tsv")
            cmd.extend(['--mutation_all', mut_file])
            self.logger.info(f"Will report mutations to {mut_file}")

        self.logger.info("🚀 MAXIMUM SPEED: Running AMRfinderPlus on %s (using ALL %d CORES)", genome_name, self.cpus)

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            hits = self._parse_amrfinder_output(output_file)
            self._create_amrfinder_html_report(genome_name, hits, output_dir)
            self._create_amrfinder_json_report(genome_name, hits, output_dir)

            if report_mutations and mut_file and os.path.exists(mut_file):
                self._create_mutation_html_report(genome_name, mut_file, output_dir)

            return {
                'genome': genome_name,
                'output_file': output_file,
                'hits': hits,
                'hit_count': len(hits),
                'mutations_file': mut_file,
                'status': 'success'
            }
        except subprocess.CalledProcessError as e:
            self.logger.error("AMRfinderPlus failed for %s: %s", genome_name, e.stderr)
            return {
                'genome': genome_name,
                'output_file': output_file,
                'hits': [],
                'hit_count': 0,
                'status': 'failed'
            }
        except Exception as e:
            self.logger.error(f"Unexpected error for {genome_name}: {e}", exc_info=True)
            return {
                'genome': genome_name,
                'output_file': output_file,
                'hits': [],
                'hit_count': 0,
                'status': 'failed',
                'error': str(e)
            }

    def _parse_amrfinder_output(self, amrfinder_file: str) -> List[Dict]:
        """Parse AMRfinderPlus output file into structured data"""
        hits = []
        try:
            with open(amrfinder_file, 'r') as f:
                lines = f.readlines()
            if not lines or len(lines) < 2:
                return hits
            headers = lines[0].strip().split('\t')
            for line_num, line in enumerate(lines[1:], 2):
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= len(headers):
                    hit = {}
                    for i, header in enumerate(headers):
                        if i < len(parts):
                            hit[header] = parts[i]
                        else:
                            hit[header] = ''

                    processed_hit = {
                        'Protein id': hit.get('Protein id', ''),
                        'Contig id': hit.get('Contig id', ''),
                        'Start': hit.get('Start', ''),
                        'Stop': hit.get('Stop', ''),
                        'Strand': hit.get('Strand', ''),
                        'Element symbol': hit.get('Element symbol', ''),
                        'Element name': hit.get('Element name', ''),
                        'Scope': hit.get('Scope', ''),
                        'Type': hit.get('Type', ''),
                        'Subtype': hit.get('Subtype', ''),
                        'Class': hit.get('Class', ''),
                        'Subclass': hit.get('Subclass', ''),
                        'Method': hit.get('Method', ''),
                        'Target length': hit.get('Target length', ''),
                        'Reference sequence length': hit.get('Reference sequence length', ''),
                        '% Coverage of reference': hit.get('% Coverage of reference', ''),
                        '% Identity to reference': hit.get('% Identity to reference', ''),
                        'Alignment length': hit.get('Alignment length', ''),
                        'Closest reference accession': hit.get('Closest reference accession', ''),
                        'Closest reference name': hit.get('Closest reference name', ''),
                        'HMM accession': hit.get('HMM accession', ''),
                        'HMM description': hit.get('HMM description', ''),
                        'protein_id': hit.get('Protein id', ''),
                        'contig_id': hit.get('Contig id', ''),
                        'start': hit.get('Start', ''),
                        'stop': hit.get('Stop', ''),
                        'strand': hit.get('Strand', ''),
                        'gene_symbol': hit.get('Element symbol', ''),
                        'sequence_name': hit.get('Element name', ''),
                        'scope': hit.get('Scope', ''),
                        'element_type': hit.get('Type', ''),
                        'element_subtype': hit.get('Subtype', ''),
                        'class': hit.get('Class', ''),
                        'subclass': hit.get('Subclass', ''),
                        'method': hit.get('Method', ''),
                        'target_length': hit.get('Target length', ''),
                        'ref_length': hit.get('Reference sequence length', ''),
                        'coverage': hit.get('% Coverage of reference', '').replace('%', ''),
                        'identity': hit.get('% Identity to reference', '').replace('%', ''),
                        'alignment_length': hit.get('Alignment length', ''),
                        'accession': hit.get('Closest reference accession', ''),
                        'closest_name': hit.get('Closest reference name', ''),
                        'hmm_id': hit.get('HMM accession', ''),
                        'hmm_description': hit.get('HMM description', ''),
                        '_original_headers': headers,
                        '_original_values': parts
                    }
                    hits.append(processed_hit)
                else:
                    self.logger.warning("Line %d has %d parts, expected %d: %s",
                                      line_num, len(parts), len(headers), line[:100] + "...")
        except Exception as e:
            self.logger.error("Error parsing %s: %s", amrfinder_file, e)
        self.logger.info("Parsed %d records from %s", len(hits), amrfinder_file)
        return hits

    def _parse_mutations_file(self, mut_file: str) -> List[Dict]:
        """Parse mutation TSV file - same as AMR parsing."""
        return self._parse_amrfinder_output(mut_file)

    def _create_mutation_html_report(self, genome_name: str, mutations_file: str, output_dir: str):
        """Create per-genome mutation HTML report with PseudoScope red/white styling."""
        mutations = self._parse_mutations_file(mutations_file)
        if not mutations:
            self.logger.info(f"No mutations found for {genome_name}, skipping mutation HTML.")
            return

        if not isinstance(mutations, list):
            self.logger.warning(f"mutations is not a list for {genome_name}: {type(mutations)}. Skipping.")
            return

        valid_mutations = []
        for idx, m in enumerate(mutations):
            if isinstance(m, dict):
                valid_mutations.append(m)
            else:
                self.logger.warning(f"Skipping non-dict mutation entry at index {idx} for {genome_name}: {repr(m)}")
        mutations = valid_mutations

        if not mutations:
            self.logger.info(f"No valid mutation records for {genome_name}, skipping mutation HTML.")
            return

        random_quote = random.choice(self.science_quotes)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PseudoScope - Mutation Report: {genome_name}</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #2c0b0e 0%, #4a1c20 50%, #6b2b2f 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{
            background: rgba(0, 0, 0, 0.6);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid #dc3545;
            overflow-x: auto;
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #ff4757;
            text-shadow: 0 0 20px rgba(255, 71, 87, 0.5);
            overflow-x: auto;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            min-height: 100px;
            display: flex;
            flex-direction: column;
            justify-content: center;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
            border: 1px solid rgba(255, 255, 255, 0.2);
            transition: opacity 0.5s ease-in-out;
        }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; color: #ffffff; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
        }}
        .report-section h2 {{
            color: #1e3a8a;
            border-bottom: 3px solid #dc3545;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
        }}
        .table-responsive {{ width: 100%; overflow-x: auto; margin: 20px 0; }}
        .mutation-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 13px;
            min-width: 1200px;
        }}
        .mutation-table th {{
            background: linear-gradient(135deg, #dc3545 0%, #a71d2a 100%);
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
        }}
        .mutation-table td {{
            padding: 10px;
            border-bottom: 1px solid #e5e7eb;
        }}
        .mutation-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        .timestamp {{ color: #fbbf24; font-weight: bold; }}
        .authorship {{ margin-top: 15px; padding: 15px; background: rgba(255, 255, 255, 0.1); border-radius: 8px; font-size: 12px; }}
    </style>
</head>
<body>
<div class="container">
    <div class="header">
        <div class="ascii-container">
            <div class="ascii-art">{self.ascii_art}</div>
        </div>
        <div class="quote-container" id="quoteContainer">
            <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
            <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
        </div>
    </div>

    <div class="report-section">
        <h2>🧬 Point Mutation Report: {genome_name}</h2>
        <p>All point mutations detected by AMRfinderPlus (including synonymous variants).</p>
        <div class="table-responsive">
            <table class="mutation-table">
                <thead>
                    <tr>
                        <th>Gene Symbol</th><th>Mutation</th><th>Class</th><th>Subclass</th>
                        <th>Contig</th><th>Start</th><th>Stop</th><th>Strand</th>
                        <th>Coverage (%)</th><th>Identity (%)</th><th>Accession</th>
                    </tr>
                </thead>
                <tbody>
"""
        for m in mutations:
            if not isinstance(m, dict):
                continue
            html += f"""
                    <tr>
                        <td>{m.get('gene_symbol', '')}</td>
                        <td>{m.get('element_name', '')}</td>
                        <td>{m.get('class', '')}</td>
                        <td>{m.get('subclass', '')}</td>
                        <td>{m.get('contig_id', '')}</td>
                        <td>{m.get('start', '')}</td>
                        <td>{m.get('stop', '')}</td>
                        <td>{m.get('strand', '')}</td>
                        <td>{m.get('coverage', '')}</td>
                        <td>{m.get('identity', '')}</td>
                        <td>{m.get('accession', '')}</td>
                    </tr>
"""
        html += f"""
                </tbody>
            </table>
        </div>
    </div>

    <div class="footer">
        <p><strong>PseudoScope</strong> - Mutation Analysis Module</p>
        <p class="timestamp">Generated: {current_time}</p>
        <div class="authorship">
            <p><strong>Technical Support & Inquiries:</strong></p>
            <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
            <p>Email: brownbeckley94@gmail.com</p>
            <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
        </div>
    </div>
</div>

<script>
    const quotes = {json.dumps(self.science_quotes)};
    const quoteContainer = document.getElementById('quoteContainer');
    const quoteText = document.getElementById('quoteText');
    const quoteAuthor = document.getElementById('quoteAuthor');
    function getRandomQuote() {{ return quotes[Math.floor(Math.random() * quotes.length)]; }}
    function displayQuote() {{
        quoteContainer.style.opacity = '0';
        setTimeout(() => {{
            const quote = getRandomQuote();
            quoteText.textContent = '"' + quote.text + '"';
            quoteAuthor.textContent = '— ' + quote.author;
            quoteContainer.style.opacity = '1';
        }}, 500);
    }}
    setInterval(displayQuote, 10000);
</script>
</body>
</html>"""
        out_file = os.path.join(output_dir, f"{genome_name}_mutations.html")
        with open(out_file, 'w') as f:
            f.write(html)
        self.logger.info(f"✓ Mutation HTML report: {out_file}")

    def create_mutation_summary(self, all_results: Dict[str, Any], output_base: str):
        """Create mutation batch summaries (TSV, HTML, JSON) across all genomes."""
        self.logger.info("Creating mutation batch summaries...")
        all_mutations = []
        genome_mutation_counts = {}
        for genome_name, result in all_results.items():
            if not isinstance(result, dict):
                self.logger.warning(f"Skipping {genome_name}: result is not a dict.")
                genome_mutation_counts[genome_name] = 0
                continue
            if 'mutations_file' in result and result['mutations_file'] and os.path.exists(result['mutations_file']):
                muts = self._parse_mutations_file(result['mutations_file'])
                valid_muts = []
                if isinstance(muts, list):
                    for m in muts:
                        if isinstance(m, dict):
                            valid_muts.append(m)
                        else:
                            self.logger.warning(f"Skipping non-dict mutation entry for {genome_name}: {repr(m)}")
                else:
                    self.logger.warning(f"mutations_file returned non-list for {genome_name}: {type(muts)}")
                if valid_muts:
                    genome_mutation_counts[genome_name] = len(valid_muts)
                    for m in valid_muts:
                        m_copy = m.copy()
                        m_copy['genome'] = genome_name
                        all_mutations.append(m_copy)
                else:
                    genome_mutation_counts[genome_name] = 0
            else:
                genome_mutation_counts[genome_name] = 0

        if not all_mutations:
            self.logger.info("No mutations found in any genome; skipping mutation summaries.")
            return

        tsv_file = os.path.join(output_base, "mutation_summary.tsv")
        with open(tsv_file, 'w') as f:
            fieldnames = ['genome', 'gene_symbol', 'element_name', 'class', 'subclass',
                          'contig_id', 'start', 'stop', 'strand', 'coverage', 'identity', 'accession']
            f.write('\t'.join(fieldnames) + '\n')
            for m in all_mutations:
                if not isinstance(m, dict):
                    continue
                row = [m.get('genome', ''),
                       m.get('gene_symbol', ''),
                       m.get('element_name', ''),
                       m.get('class', ''),
                       m.get('subclass', ''),
                       m.get('contig_id', ''),
                       m.get('start', ''),
                       m.get('stop', ''),
                       m.get('strand', ''),
                       m.get('coverage', ''),
                       m.get('identity', ''),
                       m.get('accession', '')]
                f.write('\t'.join(str(x) for x in row) + '\n')
        self.logger.info(f"✓ Mutation TSV: {tsv_file}")

        self._create_mutation_summary_html(all_mutations, genome_mutation_counts, output_base)
        self._create_mutation_json_summaries(all_mutations, genome_mutation_counts, output_base)

    def _create_mutation_summary_html(self, all_mutations: List[Dict], genome_counts: Dict[str, int], output_base: str):
        """Create mutation HTML summary across all genomes (PseudoScope styling)."""
        gene_freq = {}
        for m in all_mutations:
            if not isinstance(m, dict):
                continue
            gene = m.get('gene_symbol', 'unknown')
            mutation = m.get('element_name', '')
            key = f"{gene}_{mutation}" if mutation else gene
            if key not in gene_freq:
                gene_freq[key] = {'count': 0, 'genomes': set(), 'gene': gene, 'mutation': mutation,
                                  'class': m.get('class',''), 'subclass': m.get('subclass','')}
            gene_freq[key]['count'] += 1
            gene_freq[key]['genomes'].add(m.get('genome',''))
        for k in gene_freq:
            gene_freq[k]['genomes'] = ', '.join(sorted(gene_freq[k]['genomes']))
        sorted_freq = sorted(gene_freq.values(), key=lambda x: x['count'], reverse=True)

        random_quote = random.choice(self.science_quotes)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PseudoScope - Mutation Batch Summary</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #2c0b0e 0%, #4a1c20 50%, #6b2b2f 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{
            background: rgba(0, 0, 0, 0.6);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid #dc3545;
            overflow-x: auto;
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #ff4757;
            text-shadow: 0 0 20px rgba(255, 71, 87, 0.5);
            overflow-x: auto;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            min-height: 100px;
            display: flex;
            flex-direction: column;
            justify-content: center;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
            border: 1px solid rgba(255, 255, 255, 0.2);
            transition: opacity 0.5s ease-in-out;
        }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; color: #ffffff; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
        }}
        .report-section h2 {{
            color: #1e3a8a;
            border-bottom: 3px solid #dc3545;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
        }}
        .summary-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            font-size: 14px;
        }}
        .summary-table th {{
            background: linear-gradient(135deg, #dc3545 0%, #a71d2a 100%);
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
        }}
        .summary-table td {{
            padding: 12px;
            border-bottom: 1px solid #e5e7eb;
        }}
        .summary-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .summary-table tr:hover {{ background-color: #e0f2fe; }}
        .table-responsive {{ width: 100%; overflow-x: auto; margin: 20px 0; }}
        .sequence-cell {{ white-space: normal !important; word-wrap: break-word; max-width: 400px; min-width: 200px; }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        .timestamp {{ color: #fbbf24; font-weight: bold; }}
        .authorship {{ margin-top: 15px; padding: 15px; background: rgba(255, 255, 255, 0.1); border-radius: 8px; font-size: 12px; }}
    </style>
</head>
<body>
<div class="container">
    <div class="header">
        <div class="ascii-container">
            <div class="ascii-art">{self.ascii_art}</div>
        </div>
        <div class="quote-container" id="quoteContainer">
            <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
            <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
        </div>
    </div>

    <div class="report-section">
        <h2>🧬 Mutation Summary Across All Genomes</h2>
        <p>Total genomes with mutations: {len([c for c in genome_counts.values() if c > 0])} / {len(genome_counts)}<br>
        Total mutation events: {len(all_mutations)}</p>
    </div>

    <div class="report-section">
        <h2>📊 Mutation Frequency by Gene/Mutation</h2>
        <div class="table-responsive">
            <table class="summary-table">
                <thead>
                    <tr><th>Gene</th><th>Mutation</th><th>Count</th><th>Genomes</th><th>Class</th><th>Subclass</th></tr>
                </thead>
                <tbody>
"""
        for item in sorted_freq:
            html += f"""
                    <tr>
                        <td><strong>{item['gene']}</strong></td>
                        <td>{item['mutation']}</td>
                        <td>{item['count']}</td>
                        <td class="sequence-cell">{item['genomes']}</td>
                        <td>{item['class']}</td>
                        <td>{item['subclass']}</td>
                    </tr>
"""
        html += """
                </tbody>
            </table>
        </div>
    </div>

    <div class="footer">
        <p><strong>PseudoScope</strong> - Mutation Batch Summary Module</p>
        <p class="timestamp">Generated: """ + current_time + """</p>
        <div class="authorship">
            <p><strong>Technical Support & Inquiries:</strong></p>
            <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
            <p>Email: brownbeckley94@gmail.com</p>
            <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
        </div>
    </div>
</div>

<script>
    const quotes = """ + json.dumps(self.science_quotes) + """;
    const quoteContainer = document.getElementById('quoteContainer');
    const quoteText = document.getElementById('quoteText');
    const quoteAuthor = document.getElementById('quoteAuthor');
    function getRandomQuote() { return quotes[Math.floor(Math.random() * quotes.length)]; }
    function displayQuote() {
        quoteContainer.style.opacity = '0';
        setTimeout(() => {
            const quote = getRandomQuote();
            quoteText.textContent = '"' + quote.text + '"';
            quoteAuthor.textContent = '— ' + quote.author;
            quoteContainer.style.opacity = '1';
        }, 500);
    }
    setInterval(displayQuote, 10000);
</script>
</body>
</html>"""
        out_file = os.path.join(output_base, "mutation_summary.html")
        with open(out_file, 'w') as f:
            f.write(html)
        self.logger.info(f"✓ Mutation HTML summary: {out_file}")

    def _create_mutation_json_summaries(self, all_mutations: List[Dict], genome_counts: Dict[str, int], output_base: str):
        """Create mutation JSON summaries."""
        genome_summary = {genome: {'total_mutations': count} for genome, count in genome_counts.items()}
        gene_mutation_map = defaultdict(lambda: {'count': 0, 'genomes': set(), 'details': []})
        for m in all_mutations:
            if not isinstance(m, dict):
                continue
            gene = m.get('gene_symbol', 'unknown')
            mut_name = m.get('element_name', '')
            key = f"{gene}_{mut_name}"
            gene_mutation_map[key]['count'] += 1
            gene_mutation_map[key]['genomes'].add(m.get('genome',''))
            gene_mutation_map[key]['details'].append({
                'genome': m.get('genome'),
                'gene': gene,
                'mutation': mut_name,
                'class': m.get('class'),
                'subclass': m.get('subclass'),
                'contig': m.get('contig_id'),
                'start': m.get('start'),
                'stop': m.get('stop')
            })
        for v in gene_mutation_map.values():
            v['genomes'] = list(v['genomes'])
        master_json = {
            'metadata': {
                'tool': 'PseudoScope AMRfinderPlus Mutation Module',
                'version': self.metadata['version'],
                'analysis_date': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                'total_genomes_analyzed': len(genome_counts),
                'total_mutations_detected': len(all_mutations)
            },
            'genome_summary': genome_summary,
            'mutation_frequency': {k: {'count': v['count'], 'genomes': v['genomes']} for k, v in gene_mutation_map.items()},
            'all_mutations': all_mutations
        }
        json_file = os.path.join(output_base, "mutation_master_summary.json")
        with open(json_file, 'w') as f:
            json.dump(master_json, f, indent=2, default=str)
        self.logger.info(f"✓ Mutation master JSON: {json_file}")

    def _create_amrfinder_html_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        """Create comprehensive HTML report for AMRfinderPlus results with PseudoScope red/white styling"""
        analysis = self._analyze_pseudo_amr_results(hits)

        interactive_js = f"""
        <script>
            function searchTable(tableId, searchTerm) {{
                const table = document.getElementById(tableId);
                const rows = table.getElementsByTagName('tr');
                let visibleCount = 0;
                for (let i = 1; i < rows.length; i++) {{
                    const row = rows[i];
                    const text = row.textContent.toLowerCase();
                    if (text.includes(searchTerm.toLowerCase())) {{
                        row.style.display = '';
                        visibleCount++;
                    }} else {{
                        row.style.display = 'none';
                    }}
                }}
                const resultCounter = document.getElementById('result-counter-' + tableId);
                if (resultCounter) {{
                    resultCounter.textContent = visibleCount + ' results found';
                }}
            }}
            function exportToCSV(tableId, filename) {{
                const table = document.getElementById(tableId);
                const rows = table.getElementsByTagName('tr');
                let csv = [];
                const headerCells = rows[0].getElementsByTagName('th');
                const headerRow = [];
                for (let cell of headerCells) {{
                    headerRow.push(cell.textContent);
                }}
                csv.push(headerRow.join(','));
                for (let i = 1; i < rows.length; i++) {{
                    if (rows[i].style.display !== 'none') {{
                        const cells = rows[i].getElementsByTagName('td');
                        const row = [];
                        for (let cell of cells) {{
                            let text = cell.textContent.trim();
                            text = text.replace(/\\n/g, ' ').replace(/,/g, ';');
                            row.push(text);
                        }}
                        csv.push(row.join(','));
                    }}
                }}
                const blob = new Blob([csv.join('\\n')], {{ type: 'text/csv' }});
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = filename;
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);
            }}
            function exportToJSON(dataVar, filename) {{
                const data = window[dataVar];
                const jsonStr = JSON.stringify(data, null, 2);
                const blob = new Blob([jsonStr], {{ type: 'application/json' }});
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = filename;
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);
            }}
            function printReport() {{
                const printWindow = window.open('', '_blank');
                printWindow.document.write('<html><head><title>{genome_name} - AMR Report</title>');
                printWindow.document.write('<style>');
                printWindow.document.write(`
                    body {{ font-family: Arial, sans-serif; margin: 20px; }}
                    .no-print {{ display: none !important; }}
                    table {{ border-collapse: collapse; width: 100%; margin: 10px 0; }}
                    th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                    th {{ background-color: #f2f2f2; }}
                    .critical-row {{ background-color: #f8d7da; }}
                    .high-risk-row {{ background-color: #fff3cd; }}
                    .header {{ text-align: center; margin-bottom: 20px; }}
                `);
                printWindow.document.write('</style>');
                printWindow.document.write('</head><body>');
                const content = document.querySelector('.container').cloneNode(true);
                const noPrintElements = content.querySelectorAll('.no-print');
                noPrintElements.forEach(el => el.remove());
                printWindow.document.write(content.innerHTML);
                printWindow.document.write('</body></html>');
                printWindow.document.close();
                printWindow.print();
            }}
            function toggleAllMechanisms() {{
                const container = document.getElementById('all-mechanisms-container');
                const button = document.getElementById('toggle-mechanisms-btn');
                if (container.style.display === 'none') {{
                    container.style.display = 'block';
                    button.textContent = 'Show Less';
                }} else {{
                    container.style.display = 'none';
                    button.textContent = 'Show All Mechanisms';
                }}
            }}
            function quickSearch() {{
                const searchTerm = document.getElementById('quick-search').value.toLowerCase();
                const sections = document.querySelectorAll('.card');
                sections.forEach(section => {{
                    const text = section.textContent.toLowerCase();
                    if (text.includes(searchTerm)) {{
                        section.style.border = '2px solid #dc3545';
                        section.style.backgroundColor = 'rgba(220, 53, 69, 0.1)';
                        section.scrollIntoView({{ behavior: 'smooth', block: 'nearest' }});
                    }} else {{
                        section.style.border = '';
                        section.style.backgroundColor = '';
                    }}
                }});
            }}
            function clearSearch() {{
                document.getElementById('quick-search').value = '';
                const sections = document.querySelectorAll('.card');
                sections.forEach(section => {{
                    section.style.border = '';
                    section.style.backgroundColor = '';
                }});
            }}
            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            function rotateQuote() {{
                document.getElementById('science-quote').innerHTML = quotes[currentQuote].text;
                currentQuote = (currentQuote + 1) % quotes.length;
            }}
            document.addEventListener('DOMContentLoaded', function() {{
                rotateQuote();
                setInterval(rotateQuote, 10000);
            }});
        </script>
        """

        export_data_js = f"""
        <script>
            window.reportData = {{
                metadata: {{
                    genome: '{genome_name}',
                    date: '{self.metadata['analysis_date']}',
                    tool: '{self.metadata['tool_name']}',
                    version: '{self.metadata['version']}'
                }},
                summary: {json.dumps(analysis, indent=2)},
                hits: {json.dumps(hits, indent=2)}
            }};
        </script>
        """

        random_quote = random.choice(self.science_quotes)

        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>PseudoScope AMRfinderPlus Analysis Report</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            background: linear-gradient(135deg, #2c0b0e 0%, #4a1c20 50%, #6b2b2f 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{
            background: rgba(0, 0, 0, 0.6);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid #dc3545;
            overflow-x: auto;
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #ff4757;
            text-shadow: 0 0 20px rgba(255, 71, 87, 0.5);
            overflow-x: auto;
        }}
        .card {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            margin: 20px 0;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.2);
        }}
        .gene-table, .class-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .gene-table th, .gene-table td, .class-table th, .class-table td {{
            padding: 15px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
        }}
        .gene-table th, .class-table th {{
            background: linear-gradient(135deg, #dc3545 0%, #a71d2a 100%);
            color: white;
            font-weight: 600;
        }}
        tr:hover {{ background-color: #f8f9fa; }}
        .success {{ color: #28a745; font-weight: 600; }}
        .warning {{ color: #ffc107; font-weight: 600; }}
        .error {{ color: #dc3545; font-weight: 600; }}
        .summary-stats {{
            display: flex;
            justify-content: space-around;
            margin: 20px 0;
            flex-wrap: wrap;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #dc3545 0%, #a71d2a 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        .critical-stat-card {{
            background: linear-gradient(135deg, #8b0000 0%, #6a0000 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            color: white;
            padding: 20px;
            border-radius: 12px;
            margin: 20px 0;
            text-align: center;
            font-style: italic;
            border-left: 4px solid #dc3545;
        }}
        .footer {{
            background: rgba(0, 0, 0, 0.8);
            color: white;
            padding: 30px;
            border-radius: 12px;
            margin-top: 40px;
        }}
        .footer a {{ color: #dc3545; text-decoration: none; }}
        .footer a:hover {{ text-decoration: underline; }}
        .resistance-badge {{
            display: inline-block;
            background: #dc3545;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .critical-risk-badge {{
            display: inline-block;
            background: #8b0000;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
            font-weight: bold;
        }}
        .warning-badge {{
            display: inline-block;
            background: #ffc107;
            color: black;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .success-badge {{
            display: inline-block;
            background: #28a745;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .present {{ background-color: #d4edda; }}
        .critical-row {{ background-color: #f8d7da; font-weight: bold; border-left: 4px solid #dc3545; }}
        .high-risk-row {{ background-color: #fff3cd; border-left: 4px solid #ffc107; }}
        .interactive-controls {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            align-items: center;
        }}
        .search-box {{
            flex: 1;
            min-width: 200px;
        }}
        .search-box input {{
            width: 100%;
            padding: 8px 12px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
        }}
        .export-buttons {{
            display: flex;
            gap: 8px;
            flex-wrap: wrap;
        }}
        .export-buttons button {{
            padding: 8px 16px;
            background: #dc3545;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
            font-size: 14px;
            transition: background 0.3s;
        }}
        .export-buttons button:hover {{
            background: #c82333;
        }}
        .export-buttons button.print {{
            background: #10b981;
        }}
        .export-buttons button.print:hover {{
            background: #059669;
        }}
        .result-counter {{
            font-size: 0.9em;
            color: #666;
            font-style: italic;
        }}
        .sequence-name {{
            max-width: 300px;
            overflow-wrap: break-word;
            word-wrap: break-word;
            word-break: break-all;
        }}
        .gene-list-container {{
            max-height: 300px;
            overflow-y: auto;
            padding: 10px;
            background: #f8f9fa;
            border-radius: 5px;
            border: 1px solid #e0e0e0;
        }}
        .gene-item {{
            display: inline-block;
            margin: 2px;
            padding: 3px 8px;
            background: #e9ecef;
            border-radius: 3px;
            font-size: 0.85em;
        }}
        .critical-gene {{
            background: #f8d7da;
            color: #721c24;
            border: 1px solid #f5c6cb;
            font-weight: bold;
        }}
        .high-risk-gene {{
            background: #fff3cd;
            color: #856404;
            border: 1px solid #ffeaa7;
        }}
        .quick-search-bar {{
            background: rgba(255, 255, 255, 0.95);
            padding: 15px;
            border-radius: 8px;
            margin: 10px 0;
            display: flex;
            gap: 10px;
        }}
        .quick-search-bar input {{
            flex: 1;
            padding: 8px 12px;
            border: 2px solid #dc3545;
            border-radius: 4px;
            font-size: 14px;
        }}
        .quick-search-bar button {{
            padding: 8px 20px;
            background: #dc3545;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
        }}
        @media print {{
            .no-print, .interactive-controls, .quick-search-bar {{ display: none !important; }}
            body {{ background: white !important; color: black !important; }}
            .card {{ background: white !important; color: black !important; box-shadow: none !important; }}
            .gene-table, .class-table {{ box-shadow: none !important; }}
        }}
    </style>
    {interactive_js}
    {export_data_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.ascii_art}</div>
            </div>
            <div class="card">
                <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 PseudoScope AMRfinderPlus Analysis Report</h1>
                <p style="color: #666; font-size: 1.2em;">Comprehensive <em>Pseudomonas aeruginosa</em> Antimicrobial Resistance Analysis</p>
                <div class="quick-search-bar no-print">
                    <input type="text" id="quick-search" placeholder="🔍 Quick search across entire report...">
                    <button onclick="quickSearch()">Search</button>
                    <button onclick="clearSearch()" style="background: #6c757d;">Clear</button>
                </div>
                <div class="export-buttons no-print" style="margin-top: 15px;">
                    <button onclick="exportToJSON('reportData', '{genome_name}_amr_report.json')">📥 Export JSON</button>
                    <button onclick="printReport()" class="print">🖨️ Print Report</button>
                </div>
            </div>
        </div>
        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
"""
        if analysis['critical_risk_genes'] > 0:
            html_content += f"""
        <div class="card" style="border-left: 4px solid #dc3545; background: #f8d7da;">
            <h2 style="color: #dc3545;">🚨 CRITICAL RISK AMR GENES DETECTED</h2>
            <p><strong>{analysis['critical_risk_genes']} CRITICAL RISK antimicrobial resistance genes found:</strong></p>
            <div style="margin: 10px 0;">
                <p style="color: #721c24; font-weight: bold;">
                    ⚠️ These genes confer resistance to last-resort antibiotics and represent
                    a serious public health concern requiring immediate attention.
                </p>
"""
            for gene in analysis['critical_risk_list']:
                html_content += f'<span class="critical-risk-badge" style="font-size: 1.1em;">🚨 {gene}</span>'
            html_content += """
            </div>
        </div>
"""
        html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #dc3545; padding-bottom: 10px;">📊 P. aeruginosa AMR Summary</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total AMR Genes</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['total_genes']}</p>
                </div>
                <div class="stat-card">
                    <h3>High Risk Genes</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['high_risk_genes']}</p>
                </div>
                <div class="critical-stat-card">
                    <h3>Critical Risk</h3>
                    <p style="font-size: 2em; margin: 0;">{analysis['critical_risk_genes']}</p>
                </div>
            </div>
            <p><strong>Genome:</strong> {genome_name}</p>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>AMRfinderPlus Version:</strong> {self.metadata['amrfinder_version']}</p>
            <p><strong>Database Version:</strong> {self.metadata['database_version']}</p>
        </div>
"""
        if analysis['high_risk_genes'] > 0 and analysis['critical_risk_genes'] == 0:
            html_content += f"""
        <div class="card" style="border-left: 4px solid #ffc107;">
            <h2 style="color: #856404;">⚠️ High-Risk AMR Genes Detected</h2>
            <p><strong>{analysis['high_risk_genes']} high-risk antimicrobial resistance genes found:</strong></p>
            <div style="margin: 10px 0;">
"""
            for gene in analysis['high_risk_list']:
                html_content += f'<span class="resistance-badge">{gene}</span>'
            html_content += """
            </div>
        </div>
"""
        if any(analysis['resistance_mechanisms'].values()):
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #dc3545; padding-bottom: 10px;">🔬 Resistance Mechanism Breakdown</h2>
"""
            mechanisms = analysis['resistance_mechanisms']
            all_mechanisms_html = ""
            if mechanisms['carbapenemase']:
                all_mechanisms_html += f"""
            <div style="margin: 10px 0; padding: 10px; background: #f8d7da; border-radius: 5px;">
                <strong>Carbapenemase Genes (CRITICAL):</strong> {', '.join(mechanisms['carbapenemase'])}
            </div>
"""
            if mechanisms['esbl']:
                all_mechanisms_html += f"""
            <div style="margin: 10px 0; padding: 10px; background: #fff3cd; border-radius: 5px;">
                <strong>ESBL Genes:</strong> {', '.join(mechanisms['esbl'])}
            </div>
"""
            if mechanisms['colistin_resistance']:
                all_mechanisms_html += f"""
            <div style="margin: 10px 0; padding: 10px; background: #f8d7da; border-radius: 5px;">
                <strong>Colistin Resistance (CRITICAL):</strong> {', '.join(mechanisms['colistin_resistance'])}
            </div>
"""
            if mechanisms['aminoglycoside_resistance']:
                all_mechanisms_html += f"""
            <div style="margin: 10px 0; padding: 10px; background: #d1ecf1; border-radius: 5px;">
                <strong>Aminoglycoside Resistance:</strong> {', '.join(mechanisms['aminoglycoside_resistance'])}
            </div>
"""
            if mechanisms['fluoroquinolone_resistance']:
                all_mechanisms_html += f"""
            <div style="margin: 10px 0; padding: 10px; background: #d1ecf1; border-radius: 5px;">
                <strong>Fluoroquinolone Resistance:</strong> {', '.join(mechanisms['fluoroquinolone_resistance'])}
            </div>
"""
            if mechanisms['efflux_pumps']:
                all_mechanisms_html += f"""
            <div style="margin: 10px 0; padding: 10px; background: #e2e3e5; border-radius: 5px;">
                <strong>Efflux Pumps:</strong> {', '.join(mechanisms['efflux_pumps'])}
            </div>
"""
            if mechanisms['other_amr']:
                all_mechanisms_html += f"""
            <div style="margin: 10px 0; padding: 10px; background: #f8f9fa; border-radius: 5px;">
                <strong>Other AMR Genes:</strong> {', '.join(mechanisms['other_amr'])}
            </div>
"""
            html_content += all_mechanisms_html
            html_content += """
        </div>
"""
        if analysis['resistance_classes']:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #dc3545; padding-bottom: 10px;">🧪 Resistance Classes Detected</h2>
            <table class="class-table" id="resistance-classes-table">
                <thead>
                    <tr>
                        <th>Resistance Class</th>
                        <th>Gene Count</th>
                        <th>Genes</th>
                    </tr>
                </thead>
                <tbody>
"""
            for class_name, genes in analysis['resistance_classes'].items():
                gene_list = ", ".join(genes)
                html_content += f"""
                    <tr>
                        <td><strong>{class_name}</strong></td>
                        <td>{len(genes)}</td>
                        <td>{gene_list}</td>
                    </tr>
"""
            html_content += """
                </tbody>
            </table>
        </div>
"""
        if hits:
            html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #dc3545; padding-bottom: 10px;">🔬 Detailed AMR Genes Detected</h2>
            <div class="interactive-controls no-print">
                <div class="search-box">
                    <input type="text" id="search-detailed-amr"
                           placeholder="Search genes, classes, or sequence names..."
                           onkeyup="searchTable('detailed-amr-table', this.value)">
                </div>
                <div class="export-buttons">
                    <button onclick="exportToCSV('detailed-amr-table', '{genome_name}_amr_genes.csv')">📥 Export CSV</button>
                </div>
                <div class="result-counter" id="result-counter-detailed-amr-table">
                    {len(hits)} results found
                </div>
            </div>
            <table class="gene-table" id="detailed-amr-table">
                <thead>
                    <tr>
                        <th>Gene Symbol</th>
                        <th>Sequence Name</th>
                        <th>Class</th>
                        <th>Subclass</th>
                        <th>Coverage</th>
                        <th>Identity</th>
                        <th>Scope</th>
                    </tr>
                </thead>
                <tbody>
"""
            for hit in hits:
                row_class = "present"
                gene_symbol = hit.get('gene_symbol', '')
                if gene_symbol in analysis['critical_risk_list']:
                    row_class = "critical-row"
                elif gene_symbol in analysis['high_risk_list']:
                    row_class = "high-risk-row"
                html_content += f"""
                    <tr class="{row_class}">
                        <td><strong>{gene_symbol}</strong></td>
                        <td class="sequence-name" title="{hit.get('sequence_name', '')}">{hit.get('sequence_name', '')}</td>
                        <td>{hit.get('class', '')}</td>
                        <td>{hit.get('subclass', '')}</td>
                        <td>{hit.get('coverage', '')}%</td>
                        <td>{hit.get('identity', '')}%</td>
                        <td>{hit.get('scope', '')}</td>
                    </tr>
"""
            html_content += """
                </tbody>
            </table>
        </div>
"""
        else:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #dc3545; padding-bottom: 10px;">✅ No AMR Genes Detected</h2>
            <p>No antimicrobial resistance genes found in this P. aeruginosa genome.</p>
        </div>
"""
        html_content += f"""
        <div class="footer">
            <h3 style="color: #fff; border-bottom: 2px solid #dc3545; padding-bottom: 10px;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
            <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;">
            </p>
        </div>
    </div>
</body>
</html>
"""
        html_file = os.path.join(output_dir, f"{genome_name}_amrfinder_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info("P. aeruginosa AMRfinderPlus HTML report generated: %s", html_file)

    def _create_amrfinder_json_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        """Create JSON report for AMRfinderPlus results"""
        analysis = self._analyze_pseudo_amr_results(hits)
        json_data = {
            'metadata': {
                'genome': genome_name,
                'analysis_date': self.metadata['analysis_date'],
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'amrfinder_version': self.metadata['amrfinder_version'],
                'database_version': self.metadata['database_version']
            },
            'summary': {
                'total_genes': analysis['total_genes'],
                'high_risk_genes': analysis['high_risk_genes'],
                'critical_risk_genes': analysis['critical_risk_genes'],
                'high_risk_list': analysis['high_risk_list'],
                'critical_risk_list': analysis['critical_risk_list'],
                'resistance_classes': analysis['resistance_classes'],
                'resistance_mechanisms': analysis['resistance_mechanisms']
            },
            'hits': hits
        }
        json_file = os.path.join(output_dir, f"{genome_name}_amrfinder_report.json")
        with open(json_file, 'w') as f:
            json.dump(json_data, f, indent=2)
        self.logger.info("P. aeruginosa AMRfinderPlus JSON report generated: %s", json_file)

    def _analyze_pseudo_amr_results(self, hits: List[Dict]) -> Dict[str, Any]:
        """Analyze AMR results specifically for P. aeruginosa with enhanced risk assessment"""
        analysis = {
            'total_genes': len(hits),
            'resistance_classes': {},
            'total_classes': 0,
            'high_risk_genes': 0,
            'critical_risk_genes': 0,
            'high_risk_list': [],
            'critical_risk_list': [],
            'resistance_mechanisms': {
                'carbapenemase': [],
                'esbl': [],
                'colistin_resistance': [],
                'fluoroquinolone_resistance': [],
                'aminoglycoside_resistance': [],
                'efflux_pumps': [],
                'other_amr': []
            }
        }
        for hit in hits:
            gene_symbol = hit.get('gene_symbol', '')
            resistance_class = hit.get('class', '')
            self._categorize_pseudo_resistance_mechanism(gene_symbol, resistance_class, analysis)
            if gene_symbol in self.critical_risk_genes:
                analysis['critical_risk_genes'] += 1
                if gene_symbol not in analysis['critical_risk_list']:
                    analysis['critical_risk_list'].append(gene_symbol)
            if gene_symbol in self.high_risk_genes:
                analysis['high_risk_genes'] += 1
                if gene_symbol not in analysis['high_risk_list']:
                    analysis['high_risk_list'].append(gene_symbol)
            if resistance_class:
                if resistance_class not in analysis['resistance_classes']:
                    analysis['resistance_classes'][resistance_class] = []
                if gene_symbol not in analysis['resistance_classes'][resistance_class]:
                    analysis['resistance_classes'][resistance_class].append(gene_symbol)
        analysis['total_classes'] = len(analysis['resistance_classes'])
        return analysis

    def _categorize_pseudo_resistance_mechanism(self, gene_symbol: str, resistance_class: str, analysis: Dict[str, Any]):
        """Categorize P. aeruginosa genes by resistance mechanism"""
        gene_lower = gene_symbol.lower()
        if any(carba in gene_lower for carba in ['oxa', 'imp', 'vim', 'ndm', 'kpc', 'ges', 'spm', 'aim', 'dim']):
            if any(crit_gene in gene_lower for crit_gene in [g.lower() for g in self.critical_carbapenemases]):
                analysis['resistance_mechanisms']['carbapenemase'].append(gene_symbol)
                return
        if any(esbl in gene_lower for esbl in ['per', 'veb', 'bel', 'ges', 'tem', 'shv', 'ctx']):
            if any(crit_gene in gene_lower for crit_gene in [g.lower() for g in self.critical_esbls]):
                analysis['resistance_mechanisms']['esbl'].append(gene_symbol)
                return
        if any(mcr in gene_lower for mcr in ['mcr']) or any(pm in gene_lower for pm in ['pmr', 'pho', 'mgr', 'lpx', 'arn', 'ept', 'bas']):
            if any(crit_gene in gene_lower for crit_gene in [g.lower() for g in self.critical_colistin]):
                analysis['resistance_mechanisms']['colistin_resistance'].append(gene_symbol)
                return
        if any(ag in gene_lower for ag in ['arm', 'rmt', 'npm', 'aac', 'ant', 'aph', 'str']):
            if any(crit_gene in gene_lower for crit_gene in [g.lower() for g in self.critical_aminoglycoside]):
                analysis['resistance_mechanisms']['aminoglycoside_resistance'].append(gene_symbol)
                return
        if any(qnr in gene_lower for qnr in ['qnr', 'qep']):
            analysis['resistance_mechanisms']['fluoroquinolone_resistance'].append(gene_symbol)
            return
        if any(efflux in gene_lower for efflux in ['mex', 'opr', 'ade', 'abe', 'acr', 'emr', 'mdt']):
            analysis['resistance_mechanisms']['efflux_pumps'].append(gene_symbol)
            return
        analysis['resistance_mechanisms']['other_amr'].append(gene_symbol)

    def create_amr_summary(self, all_results: Dict[str, Any], output_base: str):
        """Create comprehensive AMR summary files and HTML reports for all P. aeruginosa samples"""
        self.logger.info("Creating P. aeruginosa AMR summary files and HTML reports...")
        summary_file = os.path.join(output_base, "pseudo_amrfinder_summary.tsv")
        with open(summary_file, 'w') as f:
            f.write("Genome\tProtein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference\t% Identity to reference\tAlignment length\tClosest reference accession\tClosest reference name\tHMM accession\tHMM description\n")
            for genome_name, result in all_results.items():
                if not isinstance(result, dict):
                    self.logger.warning(f"Skipping {genome_name} in summary: result is not a dict.")
                    continue
                for hit in result.get('hits', []):
                    if not isinstance(hit, dict):
                        continue
                    row = [
                        genome_name,
                        hit.get('Protein id', ''),
                        hit.get('Contig id', ''),
                        hit.get('Start', ''),
                        hit.get('Stop', ''),
                        hit.get('Strand', ''),
                        hit.get('Element symbol', ''),
                        hit.get('Element name', ''),
                        hit.get('Scope', ''),
                        hit.get('Type', ''),
                        hit.get('Subtype', ''),
                        hit.get('Class', ''),
                        hit.get('Subclass', ''),
                        hit.get('Method', ''),
                        hit.get('Target length', ''),
                        hit.get('Reference sequence length', ''),
                        hit.get('% Coverage of reference', ''),
                        hit.get('% Identity to reference', ''),
                        hit.get('Alignment length', ''),
                        hit.get('Closest reference accession', ''),
                        hit.get('Closest reference name', ''),
                        hit.get('HMM accession', ''),
                        hit.get('HMM description', '')
                    ]
                    f.write('\t'.join(str(x) for x in row) + '\n')
        self.logger.info("✓ P. aeruginosa AMR summary file created: %s", summary_file)

        stats_file = os.path.join(output_base, "pseudo_amrfinder_statistics_summary.tsv")
        with open(stats_file, 'w') as f:
            f.write("Genome\tTotal_AMR_Genes\tHigh_Risk_Genes\tCritical_Risk_Genes\tResistance_Classes\tGene_List\n")
            for genome_name, result in all_results.items():
                if not isinstance(result, dict):
                    self.logger.warning(f"Skipping {genome_name} in stats: result is not a dict.")
                    continue
                genes = list(set(hit.get('gene_symbol', '') for hit in result.get('hits', []) if isinstance(hit, dict) and hit.get('gene_symbol')))
                gene_list = ",".join(genes)
                high_risk_count = sum(1 for gene in genes if gene in self.high_risk_genes)
                critical_risk_count = sum(1 for gene in genes if gene in self.critical_risk_genes)
                classes = list(set(hit.get('class', '') for hit in result.get('hits', []) if isinstance(hit, dict) and hit.get('class')))
                class_list = ",".join(classes)
                f.write(f"{genome_name}\t{result.get('hit_count', 0)}\t{high_risk_count}\t{critical_risk_count}\t{class_list}\t{gene_list}\n")
        self.logger.info("✓ P. aeruginosa AMR statistics summary created: %s", stats_file)

        self.create_json_summaries(all_results, output_base)
        self._create_summary_html_report(all_results, output_base)

    def create_json_summaries(self, all_results: Dict[str, Any], output_base: str):
        """Create JSON summary files"""
        self.logger.info("Creating JSON summaries...")
        master_summary = {
            'metadata': {
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'amrfinder_version': self.metadata['amrfinder_version'],
                'database_version': self.metadata['database_version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_genomes': len(all_results)
            },
            'genome_summaries': {},
            'cross_genome_patterns': {}
        }
        all_hits_by_gene = defaultdict(lambda: {'count': 0, 'genomes': set()})
        genomes_with_critical = 0
        genomes_with_high_risk = 0
        for genome_name, result in all_results.items():
            if not isinstance(result, dict):
                self.logger.warning(f"Skipping {genome_name} in JSON: result is not a dict.")
                continue
            hits = result.get('hits', [])
            genes = [hit.get('gene_symbol', '') for hit in hits if isinstance(hit, dict) and hit.get('gene_symbol', '')]
            unique_genes = set(genes)
            critical_genes = [g for g in unique_genes if g in self.critical_risk_genes]
            high_risk_genes = [g for g in unique_genes if g in self.high_risk_genes and g not in self.critical_risk_genes]
            if critical_genes:
                genomes_with_critical += 1
            if high_risk_genes:
                genomes_with_high_risk += 1
            master_summary['genome_summaries'][genome_name] = {
                'total_hits': result.get('hit_count', 0),
                'unique_genes': len(unique_genes),
                'critical_genes': critical_genes,
                'high_risk_genes': high_risk_genes,
                'genes': list(unique_genes),
                'status': result.get('status', 'unknown')
            }
            for gene in unique_genes:
                all_hits_by_gene[gene]['count'] += 1
                all_hits_by_gene[gene]['genomes'].add(genome_name)
        cross_genome_data = {}
        for gene, data in all_hits_by_gene.items():
            cross_genome_data[gene] = {
                'frequency': data['count'],
                'genomes': list(data['genomes']),
                'risk_level': 'CRITICAL' if gene in self.critical_risk_genes else 'HIGH' if gene in self.high_risk_genes else 'STANDARD'
            }
        master_summary['cross_genome_patterns'] = {
            'total_unique_genes': len(all_hits_by_gene),
            'genomes_with_critical': genomes_with_critical,
            'genomes_with_high_risk': genomes_with_high_risk,
            'gene_frequency': cross_genome_data
        }
        master_json_file = os.path.join(output_base, "pseudo_amrfinder_master_summary.json")
        with open(master_json_file, 'w') as f:
            json.dump(master_summary, f, indent=2)
        self.logger.info("✓ Master JSON summary created: %s", master_json_file)
        for genome_name, result in all_results.items():
            if not isinstance(result, dict):
                continue
            genome_dir = os.path.join(output_base, genome_name)
            if os.path.exists(genome_dir):
                json_file = os.path.join(genome_dir, f"{genome_name}_amrfinder_summary.json")
                with open(json_file, 'w') as f:
                    json.dump({
                        'metadata': {
                            'genome': genome_name,
                            'analysis_date': self.metadata['analysis_date']
                        },
                        'summary': {
                            'total_hits': result.get('hit_count', 0),
                            'genes': list(set(hit.get('gene_symbol', '') for hit in result.get('hits', []) if isinstance(hit, dict) and hit.get('gene_symbol', '')))
                        },
                        'hits': result.get('hits', [])[:10000]
                    }, f, indent=2)

    def _create_summary_html_report(self, all_results: Dict[str, Any], output_base: str):
        """Create comprehensive HTML summary report with interactive features"""
        all_hits = []
        for genome_name, result in all_results.items():
            if not isinstance(result, dict):
                self.logger.warning(f"Skipping {genome_name} in summary HTML: result is not a dict.")
                continue
            for hit in result.get('hits', []):
                if not isinstance(hit, dict):
                    continue
                hit_with_genome = hit.copy()
                hit_with_genome['genome'] = genome_name
                all_hits.append(hit_with_genome)
        total_genomes = len(all_results)
        total_hits = len(all_hits)
        critical_genes_found = set()
        high_risk_genes_found = set()
        genomes_with_critical = 0
        genomes_with_high_risk = 0
        genes_per_genome = {}
        gene_frequency = {}
        for genome_name, result in all_results.items():
            if not isinstance(result, dict):
                continue
            genome_genes = set()
            for hit in result.get('hits', []):
                if not isinstance(hit, dict):
                    continue
                gene = hit.get('gene_symbol', '')
                if gene:
                    genome_genes.add(gene)
                    if gene not in gene_frequency:
                        gene_frequency[gene] = set()
                    gene_frequency[gene].add(genome_name)
            genes_per_genome[genome_name] = genome_genes
            has_critical = any(gene in genome_genes for gene in self.critical_risk_genes)
            has_high_risk = any(gene in genome_genes for gene in self.high_risk_genes)
            if has_critical:
                genomes_with_critical += 1
                critical_genes_found.update(genome_genes.intersection(self.critical_risk_genes))
            if has_high_risk:
                genomes_with_high_risk += 1
                high_risk_genes_found.update(genome_genes.intersection(self.high_risk_genes))

        interactive_js = f"""
        <script>
            function searchTable(tableId, searchTerm) {{
                const table = document.getElementById(tableId);
                const rows = table.getElementsByTagName('tr');
                let visibleCount = 0;
                for (let i = 1; i < rows.length; i++) {{
                    const row = rows[i];
                    const text = row.textContent.toLowerCase();
                    if (text.includes(searchTerm.toLowerCase())) {{
                        row.style.display = '';
                        visibleCount++;
                    }} else {{
                        row.style.display = 'none';
                    }}
                }}
                const resultCounter = document.getElementById('result-counter-' + tableId);
                if (resultCounter) {{
                    resultCounter.textContent = visibleCount + ' results found';
                }}
            }}
            function exportToCSV(tableId, filename) {{
                const table = document.getElementById(tableId);
                const rows = table.getElementsByTagName('tr');
                let csv = [];
                const headerCells = rows[0].getElementsByTagName('th');
                const headerRow = [];
                for (let cell of headerCells) {{
                    headerRow.push(cell.textContent);
                }}
                csv.push(headerRow.join(','));
                for (let i = 1; i < rows.length; i++) {{
                    if (rows[i].style.display !== 'none') {{
                        const cells = rows[i].getElementsByTagName('td');
                        const row = [];
                        for (let cell of cells) {{
                            let text = cell.textContent.trim();
                            text = text.replace(/\\n/g, ' ').replace(/,/g, ';');
                            row.push(text);
                        }}
                        csv.push(row.join(','));
                    }}
                }}
                const blob = new Blob([csv.join('\\n')], {{ type: 'text/csv' }});
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = filename;
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);
            }}
            function exportToJSON(dataVar, filename) {{
                const data = window[dataVar];
                const jsonStr = JSON.stringify(data, null, 2);
                const blob = new Blob([jsonStr], {{ type: 'application/json' }});
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = filename;
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);
            }}
            function printReport() {{
                const printWindow = window.open('', '_blank');
                printWindow.document.write('<html><head><title>Summary AMR Report</title>');
                printWindow.document.write('<style>');
                printWindow.document.write(`
                    body {{ font-family: Arial, sans-serif; margin: 20px; }}
                    .no-print {{ display: none !important; }}
                    table {{ border-collapse: collapse; width: 100%; margin: 10px 0; }}
                    th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                    th {{ background-color: #f2f2f2; }}
                    .critical-row {{ background-color: #f8d7da; }}
                    .high-risk-row {{ background-color: #fff3cd; }}
                    .header {{ text-align: center; margin-bottom: 20px; }}
                `);
                printWindow.document.write('</style>');
                printWindow.document.write('</head><body>');
                const content = document.querySelector('.container').cloneNode(true);
                const noPrintElements = content.querySelectorAll('.no-print');
                noPrintElements.forEach(el => el.remove());
                printWindow.document.write(content.innerHTML);
                printWindow.document.write('</body></html>');
                printWindow.document.close();
                printWindow.print();
            }}
            function quickSearch() {{
                const searchTerm = document.getElementById('quick-search').value.toLowerCase();
                const sections = document.querySelectorAll('.card');
                sections.forEach(section => {{
                    const text = section.textContent.toLowerCase();
                    if (text.includes(searchTerm)) {{
                        section.style.border = '2px solid #dc3545';
                        section.style.backgroundColor = 'rgba(220, 53, 69, 0.1)';
                        section.scrollIntoView({{ behavior: 'smooth', block: 'nearest' }});
                    }} else {{
                        section.style.border = '';
                        section.style.backgroundColor = '';
                    }}
                }});
            }}
            function clearSearch() {{
                document.getElementById('quick-search').value = '';
                const sections = document.querySelectorAll('.card');
                sections.forEach(section => {{
                    section.style.border = '';
                    section.style.backgroundColor = '';
                }});
            }}
            function filterByRisk(riskLevel) {{
                const table = document.getElementById('gene-frequency-table');
                const rows = table.getElementsByTagName('tr');
                for (let i = 1; i < rows.length; i++) {{
                    const row = rows[i];
                    const riskCell = row.cells[3];
                    if (riskLevel === 'all' || riskCell.textContent.includes(riskLevel.toUpperCase())) {{
                        row.style.display = '';
                    }} else {{
                        row.style.display = 'none';
                    }}
                }}
            }}
            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            function rotateQuote() {{
                document.getElementById('science-quote').innerHTML = quotes[currentQuote].text;
                currentQuote = (currentQuote + 1) % quotes.length;
            }}
            document.addEventListener('DOMContentLoaded', function() {{
                rotateQuote();
                setInterval(rotateQuote, 10000);
            }});
        </script>
        """

        export_data_js = f"""
        <script>
            window.summaryData = {{
                metadata: {{
                    total_genomes: {total_genomes},
                    total_hits: {total_hits},
                    date: '{self.metadata['analysis_date']}',
                    tool: '{self.metadata['tool_name']}',
                    version: '{self.metadata['version']}'
                }},
                critical_genes: {json.dumps(list(critical_genes_found))},
                high_risk_genes: {json.dumps(list(high_risk_genes_found))},
                genomes_with_critical: {genomes_with_critical},
                genomes_with_high_risk: {genomes_with_high_risk},
                gene_frequency: {json.dumps({gene: list(genomes) for gene, genomes in gene_frequency.items()})}
            }};
        </script>
        """

        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>PseudoScope AMRfinderPlus - Summary Report</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            background: linear-gradient(135deg, #2c0b0e 0%, #4a1c20 50%, #6b2b2f 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{
            background: rgba(0, 0, 0, 0.6);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid #dc3545;
            overflow-x: auto;
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #ff4757;
            text-shadow: 0 0 20px rgba(255, 71, 87, 0.5);
            overflow-x: auto;
        }}
        .card {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            margin: 20px 0;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.2);
        }}
        .gene-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
        }}
        .gene-table th, .gene-table td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
        }}
        .gene-table th {{
            background: linear-gradient(135deg, #dc3545 0%, #a71d2a 100%);
            color: white;
            font-weight: 600;
        }}
        tr:hover {{ background-color: #f8f9fa; }}
        .summary-stats {{
            display: flex;
            justify-content: space-around;
            margin: 20px 0;
            flex-wrap: wrap;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #dc3545 0%, #a71d2a 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        .critical-stat-card {{
            background: linear-gradient(135deg, #8b0000 0%, #6a0000 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            color: white;
            padding: 20px;
            border-radius: 12px;
            margin: 20px 0;
            text-align: center;
            font-style: italic;
            border-left: 4px solid #dc3545;
        }}
        .footer {{
            background: rgba(0, 0, 0, 0.8);
            color: white;
            padding: 30px;
            border-radius: 12px;
            margin-top: 40px;
        }}
        .footer a {{ color: #dc3545; text-decoration: none; }}
        .footer a:hover {{ text-decoration: underline; }}
        .resistance-badge {{
            display: inline-block;
            background: #dc3545;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .critical-resistance-badge {{
            display: inline-block;
            background: #8b0000;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
            font-weight: bold;
        }}
        .warning-badge {{
            display: inline-block;
            background: #ffc107;
            color: black;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .success-badge {{
            display: inline-block;
            background: #28a745;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .interactive-controls {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            align-items: center;
        }}
        .search-box {{ flex: 1; min-width: 200px; }}
        .search-box input {{
            width: 100%;
            padding: 8px 12px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
        }}
        .export-buttons {{
            display: flex;
            gap: 8px;
            flex-wrap: wrap;
        }}
        .export-buttons button {{
            padding: 8px 16px;
            background: #dc3545;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
            font-size: 14px;
            transition: background 0.3s;
        }}
        .export-buttons button:hover {{ background: #c82333; }}
        .export-buttons button.print {{ background: #10b981; }}
        .export-buttons button.print:hover {{ background: #059669; }}
        .risk-filter {{
            display: flex;
            gap: 5px;
            margin: 10px 0;
        }}
        .risk-filter button {{
            padding: 5px 10px;
            border: 1px solid #ddd;
            background: white;
            border-radius: 3px;
            cursor: pointer;
        }}
        .risk-filter button.active {{ background: #dc3545; color: white; }}
        .result-counter {{
            font-size: 0.9em;
            color: #666;
            font-style: italic;
        }}
        .gene-list-container {{
            padding: 10px;
            background: #f8f9fa;
            border-radius: 5px;
            border: 1px solid #e0e0e0;
        }}
        .gene-item {{
            display: inline-block;
            margin: 2px 4px 2px 0;
            padding: 2px 6px;
            background: #e9ecef;
            border-radius: 3px;
            font-size: 0.85em;
        }}
        .critical-gene {{
            background: #f8d7da;
            color: #721c24;
            border: 1px solid #f5c6cb;
            font-weight: bold;
        }}
        .high-risk-gene {{
            background: #fff3cd;
            color: #856404;
            border: 1px solid #ffeaa7;
        }}
        .sequence-cell {{
            white-space: normal !important;
            word-wrap: break-word;
            max-width: none !important;
            min-width: 300px;
        }}
        .quick-search-bar {{
            background: rgba(255, 255, 255, 0.95);
            padding: 15px;
            border-radius: 8px;
            margin: 10px 0;
            display: flex;
            gap: 10px;
        }}
        .quick-search-bar input {{
            flex: 1;
            padding: 8px 12px;
            border: 2px solid #dc3545;
            border-radius: 4px;
            font-size: 14px;
        }}
        .quick-search-bar button {{
            padding: 8px 20px;
            background: #dc3545;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
        }}
        .critical-row {{ background-color: #f8d7da; font-weight: bold; border-left: 4px solid #dc3545; }}
        .high-risk-row {{ background-color: #fff3cd; border-left: 4px solid #ffc107; }}
        .frequency-high {{ background-color: #f8d7da; font-weight: bold; border-left: 4px solid #dc3545; }}
        .frequency-medium-high {{ background-color: #ffeaa7; border-left: 4px solid #fdcb6e; }}
        .frequency-medium {{ background-color: #fff3cd; border-left: 4px solid #ffc107; }}
        .frequency-low-medium {{ background-color: #d1ecf1; border-left: 4px solid #17a2b8; }}
        .frequency-low {{ background-color: #d4edda; border-left: 4px solid #28a745; }}
        @media print {{
            .no-print, .interactive-controls, .quick-search-bar, .risk-filter {{ display: none !important; }}
            body {{ background: white !important; color: black !important; }}
            .card {{ background: white !important; color: black !important; box-shadow: none !important; }}
            .gene-table {{ box-shadow: none !important; }}
        }}
    </style>
    {interactive_js}
    {export_data_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.ascii_art}</div>
            </div>
            <div class="card">
                <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 PseudoScope AMRfinderPlus - Summary Report</h1>
                <p style="color: #666; font-size: 1.2em;">Comprehensive <em>Pseudomonas aeruginosa</em> Antimicrobial Resistance Analysis Across All Genomes</p>
                <p style="color: #666; font-size: 1.1em;">AMRfinderPlus 4.2.7 | Database: {self.metadata['database_version']}</p>
                <div class="quick-search-bar no-print">
                    <input type="text" id="quick-search"
                           placeholder="🔍 Quick search across entire report (genes, genomes, classes)...">
                    <button onclick="quickSearch()">Search</button>
                    <button onclick="clearSearch()" style="background: #6c757d;">Clear</button>
                </div>
                <div class="export-buttons no-print" style="margin-top: 15px;">
                    <button onclick="exportToJSON('summaryData', 'amr_summary_data.json')">📥 Export JSON</button>
                    <button onclick="printReport()" class="print">🖨️ Print Report</button>
                </div>
            </div>
        </div>
        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
"""
        if critical_genes_found:
            html_content += f"""
        <div class="card" style="border-left: 4px solid #dc3545; background: #f8d7da;">
            <h2 style="color: #dc3545;">🚨 CRITICAL RISK AMR GENES ACROSS ALL GENOMES</h2>
            <p><strong>{len(critical_genes_found)} unique critical risk genes found in {genomes_with_critical} genomes:</strong></p>
            <div style="margin: 10px 0;">
                <p style="color: #721c24; font-weight: bold;">
                    ⚠️ IMMEDIATE ATTENTION REQUIRED: These genes confer resistance to last-resort antibiotics
                </p>
"""
            for gene in sorted(critical_genes_found):
                html_content += f'<span class="critical-resistance-badge">🚨 {gene}</span>'
            html_content += """
            </div>
        </div>
"""
        html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #dc3545; padding-bottom: 10px;">📊 Overall Summary</h2>
            <div class="summary-stats">
                <div class="stat-card">
                    <h3>Total Genomes</h3>
                    <p style="font-size: 2em; margin: 0;">{total_genomes}</p>
                </div>
                <div class="stat-card">
                    <h3>Total AMR Genes</h3>
                    <p style="font-size: 2em; margin: 0;">{total_hits}</p>
                </div>
                <div class="critical-stat-card">
                    <h3>High-Risk Genomes</h3>
                    <p style="font-size: 2em; margin: 0;">{genomes_with_high_risk}</p>
                </div>
            </div>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>AMRfinderPlus:</strong> {self.metadata['amrfinder_version']}</p>
            <p><strong>Database:</strong> {self.metadata['database_version']}</p>
        </div>
"""
        if high_risk_genes_found and not critical_genes_found:
            html_content += f"""
        <div class="card" style="border-left: 4px solid #ffc107;">
            <h2 style="color: #856404;">⚠️ High-Risk AMR Genes Detected</h2>
            <p><strong>{len(high_risk_genes_found)} unique high-risk genes found across {genomes_with_high_risk} genomes:</strong></p>
            <div style="margin: 10px 0;">
"""
            for gene in sorted(high_risk_genes_found):
                html_content += f'<span class="resistance-badge">{gene}</span>'
            html_content += """
            </div>
        </div>
"""
        html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #dc3545; padding-bottom: 10px;">🔍 Genes by Genome</h2>
            <p style="color: #666; margin-bottom: 15px;">Complete list of AMR genes detected in each genome:</p>
            <div class="interactive-controls no-print">
                <div class="search-box">
                    <input type="text" id="search-genes-by-genome"
                           placeholder="Search genomes or genes..."
                           onkeyup="searchTable('genes-by-genome-table', this.value)">
                </div>
                <div class="export-buttons">
                    <button onclick="exportToCSV('genes-by-genome-table', 'genes_by_genome.csv')">📥 Export CSV</button>
                </div>
                <div class="result-counter" id="result-counter-genes-by-genome-table">
                    """ + str(len(genes_per_genome)) + """ results found
                </div>
            </div>
            <table class="gene-table" id="genes-by-genome-table">
                <thead>
                    <tr>
                        <th style="width: 20%;">Genome</th>
                        <th style="width: 10%;">Gene Count</th>
                        <th style="width: 70%;">Genes Detected</th>
                    </tr>
                </thead>
                <tbody>
"""
        for genome in sorted(genes_per_genome.keys()):
            genes = sorted(genes_per_genome.get(genome, set()))
            formatted_genes = ""
            for gene in genes:
                if gene in self.critical_risk_genes:
                    formatted_genes += f'<span class="gene-item critical-gene">{gene}</span>'
                elif gene in self.high_risk_genes:
                    formatted_genes += f'<span class="gene-item high-risk-gene">{gene}</span>'
                else:
                    formatted_genes += f'<span class="gene-item">{gene}</span>'
            critical_genes = [g for g in genes if g in self.critical_risk_genes]
            high_risk_genes_list = [g for g in genes if g in self.high_risk_genes and g not in self.critical_risk_genes]
            row_class = "critical-row" if critical_genes else "high-risk-row" if high_risk_genes_list else ""
            html_content += f"""
                    <tr class="{row_class}">
                        <td><strong>{genome}</strong></td>
                        <td><span style="font-weight: bold; font-size: 1.1em;">{len(genes)}</span></td>
                        <td>
                            <div class="gene-list-container">
                                {formatted_genes}
                            </div>
                        </td>
                    </tr>
"""
        html_content += """
                </tbody>
            </table>
            <div style="margin-top: 15px; font-size: 0.9em; color: #666;">
                <p><span style="display: inline-block; width: 10px; height: 10px; background: #f8d7da; margin-right: 5px;"></span> Critical Risk Genes</p>
                <p><span style="display: inline-block; width: 10px; height: 10px; background: #fff3cd; margin-right: 5px;"></span> High Risk Genes</p>
            </div>
        </div>
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #dc3545; padding-bottom: 10px;">📈 Gene Frequency</h2>
            <p style="color: #666; margin-bottom: 15px;">Frequency of each AMR gene across all genomes:</p>
            <div class="interactive-controls no-print">
                <div class="search-box">
                    <input type="text" id="search-gene-frequency"
                           placeholder="Search genes, genomes, or risk levels..."
                           onkeyup="searchTable('gene-frequency-table', this.value)">
                </div>
                <div class="risk-filter">
                    <button onclick="filterByRisk('all')" class="active">All</button>
                    <button onclick="filterByRisk('critical')">Critical</button>
                    <button onclick="filterByRisk('high')">High</button>
                    <button onclick="filterByRisk('standard')">Standard</button>
                </div>
                <div class="export-buttons">
                    <button onclick="exportToCSV('gene-frequency-table', 'gene_frequency.csv')">📥 Export CSV</button>
                </div>
                <div class="result-counter" id="result-counter-gene-frequency-table">
                    """ + str(len(gene_frequency)) + """ results found
                </div>
            </div>
            <table class="gene-table" id="gene-frequency-table">
                <thead>
                    <tr>
                        <th style="width: 20%;">Gene</th>
                        <th style="width: 15%;">Frequency</th>
                        <th style="width: 15%;">Prevalence</th>
                        <th style="width: 15%;">Risk Level</th>
                        <th style="width: 35%;">Genomes</th>
                    </tr>
                </thead>
                <tbody>
"""
        for gene, genomes in sorted(gene_frequency.items(), key=lambda x: len(x[1]), reverse=True):
            frequency = len(genomes)
            genome_list = ", ".join(sorted(genomes))
            frequency_percent = (frequency / total_genomes) * 100 if total_genomes > 0 else 0
            if gene in self.critical_risk_genes:
                risk_level = '<span class="critical-resistance-badge">CRITICAL</span>'
                risk_class = 'CRITICAL'
            elif gene in self.high_risk_genes:
                risk_level = '<span class="resistance-badge">HIGH</span>'
                risk_class = 'HIGH'
            else:
                risk_level = '<span class="success-badge">Standard</span>'
                risk_class = 'STANDARD'
            if frequency_percent >= 75:
                frequency_class = "frequency-high"
                prevalence_badge = '<span class="resistance-badge">Very High</span>'
            elif frequency_percent >= 50:
                frequency_class = "frequency-medium-high"
                prevalence_badge = '<span class="warning-badge">High</span>'
            elif frequency_percent >= 25:
                frequency_class = "frequency-medium"
                prevalence_badge = '<span class="warning-badge">Medium</span>'
            elif frequency_percent >= 10:
                frequency_class = "frequency-low-medium"
                prevalence_badge = '<span class="success-badge">Low</span>'
            else:
                frequency_class = "frequency-low"
                prevalence_badge = '<span class="success-badge">Rare</span>'
            html_content += f"""
                    <tr class="{frequency_class}" data-risk="{risk_class}">
                        <td><strong>{gene}</strong></td>
                        <td>{frequency} ({frequency_percent:.1f}%)</td>
                        <td>{prevalence_badge}</td>
                        <td>{risk_level}</td>
                        <td class="sequence-cell">{genome_list}</td>
                    </tr>
"""
        html_content += """
                </tbody>
            </table>
        </div>
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #dc3545; padding-bottom: 10px;">📁 Generated Files</h2>
            <ul style="color: #666; font-size: 1.1em;">
                <li><strong>pseudo_amrfinder_summary.tsv</strong> - Complete AMR data for all genomes</li>
                <li><strong>pseudo_amrfinder_statistics_summary.tsv</strong> - Statistical summary</li>
                <li><strong>pseudo_amrfinder_master_summary.json</strong> - Master JSON summary</li>
                <li><strong>Individual genome HTML reports</strong> - Detailed analysis per genome</li>
                <li><strong>Individual genome JSON reports</strong> - JSON data per genome</li>
                <li><strong>This summary report</strong> - Cross-genome analysis with pattern discovery</li>
            </ul>
            <div class="export-buttons no-print" style="margin-top: 20px;">
                <button onclick="window.open('pseudo_amrfinder_summary.tsv')">📄 View TSV Summary</button>
                <button onclick="window.open('pseudo_amrfinder_master_summary.json')">📄 View JSON Summary</button>
            </div>
        </div>
        <div class="footer">
            <h3 style="color: #fff; border-bottom: 2px solid #dc3545; padding-bottom: 10px;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
            <p style="margin-top: 20px; font-size: 0.9em; color: #ccc;">
            </p>
        </div>
    </div>
</body>
</html>
"""
        html_file = os.path.join(output_base, "pseudo_amrfinder_summary_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info("✓ P. aeruginosa AMRfinderPlus summary HTML report created: %s", html_file)

    def process_single_genome(self, genome_file: str, output_base: str = "pseudo_amrfinder_results",
                              min_identity: float = None, min_coverage: float = None,
                              report_mutations: bool = True) -> Dict[str, Any]:
        """Process a single P. aeruginosa genome with AMRfinderPlus"""
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        self.logger.info("=== PROCESSING P. AERUGINOSA GENOME: %s ===", genome_name)
        os.makedirs(results_dir, exist_ok=True)
        if not self.check_amrfinder_installed():
            self.logger.error("Bundled AMRfinderPlus not available!")
            return {
                'genome': genome_name,
                'hits': [],
                'hit_count': 0,
                'status': 'failed',
                'error': 'Bundled AMRfinderPlus not available'
            }
        result = self.run_amrfinder_single_genome(genome_file, results_dir,
                                                  min_identity=min_identity,
                                                  min_coverage=min_coverage,
                                                  report_mutations=report_mutations)
        status_icon = "✓" if result['status'] == 'success' else "✗"
        self.logger.info("%s %s: %d AMR hits", status_icon, genome_name, result['hit_count'])
        return result

    def process_multiple_genomes(self, genome_pattern: str, output_base: str = "pseudo_amrfinder_results",
                                  min_identity: float = None, min_coverage: float = None,
                                  report_mutations: bool = True) -> Dict[str, Any]:
        """Process multiple P. aeruginosa genomes using wildcard pattern - MAXIMUM SPEED"""
        fasta_patterns = [genome_pattern, f"{genome_pattern}.fasta", f"{genome_pattern}.fa",
                         f"{genome_pattern}.fna", f"{genome_pattern}.faa"]
        genome_files = []
        for pattern in fasta_patterns:
            genome_files.extend(glob.glob(pattern))
        genome_files = list(set(genome_files))
        if not genome_files:
            raise FileNotFoundError(f"No FASTA files found matching pattern: {genome_pattern}")
        self.logger.info("Found %d P. aeruginosa genomes: %s", len(genome_files), [Path(f).name for f in genome_files])
        os.makedirs(output_base, exist_ok=True)
        all_results = {}
        max_concurrent = max(1, min(self.cpus, len(genome_files), int(self.available_ram / 2.5)))
        self.logger.info("🚀 MAXIMUM SPEED: Using %d concurrent genome processing jobs", max_concurrent)
        self.logger.info("   Using BUNDLED AMRfinderPlus: %s", self.bundled_amrfinder)
        self.logger.info("   Using BUNDLED database: %s", self.bundled_database)
        with ThreadPoolExecutor(max_workers=max_concurrent) as executor:
            future_to_genome = {
                executor.submit(self.process_single_genome, genome, output_base,
                                min_identity, min_coverage, report_mutations): genome
                for genome in genome_files
            }
            for future in as_completed(future_to_genome):
                genome = future_to_genome[future]
                try:
                    result = future.result()
                    all_results[result['genome']] = result
                    self.logger.info("✓ COMPLETED: %s (%d AMR hits)", result['genome'], result['hit_count'])
                except Exception as e:
                    self.logger.error("✗ FAILED: %s - %s", genome, e)
                    all_results[Path(genome).stem] = {
                        'genome': Path(genome).stem,
                        'hits': [],
                        'hit_count': 0,
                        'status': 'failed'
                    }
        self.create_amr_summary(all_results, output_base)
        if report_mutations:
            self.create_mutation_summary(all_results, output_base)
        self.logger.info("=== P. AERUGINOSA AMR ANALYSIS COMPLETE ===")
        self.logger.info("Processed %d genomes", len(all_results))
        self.logger.info("Results saved to: %s", output_base)
        self.logger.info("Bundled AMRfinderPlus used: %s", self.bundled_amrfinder)
        self.logger.info("Bundled database used: %s", self.bundled_database)
        return all_results


def main():
    """Command line interface for P. aeruginosa AMR analysis with DYNAMIC DATABASE and mutation flags"""
    parser = argparse.ArgumentParser(
        description='PseudoScope AMRfinderPlus Analysis - P. aeruginosa Antimicrobial Resistance - INTERACTIVE ENHANCED VERSION with DYNAMIC DATABASE and Mutation Reporting',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python p_amrfinder.py "*.fna"
  python p_amrfinder.py "*.fna" --min-identity 0.95 --min-coverage 0.9 --skip-mutations
  python p_amrfinder.py --update-db
        """
    )
    parser.add_argument('pattern', nargs='?', help='File pattern for P. aeruginosa genomes (e.g., "*.fna")')
    parser.add_argument('--cpus', '-c', type=int, default=None,
                       help='Number of CPU cores to use (auto-detect optimal by default)')
    parser.add_argument('--output', '-o', default='pseudo_amrfinder_results',
                       help='Output directory (default: pseudo_amrfinder_results)')
    parser.add_argument('--min-identity', type=float, default=None,
                       help='Minimum identity (0..1) for hits')
    parser.add_argument('--min-coverage', type=float, default=None,
                       help='Minimum coverage of reference (0..1)')
    parser.add_argument('--skip-mutations', action='store_true',
                       help='Skip mutation reporting (mutations reported by default)')
    parser.add_argument('--update-db', action='store_true',
                       help='Update AMRfinderPlus database and exit')
    parser.add_argument('--db-version', action='store_true',
                       help='Show current database version and exit')

    args = parser.parse_args()
    executor = PseudoAMRfinderPlus(cpus=args.cpus)

    if args.update_db:
        print("Updating AMRfinderPlus database...")
        success = executor.update_database()
        sys.exit(0 if success else 1)

    if args.db_version:
        print(f"Database version: {executor.metadata['database_version']}")
        print(f"Database path: {executor.bundled_database or 'Not found'}")
        sys.exit(0)

    if not args.pattern:
        parser.error("Please provide a file pattern for genomes (or use --update-db / --db-version)")

    print("\n" + "="*80)
    print(executor.ascii_art)
    print("P. aeruginosa AMRfinderPlus Analysis - DYNAMIC DATABASE and Mutation Reporting")
    print("="*80)
    print(f"Author: Brown Beckley | Email: brownbeckley94@gmail.com")
    print(f"Affiliation: University of Ghana Medical School - Department of Medical Biochemistry")
    print("="*80)

    try:
        results = executor.process_multiple_genomes(args.pattern, args.output,
                                                    min_identity=args.min_identity,
                                                    min_coverage=args.min_coverage,
                                                    report_mutations=not args.skip_mutations)
        executor.logger.info("\n" + "="*50)
        executor.logger.info("🧬 PseudoScope AMRfinderPlus FINAL SUMMARY")
        executor.logger.info("="*50)
        total_hits = 0
        for genome, result in results.items():
            total_hits += result.get('hit_count', 0)
            executor.logger.info(f"✓ {genome}: {result.get('hit_count', 0)} AMR hits")
        executor.logger.info(f"Total genomes processed: {len(results)}")
        executor.logger.info(f"Total AMR hits: {total_hits}")
        executor.logger.info(f"Results saved to: {args.output}")
    except Exception as e:
        executor.logger.error(f"Analysis failed: {e}")
        import traceback
        executor.logger.error(traceback.format_exc())
        sys.exit(1)

if __name__ == "__main__":
    main()