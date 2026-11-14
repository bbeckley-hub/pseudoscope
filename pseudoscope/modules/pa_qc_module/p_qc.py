#!/usr/bin/env python3
"""
PseudoScope FASTA QC - Comprehensive Quality Control for Pseudomonas aeruginosa
Includes species confirmation using fastANI (P. aeruginosa panel + close relatives)
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-10-02 (Enhanced with fastANI species confirmation)
"""

import os
import sys
import glob
import json
import math
import gzip
import statistics
import tempfile
import shutil
import subprocess
from pathlib import Path
from datetime import datetime
from collections import Counter, defaultdict
import argparse
import logging
from typing import List, Dict, Any, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

# BioPython imports
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

class PseudoFASTAQC:
    """Comprehensive FASTA QC with beautiful HTML and fastANI species confirmation for P. aeruginosa"""

    def __init__(self, cpus: int = None, ref_dir: Optional[str] = None):
        # Setup logging
        self.logger = self._setup_logging()

        # ASCII Art for PseudoScope
        self.ascii_art = r"""

██████╗ ███████╗███████╗██╗   ██╗██████╗  ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██╔════╝██║   ██║██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
██████╔╝███████╗█████╗  ██║   ██║██║  ██║██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔═══╝ ╚════██║██╔══╝  ██║   ██║██║  ██║██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║     ███████║███████╗╚██████╔╝██████╔╝╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝     ╚══════╝╚══════╝ ╚═════╝ ╚═════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
"""

        # Metadata
        self.metadata = {
            "tool_name": "PseudoScope FASTA QC for Pseudomonas aeruginosa",
            "version": "1.3.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School - Department of Medical Biochemistry",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "biopython_version": self._get_biopython_version()
        }

        # PseudoScope brown palette
        self.colors = {
            'primary': "#8B4513",
            'success': '#28a745',
            'warning': '#ffc107',
            'danger': '#5C2E0B',
            'info': '#17a2b8',
            'dark': '#1f2937',
            'light': '#f8f9fa'
        }

        # Science quotes
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

        # P. aeruginosa specific thresholds
        self.thresholds = {
            'gc_normal': (60, 70),
            'short_seq': 100,
            'long_seq': 1000000,
            'max_n_run': 100,
            'max_homopolymer': 20,
            'ambiguous_critical': 5.0,
            'ambiguous_warning': 1.0,
            'genome_size_min': 5_000_000,
            'genome_size_max': 8_000_000,
            'n50_min_draft': 10000,
        }

        self.cpus = cpus or os.cpu_count() or 1

        # Species confirmation setup
        if ref_dir is None:
            ref_dir = Path(__file__).resolve().parent / "ref_db"
        self.ref_dir = Path(ref_dir)
        self.species_refs = self._load_references()
        self.fastani_available = shutil.which("fastANI") is not None
        if self.species_refs and not self.fastani_available:
            self.logger.warning("fastANI not found on PATH. Species check will be skipped for all samples.")

    def _setup_logging(self):
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)

    def _get_biopython_version(self) -> str:
        try:
            import Bio
            return Bio.__version__
        except Exception:
            return "Unknown"

    # ------------------------------------------------------------------
    # Species confirmation helpers
    # ------------------------------------------------------------------
    def _load_references(self) -> Dict[str, Path]:
        """Load reference genomes from ref_dir (expected: pseudomonas_aeruginosa_*.fna, etc.)."""
        refs: Dict[str, Path] = {}
        if not self.ref_dir.exists():
            self.logger.warning(
                f"Reference directory {self.ref_dir} not found. Species check disabled."
            )
            return refs
        for f in sorted(self.ref_dir.glob("*.fna")):
            stem = f.stem.replace('_', ' ')
            stem = ' '.join(stem.split())
            species = stem.title()
            refs[species] = f
        if refs:
            self.logger.info(
                f"Loaded {len(refs)} reference genomes for species check: {', '.join(refs.keys())}"
            )
        return refs

    @staticmethod
    def _is_p_aeruginosa(species_name: str) -> bool:
        """Strict match: name must contain both 'pseudomonas' and 'aeruginosa'."""
        s = species_name.lower()
        return ('pseudomonas' in s) and ('aeruginosa' in s)

    def _prepare_query_for_fastani(self, query_fasta: str) -> Tuple[Optional[Path], bool]:
        """
        Return (path, needs_cleanup) for fastANI input.
        If the query is gzipped, decompress to a temporary .fna file.
        """
        q = Path(query_fasta)
        if q.suffix.lower() == '.gz':
            try:
                tmp_fd, tmp_path = tempfile.mkstemp(suffix='.fna', prefix='pseudo_qc_')
                os.close(tmp_fd)
                with gzip.open(q, 'rb') as f_in, open(tmp_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                return Path(tmp_path), True
            except Exception as e:
                self.logger.warning(f"Failed to decompress {q.name}: {e}")
                return None, False
        return q, False

    def _run_species_ani(self, query_fasta: str, threads: int = 4) -> Optional[Dict]:
        """
        Run fastANI against all references; return best match, ANI, and contamination flag.
        Pass when best match is P. aeruginosa with ANI >= 95%.
        Contamination when second best is > 90%, within 5% of the best, and not P. aeruginosa.
        """
        if not self.species_refs:
            return None
        if not self.fastani_available:
            return None

        query_path, needs_cleanup = self._prepare_query_for_fastani(query_fasta)
        if query_path is None:
            return None

        results: List[Tuple[str, float]] = []
        try:
            for species, ref_path in self.species_refs.items():
                tmp_fd, tmp_out = tempfile.mkstemp(suffix='.txt', prefix='fastani_')
                os.close(tmp_fd)
                try:
                    cmd = [
                        'fastANI',
                        '-q', str(query_path),
                        '-r', str(ref_path),
                        '-o', tmp_out,
                        '-t', str(max(1, threads)),
                    ]
                    subprocess.run(cmd, capture_output=True, text=True, check=True)
                    with open(tmp_out, 'r') as fh:
                        line = fh.readline().strip()
                        if line:
                            parts = line.split('\t')
                            if len(parts) >= 3:
                                try:
                                    ani = float(parts[2])
                                    results.append((species, ani))
                                except ValueError:
                                    self.logger.warning(f"Could not parse ANI from line: {line}")
                except subprocess.CalledProcessError:
                    self.logger.warning(
                        f"fastANI failed for {Path(query_fasta).name} vs {species}"
                    )
                except FileNotFoundError:
                    self.logger.error("fastANI executable not found on PATH. Species check disabled.")
                    return None
                finally:
                    try:
                        os.remove(tmp_out)
                    except OSError:
                        pass
        finally:
            if needs_cleanup and query_path is not None:
                try:
                    os.remove(query_path)
                except OSError:
                    pass

        if not results:
            return None

        results.sort(key=lambda x: x[1], reverse=True)
        best_species, best_ani = results[0]
        second_species = results[1][0] if len(results) > 1 else None
        second_ani = results[1][1] if len(results) > 1 else None

        passed = (best_ani >= 95.0 and self._is_p_aeruginosa(best_species))

        contamination = False
        if second_ani is not None and second_ani > 90.0 and (best_ani - second_ani) < 5.0:
            if second_species is not None and not self._is_p_aeruginosa(second_species):
                contamination = True

        return {
            'best_match': best_species,
            'ani_percent': round(best_ani, 2),
            'passed': passed,
            'contamination_suspected': contamination,
            'second_best_match': second_species,
            'second_best_ani': round(second_ani, 2) if second_ani is not None else None,
            'all_matches': {sp: round(ani, 2) for sp, ani in results},
        }

    # ------------------------------------------------------------------
    # Core analysis
    # ------------------------------------------------------------------
    def analyze_file(self, fasta_file: str) -> Dict[str, Any]:
        """Comprehensive analysis of a FASTA file with P. aeruginosa specific checks."""
        try:
            self.logger.info(f"Analyzing {os.path.basename(fasta_file)}...")

            sequences = list(SeqIO.parse(fasta_file, "fasta"))
            if not sequences:
                return {
                    'filename': os.path.basename(fasta_file),
                    'filepath': fasta_file,
                    'status': 'error',
                    'error': 'No sequences found'
                }

            seq_lengths = [len(seq) for seq in sequences]
            total_length = sum(seq_lengths)
            sorted_lengths = sorted(seq_lengths, reverse=True)

            n50 = self._calculate_nx(sorted_lengths, total_length, 50)
            n75 = self._calculate_nx(sorted_lengths, total_length, 75)
            n90 = self._calculate_nx(sorted_lengths, total_length, 90)

            l50 = self._calculate_lx(sorted_lengths, total_length, 50)
            l75 = self._calculate_lx(sorted_lengths, total_length, 75)
            l90 = self._calculate_lx(sorted_lengths, total_length, 90)

            base_counts = Counter()
            gc_contents = []
            for seq in sequences:
                seq_str = str(seq.seq).upper()
                base_counts.update(seq_str)
                gc_contents.append(gc_fraction(seq_str) * 100)

            total_bases = sum(base_counts.values())

            a_count = base_counts.get('A', 0)
            t_count = base_counts.get('T', 0)
            g_count = base_counts.get('G', 0)
            c_count = base_counts.get('C', 0)

            a_percent = (a_count / total_bases) * 100 if total_bases > 0 else 0
            t_percent = (t_count / total_bases) * 100 if total_bases > 0 else 0
            g_percent = (g_count / total_bases) * 100 if total_bases > 0 else 0
            c_percent = (c_count / total_bases) * 100 if total_bases > 0 else 0
            gc_percent = g_percent + c_percent
            at_percent = a_percent + t_percent

            ambiguous_bases = sum(base_counts.get(b, 0) for b in ['N', 'Y', 'R', 'W', 'S', 'K', 'M', 'B', 'D', 'H', 'V'])
            ambiguous_percent = (ambiguous_bases / total_bases) * 100 if total_bases > 0 else 0

            sequences_with_n = sum(1 for seq in sequences if 'N' in str(seq.seq).upper())
            total_n_bases = base_counts.get('N', 0)

            max_n_run = 0
            n_runs = []
            import re
            for seq in sequences:
                seq_str = str(seq.seq).upper()
                for match in re.finditer(r'N+', seq_str):
                    run_len = len(match.group())
                    n_runs.append(run_len)
                    if run_len > max_n_run:
                        max_n_run = run_len

            homopolymers = []
            for seq in sequences:
                seq_str = str(seq.seq).upper()
                for base in ['A', 'T', 'G', 'C']:
                    for match in re.finditer(f'{base}+', seq_str):
                        if len(match.group()) > 4:
                            homopolymers.append({
                                'base': base,
                                'length': len(match.group()),
                                'position': match.start()
                            })

            max_homopolymer = max([h['length'] for h in homopolymers]) if homopolymers else 0

            seq_hashes = set()
            duplicate_sequences = 0
            for seq in sequences:
                seq_hash = hash(str(seq.seq))
                if seq_hash in seq_hashes:
                    duplicate_sequences += 1
                else:
                    seq_hashes.add(seq_hash)

            short_sequences = sum(1 for length in seq_lengths if length < self.thresholds['short_seq'])
            long_sequences = sum(1 for length in seq_lengths if length > self.thresholds['long_seq'])

            length_distribution = self._create_length_bins(seq_lengths)

            # ----- Species confirmation using fastANI -----
            species_result = self._run_species_ani(fasta_file, self.cpus)

            results = {
                'filename': os.path.basename(fasta_file),
                'filepath': fasta_file,
                'file_size_mb': os.path.getsize(fasta_file) / (1024 ** 2),
                'analysis_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'status': 'success',
                'total_sequences': len(sequences),
                'total_length': total_length,
                'total_bases': total_bases,
                'longest_sequence': max(seq_lengths),
                'shortest_sequence': min(seq_lengths),
                'mean_length': statistics.mean(seq_lengths) if seq_lengths else 0,
                'median_length': statistics.median(seq_lengths) if seq_lengths else 0,
                'n50': n50,
                'n75': n75,
                'n90': n90,
                'l50': l50,
                'l75': l75,
                'l90': l90,
                'gc_percent': gc_percent,
                'at_percent': at_percent,
                'a_percent': a_percent,
                't_percent': t_percent,
                'g_percent': g_percent,
                'c_percent': c_percent,
                'ambiguous_percent': ambiguous_percent,
                'sequences_with_n': sequences_with_n,
                'total_n_bases': total_n_bases,
                'max_n_run': max_n_run,
                'homopolymers_count': len(homopolymers),
                'max_homopolymer': max_homopolymer,
                'duplicate_sequences': duplicate_sequences,
                'short_sequences': short_sequences,
                'long_sequences': long_sequences,
                'gc_distribution': {
                    'mean': statistics.mean(gc_contents) if gc_contents else 0,
                    'median': statistics.median(gc_contents) if gc_contents else 0,
                    'min': min(gc_contents) if gc_contents else 0,
                    'max': max(gc_contents) if gc_contents else 0,
                },
                'length_distribution': length_distribution,
                'base_counts': dict(base_counts),
                'n_runs': n_runs,
                'species_check': species_result if species_result else {'error': 'Species check not available'},
                'warnings': self._generate_warnings(
                    gc_percent, ambiguous_percent, max_n_run, max_homopolymer,
                    short_sequences, long_sequences, duplicate_sequences, len(sequences),
                    total_length, n50, species_result
                )
            }

            self.logger.info(
                f"✓ {os.path.basename(fasta_file)}: {len(sequences)} sequences, "
                f"{total_length:,} bp, N50: {n50:,}, GC: {gc_percent:.1f}%"
            )
            return results

        except Exception as e:
            self.logger.error(f"✗ Error analyzing {fasta_file}: {e}")
            return {
                'filename': os.path.basename(fasta_file),
                'status': 'error',
                'error': str(e)
            }

    def _calculate_nx(self, sorted_lengths: List[int], total_length: int, x: int) -> int:
        if not sorted_lengths:
            return 0
        target = total_length * (x / 100)
        cumulative = 0
        for length in sorted_lengths:
            cumulative += length
            if cumulative >= target:
                return length
        return sorted_lengths[-1]

    def _calculate_lx(self, sorted_lengths: List[int], total_length: int, x: int) -> int:
        if not sorted_lengths:
            return 0
        target = total_length * (x / 100)
        cumulative = 0
        for i, length in enumerate(sorted_lengths, 1):
            cumulative += length
            if cumulative >= target:
                return i
        return len(sorted_lengths)

    def _create_length_bins(self, lengths: List[int]) -> Dict[str, int]:
        bins = {
            '< 100 bp': 0, '100-500 bp': 0, '500-1k bp': 0, '1k-5k bp': 0,
            '5k-10k bp': 0, '10k-50k bp': 0, '50k-100k bp': 0, '100k-500k bp': 0,
            '500k-1M bp': 0, '> 1M bp': 0
        }
        for length in lengths:
            if length < 100:
                bins['< 100 bp'] += 1
            elif length < 500:
                bins['100-500 bp'] += 1
            elif length < 1000:
                bins['500-1k bp'] += 1
            elif length < 5000:
                bins['1k-5k bp'] += 1
            elif length < 10000:
                bins['5k-10k bp'] += 1
            elif length < 50000:
                bins['10k-50k bp'] += 1
            elif length < 100000:
                bins['50k-100k bp'] += 1
            elif length < 500000:
                bins['100k-500k bp'] += 1
            elif length < 1000000:
                bins['500k-1M bp'] += 1
            else:
                bins['> 1M bp'] += 1
        return bins

    def _generate_warnings(self, gc_percent: float, ambiguous_percent: float,
                           max_n_run: int, max_homopolymer: int,
                           short_sequences: int, long_sequences: int,
                           duplicate_sequences: int, total_sequences: int,
                           total_length: int, n50: int,
                           species_result: Optional[Dict] = None) -> List[Dict]:
        """Generate warning messages including P. aeruginosa specific checks and species check."""
        warnings = []

        low, high = self.thresholds['gc_normal']
        if gc_percent < low:
            warnings.append({'level': 'warning', 'message': f'Low GC content ({gc_percent:.1f}%) - below typical P. aeruginosa range ({low}-{high}%)'})
        elif gc_percent > high:
            warnings.append({'level': 'warning', 'message': f'High GC content ({gc_percent:.1f}%) - above typical P. aeruginosa range ({low}-{high}%)'})

        genome_min = self.thresholds['genome_size_min']
        genome_max = self.thresholds['genome_size_max']
        if total_length < genome_min:
            warnings.append({'level': 'warning', 'message': f'Total assembly size ({total_length:,} bp) is smaller than expected P. aeruginosa genome (min {genome_min:,} bp)'})
        elif total_length > genome_max:
            warnings.append({'level': 'warning', 'message': f'Total assembly size ({total_length:,} bp) is larger than expected P. aeruginosa genome (max {genome_max:,} bp) - possible contamination'})

        if n50 < self.thresholds['n50_min_draft']:
            warnings.append({'level': 'warning', 'message': f'Low N50 ({n50:,} bp) - assembly may be highly fragmented (typical draft N50 > {self.thresholds["n50_min_draft"]:,} bp)'})

        if ambiguous_percent > self.thresholds['ambiguous_critical']:
            warnings.append({'level': 'danger', 'message': f'High ambiguous bases ({ambiguous_percent:.2f}%) - may indicate poor quality or contamination'})
        elif ambiguous_percent > self.thresholds['ambiguous_warning']:
            warnings.append({'level': 'warning', 'message': f'Elevated ambiguous bases ({ambiguous_percent:.2f}%)'})

        if max_n_run > 100:
            warnings.append({'level': 'danger', 'message': f'Very long N-run detected ({max_n_run} bases) - may indicate assembly gaps'})
        elif max_n_run > 10:
            warnings.append({'level': 'warning', 'message': f'Long N-run detected ({max_n_run} bases)'})

        if max_homopolymer > 20:
            warnings.append({'level': 'danger', 'message': f'Very long homopolymer ({max_homopolymer} bases) - may cause sequencing errors'})
        elif max_homopolymer > 10:
            warnings.append({'level': 'warning', 'message': f'Long homopolymer ({max_homopolymer} bases)'})

        if short_sequences > total_sequences * 0.5:
            warnings.append({'level': 'danger', 'message': f'Many short sequences ({short_sequences}) - may indicate poor assembly'})
        elif short_sequences > total_sequences * 0.1:
            warnings.append({'level': 'warning', 'message': f'Some short sequences ({short_sequences})'})

        if long_sequences > 0:
            warnings.append({'level': 'warning', 'message': f'Very long sequences detected ({long_sequences}) - may indicate contamination'})

        duplicate_percent = (duplicate_sequences / total_sequences) * 100 if total_sequences > 0 else 0
        if duplicate_percent > 10:
            warnings.append({'level': 'warning', 'message': f'Duplicate sequences detected ({duplicate_sequences}, {duplicate_percent:.1f}%)'})

        # Species confirmation warnings
        if species_result and 'error' not in species_result:
            if not species_result.get('passed', False):
                warnings.append({
                    'level': 'danger',
                    'message': (
                        f"Species mismatch: best match is {species_result.get('best_match', 'ND')} "
                        f"(ANI {species_result.get('ani_percent', 'ND')}%). "
                        f"Expected P. aeruginosa (ANI ≥ 95%)."
                    )
                })
            if species_result.get('contamination_suspected', False):
                warnings.append({
                    'level': 'danger',
                    'message': (
                        f"Contamination suspected: second best match is "
                        f"{species_result.get('second_best_match', 'unknown')} "
                        f"(ANI {species_result.get('second_best_ani', 'ND')}%), "
                        f"which is not P. aeruginosa."
                    )
                })

        return warnings

    # ------------------------------------------------------------------
    # Individual HTML report
    # ------------------------------------------------------------------
    def create_individual_html_report(self, results: Dict[str, Any], output_dir: str) -> str:
        """Create comprehensive HTML report for a single FASTA file - Brown theme."""
        filename_no_ext = Path(results['filename']).stem
        sample_dir = os.path.join(output_dir, filename_no_ext)
        os.makedirs(sample_dir, exist_ok=True)

        html_file = os.path.join(sample_dir, f"{filename_no_ext}_fasta_qc_report.html")

        js_content = f"""
        <script>
            function printReport() {{
                window.print();
            }}

            function exportToJSON() {{
                const reportData = {json.dumps(results, indent=2)};
                const dataStr = JSON.stringify(reportData, null, 2);
                const dataBlob = new Blob([dataStr], {{ type: 'application/json' }});
                const url = URL.createObjectURL(dataBlob);
                const a = document.createElement('a');
                a.href = url;
                a.download = '{filename_no_ext}_fasta_qc_report.json';
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                URL.revokeObjectURL(url);
            }}

            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            function rotateQuote() {{
                document.getElementById('science-quote').innerHTML = quotes[currentQuote];
                currentQuote = (currentQuote + 1) % quotes.length;
            }}

            document.addEventListener('DOMContentLoaded', function() {{
                rotateQuote();
                setInterval(rotateQuote, 10000);
            }});

            function sortTable(tableId, columnIndex) {{
                const table = document.getElementById(tableId);
                const tbody = table.querySelector('tbody');
                const rows = Array.from(tbody.querySelectorAll('tr'));

                rows.sort((a, b) => {{
                    const aText = a.children[columnIndex].textContent.trim();
                    const bText = b.children[columnIndex].textContent.trim();
                    const aNum = parseFloat(aText.replace(/,/g, '').replace('%', ''));
                    const bNum = parseFloat(bText.replace(/,/g, '').replace('%', ''));
                    if (!isNaN(aNum) && !isNaN(bNum)) {{
                        return aNum - bNum;
                    }}
                    return aText.localeCompare(bText);
                }});
                rows.forEach(row => tbody.appendChild(row));
            }}
        </script>
        """

        species_check = results.get('species_check', {}) or {}

        # Build species panel HTML (brown themed)
        species_html = ''
        if species_check and 'error' not in species_check:
            passed = species_check.get('passed', False)
            contam = species_check.get('contamination_suspected', False)

            pass_color = '#28a745' if passed else '#C62828'
            contam_color = '#e67e22' if contam else '#28a745'
            pass_text = '✓ Confirmed' if passed else '✗ Not Confirmed'
            contam_text = '⚠️ Suspected' if contam else '✓ None Detected'

            ani_val = species_check.get('ani_percent', 'ND')
            ani_display = f"{ani_val}%" if ani_val not in (None, 'ND') else 'ND'
            second_ani = species_check.get('second_best_ani')
            second_ani_display = f"{second_ani}%" if second_ani is not None else 'ND'

            all_matches_html = ''
            for sp, ani in species_check.get('all_matches', {}).items():
                all_matches_html += f'<li><strong>{sp}</strong> — ANI {ani}%</li>'

            species_html = f'''
        <div class="card" style="border-left: 5px solid #8B4513;">
            <h2 style="color: #5C2E0B; border-bottom: 2px solid #8B4513; padding-bottom: 10px;">🧬 Species Confirmation (fastANI)</h2>
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin-top: 15px;">
                <div style="background: linear-gradient(135deg, #A0522D 0%, #5C2E0B 100%); color: white; padding: 15px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 12px; opacity: 0.9;">Best Match</div>
                    <div style="font-size: 18px; font-weight: bold; margin: 5px 0;">{species_check.get('best_match', 'ND')}</div>
                    <div style="font-size: 12px;">ANI: {ani_display}</div>
                </div>
                <div style="background: linear-gradient(135deg, #A0522D 0%, #5C2E0B 100%); color: white; padding: 15px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 12px; opacity: 0.9;">P. aeruginosa Confirmed</div>
                    <div style="font-size: 18px; font-weight: bold; margin: 5px 0; color: {pass_color};">{pass_text}</div>
                    <div style="font-size: 12px;">Threshold ≥ 95%</div>
                </div>
                <div style="background: linear-gradient(135deg, #A0522D 0%, #5C2E0B 100%); color: white; padding: 15px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 12px; opacity: 0.9;">Contamination</div>
                    <div style="font-size: 18px; font-weight: bold; margin: 5px 0; color: {contam_color};">{contam_text}</div>
                    <div style="font-size: 12px;">Second best &gt; 90% and not P. aeruginosa</div>
                </div>
                <div style="background: linear-gradient(135deg, #A0522D 0%, #5C2E0B 100%); color: white; padding: 15px; border-radius: 8px; text-align: center;">
                    <div style="font-size: 12px; opacity: 0.9;">Second Best Match</div>
                    <div style="font-size: 18px; font-weight: bold; margin: 5px 0;">{species_check.get('second_best_match') or 'ND'}</div>
                    <div style="font-size: 12px;">ANI: {second_ani_display}</div>
                </div>
            </div>
            <details style="margin-top: 15px;">
                <summary style="cursor: pointer; color: #8B4513; font-weight: bold;">All ANI values</summary>
                <ul style="margin-top: 10px; margin-left: 20px; line-height: 1.8;">{all_matches_html}</ul>
            </details>
        </div>
'''
        elif species_check and 'error' in species_check:
            species_html = f'''
        <div class="card" style="border-left: 5px solid #8B4513;">
            <h2 style="color: #5C2E0B; border-bottom: 2px solid #8B4513; padding-bottom: 10px;">🧬 Species Confirmation (fastANI)</h2>
            <p style="color: #666;">Species check not performed: {species_check.get('error', 'unavailable')}</p>
        </div>
'''

        # Species alert banner (only if check failed)
        species_alert = ''
        if species_check and 'error' not in species_check and not species_check.get('passed', False):
            species_alert = f'''
        <div class="card" style="background: #fdecea; border-left: 5px solid #C62828; color: #7f1d1d;">
            <h2 style="color: #C62828; border-bottom: 2px solid #C62828; padding-bottom: 10px;">🚨 Species Mismatch Detected</h2>
            <p style="font-weight: bold; font-size: 1.05em;">
                Best match is <strong>{species_check.get('best_match', 'ND')}</strong>
                (ANI {species_check.get('ani_percent', 'ND')}%) — not <em>Pseudomonas aeruginosa</em>.
            </p>
            <p style="margin-top: 8px;">This genome may be a different species or a mixed culture. Review downstream analyses with caution.</p>
        </div>
'''

        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>PseudoScope FASTA QC Report - {results['filename']}</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #2A1408 0%, #5C2E0B 40%, #8B4513 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1600px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{
            background: rgba(0, 0, 0, 0.75);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(210, 180, 140, 0.4);
        
            text-align: center;
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #E8C39E;
            text-shadow: 0 0 10px rgba(232, 195, 158, 0.4);
            overflow-x: auto;
        
            display: inline-block;
            text-align: left;
        }}
        .card {{
            background: rgba(255, 255, 255, 0.96);
            color: #2c2c2c;
            padding: 25px;
            margin: 20px 0;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.2);
        }}
        .qc-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            font-size: 13px;
        }}
        .qc-table th, .qc-table td {{
            padding: 12px 8px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
            word-wrap: break-word;
        }}
        .qc-table th {{
            background: linear-gradient(135deg, #A0522D 0%, #5C2E0B 100%);
            color: white;
            font-weight: 600;
            cursor: pointer;
            position: relative;
        }}
        .qc-table th:hover {{
            background: linear-gradient(135deg, #8B4513 0%, #3D1E08 100%);
        }}
        .qc-table th.sortable::after {{
            content: "↕";
            position: absolute;
            right: 8px;
            opacity: 0.6;
        }}
        tr:hover {{ background-color: #f8f9fa; }}
        .success {{ color: #28a745; font-weight: 600; }}
        .warning {{ color: #ffc107; font-weight: 600; }}
        .danger {{ color: #5C2E0B; font-weight: 600; }}
        .summary-stats {{
            display: flex;
            justify-content: space-around;
            margin: 20px 0;
            flex-wrap: wrap;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #A0522D 0%, #5C2E0B 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 180px;
        }}
        .gc-stat-card {{
            background: linear-gradient(135deg, #28a745 0%, #1e7e34 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 180px;
        }}
        .at-stat-card {{
            background: linear-gradient(135deg, #17a2b8 0%, #117a8b 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 180px;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.08);
            color: white;
            padding: 20px;
            border-radius: 12px;
            margin: 20px 0;
            text-align: center;
            font-style: italic;
            border-left: 4px solid #8B4513;
        }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0,0,0,0.35);
            border-radius: 10px;
            font-size: 14px;
        }}
        .footer a {{ color: #E8C39E; text-decoration: none; }}
        .footer a:hover {{ text-decoration: underline; }}
        .timestamp {{ color: #E8C39E; font-weight: bold; }}
        .authorship {{
            margin-top: 15px;
            padding: 15px;
            background: rgba(255,255,255,0.08);
            border-radius: 8px;
            font-size: 12px;
        }}
        .controls {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            align-items: center;
        }}
        .export-buttons {{ display: flex; gap: 8px; flex-wrap: wrap; }}
        .export-buttons button {{
            padding: 8px 16px;
            background: #8B4513;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
            font-size: 14px;
            transition: background 0.3s;
        }}
        .export-buttons button:hover {{ background: #A0522D; }}
        .export-buttons button.print {{ background: #10b981; }}
        .export-buttons button.print:hover {{ background: #059669; }}
        .export-buttons button.json {{ background: #8b5cf6; }}
        .export-buttons button.json:hover {{ background: #7c3aed; }}
        .warning-box {{
            background: #fff3cd;
            border: 1px solid #ffc107;
            border-radius: 5px;
            padding: 15px;
            margin: 10px 0;
            color: #856404;
        }}
        .danger-box {{
            background: #f5e6d3;
            border: 1px solid #5C2E0B;
            border-radius: 5px;
            padding: 15px;
            margin: 10px 0;
            color: #3D1E08;
        }}
        .table-container {{ overflow-x: auto; margin: 20px 0; }}
        .composition-grid {{
            display: flex;
            justify-content: space-around;
            flex-wrap: wrap;
            margin: 20px 0;
        }}
        .composition-item {{ text-align: center; margin: 10px; flex: 1; min-width: 120px; }}
        .composition-value {{ font-size: 24px; font-weight: bold; }}
        .composition-label {{ font-size: 12px; color: #666; margin-top: 5px; }}
        @media (max-width: 1200px) {{
            .qc-table {{ font-size: 11px; }}
            .qc-table th, .qc-table td {{ padding: 8px 6px; }}
        }}
        @media print {{
            .controls {{ display: none !important; }}
            body {{ background: white !important; color: black !important; }}
            .card {{ background: white !important; color: black !important; box-shadow: none !important; }}
            .qc-table {{ box-shadow: none !important; }}
        }}
    </style>
    {js_content}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.ascii_art}</div>
            </div>

            <div class="card">
                <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 PseudoScope FASTA QC Report</h1>
                <p style="color: #666; font-size: 1.2em;">Comprehensive Quality Control for <em>Pseudomonas aeruginosa</em> FASTA Files</p>
                <p style="color: #666; font-size: 1.1em;"><strong>Sample:</strong> {results['filename']}</p>

                <div class="controls">
                    <div class="export-buttons">
                        <button onclick="exportToJSON()" class="json">📥 Export JSON</button>
                        <button onclick="printReport()" class="print">🖨️ Print Report</button>
                    </div>
                </div>
            </div>
        </div>

        {species_alert}

        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
"""

        html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #8B4513; padding-bottom: 10px;">📊 FASTA QC Summary</h2>
            <div class="summary-stats">
                <div class="stat-card"><h3>Total Sequences</h3><p style="font-size: 2em; margin: 0;">{results['total_sequences']:,}</p></div>
                <div class="stat-card"><h3>Total Length</h3><p style="font-size: 2em; margin: 0;">{results['total_length']:,}</p><p style="font-size: 0.9em; margin-top: 5px;">bp</p></div>
                <div class="stat-card"><h3>N50</h3><p style="font-size: 2em; margin: 0;">{results['n50']:,}</p><p style="font-size: 0.9em; margin-top: 5px;">bp</p></div>
                <div class="gc-stat-card"><h3>GC Content</h3><p style="font-size: 2em; margin: 0;">{results['gc_percent']:.1f}%</p><p style="font-size: 0.9em; margin-top: 5px;">P. aeruginosa range: {self.thresholds['gc_normal'][0]}-{self.thresholds['gc_normal'][1]}%</p></div>
                <div class="at-stat-card"><h3>AT Content</h3><p style="font-size: 2em; margin: 0;">{results['at_percent']:.1f}%</p></div>
            </div>
            <p><strong>File:</strong> {results['filename']}</p>
            <p><strong>File Size:</strong> {results['file_size_mb']:.2f} MB</p>
            <p><strong>Analysis Date:</strong> {results['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
        </div>
"""

        html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #8B4513; padding-bottom: 10px;">📈 Basic Statistics</h2>
            <div class="table-container">
                <table class="qc-table" id="basic-stats-table">
                    <thead><tr><th class="sortable" onclick="sortTable('basic-stats-table', 0)">Metric</th><th class="sortable" onclick="sortTable('basic-stats-table', 1)">Value</th><th>Description</th></tr></thead>
                    <tbody>
                        <tr><td><strong>Total Sequences</strong></td><td>{results['total_sequences']:,}</td><td>Number of sequences in the file</td></tr>
                        <tr><td><strong>Total Length</strong></td><td>{results['total_length']:,} bp</td><td>Total number of bases (including Ns)</td></tr>
                        <tr><td><strong>Total Bases</strong></td><td>{results['total_bases']:,} bp</td><td>Total bases excluding ambiguous characters</td></tr>
                        <tr><td><strong>Longest Sequence</strong></td><td>{results['longest_sequence']:,} bp</td><td>Length of the longest sequence</td></tr>
                        <tr><td><strong>Shortest Sequence</strong></td><td>{results['shortest_sequence']:,} bp</td><td>Length of the shortest sequence</td></tr>
                        <tr><td><strong>Mean Length</strong></td><td>{results['mean_length']:,.0f} bp</td><td>Average sequence length</td></tr>
                        <tr><td><strong>Median Length</strong></td><td>{results['median_length']:,} bp</td><td>Median sequence length</td></tr>
                        <tr><td><strong>N50</strong></td><td>{results['n50']:,} bp</td><td>Length for which 50% of total bases are in longer sequences</td></tr>
                        <tr><td><strong>N75</strong></td><td>{results['n75']:,} bp</td><td>Length for which 75% of total bases are in longer sequences</td></tr>
                        <tr><td><strong>N90</strong></td><td>{results['n90']:,} bp</td><td>Length for which 90% of total bases are in longer sequences</td></tr>
                        <tr><td><strong>L50</strong></td><td>{results['l50']:,}</td><td>Number of sequences that make up 50% of total length</td></tr>
                        <tr><td><strong>L75</strong></td><td>{results['l75']:,}</td><td>Number of sequences that make up 75% of total length</td></tr>
                        <tr><td><strong>L90</strong></td><td>{results['l90']:,}</td><td>Number of sequences that make up 90% of total length</td></tr>
                        <tr><td><strong>GC Content</strong></td><td>{results['gc_percent']:.1f}%</td><td>Percentage of G and C bases</td></tr>
                        <tr><td><strong>AT Content</strong></td><td>{results['at_percent']:.1f}%</td><td>Percentage of A and T bases</td></tr>
                        <tr><td><strong>Ambiguous Bases</strong></td><td>{results['ambiguous_percent']:.2f}%</td><td>Percentage of N and other ambiguous bases</td></tr>
                        <tr><td><strong>Sequences with Ns</strong></td><td>{results['sequences_with_n']:,}</td><td>Number of sequences containing N bases</td></tr>
                        <tr><td><strong>Total N Bases</strong></td><td>{results['total_n_bases']:,}</td><td>Total count of N bases</td></tr>
                        <tr><td><strong>Max N-run</strong></td><td>{results['max_n_run']:,}</td><td>Longest consecutive run of N bases</td></tr>
                        <tr><td><strong>Homopolymers</strong></td><td>{results['homopolymers_count']:,}</td><td>Number of homopolymers longer than 4 bases</td></tr>
                        <tr><td><strong>Max Homopolymer</strong></td><td>{results['max_homopolymer']:,}</td><td>Longest homopolymer run</td></tr>
                        <tr><td><strong>Duplicate Sequences</strong></td><td>{results['duplicate_sequences']:,}</td><td>Number of identical sequences</td></tr>
                        <tr><td><strong>Short Sequences (&lt;{self.thresholds['short_seq']} bp)</strong></td><td>{results['short_sequences']:,}</td><td>Sequences shorter than {self.thresholds['short_seq']} bp</td></tr>
                        <tr><td><strong>Long Sequences (&gt;{self.thresholds['long_seq']:,} bp)</strong></td><td>{results['long_sequences']:,}</td><td>Sequences longer than {self.thresholds['long_seq']:,} bp</td></tr>
                    </tbody>
                </table>
            </div>
        </div>
"""

        html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #8B4513; padding-bottom: 10px;">🧬 Nucleotide Composition</h2>
            <div class="composition-grid">
                <div class="composition-item"><div class="composition-value" style="color: #28a745;">{results['gc_percent']:.1f}%</div><div class="composition-label">GC Content</div></div>
                <div class="composition-item"><div class="composition-value" style="color: #17a2b8;">{results['at_percent']:.1f}%</div><div class="composition-label">AT Content</div></div>
                <div class="composition-item"><div class="composition-value">{results['a_percent']:.1f}%</div><div class="composition-label">Adenine (A)</div></div>
                <div class="composition-item"><div class="composition-value">{results['t_percent']:.1f}%</div><div class="composition-label">Thymine (T)</div></div>
                <div class="composition-item"><div class="composition-value">{results['g_percent']:.1f}%</div><div class="composition-label">Guanine (G)</div></div>
                <div class="composition-item"><div class="composition-value">{results['c_percent']:.1f}%</div><div class="composition-label">Cytosine (C)</div></div>
                <div class="composition-item"><div class="composition-value" style="color: #5C2E0B;">{results['ambiguous_percent']:.2f}%</div><div class="composition-label">Ambiguous Bases</div></div>
            </div>
        </div>
"""

        html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #8B4513; padding-bottom: 10px;">📊 Length Distribution</h2>
            <div class="table-container">
                <table class="qc-table" id="length-dist-table">
                    <thead><tr><th class="sortable" onclick="sortTable('length-dist-table', 0)">Length Range</th><th class="sortable" onclick="sortTable('length-dist-table', 1)">Count</th><th class="sortable" onclick="sortTable('length-dist-table', 2)">Percentage</th></tr></thead>
                    <tbody>
"""

        total_seqs = results['total_sequences']
        for length_range, count in results['length_distribution'].items():
            percentage = (count / total_seqs * 100) if total_seqs > 0 else 0
            html_content += f"""
                        <tr><td>{length_range}</td><td>{count:,}</td><td>{percentage:.1f}%</td></tr>
"""

        html_content += """
                    </tbody>
                </table>
            </div>
        </div>
"""

        # Species panel
        html_content += species_html

        # Warnings
        if results.get('warnings'):
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #8B4513; padding-bottom: 10px;">⚠️ Quality Warnings</h2>
"""
            for warning in results['warnings']:
                box_class = 'danger-box' if warning['level'] == 'danger' else 'warning-box'
                html_content += f"""
            <div class="{box_class}">
                <strong>{warning['level'].upper()}:</strong> {warning['message']}
            </div>
"""
            html_content += """
        </div>
"""

        html_content += f"""
        <div class="footer">
<p><strong>PSEUDOSCOPE</strong> - Quality Control Analysis</p>
<p class="timestamp">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
<div class="authorship">
<p><strong>Technical Support &amp; Inquiries:</strong></p>
<p>Maintainer: Brown Beckley</p>
<p>GitHub: <a href="https://github.com/bbeckley-hub" style="color:#E8C39E;">bbeckley-hub</a></p>
<p>Email: <a href="mailto:brownbeckley94@gmail.com" style="color:#E8C39E;">brownbeckley94@gmail.com</a></p>
<p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
<p style="margin-top:10px;font-size:11px;opacity:0.8;">
Powered by Biopython and fastANI
</p>
</div>
</div>
    </div>
</body>
</html>
"""

        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        self.logger.info(f"✓ HTML report generated: {html_file}")
        return html_file

    # ------------------------------------------------------------------
    # Summary report
    # ------------------------------------------------------------------
    def create_summary_report(self, all_results: List[Dict[str, Any]], output_dir: str):
        """Create summary report for multiple FASTA files."""
        successful_results = [r for r in all_results if r.get('status') == 'success']
        if not successful_results:
            self.logger.warning("No successful analyses to create summary report")
            return

        tsv_file = os.path.join(output_dir, "PseudoScope_FASTA_QC_summary.tsv")
        with open(tsv_file, 'w', encoding='utf-8') as f:
            headers = [
                'Filename', 'Total Sequences', 'Total Length', 'Total Bases',
                'GC Content (%)', 'AT Content (%)', 'N50', 'N75', 'N90',
                'Median Length', 'Mean Length', 'Longest Sequence', 'Shortest Sequence',
                'Ambiguous Bases (%)', 'Sequences with Ns', 'Max N-run',
                'Homopolymers', 'Max Homopolymer', 'Duplicate Sequences',
                'Short Sequences (<100 bp)', 'Long Sequences (>1M bp)',
                'File Size (MB)', 'Warnings',
                'Best Species', 'ANI (%)', 'P. aeruginosa Confirmed', 'Contamination Suspected'
            ]
            f.write('\t'.join(headers) + '\n')

            for result in successful_results:
                sc = result.get('species_check', {}) or {}
                best_species = sc.get('best_match', 'ND')
                ani = sc.get('ani_percent', 'ND')
                pa_confirmed = 'Yes' if sc.get('passed') else 'No'
                contam = 'Yes' if sc.get('contamination_suspected') else 'No'

                row = [
                    result['filename'],
                    str(result['total_sequences']),
                    str(result['total_length']),
                    str(result['total_bases']),
                    f"{result['gc_percent']:.2f}",
                    f"{result['at_percent']:.2f}",
                    str(result['n50']),
                    str(result['n75']),
                    str(result['n90']),
                    str(result['median_length']),
                    f"{result['mean_length']:.0f}",
                    str(result['longest_sequence']),
                    str(result['shortest_sequence']),
                    f"{result['ambiguous_percent']:.2f}",
                    str(result['sequences_with_n']),
                    str(result['max_n_run']),
                    str(result['homopolymers_count']),
                    str(result['max_homopolymer']),
                    str(result['duplicate_sequences']),
                    str(result['short_sequences']),
                    str(result['long_sequences']),
                    f"{result['file_size_mb']:.2f}",
                    str(len(result.get('warnings', []))),
                    best_species,
                    str(ani),
                    pa_confirmed,
                    contam,
                ]
                f.write('\t'.join(row) + '\n')

        json_summary = self._create_json_summary(successful_results)
        json_file = os.path.join(output_dir, "PseudoScope_FASTA_QC_summary.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(json_summary, f, indent=2, default=str)

        html_file = os.path.join(output_dir, "PseudoScope_FASTA_QC_summary.html")
        self._create_summary_html_report(successful_results, html_file, json_summary)

        self.logger.info(f"✓ TSV summary created: {tsv_file}")
        self.logger.info(f"✓ JSON summary created: {json_file}")
        self.logger.info(f"✓ HTML summary created: {html_file}")

    def _create_json_summary(self, successful_results: List[Dict[str, Any]]) -> Dict[str, Any]:
        summary_data = {
            'metadata': {
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'biopython_version': self.metadata['biopython_version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_files': len(successful_results),
                'total_sequences': sum(r['total_sequences'] for r in successful_results),
                'total_length': sum(r['total_length'] for r in successful_results),
                'species_check_enabled': bool(self.species_refs) and self.fastani_available,
                'reference_panel': list(self.species_refs.keys()),
            },
            'statistics': {
                'gc_content_range': {
                    'min': min(r['gc_percent'] for r in successful_results) if successful_results else 0,
                    'max': max(r['gc_percent'] for r in successful_results) if successful_results else 0,
                    'mean': statistics.mean(r['gc_percent'] for r in successful_results) if successful_results else 0,
                    'median': statistics.median(r['gc_percent'] for r in successful_results) if successful_results else 0,
                },
                'n50_range': {
                    'min': min(r['n50'] for r in successful_results) if successful_results else 0,
                    'max': max(r['n50'] for r in successful_results) if successful_results else 0,
                    'mean': statistics.mean(r['n50'] for r in successful_results) if successful_results else 0,
                    'median': statistics.median(r['n50'] for r in successful_results) if successful_results else 0,
                },
                'files_with_warnings': sum(1 for r in successful_results if r.get('warnings')),
                'total_warnings': sum(len(r.get('warnings', [])) for r in successful_results),
                'files_with_high_gc': sum(1 for r in successful_results if r['gc_percent'] > self.thresholds['gc_normal'][1]),
                'files_with_low_gc': sum(1 for r in successful_results if r['gc_percent'] < self.thresholds['gc_normal'][0]),
                'files_confirmed_p_aeruginosa': sum(
                    1 for r in successful_results
                    if (r.get('species_check') or {}).get('passed', False)
                ),
                'files_with_contamination_suspected': sum(
                    1 for r in successful_results
                    if (r.get('species_check') or {}).get('contamination_suspected', False)
                ),
            },
            'files': successful_results
        }
        return summary_data

    def _create_summary_html_report(self, successful_results: List[Dict[str, Any]],
                                    html_file: str, json_summary: Dict[str, Any]):
        """Create HTML summary report with scrollable table - Brown theme."""

        js_content = f"""
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
                    resultCounter.textContent = visibleCount + ' files found';
                }}
            }}

            function sortTable(tableId, columnIndex) {{
                const table = document.getElementById(tableId);
                const tbody = table.querySelector('tbody');
                const rows = Array.from(tbody.querySelectorAll('tr'));
                rows.sort((a, b) => {{
                    const aText = a.children[columnIndex].textContent.trim();
                    const bText = b.children[columnIndex].textContent.trim();
                    const aNum = parseFloat(aText.replace(/,/g, '').replace('%', ''));
                    const bNum = parseFloat(bText.replace(/,/g, '').replace('%', ''));
                    if (!isNaN(aNum) && !isNaN(bNum)) {{ return aNum - bNum; }}
                    return aText.localeCompare(bText);
                }});
                rows.forEach(row => tbody.appendChild(row));
            }}

            function exportToCSV() {{
                const rows = document.querySelectorAll('#qc-summary-table tr');
                let csv = [];
                const headerCells = rows[0].querySelectorAll('th');
                const headerRow = [];
                for (let cell of headerCells) {{ headerRow.push(cell.textContent); }}
                csv.push(headerRow.join(','));
                for (let i = 1; i < rows.length; i++) {{
                    if (rows[i].style.display !== 'none') {{
                        const cells = rows[i].querySelectorAll('td');
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
                a.download = 'PseudoScope_FASTA_QC_summary.csv';
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);
            }}

            function exportToJSON() {{
                const data = {json.dumps(json_summary, indent=2)};
                const jsonStr = JSON.stringify(data, null, 2);
                const blob = new Blob([jsonStr], {{ type: 'application/json' }});
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = 'PseudoScope_FASTA_QC_summary.json';
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);
            }}

            function printReport() {{ window.print(); }}

            function filterByWarnings(level) {{
                const table = document.getElementById('qc-summary-table');
                const rows = table.getElementsByTagName('tr');
                for (let i = 1; i < rows.length; i++) {{
                    const row = rows[i];
                    const warningsCell = row.children[22];
                    const warnings = parseInt(warningsCell.textContent);
                    if (level === 'all') {{ row.style.display = ''; }}
                    else if (level === 'with_warnings' && warnings > 0) {{ row.style.display = ''; }}
                    else if (level === 'no_warnings' && warnings === 0) {{ row.style.display = ''; }}
                    else {{ row.style.display = 'none'; }}
                }}
            }}

            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            function rotateQuote() {{
                document.getElementById('science-quote').innerHTML = quotes[currentQuote];
                currentQuote = (currentQuote + 1) % quotes.length;
            }}
            document.addEventListener('DOMContentLoaded', function() {{
                rotateQuote();
                setInterval(rotateQuote, 10000);
            }});
        </script>
        """

        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>PseudoScope FASTA QC - Summary Report</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #2A1408 0%, #5C2E0B 40%, #8B4513 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1800px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{
            background: rgba(0, 0, 0, 0.75);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(210, 180, 140, 0.4);
        
            text-align: center;
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #E8C39E;
            text-shadow: 0 0 10px rgba(232, 195, 158, 0.4);
            overflow-x: auto;
        
            display: inline-block;
            text-align: left;
        }}
        .card {{
            background: rgba(255, 255, 255, 0.96);
            color: #2c2c2c;
            padding: 25px;
            margin: 20px 0;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.2);
        }}
        .qc-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
            font-size: 12px;
        }}
        .qc-table th, .qc-table td {{
            padding: 10px 8px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
            white-space: nowrap;
        }}
        .qc-table th {{
            background: linear-gradient(135deg, #A0522D 0%, #5C2E0B 100%);
            color: white;
            font-weight: 600;
            cursor: pointer;
            position: sticky;
            top: 0;
            z-index: 10;
        }}
        .qc-table th:hover {{
            background: linear-gradient(135deg, #8B4513 0%, #3D1E08 100%);
        }}
        .qc-table th.sortable::after {{
            content: "↕";
            position: absolute;
            right: 8px;
            opacity: 0.6;
        }}
        tr:hover {{ background-color: #f8f9fa; }}
        .summary-stats {{
            display: flex;
            justify-content: space-around;
            margin: 20px 0;
            flex-wrap: wrap;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #A0522D 0%, #5C2E0B 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            margin: 10px;
            flex: 1;
            min-width: 180px;
        }}
        .seq-stat-card {{
            background: linear-gradient(135deg, #28a745 0%, #1e7e34 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            margin: 10px;
            flex: 1;
            min-width: 180px;
        }}
        .warning-stat-card {{
            background: linear-gradient(135deg, #ffc107 0%, #e0a800 100%);
            color: black;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            margin: 10px;
            flex: 1;
            min-width: 180px;
        }}
        .species-stat-card {{
            background: linear-gradient(135deg, #17a2b8 0%, #117a8b 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
            margin: 10px;
            flex: 1;
            min-width: 180px;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.08);
            color: white;
            padding: 20px;
            border-radius: 12px;
            margin: 20px 0;
            text-align: center;
            font-style: italic;
            border-left: 4px solid #8B4513;
        }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0,0,0,0.35);
            border-radius: 10px;
            font-size: 14px;
        }}
        .footer a {{ color: #E8C39E; text-decoration: none; }}
        .footer a:hover {{ text-decoration: underline; }}
        .timestamp {{ color: #E8C39E; font-weight: bold; }}
        .authorship {{
            margin-top: 15px;
            padding: 15px;
            background: rgba(255,255,255,0.08);
            border-radius: 8px;
            font-size: 12px;
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
        .filter-buttons {{ display: flex; gap: 5px; flex-wrap: wrap; }}
        .filter-buttons button {{
            padding: 6px 12px;
            background: #e9ecef;
            color: #495057;
            border: 1px solid #dee2e6;
            border-radius: 4px;
            cursor: pointer;
            font-size: 13px;
        }}
        .filter-buttons button.active {{
            background: #8B4513;
            color: white;
            border-color: #8B4513;
        }}
        .filter-buttons button:hover {{ background: #dee2e6; }}
        .export-buttons {{ display: flex; gap: 8px; flex-wrap: wrap; }}
        .export-buttons button {{
            padding: 8px 16px;
            background: #8B4513;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
            font-size: 14px;
            transition: background 0.3s;
        }}
        .export-buttons button:hover {{ background: #A0522D; }}
        .export-buttons button.print {{ background: #10b981; }}
        .export-buttons button.print:hover {{ background: #059669; }}
        .export-buttons button.json {{ background: #8b5cf6; }}
        .export-buttons button.json:hover {{ background: #7c3aed; }}
        .result-counter {{
            font-size: 0.9em;
            color: #666;
            font-style: italic;
            margin-left: auto;
        }}
        .table-container {{
            overflow-x: auto;
            margin: 20px 0;
            max-height: 600px;
            overflow-y: auto;
        }}
        .table-container::-webkit-scrollbar {{ width: 8px; height: 8px; }}
        .table-container::-webkit-scrollbar-track {{ background: #f1f1f1; border-radius: 4px; }}
        .table-container::-webkit-scrollbar-thumb {{ background: #888; border-radius: 4px; }}
        .table-container::-webkit-scrollbar-thumb:hover {{ background: #555; }}
        .warning-count {{
            display: inline-block;
            padding: 2px 8px;
            border-radius: 12px;
            font-size: 11px;
            font-weight: bold;
        }}
        .warning-count.none {{ background: #d4edda; color: #155724; }}
        .warning-count.low {{ background: #fff3cd; color: #856404; }}
        .warning-count.high {{ background: #f5e6d3; color: #3D1E08; }}
        .species-pass {{ color: #28a745; font-weight: bold; }}
        .species-fail {{ color: #C62828; font-weight: bold; }}
        .contam-warn {{ color: #e67e22; font-weight: bold; }}
        @media print {{
            .interactive-controls {{ display: none !important; }}
            body {{ background: white !important; color: black !important; }}
            .card {{ background: white !important; color: black !important; box-shadow: none !important; }}
            .qc-table {{ box-shadow: none !important; }}
            .table-container {{ max-height: none !important; overflow: visible !important; }}
        }}
    </style>
    {js_content}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.ascii_art}</div>
            </div>

            <div class="card">
                <h1 style="color: #333; margin: 0; font-size: 2.5em;">🧬 PseudoScope FASTA QC - Summary Report</h1>
                <p style="color: #666; font-size: 1.2em;">Comprehensive Quality Control + fastANI Species Confirmation for <em>Pseudomonas aeruginosa</em></p>

                <div class="interactive-controls">
                    <div class="search-box">
                        <input type="text" id="search-summary"
                               placeholder="🔍 Search filenames, species, values..."
                               onkeyup="searchTable('qc-summary-table', this.value)">
                    </div>

                    <div class="filter-buttons">
                        <button class="active" onclick="filterByWarnings('all')">All Files</button>
                        <button onclick="filterByWarnings('with_warnings')">With Warnings</button>
                        <button onclick="filterByWarnings('no_warnings')">No Warnings</button>
                    </div>

                    <div class="result-counter" id="result-counter-qc-summary-table">
                        {len(successful_results)} files found
                    </div>

                    <div class="export-buttons">
                        <button onclick="exportToCSV()">📥 Export CSV</button>
                        <button onclick="exportToJSON()" class="json">📥 Export JSON</button>
                        <button onclick="printReport()" class="print">🖨️ Print</button>
                    </div>
                </div>
            </div>
        </div>

        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
"""

        total_sequences = sum(r['total_sequences'] for r in successful_results)
        total_length = sum(r['total_length'] for r in successful_results)
        total_warnings = sum(len(r.get('warnings', [])) for r in successful_results)
        files_with_warnings = sum(1 for r in successful_results if r.get('warnings'))
        confirmed_pa = sum(1 for r in successful_results if (r.get('species_check') or {}).get('passed', False))
        contam_flagged = sum(1 for r in successful_results if (r.get('species_check') or {}).get('contamination_suspected', False))

        html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #8B4513; padding-bottom: 10px;">📊 Overall Statistics</h2>
            <div class="summary-stats">
                <div class="stat-card"><h3>Total Files</h3><p style="font-size: 2em; margin: 0;">{len(successful_results)}</p></div>
                <div class="seq-stat-card"><h3>Total Sequences</h3><p style="font-size: 2em; margin: 0;">{total_sequences:,}</p></div>
                <div class="seq-stat-card"><h3>Total Bases</h3><p style="font-size: 2em; margin: 0;">{total_length:,}</p><p style="font-size: 0.9em; margin-top: 5px;">bp</p></div>
                <div class="warning-stat-card"><h3>Quality Warnings</h3><p style="font-size: 2em; margin: 0;">{total_warnings}</p><p style="font-size: 0.9em; margin-top: 5px;">in {files_with_warnings} files</p></div>
                <div class="species-stat-card"><h3>Confirmed P. aeruginosa</h3><p style="font-size: 2em; margin: 0;">{confirmed_pa} / {len(successful_results)}</p><p style="font-size: 0.9em; margin-top: 5px;">Contamination flagged: {contam_flagged}</p></div>
            </div>
            <p><strong>Analysis Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>BioPython Version:</strong> {self.metadata['biopython_version']}</p>
            <p><strong>Species Reference Panel:</strong> {', '.join(self.species_refs.keys()) if self.species_refs else 'Disabled'}</p>
        </div>
"""

        html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #8B4513; padding-bottom: 10px;">📈 Detailed FASTA QC Results + Species Check</h2>
            <p style="color: #666; margin-bottom: 15px;">Scroll horizontally. Species column shows best match + confirmation (✓/✗) + contamination flag (⚠️). ANI (%) is identity to best match.</p>

            <div class="table-container">
                <table class="qc-table" id="qc-summary-table">
                    <thead>
                        <tr>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 0)">Filename</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 1)">Total Sequences</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 2)">Total Length</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 3)">Total Bases</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 4)">GC Content (%)</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 5)">AT Content (%)</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 6)">N50</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 7)">N75</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 8)">N90</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 9)">Median Length</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 10)">Mean Length</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 11)">Longest Sequence</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 12)">Shortest Sequence</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 13)">Ambiguous Bases (%)</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 14)">Sequences with Ns</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 15)">Max N-run</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 16)">Homopolymers</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 17)">Max Homopolymer</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 18)">Duplicate Sequences</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 19)">Short Sequences (&lt;100 bp)</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 20)">Long Sequences (&gt;1M bp)</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 21)">File Size (MB)</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 22)">Warnings</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 23)">ANI (%)</th>
                            <th class="sortable" onclick="sortTable('qc-summary-table', 24)">Species</th>
                        </tr>
                    </thead>
                    <tbody>
"""

        for result in successful_results:
            warning_count = len(result.get('warnings', []))
            warning_class = 'none' if warning_count == 0 else 'low' if warning_count < 3 else 'high'

            sc = result.get('species_check', {}) or {}
            best_species = sc.get('best_match', 'ND')
            ani = sc.get('ani_percent', 'ND')
            ani_display = f"{ani}%" if ani not in (None, 'ND') else 'ND'
            if sc.get('passed'):
                confirm_icon = '<span class="species-pass">✓</span>'
            elif 'error' in sc:
                confirm_icon = '<span style="color:#888;">?</span>'
            else:
                confirm_icon = '<span class="species-fail">✗</span>'
            contam_icon = '<span class="contam-warn">⚠️</span>' if sc.get('contamination_suspected') else ''
            species_display = f"{best_species} {confirm_icon} {contam_icon}".strip()

            html_content += f"""
                        <tr>
                            <td><strong>{result['filename']}</strong></td>
                            <td>{result['total_sequences']:,}</td>
                            <td>{result['total_length']:,}</td>
                            <td>{result['total_bases']:,}</td>
                            <td>{result['gc_percent']:.1f}</td>
                            <td>{result['at_percent']:.1f}</td>
                            <td>{result['n50']:,}</td>
                            <td>{result['n75']:,}</td>
                            <td>{result['n90']:,}</td>
                            <td>{result['median_length']:,}</td>
                            <td>{result['mean_length']:.0f}</td>
                            <td>{result['longest_sequence']:,}</td>
                            <td>{result['shortest_sequence']:,}</td>
                            <td>{result['ambiguous_percent']:.2f}</td>
                            <td>{result['sequences_with_n']:,}</td>
                            <td>{result['max_n_run']:,}</td>
                            <td>{result['homopolymers_count']:,}</td>
                            <td>{result['max_homopolymer']:,}</td>
                            <td>{result['duplicate_sequences']:,}</td>
                            <td>{result['short_sequences']:,}</td>
                            <td>{result['long_sequences']:,}</td>
                            <td>{result['file_size_mb']:.2f}</td>
                            <td><span class="warning-count {warning_class}">{warning_count}</span></td>
                            <td>{ani_display}</td>
                            <td>{species_display}</td>
                        </tr>
"""

        html_content += """
                    </tbody>
                </table>
            </div>
        </div>
"""

        html_content += f"""
        <div class="footer">
<p><strong>PSEUDOSCOPE</strong> - Quality Control Analysis</p>
<p class="timestamp">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
<div class="authorship">
<p><strong>Technical Support &amp; Inquiries:</strong></p>
<p>Maintainer: Brown Beckley</p>
<p>GitHub: <a href="https://github.com/bbeckley-hub" style="color:#E8C39E;">bbeckley-hub</a></p>
<p>Email: <a href="mailto:brownbeckley94@gmail.com" style="color:#E8C39E;">brownbeckley94@gmail.com</a></p>
<p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
<p style="margin-top:10px;font-size:11px;opacity:0.8;">
Powered by Biopython and fastANI
</p>
</div>
</div>
    </div>
</body>
</html>
"""

        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

    # ------------------------------------------------------------------
    # Driver
    # ------------------------------------------------------------------
    def process_files(self, pattern: str, output_dir: str = "fasta_qc_results") -> List[Dict[str, Any]]:
        """Process multiple FASTA files with parallel execution and comprehensive reporting."""
        print("\n" + "="*80)
        print(self.ascii_art)
        print("PseudoScope FASTA QC - Comprehensive Quality Control for Pseudomonas aeruginosa")
        print("="*80)
        print(f"Author: Brown Beckley | Email: brownbeckley94@gmail.com")
        print(f"Affiliation: University of Ghana Medical School - Department of Medical Biochemistry")
        print("="*80)
        print(f"Output directory: {output_dir}")
        print(f"Using {self.cpus} CPU cores")
        if self.species_refs and self.fastani_available:
            print(f"Species references: {', '.join(self.species_refs.keys())}")
        elif self.species_refs and not self.fastani_available:
            print("⚠️  Species check disabled: fastANI not found on PATH")
        else:
            print("⚠️  Species check disabled: no references found in ref_db/")
        print("="*80)

        fasta_extensions = ['.fna', '.fasta', '.fa', '.faa', '.fn', '.fna.gz', '.fasta.gz', '.fa.gz']
        files = []
        files.extend(glob.glob(pattern))
        if not files:
            for ext in fasta_extensions:
                files.extend(glob.glob(f"{pattern}{ext}"))
        files = list(set([f for f in files if os.path.exists(f)]))

        if not files:
            self.logger.error(f"No FASTA files found matching pattern: {pattern}")
            self.logger.info(f"Supported extensions: {', '.join(fasta_extensions)}")
            return []

        self.logger.info(f"Found {len(files)} FASTA files: {[Path(f).name for f in files]}")
        os.makedirs(output_dir, exist_ok=True)

        all_results = []
        successful_files = 0

        with ThreadPoolExecutor(max_workers=self.cpus) as executor:
            future_to_file = {executor.submit(self.analyze_file, f): f for f in files}
            for future in as_completed(future_to_file):
                file = future_to_file[future]
                try:
                    result = future.result()
                    all_results.append(result)
                    if result.get('status') == 'success':
                        successful_files += 1
                        self.create_individual_html_report(result, output_dir)
                    else:
                        self.logger.error(f"✗ {result['filename']}: {result.get('error', 'Unknown error')}")
                except Exception as e:
                    self.logger.error(f"✗ {os.path.basename(file)}: ERROR - {str(e)}")
                    all_results.append({
                        'filename': os.path.basename(file),
                        'status': 'error',
                        'error': str(e)
                    })

        if successful_files > 0:
            self.create_summary_report(all_results, output_dir)

        print("\n" + "="*80)
        print("ANALYSIS COMPLETE")
        print("="*80)
        print(f"Total files: {len(files)}")
        print(f"Successful analyses: {successful_files}")
        print(f"Failed analyses: {len(files) - successful_files}")
        print("\n📁 OUTPUT FILES:")
        print(f"   Individual HTML reports: {output_dir}/*/*_fasta_qc_report.html")
        if successful_files > 1:
            print(f"   Summary TSV: {output_dir}/PseudoScope_FASTA_QC_summary.tsv")
            print(f"   Summary JSON: {output_dir}/PseudoScope_FASTA_QC_summary.json")
            print(f"   Summary HTML: {output_dir}/PseudoScope_FASTA_QC_summary.html")
        print("\n" + "="*80)
        return all_results


def main():
    """Command line interface for FASTA QC analysis."""
    parser = argparse.ArgumentParser(
        description='PseudoScope FASTA QC - Quality Control for Pseudomonas aeruginosa with fastANI species confirmation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run on all FASTA files (uses ./ref_db next to this script by default)
  python p_qc.py "*.fna"

  # Run on specific pattern
  python p_qc.py "genomes/*.fasta"

  # Run with custom output directory
  python p_qc.py "*.fa" --output my_qc_results

  # Run with a custom reference directory
  python p_qc.py "*.fna" --ref_dir /path/to/ref_db

  # Skip the fastANI species check entirely
  python p_qc.py "*.fna" --skip-species-check

Supported FASTA extensions: .fasta, .fa, .fna, .faa, .fn (plus .gz variants)
Species check uses fastANI against all *.fna files in ref_db/.
For P. aeruginosa confirmation, ensure pseudomonas_aeruginosa_*.fna is present.
        """
    )

    parser.add_argument('pattern', help='File pattern for FASTA files (e.g., "*.fasta", "genomes/*.fna")')
    parser.add_argument('--output', '-o', default='fasta_qc_results',
                        help='Output directory (default: fasta_qc_results)')
    parser.add_argument('--cpus', '-c', type=int, default=None,
                        help='Number of CPU cores to use (default: all available)')
    parser.add_argument('--ref_dir', default=None,
                        help='Directory containing reference genomes (default: ./ref_db next to this script)')
    parser.add_argument('--skip-species-check', action='store_true',
                        help='Skip fastANI species confirmation entirely')

    args = parser.parse_args()

    ref_dir = args.ref_dir
    if args.skip_species_check:
        ref_dir = "/nonexistent/__skip__"

    try:
        qc = PseudoFASTAQC(cpus=args.cpus, ref_dir=ref_dir)
        results = qc.process_files(args.pattern, args.output)
    except Exception as e:
        print(f"FASTA QC analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()