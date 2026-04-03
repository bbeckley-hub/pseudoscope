#!/usr/bin/env python3
"""
PseudoScope - Pseudomonas aeruginosa Ultimate Reporter
Comprehensive Gene-Centric Cross-Genome Analysis with Environmental Markers
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 1.0.0 - P. aeruginosa Edition
Date: 2026-01-07
"""

import os
import sys
import json
import re
import glob
import argparse
import csv
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any, Optional
from datetime import datetime
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

from bs4 import BeautifulSoup


class PseudoHTMLParser:
    """HTML parser for PseudoScope summary reports – adapted for P. aeruginosa."""

    def __init__(self):
        # Database name mapping (from filenames to internal names)
        self.db_name_mapping = {
            'pseudo_amrfinder': 'amrfinder',
            'pseudo_card': 'card',
            'pseudo_resfinder': 'resfinder',
            'pseudo_argannot': 'argannot',
            'pseudo_megares': 'megares',
            'pseudo_bacmet2': 'bacmet2',
            'pseudo_vfdb': 'vfdb',
            'pseudo_victors': 'victors',
            'pseudo_ecoli_vf': 'ecoli_vf',
            'pseudo_plasmidfinder': 'plasmidfinder',
            'pseudo_ncbi': 'ncbi'
        }

    def normalize_sample_id(self, sample_id: str) -> str:
        """Normalize sample ID to a common format (remove extensions and paths)."""
        sample = str(sample_id)
        # Remove path if present
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        # Remove common extensions and suffixes
        extensions = ['.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv',
                      ]
        for ext in extensions:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
                break
        # Also remove any trailing .fna etc. that might remain
        if sample.endswith('.fna') or sample.endswith('.fasta') or sample.endswith('.fa'):
            sample = sample[:sample.rfind('.')]
        return sample.strip()

    def parse_html_table(self, html_content: str, table_index: int = 0) -> List[Dict[str, str]]:
        """
        Parse an HTML table into a list of dicts (rows).
        Returns empty list on failure.
        """
        try:
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if not tables or table_index >= len(tables):
                return []
            table = tables[table_index]
            rows = table.find_all('tr')
            if not rows:
                return []
            
            # Get headers from first row
            headers = []
            header_cells = rows[0].find_all(['th', 'td'])
            for th in header_cells:
                headers.append(th.get_text().strip())
            
            if not headers:
                return []
            
            # Parse data rows
            data = []
            for row in rows[1:]:
                cells = row.find_all(['td', 'th'])
                if cells:
                    row_dict = {}
                    for i, cell in enumerate(cells):
                        if i < len(headers):
                            row_dict[headers[i]] = cell.get_text().strip()
                        else:
                            # If more cells than headers, add to last header
                            last_header = headers[-1] if headers else f"col_{i}"
                            if last_header in row_dict:
                                row_dict[last_header] += " " + cell.get_text().strip()
                            else:
                                row_dict[last_header] = cell.get_text().strip()
                    data.append(row_dict)
            return data
        except Exception as e:
            return []

    def parse_mlst_summary(self, file_path: Path) -> Dict[str, Dict]:
        """Parse mlst_summary.html to get ST for each sample."""
        print(f"  🧬 Parsing MLST summary: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            rows = self.parse_html_table(html, 0)
            if not rows:
                return {}
            
            # Identify columns based on the actual format
            sample_col = None
            st_col = None
            allele_col = None
            
            if rows:
                for key in rows[0].keys():
                    key_lower = key.lower()
                    if 'sample' in key_lower:
                        sample_col = key
                    elif key_lower == 'st' or ('st' in key_lower and 'allele' not in key_lower and 'assigned' not in key_lower):
                        st_col = key
                    elif 'allele' in key_lower or 'profile' in key_lower:
                        allele_col = key
            
            # Fallbacks if columns not found
            if not sample_col and rows:
                keys = list(rows[0].keys())
                if len(keys) >= 2:
                    sample_col = keys[1]  # second column is usually Sample
            
            if not st_col and rows:
                for key in rows[0].keys():
                    if key.lower() in ['st', 'sequence type', 'type']:
                        st_col = key
                        break
                if not st_col and len(keys) >= 3:
                    st_col = keys[2]  # third column is usually ST
            
            results = {}
            for r in rows:
                sample_val = r.get(sample_col, '') if sample_col else ''
                if not sample_val and len(r) >= 2:
                    sample_val = list(r.values())[1]
                
                sample = self.normalize_sample_id(sample_val)
                if not sample:
                    continue
                
                st_val = r.get(st_col, '') if st_col else ''
                if not st_val and len(r) >= 3:
                    st_val = list(r.values())[2]
                
                # Clean ST value (remove "ST" prefix)
                st = st_val
                if st.startswith('ST'):
                    st = st[2:]
                if st == '' or st == '-' or st.lower() == 'unknown':
                    st = 'ND'
                
                allele_profile = r.get(allele_col, '') if allele_col else ''
                
                results[sample] = {
                    'ST': st,
                    'Allele_Profile': allele_profile,
                    'Scheme': 'pasteur'
                }
            return results
        except Exception as e:
            print(f"    ❌ Error parsing MLST summary: {e}")
            return {}

    def parse_past_summary(self, file_path: Path) -> Dict[str, Dict]:
        """Parse past_summary.html to get serotype for each sample."""
        print(f"  🧬 Parsing PAST summary: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            rows = self.parse_html_table(html, 0)
            if not rows:
                return {}
            
            # Identify columns
            sample_col = None
            type_col = None
            targets_col = None
            coverage_col = None
            hits_col = None
            comment_col = None
            
            if rows:
                for key in rows[0].keys():
                    k = key.lower()
                    if 'sample' in k:
                        sample_col = key
                    elif 'type' in k:
                        type_col = key
                    elif 'target' in k:
                        targets_col = key
                    elif 'cover' in k:
                        coverage_col = key
                    elif 'hit' in k:
                        hits_col = key
                    elif 'comment' in k:
                        comment_col = key
            
            if not sample_col and rows:
                keys = list(rows[0].keys())
                if len(keys) >= 2:
                    sample_col = keys[1]
            
            if not type_col and rows:
                for key in rows[0].keys():
                    if 'type' in key.lower():
                        type_col = key
                        break
                if not type_col and len(keys) >= 3:
                    type_col = keys[2]
            
            results = {}
            for r in rows:
                sample_val = r.get(sample_col, '') if sample_col else ''
                if not sample_val and len(r) >= 2:
                    sample_val = list(r.values())[1]
                
                sample = self.normalize_sample_id(sample_val)
                if not sample:
                    continue
                
                serotype = r.get(type_col, 'UNKNOWN') if type_col else 'UNKNOWN'
                if serotype == '' or serotype == '-':
                    serotype = 'UNKNOWN'
                
                results[sample] = {
                    'Serotype': serotype,
                    'Targets': r.get(targets_col, '') if targets_col else '',
                    'Coverage': r.get(coverage_col, '') if coverage_col else '',
                    'Hits': r.get(hits_col, '') if hits_col else '',
                    'Comment': r.get(comment_col, '') if comment_col else ''
                }
            return results
        except Exception as e:
            print(f"    ❌ Error parsing PAST summary: {e}")
            return {}

    def parse_qc_summary(self, file_path: Path) -> Dict[str, Dict]:
        """Parse PseudoScope_FASTA_QC_summary.html to get assembly metrics."""
        print(f"  🧬 Parsing QC summary: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            rows = self.parse_html_table(html, 0)
            if not rows:
                return {}
            
            # Identify columns based on the actual format
            filename_col = None
            seq_col = None
            length_col = None
            gc_col = None
            n50_col = None
            warnings_col = None
            
            if rows:
                for key in rows[0].keys():
                    k = key.lower()
                    if 'filename' in k:
                        filename_col = key
                    elif 'total sequences' in k:
                        seq_col = key
                    elif 'total length' in k:
                        length_col = key
                    elif 'gc content' in k:
                        gc_col = key
                    elif 'n50' in k:
                        n50_col = key
                    elif 'warnings' in k:
                        warnings_col = key
            
            if not filename_col and rows:
                keys = list(rows[0].keys())
                if keys:
                    filename_col = keys[0]  # first column is usually Filename
            
            results = {}
            for r in rows:
                sample_val = r.get(filename_col, '') if filename_col else ''
                if not sample_val and len(r) >= 1:
                    sample_val = list(r.values())[0]
                
                sample = self.normalize_sample_id(sample_val)
                if not sample:
                    continue
                
                results[sample] = {
                    'Total_Sequences': r.get(seq_col, '') if seq_col else '',
                    'Total_Length': r.get(length_col, '') if length_col else '',
                    'GC_Content': r.get(gc_col, '') if gc_col else '',
                    'N50': r.get(n50_col, '') if n50_col else '',
                    'Warnings': r.get(warnings_col, '') if warnings_col else ''
                }
            return results
        except Exception as e:
            print(f"    ❌ Error parsing QC summary: {e}")
            return {}

    def parse_gene_frequency_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        """
        Parse a gene frequency report (AMRfinder or ABRicate summary).
        Returns (genes_by_genome, gene_frequencies).
        """
        print(f"  🧬 Parsing gene report: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            soup = BeautifulSoup(html, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return {}, {}

            # Determine database name from filename
            fname_low = file_path.name.lower()
            db_name = 'unknown'
            for key, val in self.db_name_mapping.items():
                if key in fname_low:
                    db_name = val
                    break
            if 'amrfinder' in fname_low:
                db_name = 'amrfinder'

            # 1. Genes by Genome table (first table)
            genes_by_genome = {}
            rows1 = self.parse_html_table(str(tables[0]), 0)
            if rows1:
                # Find sample column
                sample_col = None
                for key in rows1[0].keys():
                    if 'genome' in key.lower() or 'sample' in key.lower():
                        sample_col = key
                        break
                if not sample_col and rows1:
                    keys = list(rows1[0].keys())
                    if keys:
                        sample_col = keys[0]  # first column is usually Genome
                
                # Find genes column
                genes_col = None
                if sample_col:
                    for key in rows1[0].keys():
                        if 'genes' in key.lower() or 'detected' in key.lower():
                            genes_col = key
                            break
                    
                    if genes_col:
                        for r in rows1:
                            sample_val = r.get(sample_col, '')
                            if not sample_val:
                                continue
                            sample = self.normalize_sample_id(sample_val)
                            if not sample:
                                continue
                            
                            gene_str = r.get(genes_col, '')
                            # Parse genes - they can be space-separated or comma-separated
                            if ',' in gene_str:
                                genes = [g.strip() for g in gene_str.split(',') if g.strip()]
                            else:
                                genes = [g.strip() for g in gene_str.split() if g.strip()]
                            
                            # Clean genes (remove any HTML artifacts)
                            genes = [re.sub(r'<[^>]+>', '', g) for g in genes]
                            genes_by_genome[sample] = genes

            # 2. Gene Frequency table (second table)
            gene_freq = {}
            rows2 = self.parse_html_table(str(tables[1]), 0)
            if rows2:
                # Find column names
                gene_col = None
                freq_col = None
                genomes_col = None
                risk_col = None
                
                if rows2:
                    for key in rows2[0].keys():
                        k = key.lower()
                        if 'gene' in k:
                            gene_col = key
                        elif 'frequency' in k:
                            freq_col = key
                        elif 'genomes' in k:
                            genomes_col = key
                        elif 'risk' in k or 'level' in k:
                            risk_col = key
                
                if gene_col:
                    for r in rows2:
                        gene_full = r.get(gene_col, '').strip()
                        if not gene_full:
                            continue
                        
                        # Clean gene name
                        gene = re.sub(r'^\([^)]+\)', '', gene_full).strip()
                        if not gene:
                            gene = gene_full
                        
                        # Get count from frequency string
                        count = 0
                        freq_str = r.get(freq_col, '0') if freq_col else '0'
                        match = re.search(r'(\d+)', freq_str)
                        if match:
                            count = int(match.group(1))
                        
                        # Get percentage from frequency string
                        percentage = 0
                        pct_match = re.search(r'\((\d+\.?\d*)%\)', freq_str)
                        if pct_match:
                            percentage = float(pct_match.group(1))
                        elif total_samples > 0:
                            percentage = (count / total_samples * 100)
                        
                        # Get genomes list
                        genomes = []
                        if genomes_col and r.get(genomes_col):
                            genomes_str = r.get(genomes_col, '')
                            # Split by comma
                            for g in genomes_str.split(','):
                                g_clean = g.strip()
                                if g_clean:
                                    genomes.append(self.normalize_sample_id(g_clean))
                        
                        # Risk level (for AMRfinder)
                        risk_level = 'Standard'
                        if risk_col and r.get(risk_col):
                            risk_level = r.get(risk_col, 'Standard')
                        
                        gene_freq[gene] = {
                            'count': count,
                            'percentage': round(percentage, 2),
                            'frequency_display': f"{count} ({percentage:.1f}%)",
                            'genomes': genomes,
                            'database': db_name,
                            'full_name': gene_full,
                            'risk_level': risk_level
                        }

            return genes_by_genome, gene_freq
        except Exception as e:
            print(f"    ❌ Error parsing {file_path.name}: {e}")
            return {}, {}


class PseudoDataAnalyzer:
    """Analyzes data for gene-centric reporting for P. aeruginosa."""

    def __init__(self):
        # ===== CRITICAL CLINICAL RESISTANCE GENES FOR P. AERUGINOSA =====
        self.critical_carbapenemases = {
            'blaOXA-2', 'blaOXA-10', 'blaOXA-14', 'blaOXA-17', 'blaOXA-19', 'blaOXA-28',
            'blaOXA-35', 'blaOXA-45', 'blaOXA-50', 'blaOXA-198', 'blaOXA-395', 'blaOXA-396',
            'blaOXA-847', 'blaOXA-904',
            'blaIMP', 'blaVIM', 'blaVIM-2', 'blaNDM', 'blaKPC', 'blaKPC-3', 'blaGES', 'blaSPM', 'blaAIM', 'blaDIM'
        }
        self.critical_esbls = {
            'blaPER', 'blaVEB', 'blaBEL', 'blaGES', 'blaTEM', 'blaSHV', 'blaCTX-M'
        }
        self.critical_ampc = {'ampC', 'blaAmpC', 'blaPDC-1', 'blaPDC-3', 'blaPDC-5', 'blaPDC-374', 'ampR', 'dacB'}
        self.critical_colistin = {
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10',
            'pmrA', 'pmrB', 'pmrB_V15I', 'phoP', 'phoQ', 'mgrB', 'lpxA', 'lpxC', 'lpxD', 'arnT', 'eptA', 'eptB', 'basS', 'basR'
        }
        self.critical_aminoglycoside = {
            'armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD', 'rmtE', 'rmtF', 'rmtG', 'rmtH', 'npmA',
            'aac(3)', 'aac(6\')', 'aac(6\')-29', 'aac(6\')-29a', 'aac(6\')-29b', 'aac(6\')-33', "aac(6')-Ib'",
            'ant(2")', 'ant(4\')', 'aph(3\')', 'aph(3\')-IIb', 'aph(6)',
            'aacC1', 'aacC2', 'aacC4', 'aacA4', 'aacA7', 'aadA1', 'aadA2', 'aadA5', 'aadA7',
            'strA', 'strB', 'aphA1', 'aphA2', 'aphA3', 'aphA6'
        }
        self.critical_efflux = {
            'mexA', 'mexB', 'mexC', 'mexD', 'mexE', 'mexF', 'mexX', 'mexY',
            'oprM', 'oprN', 'oprJ', 'oprD_V359L', 'adeA', 'adeB', 'adeC'
        }
        self.critical_virulence = {
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
            'lasI', 'lasR', 'rhlI', 'rhlR', 'pqsA', 'pqsB', 'pqsC', 'pqsD', 'pqsE', 'pqsH', 'mvfR'
        }

        # ENVIRONMENTAL & CO‑SELECTION MARKERS
        self.bacmet2_markers = {
            # Biocide resistance
            'qacA', 'qacB', 'qacC', 'qacD', 'qacE', 'qacF', 'qacG', 'qacH', 'qacI', 'qacJ',
            'qacEA1', 'qacEdelta1', 'qacG2', 'qacH2',
            'cepA', 'formA', 'formB', 'formC', 'oqxA', 'oqxB',
            # Heavy metal resistance
            'czcA', 'czcB', 'czcC', 'czcD', 'czcR', 'czcS',
            'merA', 'merB', 'merC', 'merD', 'merE', 'merF', 'merP', 'merR', 'merT',
            'arsA', 'arsB', 'arsC', 'arsD', 'arsE', 'arsF', 'arsG', 'arsH', 'arsI', 'arsJ', 'arsT',
            'copA', 'copB', 'copC', 'copD',
            'zntA', 'zntB', 'zntC', 'zntD',
            'chrA', 'chrB', 'chrC', 'chrD',
            'nikA', 'nikB', 'nikC', 'nikD', 'nikE', 'nikR',
            'cadA', 'cadB', 'cadC', 'cadD',
            'silA', 'silB', 'silC', 'silD', 'silE',
            'pbrA', 'pbrB', 'pbrC', 'pbrD', 'pbrR'
        }
        self.environmental_co_selection = {
            'soxR', 'soxS', 'marA', 'marB', 'marC', 'marR', 'robA',
            'rpoS', 'rpoH', 'oxyRkp', 'cpxR', 'baeR', 'emrAsm', 'emrBsm', 'yddg/emrE', 'lmrS', 'smeF', 'sugE',
            'phoB', 'phoR', 'phoU',
            'traA', 'traB', 'traC', 'traD', 'traE', 'traF', 'traG', 'traH', 'traI', 'traJ',
            'traK', 'traL', 'traM', 'traN', 'traO', 'traP', 'traQ', 'traR', 'traS', 'traT',
            'traU', 'traV', 'traW', 'traX', 'traY',
            'mobA', 'mobB', 'mobC', 'mobD', 'mobE', 'mobF', 'mobG', 'mobH',
            'oriT', 'oriV', 'repA', 'repB', 'repC',
            'intI1', 'intI2', 'intI3',
            'tnpA', 'tnpB', 'tnpC', 'tnpD', 'tnpE', 'tnpF',
            'istA', 'istB'
        }
        self.common_antibiotic_markers = {
            'sul1', 'sul2', 'sul3', 'sul4',
            'dfrA1', 'dfrA5', 'dfrA7', 'dfrA12', 'dfrA14', 'dfrA17', 'dfrA19', 'dfrA21',
            'dfrB1', 'dfrB2', 'dfrB3', 'dfrB4',
            'catA1', 'catA2', 'catB2', 'catB3', 'catB7', 'catB8', 'catI', 'catII', 'catIII',
            'aac', 'aad', 'ant', 'aph',
            'tet', 'tetR', 'tetA', 'tetB', 'tetC', 'tetD', 'tetE', 'tetG', 'tetH',
            'tetJ', 'tetK', 'tetL', 'tetM', 'tetO', 'tetQ', 'tetS', 'tetW', 'tetX',
            'erm', 'ere', 'mef', 'msr',
            'mdt', 'emr', 'acr', 'tolC',
            'blaTEM', 'blaSHV', 'blaCTX-M', 'blaOXA', 'blaCARB-2'
        }
        self.victors_markers = {
            'fimA', 'fimB', 'fimC', 'fimD', 'fimE', 'fimF', 'fimG', 'fimH',
            'papA', 'papB', 'papC', 'papD', 'papE', 'papF', 'papG', 'papH', 'papI', 'papJ', 'papK',
            'afa', 'dra',
            'hlyA', 'hlyB', 'hlyC', 'hlyD',
            'cnf1', 'cnf2',
            'sat', 'astA', 'stx1', 'stx2',
            'iss', 'traT', 'kpsM', 'kpsT',
            'ibeA', 'ibeB',
            'iutA', 'iroN', 'fyuA', 'irp1', 'irp2',
            'chu', 'usp', 'vat', 'pic', 'sigA', 'tia'
        }

        # ASCII Art for reports
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
            {"text": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie"},
            {"text": "Research is what I'm doing when I don't know what I'm doing.", "author": "Wernher von Braun"},
            {"text": "The science of today is the technology of tomorrow.", "author": "Edward Teller"},
        ]

    def categorize_gene(self, gene: str) -> str:
        """Categorize a gene into a functional class."""
        gene_lower = gene.lower()
        
        # Critical clinical
        if any(carb in gene_lower for carb in [g.lower() for g in self.critical_carbapenemases]):
            return 'Carbapenemase'
        if any(esbl in gene_lower for esbl in [g.lower() for g in self.critical_esbls]):
            return 'ESBL'
        if any(ampc in gene_lower for ampc in [g.lower() for g in self.critical_ampc]):
            return 'AmpC'
        if any(col in gene_lower for col in [g.lower() for g in self.critical_colistin]):
            return 'Colistin Resistance'
        if any(ag in gene_lower for ag in [g.lower() for g in self.critical_aminoglycoside]):
            return 'Aminoglycoside Resistance'
        if any(eff in gene_lower for eff in [g.lower() for g in self.critical_efflux]):
            return 'Efflux Pump'
        if any(vir in gene_lower for vir in [g.lower() for g in self.critical_virulence]):
            return 'Virulence Factor'
        
        # Environmental
        if any(bac in gene_lower for bac in [g.lower() for g in self.bacmet2_markers]):
            if any(bio in gene_lower for bio in ['qac', 'cep', 'form', 'oqx']):
                return 'Biocide Resistance (BACMET2)'
            elif any(metal in gene_lower for metal in ['czc', 'mer', 'ars', 'cop', 'znt', 'chr', 'nik', 'cad', 'sil', 'pbr']):
                return 'Metal Resistance (BACMET2)'
            else:
                return 'BACMET2 Other'
        
        if any(env in gene_lower for env in [g.lower() for g in self.environmental_co_selection]):
            if any(mob in gene_lower for mob in ['tra', 'mob', 'rep', 'ori']):
                return 'Plasmid Transfer'
            elif any(stress in gene_lower for stress in ['sox', 'mar', 'rob', 'rpo']):
                return 'Stress Response'
            elif any(mobile in gene_lower for mobile in ['int', 'tnp', 'ist']):
                return 'Mobile Genetic Elements'
            else:
                return 'Environmental Co-Selection'
        
        if any(env_abx in gene_lower for env_abx in [g.lower() for g in self.common_antibiotic_markers]):
            if 'sul' in gene_lower:
                return 'Sulfonamide Resistance'
            elif 'dfr' in gene_lower:
                return 'Trimethoprim Resistance'
            elif 'cat' in gene_lower:
                return 'Chloramphenicol Resistance'
            elif any(ag in gene_lower for ag in ['aac', 'aad', 'ant', 'aph']):
                return 'Aminoglycoside Resistance'
            elif 'tet' in gene_lower and 'tet(x)' not in gene_lower:
                return 'Tetracycline Resistance'
            elif any(ml in gene_lower for ml in ['erm', 'mef', 'msr']):
                return 'Macrolide Resistance'
            elif 'bla' in gene_lower:
                return 'Beta-lactamase'
            else:
                return 'Antibiotic Resistance'
        
        if any(vic in gene_lower for vic in [g.lower() for g in self.victors_markers]):
            if any(adh in gene_lower for adh in ['fim', 'pap', 'afa', 'dra']):
                return 'VICTORS Adhesion'
            elif any(tox in gene_lower for tox in ['hly', 'cnf', 'sat', 'ast', 'stx']):
                return 'VICTORS Toxins'
            elif any(iron in gene_lower for iron in ['iut', 'iro', 'fyu', 'irp', 'chu']):
                return 'VICTORS Iron Acquisition'
            else:
                return 'VICTORS Virulence'
        
        # Default
        if any(abx in gene_lower for abx in ['sul', 'dfr', 'cat', 'aac', 'aad', 'ant', 'aph', 'tet', 'erm']):
            return 'Antibiotic Resistance'
        elif any(vir in gene_lower for vir in ['tox', 'hly', 'cnf', 'iss', 'fim', 'pap']):
            return 'Virulence Factors'
        else:
            return 'Other'

    def create_gene_centric_tables(self, integrated_data: Dict[str, Any], total_samples: int) -> Dict[str, Any]:
        """Build gene-centric tables from parsed data."""
        gene_centric = {
            'amr_databases': {},
            'virulence_databases': {},
            'plasmid_databases': {},
            'environmental_databases': {},
            'combined_gene_frequencies': [],
            'gene_categories': defaultdict(list),
            'database_stats': {},
            'environmental_summary': defaultdict(dict)
        }

        # Process each database that has gene frequencies
        all_gene_freqs = integrated_data.get('gene_frequencies', {})
        
        # Handle amrfinder separately
        if 'amrfinder' in all_gene_freqs:
            db_name = 'amrfinder'
            freq_dict = all_gene_freqs['amrfinder']
            if freq_dict:
                gene_list = []
                for gene, data in freq_dict.items():
                    category = self.categorize_gene(gene)
                    gene_list.append({
                        'gene': gene,
                        'category': category,
                        'database': db_name,
                        'count': data['count'],
                        'percentage': data['percentage'],
                        'frequency_display': data['frequency_display'],
                        'risk_level': data.get('risk_level', 'Standard'),
                        'genomes': data['genomes'],
                        'full_name': data.get('full_name', gene)
                    })
                if gene_list:
                    gene_list.sort(key=lambda x: x['count'], reverse=True)
                    gene_centric['amr_databases'][db_name] = gene_list
                    for gd in gene_list:
                        gene_centric['gene_categories'][gd['category']].append(gd)

        # Handle abricate databases
        if 'abricate' in all_gene_freqs:
            for sub_db, sub_freq in all_gene_freqs['abricate'].items():
                if not sub_freq:
                    continue
                gene_list = []
                for gene, data in sub_freq.items():
                    category = self.categorize_gene(gene)
                    gene_list.append({
                        'gene': gene,
                        'category': category,
                        'database': sub_db,
                        'count': data['count'],
                        'percentage': data['percentage'],
                        'frequency_display': data['frequency_display'],
                        'risk_level': data.get('risk_level', 'Standard'),
                        'genomes': data['genomes'],
                        'full_name': data.get('full_name', gene)
                    })
                if not gene_list:
                    continue
                gene_list.sort(key=lambda x: x['count'], reverse=True)

                # Categorize database type
                if sub_db in ['vfdb', 'victors', 'ecoli_vf']:
                    gene_centric['virulence_databases'][sub_db] = gene_list
                elif sub_db == 'plasmidfinder':
                    gene_centric['plasmid_databases'][sub_db] = gene_list
                elif sub_db == 'bacmet2':
                    gene_centric['environmental_databases'][sub_db] = gene_list
                else:
                    gene_centric['amr_databases'][sub_db] = gene_list

                for gd in gene_list:
                    gene_centric['gene_categories'][gd['category']].append(gd)

        # Build environmental summary by category
        self._build_environmental_summary(gene_centric, total_samples)

        # Combined gene frequencies
        combined = []
        for db_type in ['amr_databases', 'virulence_databases', 'plasmid_databases', 'environmental_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                combined.extend(genes)
        combined.sort(key=lambda x: x['count'], reverse=True)
        gene_centric['combined_gene_frequencies'] = combined

        # Sort categories
        for cat in gene_centric['gene_categories']:
            gene_centric['gene_categories'][cat].sort(key=lambda x: x['count'], reverse=True)

        # Database stats
        for db_type in ['amr_databases', 'virulence_databases', 'plasmid_databases', 'environmental_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                total_genes = len(genes)
                total_occ = sum(g['count'] for g in genes)
                critical_count = sum(1 for g in genes if g['category'] in
                                     ['Carbapenemase', 'ESBL', 'AmpC', 'Colistin Resistance',
                                      'Aminoglycoside Resistance', 'Efflux Pump'])
                bacmet2_count = sum(1 for g in genes if 'BACMET2' in g['category'])
                env_abx_count = sum(1 for g in genes if g['category'] in
                                     ['Sulfonamide Resistance', 'Trimethoprim Resistance',
                                      'Chloramphenicol Resistance', 'Tetracycline Resistance',
                                      'Macrolide Resistance', 'Antibiotic Resistance'])
                gene_centric['database_stats'][db_name] = {
                    'total_genes': total_genes,
                    'total_occurrences': total_occ,
                    'critical_genes': critical_count,
                    'bacmet2_genes': bacmet2_count,
                    'environmental_antibiotic_genes': env_abx_count,
                    'coverage': f"{(total_occ / (total_samples * max(total_genes, 1))) * 100:.1f}%" if total_genes else "0%"
                }
        return gene_centric

    def _build_environmental_summary(self, gene_centric: Dict, total_samples: int):
        """Group environmental markers into summary categories."""
        env_categories = {
            'BACMET2 - Biocide Resistance': [],
            'BACMET2 - Metal Resistance': [],
            'Environmental Co-Selection': [],
            'Environmental Antibiotic Resistance': [],
            'Mobile Genetic Elements': [],
            'VICTORS Virulence Factors': []
        }
        for db_type in ['amr_databases', 'virulence_databases', 'environmental_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for g in genes:
                    cat = g['category']
                    if cat == 'Biocide Resistance (BACMET2)':
                        env_categories['BACMET2 - Biocide Resistance'].append(g)
                    elif cat == 'Metal Resistance (BACMET2)':
                        env_categories['BACMET2 - Metal Resistance'].append(g)
                    elif cat in ['Environmental Co-Selection', 'Stress Response']:
                        env_categories['Environmental Co-Selection'].append(g)
                    elif cat in ['Plasmid Transfer', 'Mobile Genetic Elements']:
                        env_categories['Mobile Genetic Elements'].append(g)
                    elif cat in ['Sulfonamide Resistance', 'Trimethoprim Resistance',
                                 'Chloramphenicol Resistance', 'Tetracycline Resistance',
                                 'Macrolide Resistance', 'Antibiotic Resistance']:
                        env_categories['Environmental Antibiotic Resistance'].append(g)
                    elif 'VICTORS' in cat:
                        env_categories['VICTORS Virulence Factors'].append(g)
        for cat_name, genes in env_categories.items():
            if genes:
                genes.sort(key=lambda x: x['count'], reverse=True)
                gene_centric['environmental_summary'][cat_name] = {
                    'total_genes': len(genes),
                    'total_occurrences': sum(g['count'] for g in genes),
                    'genes': genes[:10000],  # keep all
                    'percentage_of_samples': (sum(g['count'] for g in genes) / (total_samples * len(genes))) * 100 if genes else 0
                }

    def create_cross_genome_patterns(self, integrated_data: Dict[str, Any], total_samples: int) -> Dict[str, Any]:
        """Identify patterns like high-risk combinations, ST distributions, etc."""
        patterns = {
            'st_distribution': Counter(),
            'serotype_distribution': Counter(),
            'gene_cooccurrence': {},
            'high_risk_combinations': [],
            'mdr_patterns': [],
            'carbapenemase_patterns': {},
            'environmental_patterns': {},
            'database_coverage': {},
            'co_selection_risk': []  # new: co-selection risk per sample
        }
        samples_data = integrated_data.get('samples', {})
        gene_centric = integrated_data.get('gene_centric', {})

        # Collect genes per sample
        sample_genes = defaultdict(set)
        for db_type in ['amr_databases', 'virulence_databases', 'plasmid_databases', 'environmental_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                for gene_data in genes:
                    for genome in gene_data['genomes']:
                        sample_genes[genome].add(gene_data['gene'])

        # Analyze each sample
        for sample, data in samples_data.items():
            st = data.get('mlst', {}).get('ST', 'ND')
            serotype = data.get('past', {}).get('Serotype', 'ND')
            if st != 'ND':
                patterns['st_distribution'][st] += 1
            if serotype != 'ND':
                patterns['serotype_distribution'][serotype] += 1

            genes = list(sample_genes.get(sample, set()))
            # Co-occurrence
            for i, g1 in enumerate(genes):
                if g1 not in patterns['gene_cooccurrence']:
                    patterns['gene_cooccurrence'][g1] = Counter()
                for g2 in genes[i+1:]:
                    patterns['gene_cooccurrence'][g1][g2] += 1

            # Categorize resistance types
            carb = [g for g in genes if self.categorize_gene(g) == 'Carbapenemase']
            esbl = [g for g in genes if self.categorize_gene(g) == 'ESBL']
            col = [g for g in genes if self.categorize_gene(g) == 'Colistin Resistance']
            amino = [g for g in genes if self.categorize_gene(g) == 'Aminoglycoside Resistance']
            efflux = [g for g in genes if self.categorize_gene(g) == 'Efflux Pump']
            env = [g for g in genes if self.categorize_gene(g) in 
                   ['Biocide Resistance (BACMET2)', 'Metal Resistance (BACMET2)', 'BACMET2 Other',
                    'Environmental Co-Selection', 'Stress Response', 'Plasmid Transfer',
                    'Mobile Genetic Elements', 'Sulfonamide Resistance', 'Trimethoprim Resistance',
                    'Chloramphenicol Resistance', 'Tetracycline Resistance', 'Macrolide Resistance',
                    'Antibiotic Resistance', 'VICTORS Adhesion', 'VICTORS Toxins',
                    'VICTORS Iron Acquisition', 'VICTORS Virulence']]

            if carb:
                key = '|'.join(sorted(set(carb)))  # Use string key instead of tuple
                if key not in patterns['carbapenemase_patterns']:
                    patterns['carbapenemase_patterns'][key] = []
                patterns['carbapenemase_patterns'][key].append(sample)
            
            if env:
                env_key = '|'.join(sorted(set(env)))  # Use string key instead of tuple
                if env_key not in patterns['environmental_patterns']:
                    patterns['environmental_patterns'][env_key] = []
                patterns['environmental_patterns'][env_key].append(sample)

            # High-risk: carbapenemase + colistin resistance
            if carb and col:
                patterns['high_risk_combinations'].append({
                    'sample': sample,
                    'st': st,
                    'serotype': serotype,
                    'carbapenemases': carb,
                    'colistin_resistance': col
                })

            # MDR: 3+ classes among carb, esbl, col, amino, efflux
            classes_present = sum([bool(carb), bool(esbl), bool(col), bool(amino), bool(efflux)])
            if classes_present >= 3:
                patterns['mdr_patterns'].append({
                    'sample': sample,
                    'st': st,
                    'serotype': serotype,
                    'classes': classes_present,
                    'carbapenemases': carb,
                    'esbls': esbl,
                    'colistin_resistance': col,
                    'aminoglycoside_resistance': amino,
                    'efflux_pumps': efflux,
                    'environmental_markers': env
                })

            # Co-selection risk: count of biocide/metal markers + presence of AMR genes
            biocide_metal = [g for g in genes if 'Biocide Resistance' in self.categorize_gene(g) or 'Metal Resistance' in self.categorize_gene(g)]
            amr_count = len(carb) + len(esbl) + len(amino) + len(efflux)
            # Define risk categories
            if len(biocide_metal) >= 3 and amr_count >= 5:
                risk = 'High'
            elif len(biocide_metal) >= 2 and amr_count >= 3:
                risk = 'Medium'
            else:
                risk = 'Low'
            patterns['co_selection_risk'].append({
                'sample': sample,
                'st': st,
                'serotype': serotype,
                'biocide_metal_markers': len(biocide_metal),
                'biocide_metal_genes': biocide_metal,
                'amr_genes_count': amr_count,
                'risk_category': risk
            })

        # Convert Counters to dict for serialization
        patterns['st_distribution'] = dict(patterns['st_distribution'])
        patterns['serotype_distribution'] = dict(patterns['serotype_distribution'])
        
        # Convert gene_cooccurrence Counters to dict
        new_gene_cooccurrence = {}
        for g1, counters in patterns['gene_cooccurrence'].items():
            new_gene_cooccurrence[g1] = dict(counters)
        patterns['gene_cooccurrence'] = new_gene_cooccurrence

        # Database coverage
        all_samples_set = set(samples_data.keys())
        for db_type in ['amr_databases', 'virulence_databases', 'plasmid_databases', 'environmental_databases']:
            for db_name, genes in gene_centric.get(db_type, {}).items():
                samples_with_hits = set()
                for g in genes:
                    samples_with_hits.update(g['genomes'])
                coverage = len(samples_with_hits) / len(all_samples_set) * 100 if all_samples_set else 0
                patterns['database_coverage'][db_name] = {
                    'samples_with_hits': len(samples_with_hits),
                    'total_samples': len(all_samples_set),
                    'coverage_percentage': round(coverage, 2),
                    'coverage_display': f"{len(samples_with_hits)} ({coverage:.1f}%)"
                }
        return patterns


class PseudoHTMLGenerator:
    """Generates comprehensive HTML reports for P. aeruginosa with template style."""
    
    def __init__(self, analyzer: PseudoDataAnalyzer):
        self.analyzer = analyzer
        self.tab_colors = {
            'summary': "#4CAF50",
            'samples': '#2196F3',
            'mlst': '#FF9800',
            'past': '#9C27B0',
            'amr': '#F44336',
            'virulence': '#E91E63',
            'environmental': '#795548',
            'categories': '#009688',
            'patterns': '#FF5722',
            'export': '#3F51B5',
            'databases': '#607D8B',
            'qc': '#17a2b8'
        }
    
    def generate_main_report(self, integrated_data: Dict[str, Any], output_dir: Path) -> str:
        """Generate the ultimate HTML report."""
        print("\n🎨 Generating PseudoScope Ultimate HTML report...")
        
        samples_data = integrated_data.get('samples', {})
        patterns = integrated_data.get('patterns', {})
        gene_centric = integrated_data.get('gene_centric', {})
        metadata = integrated_data.get('metadata', {})
        qc_data = integrated_data.get('qc_data', {})
        total_samples = len(samples_data)
        
        # Compute additional patterns for novelty
        st_serotype_patterns = self._compute_st_serotype_patterns(samples_data)
        carb_patterns = patterns.get('carbapenemase_patterns', {})
        co_risk = patterns.get('co_selection_risk', [])
        
        html = self._create_ultimate_html(
            metadata=metadata,
            samples_data=samples_data,
            patterns=patterns,
            gene_centric=gene_centric,
            qc_data=qc_data,
            total_samples=total_samples,
            st_serotype_patterns=st_serotype_patterns,
            carb_patterns=carb_patterns,
            co_risk=co_risk
        )
        
        output_file = output_dir / "pseudoscope_ultimate_report.html"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        
        print(f"    ✅ HTML report saved: {output_file}")
        return str(output_file)
    
    def _compute_st_serotype_patterns(self, samples_data: Dict) -> Dict[str, List[str]]:
        """Compute ST-Serotype combinations."""
        patterns = defaultdict(list)
        for sample, data in samples_data.items():
            st = data.get('mlst', {}).get('ST', 'ND')
            sero = data.get('past', {}).get('Serotype', 'ND')
            if st != 'ND' and sero != 'ND' and sero != 'UNKNOWN':
                key = f"ST{st}:{sero}"
                patterns[key].append(sample)
        return dict(sorted(patterns.items(), key=lambda x: len(x[1]), reverse=True))
    
    def _create_ultimate_html(self, **kwargs) -> str:
        metadata = kwargs['metadata']
        samples_data = kwargs['samples_data']
        patterns = kwargs['patterns']
        gene_centric = kwargs['gene_centric']
        qc_data = kwargs['qc_data']
        total_samples = kwargs['total_samples']
        st_serotype_patterns = kwargs['st_serotype_patterns']
        carb_patterns = kwargs['carb_patterns']
        co_risk = kwargs['co_risk']
        
        # Pre‑compute timestamp to avoid f‑string issues
        timestamp = metadata.get('analysis_date', datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        
        # Calculate statistics
        total_amr_genes = sum(len(db) for db in gene_centric.get('amr_databases', {}).values())
        total_virulence_genes = sum(len(db) for db in gene_centric.get('virulence_databases', {}).values())
        total_environmental_genes = sum(len(db) for db in gene_centric.get('environmental_databases', {}).values())
        high_risk_count = len(patterns.get('high_risk_combinations', []))
        mdr_count = len(patterns.get('mdr_patterns', []))
        
        # Count carbapenemases
        carbapenemase_count = 0
        for db_name, genes in gene_centric.get('amr_databases', {}).items():
            for gene_data in genes:
                if gene_data['category'] == 'Carbapenemase':
                    carbapenemase_count += 1
        
        # Count environmental markers
        environmental_marker_count = len(gene_centric.get('environmental_summary', {}))
        
        # Get a random quote
        import random
        quote = random.choice(self.analyzer.science_quotes)
        
        # Build the HTML
        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PseudoScope Ultimate Report - Pseudomonas aeruginosa</title>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            background: linear-gradient(135deg, #2c0b0e 0%, #4a1c20 50%, #6b2b2f 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #fff;
            padding: 20px;
            min-height: 100vh;
        }}
        
        .container {{
            max-width: 1800px;
            margin: 0 auto;
        }}
        
        .header {{
            text-align: center;
            margin-bottom: 30px;
        }}
        
        .ascii-container {{
            background: rgba(0,0,0,0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            border: 2px solid #dc3545;
            overflow-x: auto;
            text-align: center;
        }}
        
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #dc3545;
            text-shadow: 0 0 10px rgba(220,53,69,0.5);
            display: inline-block;
            text-align: left;
        }}
        
        .quote-container {{
            background: rgba(255,255,255,0.1);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            border-left: 4px solid #dc3545;
        }}
        
        .quote-text {{
            font-size: 18px;
            font-style: italic;
            margin-bottom: 5px;
        }}
        
        .quote-author {{
            font-size: 14px;
            color: #fbbf24;
        }}
        
        .dashboard-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-bottom: 25px;
        }}
        
        .dashboard-card {{
            background: rgba(255,255,255,0.1);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            text-align: center;
            border: 1px solid rgba(255,255,255,0.2);
            transition: all 0.3s ease;
            cursor: pointer;
        }}
        
        .dashboard-card:hover {{
            transform: translateY(-5px);
            background: rgba(220,53,69,0.2);
            border-color: #dc3545;
        }}
        
        .card-value {{
            font-size: 36px;
            font-weight: bold;
            color: #dc3545;
            margin: 10px 0;
        }}
        
        .card-label {{
            font-size: 14px;
            opacity: 0.9;
        }}
        
        .tab-navigation {{
            display: flex;
            gap: 5px;
            margin-bottom: 20px;
            flex-wrap: wrap;
            background: rgba(255,255,255,0.1);
            padding: 15px;
            border-radius: 10px;
            backdrop-filter: blur(10px);
        }}
        
        .tab-button {{
            padding: 10px 20px;
            background: transparent;
            border: 1px solid rgba(255,255,255,0.2);
            border-radius: 5px;
            color: white;
            cursor: pointer;
            font-size: 14px;
            transition: all 0.3s ease;
            display: flex;
            align-items: center;
            gap: 8px;
        }}
        
        .tab-button:hover {{
            background: rgba(220,53,69,0.3);
            border-color: #dc3545;
        }}
        
        .tab-button.active {{
            background: #dc3545;
            border-color: #dc3545;
        }}
        
        .tab-content {{
            display: none;
            background: rgba(255,255,255,0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
        }}
        
        .tab-content.active {{
            display: block;
        }}
        
        .section-header {{
            color: #a71d2a;
            border-bottom: 3px solid #dc3545;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}
        
        .search-box {{
            width: 100%;
            padding: 10px;
            margin-bottom: 20px;
            border: 2px solid #e0e0e0;
            border-radius: 5px;
            font-size: 14px;
        }}
        
        .search-box:focus {{
            outline: none;
            border-color: #dc3545;
        }}
        
        .table-container {{
            overflow-x: auto;
            max-height: 500px;
            overflow-y: auto;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            margin-bottom: 20px;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
            font-size: 14px;
        }}
        
        th {{
            background: #dc3545;
            color: white;
            padding: 12px;
            position: sticky;
            top: 0;
            z-index: 10;
            white-space: nowrap;
        }}
        
        td {{
            padding: 10px 12px;
            border-bottom: 1px solid #e0e0e0;
            vertical-align: top;
            word-break: break-word;
        }}
        
        tr:hover {{
            background-color: #f8f9fa;
        }}
        
        .badge {{
            display: inline-block;
            padding: 3px 8px;
            border-radius: 12px;
            font-size: 11px;
            font-weight: 600;
            white-space: nowrap;
        }}
        
        .badge-critical {{ background: #dc3545; color: white; }}
        .badge-high {{ background: #ffc107; color: black; }}
        .badge-medium {{ background: #17a2b8; color: white; }}
        .badge-low {{ background: #28a745; color: white; }}
        
        .genome-list {{
            max-height: 150px;
            overflow-y: auto;
            padding: 5px;
            background: #f8f9fa;
            border-radius: 5px;
        }}
        
        .genome-tag {{
            display: inline-block;
            background: #e0f2f1;
            color: #00695c;
            padding: 2px 8px;
            border-radius: 12px;
            font-size: 11px;
            margin: 2px;
            white-space: nowrap;
        }}
        
        .alert-box {{
            padding: 15px;
            border-radius: 8px;
            margin: 20px 0;
            display: flex;
            align-items: center;
            gap: 15px;
        }}
        
        .alert-danger {{ background: #f8d7da; color: #721c24; border-left: 4px solid #dc3545; }}
        .alert-warning {{ background: #fff3cd; color: #856404; border-left: 4px solid #ffc107; }}
        .alert-info {{ background: #d1ecf1; color: #0c5460; border-left: 4px solid #17a2b8; }}
        .alert-success {{ background: #d4edda; color: #155724; border-left: 4px solid #28a745; }}
        
        .action-buttons {{
            display: flex;
            gap: 10px;
            margin: 20px 0;
            flex-wrap: wrap;
        }}
        
        .action-btn {{
            padding: 8px 16px;
            border: none;
            border-radius: 5px;
            cursor: pointer;
            font-weight: 600;
            display: inline-flex;
            align-items: center;
            gap: 8px;
            transition: all 0.3s ease;
        }}
        
        .btn-primary {{ background: #dc3545; color: white; }}
        .btn-primary:hover {{ background: #c82333; }}
        .btn-success {{ background: #28a745; color: white; }}
        .btn-success:hover {{ background: #218838; }}
        
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0,0,0,0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        
        .footer a {{
            color: #dc3545;
            text-decoration: none;
        }}
        
        .footer a:hover {{
            text-decoration: underline;
        }}
        
        .timestamp {{
            color: #fbbf24;
        }}
        
        @media print {{
            .tab-navigation, .action-buttons, .search-box {{ display: none; }}
            body {{ background: white; color: black; }}
            .tab-content {{ display: block; }}
        }}
    </style>
    <script>
        function switchTab(tabName, btn) {{
            // Hide all tabs
            document.querySelectorAll('.tab-content').forEach(tab => {{
                tab.classList.remove('active');
            }});
            
            // Remove active class from all buttons
            document.querySelectorAll('.tab-button').forEach(button => {{
                button.classList.remove('active');
            }});
            
            // Show selected tab
            document.getElementById(tabName + '-tab').classList.add('active');
            
            // Activate clicked button
            btn.classList.add('active');
        }}
        
        function searchTable(tableId, inputId) {{
            const input = document.getElementById(inputId);
            const filter = input.value.toLowerCase();
            const table = document.getElementById(tableId);
            const rows = table.getElementsByTagName('tr');
            
            for (let i = 1; i < rows.length; i++) {{
                const row = rows[i];
                const text = row.textContent.toLowerCase();
                row.style.display = text.includes(filter) ? '' : 'none';
            }}
        }}
        
        function exportTableToCSV(tableId, filename) {{
            const table = document.getElementById(tableId);
            const rows = table.querySelectorAll('tr');
            const csv = [];
            
            for (let i = 0; i < rows.length; i++) {{
                const row = [], cols = rows[i].querySelectorAll('td, th');
                for (let j = 0; j < cols.length; j++) {{
                    row.push('"' + cols[j].innerText.replace(/"/g, '""') + '"');
                }}
                csv.push(row.join(','));
            }}
            
            const csvContent = csv.join('\\n');
            const blob = new Blob([csvContent], {{ type: 'text/csv' }});
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = filename;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
        }}
        
        function rotateQuote() {{
            const quotes = {json.dumps(self.analyzer.science_quotes)};
            const idx = Math.floor(Math.random() * quotes.length);
            document.getElementById('quote-text').innerHTML = quotes[idx].text;
            document.getElementById('quote-author').innerHTML = '— ' + quotes[idx].author;
        }}
        
        setInterval(rotateQuote, 10000);
    </script>
</head>
<body>
    <div class="container">
        <!-- Header -->
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.analyzer.ascii_art}</div>
            </div>
            
            <div class="quote-container">
                <div id="quote-text" class="quote-text">"{quote['text']}"</div>
                <div id="quote-author" class="quote-author">— {quote['author']}</div>
            </div>
        </div>
        
        <!-- Dashboard Cards -->
        <div class="dashboard-grid">
            <div class="dashboard-card" onclick="switchTab('summary', this)">
                <div class="card-value">{total_samples}</div>
                <div class="card-label">Samples Analyzed</div>
            </div>
            <div class="dashboard-card" onclick="switchTab('mlst', this)">
                <div class="card-value">{len(patterns.get('st_distribution', {}))}</div>
                <div class="card-label">Unique STs</div>
            </div>
            <div class="dashboard-card" onclick="switchTab('past', this)">
                <div class="card-value">{len(patterns.get('serotype_distribution', {}))}</div>
                <div class="card-label">Serotypes</div>
            </div>
            <div class="dashboard-card" onclick="switchTab('amr', this)">
                <div class="card-value">{total_amr_genes}</div>
                <div class="card-label">AMR Genes</div>
            </div>
            <div class="dashboard-card" onclick="switchTab('environmental', this)">
                <div class="card-value">{environmental_marker_count}</div>
                <div class="card-label">Environmental Markers</div>
            </div>
            <div class="dashboard-card" onclick="switchTab('qc', this)">
                <div class="card-value">{len(qc_data)}</div>
                <div class="card-label">QC Data Available</div>
            </div>
        </div>
        
        <!-- Tab Navigation -->
        <div class="tab-navigation">
            <button class="tab-button active" onclick="switchTab('summary', this)"><i class="fas fa-chart-pie"></i> Summary</button>
            <button class="tab-button" onclick="switchTab('samples', this)"><i class="fas fa-list"></i> Samples</button>
            <button class="tab-button" onclick="switchTab('mlst', this)"><i class="fas fa-dna"></i> MLST</button>
            <button class="tab-button" onclick="switchTab('past', this)"><i class="fas fa-tag"></i> Serotypes</button>
            <button class="tab-button" onclick="switchTab('amr', this)"><i class="fas fa-biohazard"></i> AMR Genes</button>
            <button class="tab-button" onclick="switchTab('virulence', this)"><i class="fas fa-virus"></i> Virulence</button>
            <button class="tab-button" onclick="switchTab('environmental', this)"><i class="fas fa-globe"></i> Environmental</button>
            <button class="tab-button" onclick="switchTab('categories', this)"><i class="fas fa-tags"></i> Categories</button>
            <button class="tab-button" onclick="switchTab('patterns', this)"><i class="fas fa-project-diagram"></i> Patterns</button>
            <button class="tab-button" onclick="switchTab('qc', this)"><i class="fas fa-chart-line"></i> QC Metrics</button>
            <button class="tab-button" onclick="switchTab('databases', this)"><i class="fas fa-database"></i> Databases</button>
        </div>
        
        <!-- Summary Tab -->
        <div id="summary-tab" class="tab-content active">
            <h2 class="section-header">
                <i class="fas fa-chart-pie"></i> Executive Summary
            </h2>
            
            <div class="alert-box alert-info">
                <i class="fas fa-info-circle fa-2x"></i>
                <div>
                    <h3>P. aeruginosa Analysis Overview</h3>
                    <p>This tab provides a high-level summary of the dataset. The cards above give quick counts, and the table below shows key metrics including the number of carbapenemase genes, MDR patterns, and environmental markers. Use the tabs above to explore detailed results for each analysis type.</p>
                </div>
            </div>
            
            {f'<div class="alert-box alert-danger"><i class="fas fa-exclamation-triangle fa-2x"></i><div><h3>⚠️ Critical Alert</h3><p><strong>{carbapenemase_count} carbapenemase genes</strong> detected. Carbapenem-resistant P. aeruginosa is a critical public health threat.</p></div></div>' if carbapenemase_count > 0 else ''}
            
            {f'<div class="alert-box alert-warning"><i class="fas fa-exclamation-triangle fa-2x"></i><div><h3>⚠️ High-Risk Combinations</h3><p><strong>{high_risk_count} samples</strong> contain carbapenemase + colistin resistance genes.</p></div></div>' if high_risk_count > 0 else ''}
            
            {f'<div class="alert-box alert-warning"><i class="fas fa-radiation fa-2x"></i><div><h3>⚠️ MDR Detected</h3><p><strong>{mdr_count} samples</strong> show resistance to 3+ antibiotic classes.</p></div></div>' if mdr_count > 0 else ''}
            
            <h3>Key Statistics</h3>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Metric</th>
                            <th>Value</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr><td>Total Samples</td><td><strong>{total_samples}</strong></td></tr>
                        <tr><td>Unique STs</td><td><strong>{len(patterns.get('st_distribution', {}))}</strong></td></tr>
                        <tr><td>Unique Serotypes</td><td><strong>{len(patterns.get('serotype_distribution', {}))}</strong></td></tr>
                        <tr><td>Total AMR Genes</td><td><strong>{total_amr_genes}</strong></td></tr>
                        <tr><td>Carbapenemase Genes</td><td><span class="badge badge-critical">{carbapenemase_count}</span></td></tr>
                        <tr><td>Virulence Genes</td><td><strong>{total_virulence_genes}</strong></td></tr>
                        <tr><td>Environmental Markers</td><td><strong>{total_environmental_genes}</strong></td></tr>
                        <tr><td>High-Risk Combinations</td><td><span class="badge badge-critical">{high_risk_count}</span></td></tr>
                        <tr><td>MDR Patterns</td><td><span class="badge badge-high">{mdr_count}</span></td></tr>
                    </tbody>
                </table>
            </div>
        </div>
        
        <!-- Samples Tab -->
        <div id="samples-tab" class="tab-content">
            <h2 class="section-header">
                <i class="fas fa-list"></i> Sample Overview
            </h2>
            
            <div class="alert-box alert-info">
                <i class="fas fa-info-circle fa-2x"></i>
                <div>
                    <h4>About This Tab</h4>
                    <p>This table lists each sample with its MLST sequence type (ST), serotype, and the total number of genes detected (AMR + other). Use the search box to filter samples, and the Export button to download the data.</p>
                </div>
            </div>
            
            <input type="text" id="search-samples" class="search-box" 
                   placeholder="Search samples..." onkeyup="searchTable('samples-table', 'search-samples')">
            
            <div class="action-buttons">
                <button class="action-btn btn-primary" onclick="exportTableToCSV('samples-table', 'samples.csv')">
                    <i class="fas fa-download"></i> Export CSV
                </button>
            </div>
            
            <div class="table-container">
                <table id="samples-table">
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>ST</th>
                            <th>Serotype</th>
                            <th>Total Genes</th>
                        </tr>
                    </thead>
                    <tbody>
        '''
        
        for sample, data in sorted(samples_data.items()):
            amr_count = len(data.get('amr_genes', []))
            other_count = len(data.get('other_genes', []))
            total = amr_count + other_count
            
            html += f'''
                        <tr>
                            <td><strong>{sample}</strong></td>
                            <td>{data['mlst'].get('ST', 'ND')}</td>
                            <td>{data['past'].get('Serotype', 'ND')}</td>
                            <td>{total}</td>
                        </tr>
            '''
        
        html += '''
                    </tbody>
                </table>
            </div>
        </div>
        
        <!-- MLST Tab -->
        <div id="mlst-tab" class="tab-content">
            <h2 class="section-header">
                <i class="fas fa-dna"></i> MLST Analysis
            </h2>
            
            <div class="alert-box alert-info">
                <i class="fas fa-info-circle fa-2x"></i>
                <div>
                    <h4>About This Tab</h4>
                    <p>Distribution of sequence types (ST) across samples. Each row shows the ST, its frequency (count and percentage), and the list of samples carrying that ST.</p>
                </div>
            </div>
            
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>ST</th>
                            <th>Count</th>
                            <th>Percentage</th>
                            <th>Samples</th>
                        </tr>
                    </thead>
                    <tbody>
        '''
        
        st_dist = patterns.get('st_distribution', {})
        total_st = sum(st_dist.values())
        for st, count in sorted(st_dist.items(), key=lambda x: x[1], reverse=True):
            if st == 'ND':
                continue
            percentage = (count / total_st * 100) if total_st else 0
            samples = []
            for sample, data in samples_data.items():
                if data.get('mlst', {}).get('ST') == st:
                    samples.append(sample)
            sample_list = ', '.join(samples)
            html += f'''
                        <tr>
                            <td><strong>ST{st}</strong></td>
                            <td>{count}</td>
                            <td>{percentage:.1f}%</td>
                            <td>{sample_list}</td>
                        </tr>
            '''
        
        html += '''
                    </tbody>
                </table>
            </div>
        </div>
        
        <!-- PAST Tab -->
        <div id="past-tab" class="tab-content">
            <h2 class="section-header">
                <i class="fas fa-tag"></i> Serotype Analysis
            </h2>
            
            <div class="alert-box alert-info">
                <i class="fas fa-info-circle fa-2x"></i>
                <div>
                    <h4>About This Tab</h4>
                    <p>Distribution of O‑serotypes across samples. Each row shows the serotype, its frequency (count and percentage), and the list of samples carrying that serotype.</p>
                </div>
            </div>
            
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Serotype</th>
                            <th>Count</th>
                            <th>Percentage</th>
                            <th>Samples</th>
                        </tr>
                    </thead>
                    <tbody>
        '''
        
        sero_dist = patterns.get('serotype_distribution', {})
        total_sero = sum(sero_dist.values())
        for sero, count in sorted(sero_dist.items(), key=lambda x: x[1], reverse=True):
            if sero == 'ND' or sero == 'UNKNOWN':
                continue
            percentage = (count / total_sero * 100) if total_sero else 0
            samples = []
            for sample, data in samples_data.items():
                if data.get('past', {}).get('Serotype') == sero:
                    samples.append(sample)
            sample_list = ', '.join(samples)
            html += f'''
                        <tr>
                            <td><strong>{sero}</strong></td>
                            <td>{count}</td>
                            <td>{percentage:.1f}%</td>
                            <td>{sample_list}</td>
                        </tr>
            '''
        
        html += '''
                    </tbody>
                </table>
            </div>
        </div>
        
        <!-- AMR Genes Tab -->
        <div id="amr-tab" class="tab-content">
            <h2 class="section-header">
                <i class="fas fa-biohazard"></i> Antimicrobial Resistance Genes
            </h2>
            
            <div class="alert-box alert-info">
                <i class="fas fa-info-circle fa-2x"></i>
                <div>
                    <h4>Gene‑Centric View</h4>
                    <p>This table lists all detected AMR genes from AMRfinder and ABRicate databases. Each gene is shown with its functional category, database source, frequency (count and percentage), risk level, and the genomes that contain it. The genome list is scrollable – hover or click to see all. Use the search box to filter by gene name, category, or database.</p>
                </div>
            </div>
            
            <input type="text" id="search-amr" class="search-box" 
                   placeholder="Search genes..." onkeyup="searchTable('amr-table', 'search-amr')">
            
            <div class="action-buttons">
                <button class="action-btn btn-primary" onclick="exportTableToCSV('amr-table', 'amr_genes.csv')">
                    <i class="fas fa-download"></i> Export CSV
                </button>
            </div>
            
            <div class="table-container">
                <table id="amr-table">
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>Category</th>
                            <th>Database</th>
                            <th>Frequency</th>
                            <th>Risk Level</th>
                            <th>Genomes</th>
                        </tr>
                    </thead>
                    <tbody>
        '''
        
        all_amr_genes = []
        for db_name, genes in gene_centric.get('amr_databases', {}).items():
            all_amr_genes.extend(genes)
        all_amr_genes.sort(key=lambda x: x['count'], reverse=True)
        
        for gene_data in all_amr_genes:
            risk_class = 'badge-critical' if gene_data['risk_level'] in ['CRITICAL', 'HIGH'] else 'badge-low'
            genome_tags = ''.join([f'<span class="genome-tag">{g}</span>' for g in gene_data['genomes']])
            html += f'''
                        <tr>
                            <td><strong>{gene_data['gene']}</strong></td>
                            <td>{gene_data['category']}</td>
                            <td>{gene_data['database']}</td>
                            <td>{gene_data['frequency_display']}</td>
                            <td><span class="badge {risk_class}">{gene_data['risk_level']}</span></td>
                            <td><div class="genome-list">{genome_tags}</div></td>
                        </tr>
            '''
        
        html += '''
                    </tbody>
                </table>
            </div>
        </div>
        
        <!-- Virulence Genes Tab -->
        <div id="virulence-tab" class="tab-content">
            <h2 class="section-header">
                <i class="fas fa-virus"></i> Virulence Genes
            </h2>
            
            <div class="alert-box alert-info">
                <i class="fas fa-info-circle fa-2x"></i>
                <div>
                    <h4>About This Tab</h4>
                    <p>Virulence factors detected across samples, with their frequencies and the genomes carrying them. The list is scrollable, and you can search by gene name or category.</p>
                </div>
            </div>
            
            <input type="text" id="search-virulence" class="search-box" 
                   placeholder="Search genes..." onkeyup="searchTable('virulence-table', 'search-virulence')">
            
            <div class="action-buttons">
                <button class="action-btn btn-primary" onclick="exportTableToCSV('virulence-table', 'virulence_genes.csv')">
                    <i class="fas fa-download"></i> Export CSV
                </button>
            </div>
            
            <div class="table-container">
                <table id="virulence-table">
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>Category</th>
                            <th>Database</th>
                            <th>Frequency</th>
                            <th>Genomes</th>
                        </tr>
                    </thead>
                    <tbody>
        '''
        
        all_vir_genes = []
        for db_name, genes in gene_centric.get('virulence_databases', {}).items():
            all_vir_genes.extend(genes)
        all_vir_genes.sort(key=lambda x: x['count'], reverse=True)
        
        for gene_data in all_vir_genes:
            genome_tags = ''.join([f'<span class="genome-tag">{g}</span>' for g in gene_data['genomes']])
            html += f'''
                        <tr>
                            <td><strong>{gene_data['gene']}</strong></td>
                            <td>{gene_data['category']}</td>
                            <td>{gene_data['database']}</td>
                            <td>{gene_data['frequency_display']}</td>
                            <td><div class="genome-list">{genome_tags}</div></td>
                        </tr>
            '''
        
        html += '''
                    </tbody>
                </table>
            </div>
        </div>
        
        <!-- Environmental Tab -->
        <div id="environmental-tab" class="tab-content">
            <h2 class="section-header">
                <i class="fas fa-globe"></i> Environmental Markers
            </h2>
            
            <div class="alert-box alert-info">
                <i class="fas fa-info-circle fa-2x"></i>
                <div>
                    <h3>Environmental Co-Selection Analysis</h3>
                    <p>These markers indicate potential for co-selection of antibiotic resistance through heavy metal and biocide exposure in hospital environments. Each category groups related genes, with their frequencies and the genomes carrying them.</p>
                </div>
            </div>
        '''
        
        env_summary = gene_centric.get('environmental_summary', {})
        for category, summary in env_summary.items():
            genes = summary.get('genes', [])
            if not genes:
                continue
            html += f'''
            <h3>{category} ({len(genes)} genes, {summary.get('total_occurrences', 0)} occurrences)</h3>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>Database</th>
                            <th>Frequency</th>
                            <th>Genomes</th>
                        </tr>
                    </thead>
                    <tbody>
            '''
            for gene_data in genes:
                genome_tags = ''.join([f'<span class="genome-tag">{g}</span>' for g in gene_data['genomes']])
                html += f'''
                        <tr>
                            <td><strong>{gene_data['gene']}</strong></td>
                            <td>{gene_data['database']}</td>
                            <td>{gene_data['frequency_display']}</td>
                            <td><div class="genome-list">{genome_tags}</div></td>
                        </tr>
                '''
            html += '''
                    </tbody>
                </table>
            </div>
            '''
        
        html += '''
        </div>
        
        <!-- Categories Tab -->
        <div id="categories-tab" class="tab-content">
            <h2 class="section-header">
                <i class="fas fa-tags"></i> Gene Categories
            </h2>
            
            <div class="alert-box alert-info">
                <i class="fas fa-info-circle fa-2x"></i>
                <div>
                    <h4>About This Tab</h4>
                    <p>Genes grouped by functional category. The overview table shows counts per category; below are detailed lists for the top five categories.</p>
                </div>
            </div>
            
            <h3>Category Overview</h3>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Category</th>
                            <th>Unique Genes</th>
                            <th>Total Occurrences</th>
                        </tr>
                    </thead>
                    <tbody>
        '''
        
        categories = gene_centric.get('gene_categories', {})
        for category, genes in sorted(categories.items(), key=lambda x: len(x[1]), reverse=True):
            unique_genes = len(set(g['gene'] for g in genes))
            total_occ = sum(g['count'] for g in genes)
            html += f'''
                        <tr>
                            <td><strong>{category}</strong></td>
                            <td>{unique_genes}</td>
                            <td>{total_occ}</td>
                        </tr>
            '''
        
        html += '''
                    </tbody>
                </table>
            </div>
        '''
        
        # Detailed category lists
        for category, genes in sorted(categories.items(), key=lambda x: len(x[1]), reverse=True)[:5]:
            if not genes:
                continue
            genes.sort(key=lambda x: x['count'], reverse=True)
            html += f'''
            <h3 style="margin-top: 30px;">{category} Details (Top 20)</h3>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>Database</th>
                            <th>Frequency</th>
                            <th>Genomes</th>
                        </tr>
                    </thead>
                    <tbody>
            '''
            for gene_data in genes[:20]:
                genome_tags = ''.join([f'<span class="genome-tag">{g}</span>' for g in gene_data['genomes']])
                html += f'''
                        <tr>
                            <td><strong>{gene_data['gene']}</strong></td>
                            <td>{gene_data['database']}</td>
                            <td>{gene_data['frequency_display']}</td>
                            <td><div class="genome-list">{genome_tags}</div></td>
                        </tr>
                '''
            html += '''
                    </tbody>
                </table>
            </div>
            '''
        
        html += '''
        </div>
        
        <!-- Patterns Tab -->
        <div id="patterns-tab" class="tab-content">
            <h2 class="section-header">
                <i class="fas fa-project-diagram"></i> Pattern Discovery
            </h2>
            
            <div class="alert-box alert-info">
                <i class="fas fa-info-circle fa-2x"></i>
                <div>
                    <h4>Novel AMR Insights</h4>
                    <p>This tab reveals important associations: ST‑serotype combinations (lineage tracking), carbapenemase gene co‑occurrence patterns, and a novel co‑selection risk assessment that quantifies the potential for environmental maintenance of resistance.</p>
                </div>
            </div>
        '''
        
        # ST-Serotype Patterns
        if st_serotype_patterns:
            total_st_sero = len(st_serotype_patterns)
            html += '''
            <h3>ST-Serotype Combinations</h3>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>ST-Serotype</th>
                            <th>Frequency</th>
                            <th>Percentage</th>
                            <th>Samples</th>
                        </tr>
                    </thead>
                    <tbody>
            '''
            for combo, samples in st_serotype_patterns.items():
                sample_list = ', '.join(samples)
                percentage = (len(samples) / total_samples * 100) if total_samples else 0
                html += f'''
                        <tr>
                            <td><strong>{combo}</strong></td>
                            <td>{len(samples)}</td>
                            <td>{percentage:.1f}%</td>
                            <td>{sample_list}</td>
                        </tr>
                '''
            html += '''
                    </tbody>
                </table>
            </div>
            '''
        
        # Carbapenemase Patterns
        if carb_patterns:
            html += '''
            <h3>Carbapenemase Combinations</h3>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Carbapenemase Combination</th>
                            <th>Frequency</th>
                            <th>Samples</th>
                        </tr>
                    </thead>
                    <tbody>
            '''
            for genes_str, samples in sorted(carb_patterns.items(), key=lambda x: len(x[1]), reverse=True):
                genes_display = genes_str.replace('|', ' + ')
                sample_list = ', '.join(samples)
                html += f'''
                        <tr>
                            <td><strong>{genes_display}</strong></td>
                            <td>{len(samples)}</td>
                            <td>{sample_list}</td>
                        </tr>
                '''
            html += '''
                    </tbody>
                </table>
            </div>
            '''
        
        # Co-selection Risk Assessment (NOVEL)
        if co_risk:
            html += '''
            <h3>Co-selection Risk Assessment</h3>
            <div class="alert-box alert-info">
                <i class="fas fa-info-circle fa-2x"></i>
                <div>
                    <h4>Fighting AMR with Co-selection Risk</h4>
                    <p>This table quantifies the co-selection potential: high numbers of biocide/metal resistance markers together with many AMR genes indicate strains at risk of being maintained under environmental pressure (disinfectants, heavy metals). High-risk isolates may persist in hospital settings and warrant enhanced infection control.</p>
                </div>
            </div>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>ST</th>
                            <th>Serotype</th>
                            <th>Biocide/Metal Markers</th>
                            <th>AMR Genes Count</th>
                            <th>Risk Category</th>
                            <th>Biocide/Metal Genes</th>
                        </tr>
                    </thead>
                    <tbody>
            '''
            for risk in sorted(co_risk, key=lambda x: (x['risk_category'] != 'High', x['risk_category'] != 'Medium'), reverse=True):
                badge_class = 'badge-critical' if risk['risk_category'] == 'High' else 'badge-medium' if risk['risk_category'] == 'Medium' else 'badge-low'
                genes_str = ', '.join(risk['biocide_metal_genes'])
                html += f'''
                        <tr>
                            <td><strong>{risk['sample']}</strong></td>
                            <td>ST{risk['st']}</td>
                            <td>{risk['serotype']}</td>
                            <td>{risk['biocide_metal_markers']}</td>
                            <td>{risk['amr_genes_count']}</td>
                            <td><span class="badge {badge_class}">{risk['risk_category']}</span></td>
                            <td>{genes_str}</td>
                        </tr>
                '''
            html += '''
                    </tbody>
                </table>
            </div>
            '''
        
        # High-Risk Combinations
        if patterns.get('high_risk_combinations'):
            html += '''
            <h3>High-Risk Combinations (Carbapenemase + Colistin)</h3>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>ST</th>
                            <th>Serotype</th>
                            <th>Carbapenemases</th>
                            <th>Colistin Resistance</th>
                        </tr>
                    </thead>
                    <tbody>
            '''
            for hr in patterns['high_risk_combinations']:
                carb_str = ', '.join(hr['carbapenemases'])
                col_str = ', '.join(hr['colistin_resistance'])
                html += f'''
                        <tr>
                            <td><strong>{hr['sample']}</strong></td>
                            <td>ST{hr['st']}</td>
                            <td>{hr['serotype']}</td>
                            <td>{carb_str}</td>
                            <td>{col_str}</td>
                        </tr>
                '''
            html += '''
                    </tbody>
                </table>
            </div>
            '''
        
        html += '''
        </div>
        
        <!-- QC Tab -->
        <div id="qc-tab" class="tab-content">
            <h2 class="section-header">
                <i class="fas fa-chart-line"></i> Quality Control Metrics
            </h2>
            
            <div class="alert-box alert-info">
                <i class="fas fa-info-circle fa-2x"></i>
                <div>
                    <h4>Assembly Statistics</h4>
                    <p>This table displays quality metrics from FASTA QC analysis: genome size (total length), GC content, N50 (contiguity), number of contigs, and any quality warnings. These metrics help assess assembly completeness and potential contamination.</p>
                </div>
            </div>
            
            <input type="text" id="search-qc" class="search-box" 
                   placeholder="Search samples..." onkeyup="searchTable('qc-table', 'search-qc')">
            
            <div class="action-buttons">
                <button class="action-btn btn-primary" onclick="exportTableToCSV('qc-table', 'qc_metrics.csv')">
                    <i class="fas fa-download"></i> Export CSV
                </button>
            </div>
            
            <div class="table-container">
                <table id="qc-table">
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>Genome Size (bp)</th>
                            <th>GC%</th>
                            <th>N50 (bp)</th>
                            <th>Contigs</th>
                            <th>Warnings</th>
                        </tr>
                    </thead>
                    <tbody>
        '''
        
        for sample, data in sorted(qc_data.items()):
            html += f'''
                        <tr>
                            <td><strong>{sample}</strong></td>
                            <td>{data.get('Total_Length', 'N/A')}</td>
                            <td>{data.get('GC_Content', 'N/A')}</td>
                            <td>{data.get('N50', 'N/A')}</td>
                            <td>{data.get('Total_Sequences', 'N/A')}</td>
                            <td>{data.get('Warnings', 'None')}</td>
                        </tr>
            '''
        
        html += '''
                    </tbody>
                </table>
            </div>
        </div>
        
        <!-- Databases Tab -->
        <div id="databases-tab" class="tab-content">
            <h2 class="section-header">
                <i class="fas fa-database"></i> Database Coverage
            </h2>
            
            <div class="alert-box alert-info">
                <i class="fas fa-info-circle fa-2x"></i>
                <div>
                    <h4>About This Tab</h4>
                    <p>Shows the percentage of samples with hits in each database, along with unique gene counts and critical gene tallies. This helps evaluate database sensitivity.</p>
                </div>
            </div>
            
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Database</th>
                            <th>Coverage</th>
                            <th>Unique Genes</th>
                            <th>Total Hits</th>
                            <th>Critical Genes</th>
                        </tr>
                    </thead>
                    <tbody>
        '''
        
        db_coverage = patterns.get('database_coverage', {})
        db_stats = gene_centric.get('database_stats', {})
        for db_name, coverage in sorted(db_coverage.items(), key=lambda x: x[1]['coverage_percentage'], reverse=True):
            stats = db_stats.get(db_name, {})
            cov_class = 'badge-low' if coverage['coverage_percentage'] >= 80 else 'badge-medium' if coverage['coverage_percentage'] >= 50 else 'badge-high'
            html += f'''
                        <tr>
                            <td><strong>{db_name.upper()}</strong></td>
                            <td><span class="badge {cov_class}">{coverage['coverage_display']}</span></td>
                            <td>{stats.get('total_genes', 0)}</td>
                            <td>{stats.get('total_occurrences', 0)}</td>
                            <td><span class="badge badge-critical">{stats.get('critical_genes', 0)}</span></td>
                        </tr>
            '''
        
        html += '''
                    </tbody>
                </table>
            </div>
        </div>
        
        <!-- Footer -->
        <div class="footer">
            <p><strong>PSEUDOSCOPE</strong> — P. aeruginosa Ultimate Reporter v1.0.0</p>
            <p>Author: <a href="mailto:brownbeckley94@gmail.com">Brown Beckley</a> | GitHub: <a href="https://github.com/bbeckley-hub" target="_blank">bbeckley-hub</a></p>
            <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
        </div>
    </div>
</body>
</html>
'''
        
        return html


class PseudoUltimateReporter:
    """Master class that orchestrates parsing, analysis, and report generation."""

    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "PSEUDO_ULTIMATE_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = PseudoHTMLParser()
        self.analyzer = PseudoDataAnalyzer()
        self.html_generator = PseudoHTMLGenerator(self.analyzer)
        self.metadata = {
            "tool_name": "PseudoScope Ultimate Reporter",
            "version": "1.0.0",
            "author": "Brown Beckley <brownbeckley94@gmail.com>",
            "affiliation": "University of Ghana Medical School",
            "pathogen": "Pseudomonas aeruginosa",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir)
        }

    def find_files(self) -> Dict[str, List[Path]]:
        """Categorize all summary HTML files."""
        print("🔍 Searching for summary HTML files...")
        files = {
            'mlst': [],
            'past': [],
            'qc': [],
            'amrfinder': [],
            'abricate': defaultdict(list)
        }
        all_html = list(self.input_dir.glob("*.html"))
        for f in all_html:
            name = f.name.lower()
            if 'mlst' in name:
                files['mlst'].append(f)
            elif 'past' in name:
                files['past'].append(f)
            elif 'fasta_qc' in name or 'qc' in name:
                files['qc'].append(f)
            elif 'amrfinder' in name:
                files['amrfinder'].append(f)
            else:
                for key in self.parser.db_name_mapping:
                    if key in name:
                        files['abricate'][key].append(f)
                        break
        # Print summary
        print(f"  ✅ MLST: {len(files['mlst'])} file(s)")
        print(f"  ✅ PAST: {len(files['past'])} file(s)")
        print(f"  ✅ QC: {len(files['qc'])} file(s)")
        print(f"  ✅ AMRfinder: {len(files['amrfinder'])} file(s)")
        abricate_count = sum(len(v) for v in files['abricate'].values())
        print(f"  ✅ ABRicate databases: {abricate_count} file(s)")
        return files

    def integrate_data(self, files: Dict) -> Dict[str, Any]:
        """Parse all files and build integrated data structure."""
        print("\n🔗 Integrating data from all reports...")
        integrated = {
            'metadata': self.metadata,
            'samples': {},
            'gene_frequencies': {},
            'patterns': {},
            'gene_centric': {},
            'qc_data': {}
        }

        # Parse MLST
        mlst_data = {}
        if files['mlst']:
            mlst_data = self.parser.parse_mlst_summary(files['mlst'][0])
            print(f"    ✅ MLST parsed: {len(mlst_data)} samples")

        # Parse PAST
        past_data = {}
        if files['past']:
            past_data = self.parser.parse_past_summary(files['past'][0])
            print(f"    ✅ PAST parsed: {len(past_data)} samples")

        # Parse QC
        qc_data = {}
        if files['qc']:
            qc_data = self.parser.parse_qc_summary(files['qc'][0])
            print(f"    ✅ QC parsed: {len(qc_data)} samples")

        # Parse AMRfinder
        amr_by_genome, amr_freq = {}, {}
        if files['amrfinder']:
            amr_by_genome, amr_freq = self.parser.parse_gene_frequency_report(files['amrfinder'][0], total_samples=0)
            print(f"    ✅ AMRfinder parsed: {len(amr_by_genome)} samples, {len(amr_freq)} genes")
            integrated['gene_frequencies']['amrfinder'] = amr_freq

        # Parse ABRicate databases
        abricate_by_genome = defaultdict(dict)
        abricate_freq = {}
        for db_key, flist in files['abricate'].items():
            if flist:
                db_name = self.parser.db_name_mapping.get(db_key, db_key)
                gbg, gf = self.parser.parse_gene_frequency_report(flist[0], total_samples=0)
                if gbg:
                    abricate_by_genome[db_name] = gbg
                if gf:
                    abricate_freq[db_name] = gf
        if abricate_freq:
            integrated['gene_frequencies']['abricate'] = abricate_freq
        print(f"    ✅ ABRicate databases parsed: {len(abricate_freq)} databases")

        # Collect all unique sample names
        all_samples = set()
        for d in [mlst_data, past_data, qc_data, amr_by_genome]:
            all_samples.update(d.keys())
        for gbg in abricate_by_genome.values():
            all_samples.update(gbg.keys())
        all_samples = sorted(all_samples)
        total_samples = len(all_samples)
        print(f"\n📊 Total unique samples: {total_samples}")

        # Update percentages in all gene frequency data now that we know total_samples
        for db_name, freq_data in integrated['gene_frequencies'].items():
            if db_name == 'abricate':
                for sub_db, sub_freq in freq_data.items():
                    for gene, data in sub_freq.items():
                        data['percentage'] = (data['count'] / total_samples * 100) if total_samples else 0
                        data['frequency_display'] = f"{data['count']} ({data['percentage']:.1f}%)"
            else:
                for gene, data in freq_data.items():
                    data['percentage'] = (data['count'] / total_samples * 100) if total_samples else 0
                    data['frequency_display'] = f"{data['count']} ({data['percentage']:.1f}%)"

        # Build integrated sample records
        for sample in all_samples:
            sample_rec = {
                'mlst': mlst_data.get(sample, {'ST': 'ND', 'Allele_Profile': ''}),
                'past': past_data.get(sample, {'Serotype': 'ND', 'Targets': '', 'Coverage': '', 'Hits': ''}),
                'qc': qc_data.get(sample, {}),
                'amr_genes': amr_by_genome.get(sample, []),
                'other_genes': []
            }
            # Collect genes from ABRicate databases
            for db_name, gbg in abricate_by_genome.items():
                genes = gbg.get(sample, [])
                if genes:
                    sample_rec['other_genes'].extend(genes)
            integrated['samples'][sample] = sample_rec

        # Store QC data separately for summary
        integrated['qc_data'] = qc_data

        # Gene-centric analysis
        print("\n🧠 Performing gene-centric analysis...")
        integrated['gene_centric'] = self.analyzer.create_gene_centric_tables(integrated, total_samples)

        # Cross-genome patterns
        integrated['patterns'] = self.analyzer.create_cross_genome_patterns(integrated, total_samples)

        return integrated

    def generate_reports(self, integrated_data: Dict[str, Any]):
        """Generate JSON and CSV reports."""
        # HTML
        self.html_generator.generate_main_report(integrated_data, self.output_dir)
        
        # JSON
        json_file = self.output_dir / "pseudoscope_ultimate_data.json"
        with open(json_file, 'w', encoding='utf-8') as f:
            # Make a serializable copy
            def make_serializable(obj):
                if isinstance(obj, (Counter, defaultdict)):
                    return dict(obj)
                if isinstance(obj, set):
                    return list(obj)
                if isinstance(obj, Path):
                    return str(obj)
                if isinstance(obj, datetime):
                    return obj.isoformat()
                # Handle any other non-serializable types
                try:
                    json.dumps(obj)
                    return obj
                except:
                    return str(obj)
            
            serializable = json.loads(json.dumps(integrated_data, default=make_serializable))
            json.dump(serializable, f, indent=2)
        print(f"    ✅ JSON data saved: {json_file}")

        # CSV summaries
        self._generate_csv_reports(integrated_data)

    def _generate_csv_reports(self, data: Dict[str, Any]):
        """Create CSV files for each major analysis using built-in csv module."""
        out = self.output_dir

        # Samples overview
        with open(out / "pseudoscope_samples.csv", 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['Sample', 'ST', 'Serotype', 'Total Genes'])
            for sample, rec in data['samples'].items():
                amr_count = len(rec['amr_genes'])
                other_count = len(rec.get('other_genes', []))
                total = amr_count + other_count
                
                w.writerow([sample, rec['mlst']['ST'], rec['past']['Serotype'], total])

        # AMR genes (gene-centric)
        rows = []
        for db_name, genes in data['gene_centric'].get('amr_databases', {}).items():
            for g in genes:
                rows.append([
                    g['gene'], g['category'], g['database'],
                    g['frequency_display'], g['count'], g['percentage'],
                    g.get('risk_level', ''), ';'.join(g['genomes'])
                ])
        if rows:
            with open(out / "pseudoscope_amr_genes.csv", 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(['Gene', 'Category', 'Database', 'Frequency', 'Count', 'Percentage', 'Risk Level', 'Genomes'])
                w.writerows(rows)

        # Virulence genes
        rows = []
        for db_name, genes in data['gene_centric'].get('virulence_databases', {}).items():
            for g in genes:
                rows.append([
                    g['gene'], g['category'], g['database'],
                    g['frequency_display'], g['count'], g['percentage'],
                    ';'.join(g['genomes'])
                ])
        if rows:
            with open(out / "pseudoscope_virulence_genes.csv", 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(['Gene', 'Category', 'Database', 'Frequency', 'Count', 'Percentage', 'Genomes'])
                w.writerows(rows)

        # Environmental markers
        rows = []
        for cat, summary in data['gene_centric'].get('environmental_summary', {}).items():
            for g in summary.get('genes', []):
                rows.append([
                    cat, g['gene'], g['database'],
                    g['frequency_display'], g['count'], g['percentage'],
                    ';'.join(g['genomes'])
                ])
        if rows:
            with open(out / "pseudoscope_environmental_markers.csv", 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(['Category', 'Gene', 'Database', 'Frequency', 'Count', 'Percentage', 'Genomes'])
                w.writerows(rows)

        # High-risk combinations
        if data['patterns'].get('high_risk_combinations'):
            rows = []
            for hr in data['patterns']['high_risk_combinations']:
                rows.append([
                    hr['sample'], hr['st'], hr['serotype'],
                    ';'.join(hr['carbapenemases']),
                    ';'.join(hr['colistin_resistance'])
                ])
            with open(out / "pseudoscope_high_risk.csv", 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(['Sample', 'ST', 'Serotype', 'Carbapenemases', 'Colistin Resistance'])
                w.writerows(rows)

        # MDR patterns
        if data['patterns'].get('mdr_patterns'):
            rows = []
            for mdr in data['patterns']['mdr_patterns']:
                rows.append([
                    mdr['sample'], mdr['st'], mdr['serotype'], mdr['classes'],
                    ';'.join(mdr['carbapenemases']),
                    ';'.join(mdr['esbls']),
                    ';'.join(mdr['colistin_resistance']),
                    ';'.join(mdr['aminoglycoside_resistance']),
                    ';'.join(mdr['efflux_pumps']),
                    ';'.join(mdr['environmental_markers'])
                ])
            with open(out / "pseudoscope_mdr_patterns.csv", 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(['Sample', 'ST', 'Serotype', 'Classes', 'Carbapenemases', 'ESBLs',
                           'Colistin Resistance', 'Aminoglycoside Resistance', 'Efflux Pumps', 'Environmental Markers'])
                w.writerows(rows)

        # Carbapenemase patterns
        if data['patterns'].get('carbapenemase_patterns'):
            rows = []
            for genes_str, samples in data['patterns']['carbapenemase_patterns'].items():
                rows.append([
                    genes_str.replace('|', ', '),
                    len(samples),
                    ';'.join(samples)
                ])
            with open(out / "pseudoscope_carbapenemase_patterns.csv", 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(['Carbapenemase Combination', 'Frequency', 'Samples'])
                w.writerows(rows)

        # Co-selection risk
        if data['patterns'].get('co_selection_risk'):
            rows = []
            for risk in data['patterns']['co_selection_risk']:
                rows.append([
                    risk['sample'], risk['st'], risk['serotype'],
                    risk['biocide_metal_markers'], risk['amr_genes_count'],
                    risk['risk_category'], ';'.join(risk['biocide_metal_genes'])
                ])
            with open(out / "pseudoscope_co_selection_risk.csv", 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(['Sample', 'ST', 'Serotype', 'Biocide/Metal Markers', 'AMR Count', 'Risk Category', 'Biocide/Metal Genes'])
                w.writerows(rows)

        # Database coverage
        if data['patterns'].get('database_coverage'):
            rows = []
            for db, cov in data['patterns']['database_coverage'].items():
                rows.append([
                    db, cov['samples_with_hits'], cov['total_samples'],
                    cov['coverage_percentage'], cov['coverage_display']
                ])
            with open(out / "pseudoscope_database_coverage.csv", 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(['Database', 'Samples with Hits', 'Total Samples', 'Coverage %', 'Coverage'])
                w.writerows(rows)

        print("    ✅ CSV reports generated.")

    def run(self):
        """Execute the full analysis pipeline."""
        print("=" * 80)
        print("🧬 PseudoScope - Pseudomonas aeruginosa Ultimate Reporter v1.0.0")
        print("=" * 80)
        print(f"📁 Input directory: {self.input_dir}")
        print("=" * 80)

        files = self.find_files()
        if not any(files.values()):
            print("❌ No summary HTML files found!")
            return False

        integrated = self.integrate_data(files)
        if not integrated['samples']:
            print("❌ No sample data could be extracted.")
            return False

        self.generate_reports(integrated)

        total = len(integrated['samples'])
        carb_count = len(integrated['patterns'].get('carbapenemase_patterns', {}))
        high_risk = len(integrated['patterns'].get('high_risk_combinations', []))
        mdr = len(integrated['patterns'].get('mdr_patterns', []))
        env_count = len(integrated['gene_centric'].get('environmental_summary', {}))

        print("\n" + "=" * 80)
        print("✅ ULTIMATE ANALYSIS COMPLETE!")
        print("=" * 80)
        print(f"📁 Output directory: {self.output_dir}")
        print(f"\n📊 Summary:")
        print(f"   • Samples analyzed: {total}")
        print(f"   • Carbapenemase patterns: {carb_count}")
        print(f"   • High-risk combinations: {high_risk}")
        print(f"   • MDR patterns: {mdr}")
        print(f"   • Environmental marker categories: {env_count}")
        print(f"\n📄 Reports generated:")
        print(f"   • pseudoscope_ultimate_report.html (interactive)")
        print(f"   • pseudoscope_ultimate_data.json (complete data)")
        print(f"   • pseudoscope_samples.csv")
        print(f"   • pseudoscope_amr_genes.csv")
        print(f"   • pseudoscope_virulence_genes.csv")
        print(f"   • pseudoscope_environmental_markers.csv")
        print(f"   • pseudoscope_high_risk.csv")
        print(f"   • pseudoscope_mdr_patterns.csv")
        print(f"   • pseudoscope_carbapenemase_patterns.csv")
        print(f"   • pseudoscope_co_selection_risk.csv (NOVEL)")
        print(f"   • pseudoscope_database_coverage.csv")
        print("\n🏫 University of Ghana Medical School")
        print("📧 brownbeckley94@gmail.com")
        print("=" * 80)
        return True


def main():
    parser = argparse.ArgumentParser(
        description='PseudoScope Ultimate Reporter for Pseudomonas aeruginosa',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python p_ultimate.py -i /path/to/summary/module/directory
  python p_ultimate.py -i ./summary_module -o ./my_reports
        """
    )
    parser.add_argument('-i', '--input-dir', required=True,
                        help='Directory containing summary HTML files')
    parser.add_argument('-o', '--output-dir',
                        help='Custom output directory (default: input_dir/PSEUDO_ULTIMATE_REPORTS)')
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        print(f"❌ Input directory not found: {input_dir}")
        sys.exit(1)

    reporter = PseudoUltimateReporter(input_dir)
    if args.output_dir:
        reporter.output_dir = Path(args.output_dir)
        reporter.output_dir.mkdir(parents=True, exist_ok=True)

    success = reporter.run()
    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()