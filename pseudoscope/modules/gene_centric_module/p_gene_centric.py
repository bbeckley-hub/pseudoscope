#!/usr/bin/env python3
"""
PSEUDOSCOPE ULTIMATE REPORTER - GENE-CENTRIC P. AERUGINOSA ANALYSIS
===================================================================
Version 3.0.0
Gene-centric cross-genome analysis with dynamic grouping by typing.
Gene-centric tabs (AMR, Virulence, Bacmet2, Plasmids, Mutations, MGE, Carbapenemase & OXA)
with dropdown grouping, filter buttons, gene-family info boxes, and rich educational content.

New in v3.0.0:
  - Dedicated Carbapenemase & OXA tab with literature-accurate family classification
  - Mutation tab simplified (Subclass column removed)
  - Footer: GitHub link + star call

New in v3.0.0:
  - ExoU/ExoS Pathotype tab
  - Mobile Genetic Elements tab (mobileOG-db summary)
  - Compare Isolates tab (pairwise + cluster detection)
  - Virulence gene expander in Sample Overview
  - fastANI species-confirmation columns in FASTA QC
  - Redesigned Summary tab (hero + KPI cards + donut + horizontal bar charts + tour cards)

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 3.0.0License: MIT
===================================================================
"""

import os
import sys
import json
import re
import io
import math
import glob
import argparse
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
from datetime import datetime
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

from bs4 import BeautifulSoup


def esc(v) -> str:
    """Escape a value for safe HTML embedding."""
    if v is None:
        return ""
    return (str(v).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


# =============================================================================
# PARSER
# =============================================================================
class PseudomonasHTMLParser:
    """Parses PseudoScope module outputs into structured dictionaries."""

    def __init__(self):
        self.abricate_databases = [
            'ncbi', 'card', 'resfinder', 'vfdb', 'argannot',
            'plasmidfinder', 'megares', 'ecoli_vf', 'bacmet2'
        ]

    def normalize_sample_id(self, sample_id: str) -> str:
        sample = str(sample_id)
        for ext in ('.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv'):
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        return sample.strip()

    def parse_html_table(self, html_content: str, table_index: int = 0) -> pd.DataFrame:
        try:
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if not tables or table_index >= len(tables):
                return pd.DataFrame()
            table = tables[table_index]
            rows = table.find_all('tr')
            if not rows:
                return pd.DataFrame()
            headers = [th.get_text().strip() for th in rows[0].find_all(['th', 'td'])]
            data = []
            for row in rows[1:]:
                cols = row.find_all(['td', 'th'])
                if cols and len(cols) == len(headers):
                    data.append([c.get_text().strip() for c in cols])
            return pd.DataFrame(data, columns=headers) if data else pd.DataFrame()
        except Exception as e:
            print(f"  ⚠️ Table parsing error: {e}")
            return pd.DataFrame()

    def parse_qc_report(self, file_path: Path) -> Dict[str, Dict]:
        print(f"  🧬 Parsing FASTA QC: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            df = self.parse_html_table(html, 0)
            if df.empty:
                return {}
            sample_col = None
            for col in df.columns:
                if 'filename' in col.lower() or 'sample' in col.lower() or col == df.columns[0]:
                    sample_col = col
                    break
            if not sample_col:
                return {}
            results = {}
            for _, row in df.iterrows():
                sample_raw = row[sample_col]
                if not sample_raw:
                    continue
                sample = self.normalize_sample_id(sample_raw)
                qc_data = {}
                for col in df.columns:
                    if col == sample_col:
                        continue
                    val = row[col]
                    if pd.isna(val) or val == '' or val == 'ND':
                        qc_data[col] = 'ND'
                    else:
                        cleaned = str(val).replace('%', '').replace(',', '').strip()
                        try:
                            qc_data[col] = float(cleaned)
                        except Exception:
                            qc_data[col] = str(val)
                results[sample] = qc_data
            print(f"    ✓ Parsed {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing QC: {e}")
            return {}

    def parse_mlst_report(self, file_path: Path) -> Dict[str, Dict]:
        print(f"  🧬 Parsing MLST: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            mlst_table = None
            for table in tables:
                if table.find(string=re.compile(r'Sample|ST|Allele', re.I)):
                    mlst_table = table
                    break
            if not mlst_table:
                mlst_table = soup.find('table')
            if not mlst_table:
                return {}
            rows = mlst_table.find_all('tr')
            if len(rows) < 2:
                return {}
            headers = [cell.get_text().strip() for cell in rows[0].find_all(['th', 'td'])]
            data = []
            for row in rows[1:]:
                cols = row.find_all(['td', 'th'])
                if cols:
                    row_data = [col.get_text().strip() for col in cols]
                    if len(row_data) >= 2:
                        data.append(row_data)
            if not data:
                return {}
            df = pd.DataFrame(data)
            if len(df.columns) > len(headers):
                df = df.iloc[:, :len(headers)]
            df.columns = headers[:len(df.columns)]
            df.columns = [col.strip() for col in df.columns]
            if 'Sample' in df.columns:
                df['normalized_sample'] = df['Sample'].apply(self.normalize_sample_id)
            elif any('sample' in col.lower() for col in df.columns):
                sample_col = [col for col in df.columns if 'sample' in col.lower()][0]
                df['normalized_sample'] = df[sample_col].apply(self.normalize_sample_id)
            else:
                df['normalized_sample'] = df.iloc[:, 0].apply(self.normalize_sample_id)
            results = {}
            for _, row in df.iterrows():
                sample = row['normalized_sample']
                st = 'ND'
                if 'ST' in df.columns:
                    st_raw = str(row['ST']) if pd.notna(row['ST']) else 'ND'
                    st = self._clean_st(st_raw)
                elif 'st' in [col.lower() for col in df.columns]:
                    st_col = [col for col in df.columns if col.lower() == 'st'][0]
                    st_raw = str(row[st_col]) if pd.notna(row[st_col]) else 'ND'
                    st = self._clean_st(st_raw)
                allele_profile = 'ND'
                if 'Allele Profile' in df.columns:
                    allele_profile = str(row['Allele Profile']) if pd.notna(row['Allele Profile']) else 'ND'
                results[sample] = {'ST': st, 'Allele_Profile': allele_profile}
            print(f"    ✓ Found {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing MLST: {e}")
            return {}

    def _clean_st(self, st_str: str) -> str:
        if not st_str or st_str == 'ND':
            return 'ND'
        st_str = str(st_str).strip()
        if not st_str.upper().startswith('ST'):
            st_str = 'ST' + st_str
        return st_str

    def parse_past_report(self, file_path: Path) -> Dict[str, Dict]:
        print(f"  🧬 Parsing PAST serotype: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            df = self.parse_html_table(html_content, 0)
            if df.empty:
                return {}
            df.columns = [col.strip() for col in df.columns]
            sample_col = None
            for col in df.columns:
                if 'sample' in col.lower() or 'id' in col.lower():
                    sample_col = col
                    break
            if not sample_col:
                sample_col = df.columns[0]
            df['normalized_sample'] = df[sample_col].apply(self.normalize_sample_id)
            results = {}
            for _, row in df.iterrows():
                sample = row['normalized_sample']
                o_type = 'ND'
                if 'Type' in df.columns:
                    o_type = str(row['Type']) if pd.notna(row['Type']) else 'ND'
                elif 'O-type' in df.columns:
                    o_type = str(row['O-type']) if pd.notna(row['O-type']) else 'ND'
                results[sample] = {'O_Type': o_type}
            print(f"    ✓ Found {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing PAST: {e}")
            return {}

    def parse_amrfinder_report(self, file_path: Path) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        print(f"  🧬 Parsing AMRfinder: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return {}, {}
            gene_frequencies = {}
            df2 = self.parse_html_table(str(tables[1]), 0)
            if not df2.empty and 'Gene' in df2.columns:
                for _, row in df2.iterrows():
                    gene = str(row['Gene']).strip()
                    frequency = str(row.get('Frequency', '0')).strip()
                    genomes = []
                    if 'Genomes' in df2.columns and pd.notna(row.get('Genomes')):
                        genomes_str = str(row['Genomes'])
                        genomes = [self.normalize_sample_id(g.strip())
                                   for g in genomes_str.split(',') if g.strip()]
                    count = 0
                    match = re.search(r'(\d+)', frequency)
                    if match:
                        count = int(match.group(1))
                    gene_frequencies[gene] = {
                        'frequency': frequency, 'count': count,
                        'genomes': genomes, 'database': 'amrfinder'
                    }
            genes_by_genome = {}
            df1 = self.parse_html_table(str(tables[0]), 0)
            if not df1.empty:
                sample_col = None
                for col in df1.columns:
                    if 'genome' in col.lower() or 'sample' in col.lower():
                        sample_col = col
                        break
                if sample_col:
                    for _, row in df1.iterrows():
                        sample = self.normalize_sample_id(row[sample_col])
                        genes = []
                        genes_col = None
                        for col in df1.columns:
                            if 'genes' in col.lower() or 'detected' in col.lower():
                                genes_col = col
                                break
                        if genes_col and pd.notna(row.get(genes_col)):
                            gene_str = str(row[genes_col])
                            genes = [g.strip() for g in gene_str.split() if g.strip()]
                        genes_by_genome[sample] = genes
            print(f"    ✓ Found {len(genes_by_genome)} samples, {len(gene_frequencies)} genes")
            return genes_by_genome, gene_frequencies
        except Exception as e:
            print(f"    ❌ Error parsing AMRfinder: {e}")
            return {}, {}

    def parse_abricate_report(self, file_path: Path) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        print(f"  🧬 Parsing ABRicate: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return {}, {}
            db_name = 'unknown'
            for db in self.abricate_databases:
                if db in file_path.name.lower():
                    db_name = db
                    break
            gene_frequencies = {}
            df2 = self.parse_html_table(str(tables[1]), 0)
            if not df2.empty and 'Gene' in df2.columns:
                for _, row in df2.iterrows():
                    gene = str(row['Gene']).strip()
                    frequency = str(row.get('Frequency', '0')).strip()
                    genomes = []
                    if 'Genomes' in df2.columns and pd.notna(row.get('Genomes')):
                        genomes_str = str(row['Genomes'])
                        genomes = [self.normalize_sample_id(g.strip())
                                   for g in genomes_str.split(',') if g.strip()]
                    count = 0
                    match = re.search(r'(\d+)', frequency)
                    if match:
                        count = int(match.group(1))
                    gene_frequencies[gene] = {
                        'frequency': frequency, 'count': count,
                        'genomes': genomes, 'database': db_name
                    }
            genes_by_genome = {}
            df1 = self.parse_html_table(str(tables[0]), 0)
            if not df1.empty:
                sample_col = None
                for col in df1.columns:
                    if 'genome' in col.lower() or 'sample' in col.lower():
                        sample_col = col
                        break
                if sample_col:
                    for _, row in df1.iterrows():
                        sample = self.normalize_sample_id(row[sample_col])
                        genes = []
                        genes_col = None
                        for col in df1.columns:
                            if 'genes' in col.lower() or 'detected' in col.lower():
                                genes_col = col
                                break
                        if genes_col and pd.notna(row.get(genes_col)):
                            gene_str = str(row[genes_col])
                            genes = [g.strip() for g in gene_str.split(',') if g.strip()]
                        genes_by_genome[sample] = genes
            print(f"    ✓ {db_name}: {len(genes_by_genome)} samples, {len(gene_frequencies)} genes")
            return genes_by_genome, gene_frequencies
        except Exception as e:
            print(f"    ❌ Error parsing ABRicate report: {e}")
            return {}, {}

    def parse_mutation_summary_html(self, file_path: Path) -> Dict[str, Any]:
        print(f"  🧬 Parsing mutation summary: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            mutation_table = None
            for table in soup.find_all('table'):
                if (table.find(string=re.compile(r'Gene', re.I)) and
                        table.find(string=re.compile(r'Mutation', re.I))):
                    mutation_table = table
                    break
            if not mutation_table:
                return {}
            header_row = None
            thead = mutation_table.find('thead')
            if thead:
                header_row = thead.find('tr')
            if not header_row:
                header_row = mutation_table.find('tr')
            if not header_row:
                return {}
            headers = [c.get_text().strip()
                       for c in header_row.find_all(['th', 'td'])
                       if c.get_text().strip()]
            col_idx = {}
            for idx, h in enumerate(headers):
                hl = h.lower()
                if 'gene' in hl: col_idx['gene'] = idx
                elif 'mutation' in hl: col_idx['mutation'] = idx
                elif 'count' in hl: col_idx['count'] = idx
                elif 'genome' in hl: col_idx['genomes'] = idx
                elif 'class' in hl: col_idx['class'] = idx
            for req in ['gene', 'mutation', 'count', 'genomes']:
                if req not in col_idx:
                    return {}
            tbody = mutation_table.find('tbody')
            rows = tbody.find_all('tr') if tbody else mutation_table.find_all('tr')[1:]
            mutations_list, genome_counts = [], defaultdict(int)
            for row in rows:
                cells = row.find_all('td')
                if len(cells) <= max(col_idx.values()):
                    continue
                gene = cells[col_idx['gene']].get_text().strip()
                mutation = cells[col_idx['mutation']].get_text().strip()
                count_str = cells[col_idx['count']].get_text().strip()
                m = re.search(r'(\d+)', count_str)
                count = int(m.group(1)) if m else 0
                genomes = [g.strip() for g in
                           cells[col_idx['genomes']].get_text().strip().split(',') if g.strip()]
                if not genomes:
                    continue
                for g in genomes:
                    genome_counts[g] += 1
                class_name = cells[col_idx['class']].get_text().strip() if 'class' in col_idx else ''
                mutations_list.append({
                    'gene': gene, 'mutation': mutation,
                    'class': class_name, 'subclass': '',
                    'count': count, 'genomes': genomes
                })
            mutations_list.sort(key=lambda x: x['count'], reverse=True)
            print(f"    ✓ {len(mutations_list)} mutations across {len(genome_counts)} genomes")
            return {
                'mutations': mutations_list,
                'genome_mutation_counts': dict(genome_counts)
            }
        except Exception as e:
            print(f"    ❌ Error parsing mutation summary: {e}")
            return {}

    def parse_mge_summary_html(self, file_path: Path) -> Dict[str, Any]:
        print(f"  🧬 Parsing MGE summary: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')

            stats = {}
            for label in ('Total Samples', 'With mobileOG hits', 'Without mobileOG hits',
                          'Total mobileOG hits', 'Key MGE-associated genes', 'Runtime'):
                for tag in soup.find_all(string=re.compile(re.escape(label), re.I)):
                    parent = tag.find_parent()
                    if not parent:
                        continue
                    text_block = parent.get_text(" ", strip=True)
                    m = re.search(rf'{re.escape(label)}\s*[:]?\s*([\d,\.]+s?)\b',
                                  text_block, re.I)
                    if m:
                        stats[label] = m.group(1)
                        break
                    for sib in parent.find_all_next(limit=3):
                        s = sib.get_text(strip=True)
                        if re.fullmatch(r'[\d,\.]+s?', s):
                            stats[label] = s
                            break
                    if label in stats:
                        break

            per_sample = {}
            mge_table = None
            for table in soup.find_all('table'):
                header_text = ' '.join(th.get_text().strip() for th in table.find_all('th'))
                if 'Sample' in header_text and 'mobileOG' in header_text:
                    mge_table = table
                    break
            if mge_table:
                try:
                    df = pd.read_html(io.StringIO(str(mge_table)))[0]
                    df.columns = [str(c).strip() for c in df.columns]
                    for _, row in df.iterrows():
                        sample = self.normalize_sample_id(row.get('Sample', ''))
                        if not sample:
                            continue
                        entry = {}
                        for col in df.columns:
                            if col == 'Sample':
                                continue
                            v = row.get(col, 0)
                            try:
                                entry[col] = int(str(v).replace(',', '').strip())
                            except Exception:
                                entry[col] = v
                        per_sample[sample] = entry
                except Exception as e:
                    print(f"    Error parsing MGE table: {e}")

            note = ""
            for tag in soup.find_all(string=re.compile(r'Interpretation note', re.I)):
                parent = tag.find_parent()
                if parent:
                    raw = parent.get_text(" ", strip=True)
                    raw = re.sub(r'^Interpretation note\s*:?\s*', '', raw, flags=re.I)
                    if len(raw) > 40:
                        note = raw
                        break

            print(f"    ✓ MGE: {len(stats)} cards, {len(per_sample)} samples")
            return {'stats': stats, 'per_sample': per_sample, 'note': note}
        except Exception as e:
            print(f"    ❌ Error parsing MGE summary: {e}")
            return {'stats': {}, 'per_sample': {}, 'note': ''}


# =============================================================================
# ANALYZER
# =============================================================================
class PseudomonasDataAnalyzer:
    """Cross-genome patterns, gene-centric tables, T3SS/carbapenemase derivation.
    Hardened to handle disparate gene naming conventions across CARD, NCBI, and ResFinder.
    """

    def __init__(self):
        # 1. True acquired carbapenemases
        self.carbapenemase_families = {
            'KPC': ['kpc'],
            'NDM': ['ndm'],
            'VIM': ['vim'],
            'IMP': ['imp'],
            'GES': ['ges'],
            'SPM': ['spm'],
            'OXA-23-like': ['oxa-23'],
            'OXA-24/40-like': ['oxa-24', 'oxa-40'],
            'OXA-48-like': ['oxa-48', 'oxa-181', 'oxa-232', 'oxa-204', 'oxa-244'],
            'OXA-58-like': ['oxa-58'],
        }

        # 2. Class D non-carbapenemases / ESBL variants
        self.oxa_class_d_esbl_families = {
            'OXA-2-like': ['oxa-2', 'oxa-3', 'oxa-15', 'oxa-20',
                           'oxa-21', 'oxa-32', 'oxa-34', 'oxa-45'],
            'OXA-10-like': ['oxa-10', 'oxa-11', 'oxa-13', 'oxa-14',
                            'oxa-16', 'oxa-17', 'oxa-19', 'oxa-28',
                            'oxa-35', 'oxa-74', 'oxa-101', 'oxa-128'],
        }

        # 3. Intrinsic chromosomal variants
        self.intrinsic_oxa_families = {
            'OXA-50-like': ['oxa-50', 'oxa-395', 'oxa-396'],
        }

        # 4. Non-OXA ESBLs & other β-lactamase classes
        self.other_esbl_families = {
            'PER': ['per'], 'VEB': ['veb'], 'BEL': ['bel'],
            'CTX-M': ['ctx-m'], 'TEM': ['tem'], 'SHV': ['shv'],
        }

        # 5. Colistin resistance
        self.colistin_resistance = {f'mcr-{i}' for i in range(1, 11)}

        # 6. 16S rRNA methyltransferases
        self.methyltransferases = {
            'arma', 'rmta', 'rmtb', 'rmtc', 'rmtd', 'rmte', 'rmtf', 'rmtg', 'rmth',
        }

        # ---- Prefix patterns for fallback "unclassified variant" detection ----
    
        self.unclassified_prefixes = {
            'oxa': 'OXA variant (unclassified)',
            'kpc': 'KPC variant (unclassified)',
            'ndm': 'NDM variant (unclassified)',
            'vim': 'VIM variant (unclassified)',
            'imp': 'IMP variant (unclassified)',
            'ges': 'GES variant (unclassified)',
            'spm': 'SPM variant (unclassified)',
            'per': 'PER variant (unclassified)',
            'veb': 'VEB variant (unclassified)',
            'bel': 'BEL variant (unclassified)',
            'ctx-m': 'CTX-M variant (unclassified)',
            'tem': 'TEM variant (unclassified)',
            'shv': 'SHV variant (unclassified)',
        }

        # Auto-compiled critical AMR gene set (used to ⚠️-mark critical rows)
        self.critical_amr_genes = set()
        for family_dict in (self.carbapenemase_families,
                            self.oxa_class_d_esbl_families,
                            self.intrinsic_oxa_families,
                            self.other_esbl_families):
            for sublist in family_dict.values():
                self.critical_amr_genes.update(sublist)
        self.critical_amr_genes.update(self.colistin_resistance)
        self.critical_amr_genes.update(self.methyltransferases)

        # Virulence signatures
        self.critical_virulence_genes = {
            'exou', 'exos', 'exot', 'exoy', 'exoa', 'toxa',
            'plch', 'plcn', 'lasb', 'lasa', 'apra',
        }
        self.t3ss_effectors = {'exou', 'exos', 'exot', 'exoy'}

    # ------------------------------------------------------------------
    # Gene-name normalization
    # ------------------------------------------------------------------
    def normalize_gene_name(self, raw_name: str) -> str:
        """Standardize gene tokens across NCBI, CARD, ResFinder."""
        name = str(raw_name).lower().strip()
        name = name.replace("beta-lactamase ", "")
        if name.startswith("bla"):
            name = name[3:]
        tokens = re.split(r'[_/]', name)
        core_token = tokens[0] if tokens else name
        return core_token

    def match_gene(self, input_gene: str) -> str:
        """Return the clean matched family token, or empty string if none."""
        clean = self.normalize_gene_name(input_gene)
        if clean in self.critical_amr_genes:
            return clean
        for baseline in self.critical_amr_genes:
            if clean.startswith(baseline + "-") or clean == baseline:
                return baseline
        return ""

    def classify_carbapenemase_family(self, raw_gene: str) -> Tuple[str, str]:
        """Return (family_name, category) or ('', '') if no match.
        Category ∈ {'acquired_carbapenemase', 'oxa_esbl', 'intrinsic_oxa',
                    'other_esbl', 'colistin', 'methyltransferase',
                    'unclassified', 'ambiguous_class_level'}.
        """
        clean = self.normalize_gene_name(raw_gene)
        if not clean:
            return '', ''

        # 1. Curated acquired carbapenemases
        for family, patterns in self.carbapenemase_families.items():
            if any(clean == p or clean.startswith(p + '-') or clean.startswith(p) for p in patterns):
                return family, 'acquired_carbapenemase'

        # 2. Curated Class D ESBL variants
        for family, patterns in self.oxa_class_d_esbl_families.items():
            if any(clean == p or clean.startswith(p + '-') or clean.startswith(p) for p in patterns):
                return family, 'oxa_esbl'

        # 3. Curated intrinsic OXA-50-like
        for family, patterns in self.intrinsic_oxa_families.items():
            if any(clean == p or clean.startswith(p + '-') or clean.startswith(p) for p in patterns):
                return family, 'intrinsic_oxa'

        # 4. Non-OXA ESBLs
        for family, patterns in self.other_esbl_families.items():
            if any(clean == p or clean.startswith(p + '-') or clean.startswith(p) for p in patterns):
                return family, 'other_esbl'

        # 5. Colistin
        for pattern in self.colistin_resistance:
            if clean == pattern or clean.startswith(pattern):
                return pattern, 'colistin'

        # 6. Methyltransferases
        for pattern in self.methyltransferases:
            if clean == pattern or clean.startswith(pattern):
                return pattern, 'methyltransferase'

        # 7. Fallback 1 — ambiguous class-level hit (e.g., bare 'oxa' from MEGARES)
        #    Only fires when the token is exactly the family prefix with no numeric suffix.
        for prefix in self.unclassified_prefixes:
            if clean == prefix:
                return f'{prefix.upper()} (class-level)', 'ambiguous_class_level'

        # 8. Fallback 2 — unclassified variant (e.g., oxa-904, oxa-1018, kpc-99)
        #    Detected by prefix + numeric suffix pattern. Preserved for manual review.
        m = re.match(r'^([a-z]+(?:-[a-z]+)?)-(\d+)$', clean)
        if m:
            prefix = m.group(1)
            if prefix in self.unclassified_prefixes:
                return f'{self.unclassified_prefixes[prefix]} — {clean}', 'unclassified'

        return '', ''

    # ------------------------------------------------------------------
    # Derivations
    # ------------------------------------------------------------------
    @staticmethod
    def derive_t3ss_pathotype(all_genes: List[str]) -> str:
        genes_lower = [str(g).lower() for g in all_genes]
        has_exou = any('exou' in g for g in genes_lower)
        has_exos = any('exos' in g for g in genes_lower)
        if has_exou and has_exos:
            return 'ExoU+ExoS'
        if has_exou:
            return 'ExoU'
        if has_exos:
            return 'ExoS'
        return 'Neither'

    # ------------------------------------------------------------------
    # Gene-centric tables
    # ------------------------------------------------------------------
    def create_gene_centric_tables(self, integrated_data: Dict[str, Any]) -> Dict[str, Any]:
        gene_centric = {
            'amr_databases': {},
            'virulence_databases': {},
            'plasmid_databases': {},
            'bacmet_databases': {},
            'combined_gene_frequencies': []
        }
        if 'amrfinder' in integrated_data.get('gene_frequencies', {}):
            amr_data = integrated_data['gene_frequencies']['amrfinder']
            gene_list = [{'gene': g, 'database': 'AMRfinder',
                          'frequency': d.get('frequency', '0'),
                          'count': d.get('count', 0),
                          'genomes': d.get('genomes', [])}
                         for g, d in amr_data.items()]
            if gene_list:
                gene_centric['amr_databases']['amrfinder'] = sorted(
                    gene_list, key=lambda x: x['count'], reverse=True)

        if 'abricate' in integrated_data.get('gene_frequencies', {}):
            for db_name, db_genes in integrated_data['gene_frequencies']['abricate'].items():
                gene_list = [{'gene': g, 'database': db_name.upper(),
                              'frequency': d.get('frequency', '0'),
                              'count': d.get('count', 0),
                              'genomes': d.get('genomes', [])}
                             for g, d in db_genes.items()]
                if not gene_list:
                    continue
                gene_list.sort(key=lambda x: x['count'], reverse=True)
                if db_name in ('vfdb', 'ecoli_vf', 'pa_vf', 'pseudomonas_vf'):
                    gene_centric['virulence_databases'][db_name] = gene_list
                elif db_name == 'plasmidfinder':
                    gene_centric['plasmid_databases'][db_name] = gene_list
                elif db_name == 'bacmet2':
                    gene_centric['bacmet_databases'][db_name] = gene_list
                else:
                    gene_centric['amr_databases'][db_name] = gene_list

        all_genes = []
        for db_type in ('amr_databases', 'virulence_databases',
                        'plasmid_databases', 'bacmet_databases'):
            for genes in gene_centric.get(db_type, {}).values():
                all_genes.extend(genes)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        gene_centric['combined_gene_frequencies'] = all_genes
        return gene_centric

    # ------------------------------------------------------------------
    # Cross-genome patterns
    # ------------------------------------------------------------------
    def create_cross_genome_patterns(self, integrated_data: Dict[str, Any]) -> Dict[str, Any]:
        samples_data = integrated_data.get('samples', {})

        P = {
            'mlst_distribution': Counter(),
            'o_type_distribution': Counter(),
            't3ss_distribution': Counter(),
            'st_to_samples': defaultdict(list),
            'o_type_to_samples': defaultdict(list),
            't3ss_to_samples': defaultdict(list),
            'gene_cooccurrence': defaultdict(Counter),
            'high_risk_combinations': [],
            # Three-level structure: sample -> category -> family -> {raw_gene: [databases]}
            'carbapenemase_by_sample': defaultdict(
                lambda: defaultdict(lambda: defaultdict(dict))
            ),
            'carbapenemase_family_distribution': Counter(),
            'carbapenemase_category_distribution': Counter(),
        }

        COMBO_FIELDS = {
            'mlst_otype':       ['MLST', 'O_Type'],
            'mlst_t3ss':        ['MLST', 'T3SS'],
            'otype_t3ss':       ['O_Type', 'T3SS'],
            'mlst_otype_t3ss':  ['MLST', 'O_Type', 'T3SS'],
        }
        combo_dicts = {k: defaultdict(list) for k in COMBO_FIELDS}

        sample_genes = defaultdict(list)
        gene_centric = integrated_data.get('gene_centric', {})
        for db_type in ('amr_databases', 'virulence_databases',
                        'plasmid_databases', 'bacmet_databases'):
            for genes in gene_centric.get(db_type, {}).values():
                for g in genes:
                    for genome in g['genomes']:
                        if g['gene'] not in sample_genes[genome]:
                            sample_genes[genome].append(g['gene'])

        # ---- Build per-sample carbapenemase classification from ALL AMR sources ----
        for sample, data in samples_data.items():
            seen_hits_per_sample = set()

            sources = []
            for raw in data.get('amrfinder', {}).get('all_genes', []):
                sources.append((raw, 'AMRfinder'))
            for db_name, db_genes in data.get('abricate_databases', {}).items():
                for raw in db_genes:
                    sources.append((raw, db_name.upper()))

            for raw_g, db_name in sources:
                result = self.classify_carbapenemase_family(raw_g)
                # Defensive: support both 2-tuple and legacy string return
                if isinstance(result, tuple) and len(result) == 2:
                    family, category = result
                else:
                    family, category = str(result), ''
                if not family or not category:
                    continue

                dedup_key = (family, category, raw_g)
                if dedup_key in seen_hits_per_sample:
                    continue
                seen_hits_per_sample.add(dedup_key)

                bucket = P['carbapenemase_by_sample'][sample][category][family]
                if raw_g not in bucket:
                    bucket[raw_g] = []
                if db_name not in bucket[raw_g]:
                    bucket[raw_g].append(db_name)

                if category == 'acquired_carbapenemase':
                    P['carbapenemase_family_distribution'][family] += 1
                P['carbapenemase_category_distribution'][category] += 1

        # ---- Per-sample typing / combos / high-risk ----
        for sample, data in samples_data.items():
            typing = data.get('typing', {})
            values = {
                'MLST':   typing.get('MLST', 'ND'),
                'O_Type': typing.get('O_Type', 'ND'),
                'T3SS':   typing.get('T3SS', 'ND'),
            }

            def ok(v):
                return v and v not in ('Not Assigned', 'ND', '', 'nan')

            if ok(values['MLST']):
                P['mlst_distribution'][values['MLST']] += 1
                P['st_to_samples'][values['MLST']].append(sample)
            if ok(values['O_Type']):
                P['o_type_distribution'][values['O_Type']] += 1
                P['o_type_to_samples'][values['O_Type']].append(sample)
            if ok(values['T3SS']):
                P['t3ss_distribution'][values['T3SS']] += 1
                P['t3ss_to_samples'][values['T3SS']].append(sample)

            for key, fields in COMBO_FIELDS.items():
                if all(ok(values[f]) for f in fields):
                    label = ' - '.join(values[f] for f in fields)
                    combo_dicts[key][label].append(sample)

            amr_genes = data.get('amrfinder', {}).get('all_genes', [])
            vir_genes = []
            for vir_db in ('vfdb', 'ecoli_vf', 'pa_vf', 'pseudomonas_vf'):
                vir_genes.extend(data.get('abricate_databases', {}).get(vir_db, []))
            critical_amr = [g for g in amr_genes if self.match_gene(g)]
            critical_vir = [g for g in vir_genes
                            if any(c in str(g).lower() for c in self.critical_virulence_genes)]
            if critical_amr and critical_vir:
                P['high_risk_combinations'].append({
                    'sample': sample,
                    'mlst': values['MLST'],
                    'o_type': values['O_Type'],
                    't3ss': values['T3SS'],
                    'critical_amr_genes': critical_amr,
                    'critical_virulence_genes': critical_vir,
                })

            genes = sample_genes.get(sample, [])
            for i, g1 in enumerate(genes):
                for g2 in genes[i + 1:]:
                    P['gene_cooccurrence'][g1][g2] += 1

        for key in combo_dicts:
            P[key] = dict(combo_dicts[key])
        return P

    # ------------------------------------------------------------------
    # MGE aggregation
    # ------------------------------------------------------------------
    def compute_mge_stats(self, mge_data: Dict[str, Any]) -> Dict[str, Any]:
        per_sample = mge_data.get('per_sample', {})
        if not per_sample:
            return {'category_totals': [], 'grand_total': 0, 'n_samples': 0}
        all_columns = set()
        for s in per_sample.values():
            all_columns.update(s.keys())
        numeric_cols = [c for c in all_columns
                        if all(isinstance(s.get(c, 0), (int, float))
                               for s in per_sample.values())]
        category_totals = []
        n = len(per_sample)
        grand_total = 0
        for col in numeric_cols:
            vals = {s: per_sample[s].get(col, 0) for s in per_sample}
            total = sum(vals.values())
            grand_total += total
            max_sample = max(vals, key=vals.get) if vals else 'ND'
            category_totals.append({
                'category': col,
                'total': total,
                'mean': round(total / n, 2) if n else 0,
                'max': vals[max_sample] if vals else 0,
                'max_sample': max_sample,
            })
        category_totals.sort(key=lambda x: x['total'], reverse=True)
        for c in category_totals:
            c['pct_of_all'] = round((c['total'] / grand_total) * 100, 2) if grand_total else 0
        return {
            'category_totals': category_totals,
            'grand_total': grand_total,
            'n_samples': n,
        }


# =============================================================================
# HTML GENERATOR
# =============================================================================
class PseudomonasHTMLGenerator:
    """Builds the interactive multi-tab HTML report."""

    def __init__(self, analyzer: PseudomonasDataAnalyzer):
        self.analyzer = analyzer
        self._current_samples_data = {}
        self.tab_colors = {
            'summary': '#8B4513', 'samples': '#2E7D32', 'qc': '#00695C',
            'mlst': '#1565C0', 'serotype': '#6A1B9A', 't3ss': '#D32F2F',
            'carb': '#B71C1C', 'amr': '#C62828', 'virulence': '#E65100',
            'plasmid': '#3F51B5', 'bacmet': '#607D8B', 'mutation': '#00BCD4',
            'mge': '#16A085', 'patterns': '#AD1457', 'compare': '#0d9488',
            'citation': '#00838F', 'guide': '#4527A0',
            'calltoaction': '#2E7D32', 'export': '#37474F',
        }

    # ------------------------------------------------------------------
    # Reusable helpers
    # ------------------------------------------------------------------
    def _credit_bar(self, color: str, icon: str, title: str, body: str) -> str:
        return f'''
        <div class="scientific-note" style="background: linear-gradient(135deg, #fdf6ec 0%, #f5e6d3 100%); border-left: 6px solid {color}; margin-bottom: 20px; padding: 15px; border-radius: 8px;">
            <div style="display: flex; align-items: center; gap: 12px; flex-wrap: wrap;">
                <span style="font-size: 1.4em;">{icon}</span>
                <div>
                    <strong style="font-size: 1.1em; color: {color};">{title}</strong><br>
                    <span style="font-size: 0.95em; color: #333;">{body}</span>
                </div>
            </div>
        </div>'''

    def _alert(self, kind: str, icon: str, html_body: str) -> str:
        return f'''
        <div class="alert-box alert-{kind}">
            <i class="fas {icon} fa-2x"></i>
            <div>{html_body}</div>
        </div>'''

    def _stat_card(self, value, label: str, color: str = '#8B4513', icon: str = '') -> str:
        return f'''
        <div class="stat-card" style="background: linear-gradient(135deg, {color} 0%, {color}dd 100%);">
            {f'<i class="fas {icon} fa-2x" style="opacity:0.9; margin-bottom:8px;"></i>' if icon else ''}
            <div class="stat-value">{value}</div>
            <div class="stat-label">{label}</div>
        </div>'''

    def _filter_buttons(self, table_id: str, buttons: list) -> str:
        html = f'''<div class="action-buttons">
            <button class="action-btn btn-primary"
                onclick="exportTableToCSV('{table_id}', '{table_id}.csv')">
                <i class="fas fa-download"></i> Export</button>'''
        for btn in buttons:
            label, search_val = btn[0], btn[1]
            css_class = btn[2] if len(btn) > 2 else 'btn-info'
            icon = btn[3] if len(btn) > 3 else 'fa-filter'
            html += (f'''<button class="action-btn {css_class}"
                onclick="document.getElementById('search-{table_id}').value='{search_val}';'''
                     f'''searchTable('{table_id}','search-{table_id}')">
                <i class="fas {icon}"></i> {label}</button>''')
        html += (f'''<button class="action-btn btn-light"
            onclick="document.getElementById('search-{table_id}').value='';'''
                 f'''searchTable('{table_id}','search-{table_id}')">
            <i class="fas fa-sync"></i> Clear</button></div>''')
        return html

    def _gene_family_info(self, border_color: str, title: str, items: list) -> str:
        html = (f'<div style="margin:10px 0 20px 0;background:#fdf6ec;padding:15px;'
                f'border-radius:8px;font-size:.9em;border-left:4px solid {border_color};">'
                f'<strong><i class="fas fa-info-circle"></i> {title}</strong><br>')
        for name, desc in items:
            html += f'• <strong>{name}</strong> – {desc}<br>'
        html += '</div>'
        return html

    def _database_cards(self, gene_dict: dict, db_labels: dict = None) -> str:
        db_labels = db_labels or {}
        html = ('<h3 style="margin-top:30px;"><i class="fas fa-database"></i> '
                'Database Summary</h3>'
                '<div style="display:grid;grid-template-columns:'
                'repeat(auto-fit,minmax(300px,1fr));gap:20px;margin:20px 0;">')
        for db, genes in gene_dict.items():
            label = db_labels.get(db, db.upper() if db != 'amrfinder' else 'AMRfinder')
            top = ', '.join(f"{g['gene']} ({g['count']})" for g in genes[:3])
            total = sum(g['count'] for g in genes)
            html += f'''<div class="database-section">
                <h4>{label}</h4>
                <p><strong>{len(genes)} unique genes</strong>
                (Total occurrences: {total})</p>
                <p>Top genes: {top}</p></div>'''
        html += '</div>'
        return html

    def _grouping_dropdown(self, table_id: str) -> str:
        return f'''
        <div class="grouping-controls">
            <strong><i class="fas fa-layer-group"></i> Group genomes by:</strong>
            <select class="group-select" onchange="groupGenomesByTyping('{table_id}', this.value)">
                <option value="">— None (flat list) —</option>
                <optgroup label="Single field">
                    <option value="MLST">MLST (ST)</option>
                    <option value="Otype">O-type</option>
                    <option value="T3SS">ExoU/ExoS Pathotype</option>
                </optgroup>
                <optgroup label="Two fields">
                    <option value="MLST-Otype">MLST + O-type</option>
                    <option value="MLST-T3SS">MLST + T3SS</option>
                    <option value="Otype-T3SS">O-type + T3SS</option>
                </optgroup>
                <optgroup label="Three fields">
                    <option value="MLST-Otype-T3SS">MLST + O-type + T3SS</option>
                </optgroup>
            </select>
            <button class="group-btn" onclick="resetGenomeList('{table_id}')"><i class="fas fa-undo"></i> Reset</button>
        </div>'''

    def _search_pair(self, table_id: str, placeholder: str) -> str:
        return f'''
        <input type="text" class="search-box" id="search-{table_id}"
               onkeyup="searchTable('{table_id}', 'search-{table_id}')"
               placeholder="🔍 {placeholder}">
        <input type="text" class="search-box" id="highlight-{table_id}"
               onkeyup="highlightGenome('{table_id}', 'highlight-{table_id}')"
               placeholder="🔍 Highlight genomes containing specific text">'''

    def _combo_table(self, table_id: str, title: str, combo_dict: dict,
                     sample_key: str = 'Samples') -> str:
        if not combo_dict:
            return ''
        rows = ''
        for combo, samples in sorted(combo_dict.items(),
                                     key=lambda x: len(x[1]), reverse=True):
            tags = ''.join(f'<span class="genome-tag">{esc(s)}</span>' for s in samples)
            rows += (f'<tr><td><strong>{esc(combo)}</strong></td>'
                     f'<td>{len(samples)}</td>'
                     f'<td><div class="genome-list">{tags}</div></td></tr>')
        return f'''
        <h3>🔗 {title}</h3>
        {self._search_pair(table_id, f'Search {title}...')}
        <div class="master-scrollable-container">
            <table id="{table_id}" class="data-table">
                <thead><tr>
                    <th data-sort="string">{title}</th>
                    <th data-sort="number">Count</th>
                    <th data-sort="string">{sample_key}</th>
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    def _gene_table(self, table_id: str, title: str, gene_list: list,
                    total_samples: int, critical_set: set = None) -> str:
        critical_set = critical_set or set()
        rows = ''
        for g in gene_list:
            norm = self.analyzer.normalize_gene_name(g['gene'])
            is_crit = any(c in norm for c in critical_set) or norm in critical_set
            gene_display = f"<strong>{esc(g['gene'])}</strong>" + (" ⚠️" if is_crit else "")
            freq = g.get('frequency', str(g['count']))
            pct = '0%'
            if '(' in freq:
                pct = freq.split('(')[-1].replace(')', '').strip()
            elif g['count'] > 0 and total_samples > 0:
                pct = f"{(g['count'] / total_samples) * 100:.1f}%"
            tags = ''.join(f'<span class="genome-tag">{esc(x)}</span>' for x in g.get('genomes', []))
            rows += (f'<tr><td>{gene_display}</td><td>{esc(g["database"])}</td>'
                     f'<td><strong>{g["count"]}</strong></td><td>{pct}</td>'
                     f'<td><div class="genome-list">{tags}</div></td></tr>')
        return f'''
        <h3><i class="fas fa-list"></i> {title}</h3>
        {self._grouping_dropdown(table_id)}
        {self._search_pair(table_id, f'Search {title}...')}
        <div class="master-scrollable-container">
            <table id="{table_id}" class="data-table">
                <thead><tr>
                    <th data-sort="string">Gene</th>
                    <th data-sort="string">Database</th>
                    <th data-sort="number">Count</th>
                    <th data-sort="number">Percentage</th>
                    <th data-sort="string">Genomes (scrollable, groupable)</th>
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    def _distribution_table(self, table_id: str, dist: Counter,
                            label: str, extra_headers: list = None) -> str:
        extra_headers = extra_headers or []
        total = sum(dist.values())
        rows = ''
        for val, cnt in dist.most_common():
            if val in ('Not Assigned', 'ND', ''):
                continue
            pct = (cnt / total * 100) if total else 0
            extras = ''.join('<td>—</td>' for _ in extra_headers)
            rows += f'<tr><td><strong>{esc(val)}</strong></td><td>{cnt}</td><td>{pct:.1f}%</td>{extras}</tr>'
        header_extra = ''.join(f'<th data-sort="string">{h}</th>' for h in extra_headers)
        return f'''
        <div class="scrollable-table">
            <table id="{table_id}" class="data-table">
                <thead><tr>
                    <th data-sort="string">{label}</th>
                    <th data-sort="number">Count</th>
                    <th data-sort="number">Percentage</th>
                    {header_extra}
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    def _donut_chart(self, segments, size=210, thickness=34,
                     center_text='', center_sub=''):
        total = sum(v for _, v, _ in segments) or 1
        radius = (size - thickness) / 2
        circumference = 2 * math.pi * radius
        cx = cy = size / 2
        offset_acc = 0.0
        segments_svg = ''
        for label, value, color in segments:
            if value == 0:
                continue
            pct = value / total
            dash = pct * circumference
            gap = circumference - dash
            segments_svg += (
                f'<circle cx="{cx}" cy="{cy}" r="{radius}" fill="none" '
                f'stroke="{color}" stroke-width="{thickness}" '
                f'stroke-dasharray="{dash:.3f} {gap:.3f}" '
                f'stroke-dashoffset="{-offset_acc:.3f}" '
                f'transform="rotate(-90 {cx} {cy})"/>'
            )
            offset_acc += dash
        center = ''
        if center_text:
            center = (f'<text x="{cx}" y="{cy + 4}" text-anchor="middle" '
                      f'font-size="30" font-weight="800" fill="#5C2E0B">{center_text}</text>')
            if center_sub:
                center += (f'<text x="{cx}" y="{cy + 26}" text-anchor="middle" '
                           f'font-size="11" fill="#8B4513" font-weight="600">'
                           f'{center_sub}</text>')
        return (f'<svg viewBox="0 0 {size} {size}" width="{size}" height="{size}" '
                f'class="donut-svg" xmlns="http://www.w3.org/2000/svg">'
                f'<circle cx="{cx}" cy="{cy}" r="{radius}" fill="none" '
                f'stroke="#f5e6d3" stroke-width="{thickness}"/>'
                f'{segments_svg}{center}</svg>')

    def _horizontal_bars(self, items, total, show_pct_of_total=True):
        if not items:
            return '<p style="color:#888;font-style:italic;">No data available</p>'
        max_val = max(v for _, v, _ in items) or 1
        rows = ''
        for label, value, color in items:
            pct_of_max = (value / max_val) * 100
            pct_of_total = (value / total * 100) if total else 0
            pct_html = (f' <span class="hbar-pct">({pct_of_total:.0f}%)</span>'
                        if show_pct_of_total else '')
            rows += (f'<div class="hbar-row">'
                     f'<div class="hbar-label" title="{esc(label)}">{esc(label)}</div>'
                     f'<div class="hbar-track">'
                     f'<div class="hbar-fill" style="width:{pct_of_max:.1f}%;'
                     f'background:linear-gradient(90deg,{color}cc,{color});"></div>'
                     f'</div>'
                     f'<div class="hbar-value">{value}{pct_html}</div>'
                     f'</div>')
        return f'<div class="hbar-container">{rows}</div>'

    # ------------------------------------------------------------------
    # CSS
    # ------------------------------------------------------------------
    def _get_css(self) -> str:
        return """
<style>
:root{
  --brown:#8B4513;--brown-light:#A0522D;--brown-dark:#5C2E0B;
  --summary-color:#8B4513;--samples-color:#2E7D32;--qc-color:#00695C;
  --mlst-color:#1565C0;--serotype-color:#6A1B9A;--t3ss-color:#D32F2F;
  --carb-color:#B71C1C;--amr-color:#C62828;--virulence-color:#E65100;
  --plasmid-color:#3F51B5;--bacmet-color:#607D8B;--mutation-color:#00BCD4;
  --mge-color:#16A085;--patterns-color:#AD1457;--compare-color:#0d9488;
  --citation-color:#00838F;--guide-color:#4527A0;
  --calltoaction-color:#2E7D32;--export-color:#37474F;
}
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:'Segoe UI',Tahoma,Geneva,Verdana,sans-serif;line-height:1.6;color:#333;background:#f5f5f5;min-width:1200px}
.container{max-width:none;margin:0 auto;padding:20px;width:100%;overflow-x:auto}
.main-header{background:linear-gradient(135deg,#2A1408 0%,#5C2E0B 40%,#8B4513 100%);color:white;padding:30px;border-radius:15px;margin-bottom:30px;text-align:center;box-shadow:0 10px 30px rgba(0,0,0,.2)}
.main-header h1{font-size:2.8em;margin-bottom:10px;color:white}
.metadata-bar{background:rgba(255,255,255,.1);padding:15px;border-radius:10px;margin:20px 0;display:flex;justify-content:space-around;flex-wrap:wrap;gap:15px;backdrop-filter:blur(10px)}
.metadata-item{display:flex;align-items:center;gap:8px;font-size:.95em}
.dashboard-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:20px;margin-bottom:30px}
.dashboard-card{background:white;padding:20px;border-radius:12px;box-shadow:0 5px 20px rgba(0,0,0,.1);text-align:center;transition:all .3s ease;cursor:pointer;border-left:5px solid;position:relative;overflow:hidden}
.dashboard-card:hover{transform:translateY(-10px);box-shadow:0 15px 30px rgba(0,0,0,.2)}
.dashboard-card .card-number{font-size:2.5em;font-weight:bold;margin:10px 0;background:linear-gradient(90deg,#5C2E0B,#8B4513);-webkit-background-clip:text;-webkit-text-fill-color:transparent}
.dashboard-card .card-label{font-size:.9em;color:#555;font-weight:600}
.tab-navigation{display:flex;gap:5px;margin-bottom:20px;flex-wrap:wrap;background:white;padding:15px;border-radius:12px;box-shadow:0 5px 20px rgba(0,0,0,.1);position:sticky;top:10px;z-index:100}
.tab-button{padding:10px 16px;background:#f5f5f5;border:none;border-radius:8px;cursor:pointer;font-weight:600;color:#666;transition:all .3s ease;display:flex;align-items:center;gap:6px;font-size:.88em}
.tab-button.active{color:white}
.tab-button.summary.active{background:var(--summary-color)}
.tab-button.samples.active{background:var(--samples-color)}
.tab-button.qc.active{background:var(--qc-color)}
.tab-button.mlst.active{background:var(--mlst-color)}
.tab-button.serotype.active{background:var(--serotype-color)}
.tab-button.t3ss.active{background:var(--t3ss-color)}
.tab-button.carb.active{background:var(--carb-color)}
.tab-button.amr.active{background:var(--amr-color)}
.tab-button.virulence.active{background:var(--virulence-color)}
.tab-button.plasmid.active{background:var(--plasmid-color)}
.tab-button.bacmet.active{background:var(--bacmet-color)}
.tab-button.mutation.active{background:var(--mutation-color)}
.tab-button.mge.active{background:var(--mge-color)}
.tab-button.patterns.active{background:var(--patterns-color)}
.tab-button.compare.active{background:var(--compare-color)}
.tab-button.citation.active{background:var(--citation-color)}
.tab-button.guide.active{background:var(--guide-color)}
.tab-button.calltoaction.active{background:var(--calltoaction-color)}
.tab-button.export.active{background:var(--export-color)}
.tab-content{display:none;background:white;padding:30px;border-radius:15px;box-shadow:0 10px 30px rgba(0,0,0,.1);margin-bottom:30px;animation:fadeIn .5s ease;width:100%;overflow-x:auto}
.tab-content.active{display:block}
@keyframes fadeIn{from{opacity:0;transform:translateY(20px)}to{opacity:1;transform:translateY(0)}}
.section-header{color:#2c3e50;margin-bottom:25px;padding-bottom:15px;border-bottom:3px solid;font-size:1.8em;display:flex;align-items:center;justify-content:space-between}
.section-header span{display:flex;align-items:center;gap:12px}
.data-table{width:100%;border-collapse:collapse;margin:20px 0;font-size:.95em;box-shadow:0 2px 10px rgba(0,0,0,.1);border-radius:8px;overflow:hidden}
.data-table th{background:#5C2E0B;color:white;padding:14px;text-align:left;font-weight:600;position:sticky;top:0;white-space:nowrap;cursor:pointer}
.data-table th:hover{background:#3D1E08}
.data-table td{padding:12px;border-bottom:1px solid #e0e0e0;vertical-align:top}
.data-table tr:hover td{background:#fdf6ec}
.scrollable-table{max-height:600px;overflow-y:auto;border:1px solid #e0e0e0;border-radius:8px;margin:20px 0;width:100%}
.master-scrollable-container{width:100%;overflow-x:auto;border:1px solid #e0e0e0;border-radius:8px;margin:20px 0}
.genome-list{display:flex;flex-wrap:wrap;gap:6px;padding:8px;background:#fdf6ec;border-radius:5px;width:100%;max-height:280px;overflow:auto;box-sizing:border-box}
.genome-group{margin-bottom:10px;width:100%;box-sizing:border-box}
.genome-group-header{font-weight:bold;background:#E8C39E;padding:5px 10px;border-radius:4px;margin:6px 0;font-size:.9em;display:inline-block;max-width:100%;word-break:break-word}
.genome-group-tags{display:flex;flex-wrap:wrap;gap:6px;margin-left:10px;max-width:calc(100% - 10px)}
.genome-tag{display:inline-block;background:#e6ffe6;color:#5C2E0B;padding:4px 12px;border-radius:12px;font-size:.9em;border:1px solid #E8C39E;white-space:nowrap;margin:2px;max-width:100%}
.genome-tag.highlight{background:#ffff99 !important;color:#000 !important;border:1px solid #ffc107}
.search-box{width:100%;padding:12px;margin-bottom:12px;border:2px solid #e0e0e0;border-radius:8px;font-size:1em;transition:all .3s ease}
.search-box:focus{outline:none;border-color:#8B4513;box-shadow:0 0 0 3px rgba(139,69,19,.1)}
.badge{display:inline-block;padding:5px 15px;border-radius:20px;font-size:.85em;font-weight:600;margin:2px}
.badge-critical{background:#8B0000;color:white}
.badge-exou{background:#C62828;color:white}
.badge-exos{background:#1565C0;color:white}
.alert-box{padding:20px;border-radius:10px;margin:20px 0;display:flex;align-items:flex-start;gap:20px;border-left:5px solid}
.alert-success{background:#d4edda;color:#155724;border-left-color:#28a745}
.alert-warning{background:#fff3cd;color:#856404;border-left-color:#ffc107}
.alert-danger{background:#f5e6d3;color:#3D1E08;border-left-color:#8B4513}
.alert-info{background:#d1ecf1;color:#0c5460;border-left-color:#17a2b8}
.action-buttons{display:flex;gap:10px;margin:20px 0;flex-wrap:wrap}
.action-btn{padding:10px 20px;border:none;border-radius:8px;cursor:pointer;font-weight:600;display:flex;align-items:center;gap:8px;transition:all .3s ease;text-decoration:none}
.action-btn:hover{transform:translateY(-2px);box-shadow:0 5px 15px rgba(0,0,0,.2)}
.btn-primary{background:#8B4513;color:white}
.btn-success{background:#28a745;color:white}
.btn-danger{background:#C62828;color:white}
.btn-warning{background:#ffc107;color:black}
.btn-info{background:#17a2b8;color:white}
.btn-secondary{background:#6c757d;color:white}
.btn-light{background:#f8f9fa;color:#212529;border:1px solid #dee2e6}
.database-section{margin:30px 0;padding:25px;border-radius:12px;background:#fdf6ec;box-shadow:0 3px 15px rgba(0,0,0,.08)}
.print-section-btn{background:#8B4513;color:white;border:none;border-radius:5px;padding:8px 15px;cursor:pointer;display:flex;align-items:center;gap:5px;font-size:.9em}
.print-section-btn:hover{background:#A0522D}
.footer{text-align:center;padding:30px;color:white;margin-top:40px;border-radius:15px;background:linear-gradient(135deg,#2A1408 0%,#5C2E0B 100%)}
.footer a{color:#E8C39E;text-decoration:none}
.footer a:hover{text-decoration:underline}
.footer .star-cta{display:inline-block;background:linear-gradient(135deg,#ffc107,#ff9800);color:#2A1408;font-weight:800;padding:10px 22px;border-radius:24px;margin-top:12px;text-decoration:none;box-shadow:0 3px 10px rgba(0,0,0,.25);transition:transform .2s}
.footer .star-cta:hover{transform:translateY(-2px);text-decoration:none}
.sort-icon{margin-left:5px;font-size:.8em;opacity:.6}
.grouping-controls{background:#fdf6ec;padding:12px;border-radius:8px;margin:15px 0;display:flex;flex-wrap:wrap;gap:10px;align-items:center;border-left:4px solid #8B4513}
.grouping-controls strong{font-weight:bold;margin-right:5px;color:#5C2E0B}
.group-select{padding:8px 14px;border:2px solid #8B4513;border-radius:6px;background:white;color:#333;font-size:.9em;font-weight:600;cursor:pointer;min-width:240px}
.group-select:focus{outline:none;box-shadow:0 0 0 3px rgba(139,69,19,.15)}
.group-btn{background:white;border:1px solid #8B4513;color:#8B4513;padding:6px 12px;border-radius:20px;cursor:pointer;font-size:.85em;transition:all .2s}
.group-btn:hover{background:#8B4513;color:white}
.stat-card{color:white;padding:18px;border-radius:10px;text-align:center;box-shadow:0 4px 15px rgba(0,0,0,.15);transition:transform .2s}
.stat-card:hover{transform:translateY(-3px)}
.stat-card .stat-value{font-size:1.9em;font-weight:bold;margin-bottom:4px}
.stat-card .stat-label{font-size:.85em;opacity:.95;text-transform:uppercase;letter-spacing:.5px}
.stats-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:15px;margin:20px 0}
.typing-badge{display:inline-block;padding:3px 10px;border-radius:12px;font-size:.8em;font-weight:600;background:#e0e0e0;color:#333;border:1px solid #ccc}
.typing-badge.t3ss-exou{background:#C62828;color:white;border-color:#8B0000}
.typing-badge.t3ss-exos{background:#1565C0;color:white;border-color:#0D47A1}
.typing-badge.t3ss-both{background:#6A1B9A;color:white;border-color:#4A148C}
.typing-badge.t3ss-neither{background:#6b7280;color:white;border-color:#4B5563}
.vir-details{margin-top:8px}
.vir-details summary{cursor:pointer;color:#8B4513;font-weight:600;font-size:.85em;list-style:none}
.vir-details summary::-webkit-details-marker{display:none}
.vir-details summary::before{content:'▸ ';transition:transform .2s}
.vir-details[open] summary::before{content:'▾ '}
.vir-details .vir-gene-tags{margin-top:6px;display:flex;flex-wrap:wrap;gap:4px}
.accordion{margin:20px 0}
.accordion-item{background:#fdf6ec;border:1px solid #E8C39E;margin-bottom:10px;border-radius:8px;overflow:hidden}
.accordion-header{background:#f5e6d3;padding:12px 20px;cursor:pointer;font-weight:bold;color:#5C2E0B;display:flex;justify-content:space-between;align-items:center}
.accordion-header:hover{background:#E8C39E}
.accordion-content{padding:15px 20px;border-top:1px solid #E8C39E;background:white}
.copy-btn{background:#8B4513;color:white;border:none;padding:4px 12px;border-radius:20px;cursor:pointer;font-size:.8em;margin-left:10px}
.copy-btn:hover{background:#A0522D}
@media print{body *{visibility:hidden}.tab-content.active,.tab-content.active *{visibility:visible}.tab-content.active{position:absolute;left:0;top:0;width:100%;padding:20px;box-shadow:none;border-radius:0}.print-section-btn,.tab-navigation,.dashboard-grid,.search-box,.action-buttons,.grouping-controls{display:none !important}}

/* ========== Summary Visuals ========== */
.kpi-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:16px;margin:24px 0}
.kpi-card{color:white;padding:22px 18px;border-radius:14px;text-align:center;cursor:pointer;transition:all .3s ease;box-shadow:0 4px 14px rgba(0,0,0,.15);position:relative;overflow:hidden;border:none}
.kpi-card::after{content:"";position:absolute;top:-50%;right:-50%;width:100%;height:100%;background:radial-gradient(circle,rgba(255,255,255,.15),transparent 70%);pointer-events:none}
.kpi-card:hover{transform:translateY(-6px);box-shadow:0 12px 28px rgba(0,0,0,.22)}
.kpi-card .kpi-icon{font-size:1.6em;opacity:.95;margin-bottom:10px;display:block}
.kpi-card .kpi-value{font-size:2.4em;font-weight:800;line-height:1;margin-bottom:6px;font-variant-numeric:tabular-nums}
.kpi-card .kpi-label{font-size:.82em;font-weight:700;text-transform:uppercase;letter-spacing:.6px;opacity:.95}
.kpi-card .kpi-sub{font-size:.75em;opacity:.85;margin-top:6px;font-style:italic}
.visual-grid{display:grid;grid-template-columns:minmax(280px,340px) 1fr;gap:24px;margin:24px 0;align-items:start}
@media (max-width:900px){.visual-grid{grid-template-columns:1fr}}
.visual-card-big{background:white;border-radius:14px;padding:24px;border:1px solid #E8C39E;box-shadow:0 3px 12px rgba(139,69,19,.08)}
.visual-card-big h3{color:#5C2E0B;margin-bottom:16px;font-size:1.05em;display:flex;align-items:center;gap:10px;border-bottom:2px solid #E8C39E;padding-bottom:10px}
.donut-wrapper{display:flex;flex-direction:column;align-items:center;gap:18px}
.legend-list{display:flex;flex-direction:column;gap:8px;width:100%}
.legend-item{display:flex;align-items:center;gap:10px;font-size:.9em;color:#5C2E0B;padding:4px 0}
.legend-swatch{width:14px;height:14px;border-radius:4px;flex-shrink:0;box-shadow:0 1px 3px rgba(0,0,0,.15)}
.hbar-container{display:flex;flex-direction:column;gap:11px}
.hbar-row{display:grid;grid-template-columns:150px 1fr 100px;gap:12px;align-items:center;font-size:.88em}
@media (max-width:700px){.hbar-row{grid-template-columns:110px 1fr 90px}}
.hbar-label{font-weight:700;color:#5C2E0B;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;font-family:'Consolas','Monaco',monospace;font-size:.95em}
.hbar-track{background:#f5e6d3;border-radius:8px;height:22px;overflow:hidden;position:relative;box-shadow:inset 0 1px 3px rgba(139,69,19,.1)}
.hbar-fill{height:100%;border-radius:8px;transition:width 1s cubic-bezier(.4,0,.2,1)}
.hbar-value{text-align:right;font-weight:700;color:#5C2E0B;font-variant-numeric:tabular-nums;font-size:.92em}
.hbar-pct{font-weight:400;color:#8B4513;font-size:.85em;margin-left:2px}
.tour-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:14px;margin:20px 0}
.tour-card{display:flex;gap:14px;padding:16px;background:white;border-radius:10px;border-left:5px solid;box-shadow:0 2px 8px rgba(0,0,0,.06);cursor:pointer;transition:all .2s}
.tour-card:hover{transform:translateX(5px);box-shadow:0 6px 18px rgba(139,69,19,.15)}
.tour-icon{width:44px;height:44px;border-radius:10px;display:flex;align-items:center;justify-content:center;color:white;font-size:1.15em;flex-shrink:0;box-shadow:0 2px 6px rgba(0,0,0,.15)}
.tour-body{display:flex;flex-direction:column;gap:4px;min-width:0}
.tour-title{font-weight:700;color:#5C2E0B;font-size:.95em}
.tour-desc{font-size:.82em;color:#666;line-height:1.45}
.collapsible-section{background:#fdf6ec;border-radius:12px;overflow:hidden;margin:20px 0;border:1px solid #E8C39E;box-shadow:0 2px 8px rgba(139,69,19,.06)}
.collapsible-header{padding:16px 20px;cursor:pointer;font-weight:700;color:#5C2E0B;display:flex;justify-content:space-between;align-items:center;user-select:none;background:#fdf6ec;transition:background .2s}
.collapsible-header:hover{background:#f5e6d3}
.collapsible-toggle{transition:transform .3s ease}
.collapsible-section.open .collapsible-toggle{transform:rotate(180deg)}
.collapsible-content{padding:0 20px;background:white;max-height:0;overflow:hidden;transition:max-height .4s ease,padding .4s ease}
.collapsible-section.open .collapsible-content{max-height:3000px;padding:20px}
.summary-hero{background:linear-gradient(135deg,#2A1408 0%,#5C2E0B 50%,#8B4513 100%);color:white;padding:28px 32px;border-radius:16px;margin-bottom:24px;box-shadow:0 8px 24px rgba(0,0,0,.2)}
.summary-hero h2{font-size:1.5em;margin-bottom:8px;color:white}
.summary-hero p{opacity:.92;font-size:.95em;line-height:1.55}
.summary-hero .hero-badge{display:inline-block;background:rgba(255,255,255,.15);padding:4px 12px;border-radius:20px;font-size:.8em;font-weight:600;margin-right:8px;margin-top:12px}

/* Carbapenemase categories */
.carb-section{background:white;border-radius:12px;padding:20px;margin:20px 0;box-shadow:0 3px 12px rgba(0,0,0,.08);border-left:5px solid}
.carb-section.cat-acquired{border-left-color:#B71C1C}
.carb-section.cat-oxa-esbl{border-left-color:#E65100}
.carb-section.cat-intrinsic{border-left-color:#607D8B}
.carb-section.cat-other-esbl{border-left-color:#F57C00}
.carb-section.cat-colistin{border-left-color:#6A1B9A}
.carb-section.cat-methyl{border-left-color:#00838F}
.carb-section h3{margin-bottom:12px;font-size:1.05em;display:flex;align-items:center;gap:10px;padding-bottom:8px;border-bottom:2px solid #E8C39E}
.cat-acquired h3{color:#B71C1C}
.cat-oxa-esbl h3{color:#E65100}
.cat-intrinsic h3{color:#607D8B}
.cat-other-esbl h3{color:#F57C00}
.cat-colistin h3{color:#6A1B9A}
.cat-methyl h3{color:#00838F}

/* Compare tab */
.mode-toggle{display:flex;gap:8px;margin-bottom:18px}
.mode-btn{padding:10px 22px;border:2px solid #E8C39E;background:white;color:#5C2E0B;border-radius:8px;font-weight:600;cursor:pointer;font-size:.95em;display:flex;align-items:center;gap:8px;transition:all .2s}
.mode-btn:hover{border-color:#8B4513;color:#8B4513}
.mode-btn.active{background:#8B4513;border-color:#8B4513;color:white}
.compare-mode{display:none}
.compare-mode.active{display:block;animation:fadeIn .3s ease}
.compare-visuals-grid{display:grid;grid-template-columns:180px 1fr;gap:16px;margin:20px 0;align-items:center}
.visual-card{background:white;border-radius:12px;padding:20px;box-shadow:0 2px 8px rgba(0,0,0,.06);border:1px solid #E8C39E}
.gauge-card{display:flex;justify-content:center;align-items:center;padding:12px}
.similarity-gauge{width:160px;height:160px;display:block}
.metrics-card{display:flex;flex-direction:column;gap:18px}
.metric-block{display:flex;flex-direction:column;gap:6px}
.metric-label{font-size:.85em;font-weight:700;color:#5C2E0B;text-transform:uppercase;letter-spacing:.5px}
.metric-bar-track{background:#f5e6d3;height:22px;border-radius:11px;overflow:hidden;position:relative}
.metric-bar-fill{height:100%;border-radius:11px;transition:width .8s ease;min-width:4px}
.metric-value{font-size:.88em;color:#5C2E0B;font-weight:600}
.typing-pill{display:inline-block;padding:5px 12px;border-radius:12px;font-size:.85em;font-weight:600;border:1.5px solid;font-family:'Consolas',monospace}
.gene-proportion-bar{display:flex;height:12px;border-radius:6px;overflow:hidden;margin-bottom:14px;background:#f5e6d3}
.gene-prop-segment{height:100%;transition:width .6s ease;min-width:2px}
.gene-prop-segment.shared{background:#10b981}
.gene-prop-segment.only-a{background:#3b82f6}
.gene-prop-segment.only-b{background:#ec4899}
.cluster-card{background:white;border-radius:10px;padding:16px 20px;margin-bottom:12px;border-left:5px solid;box-shadow:0 2px 8px rgba(0,0,0,.06);transition:transform .15s}
.cluster-card:hover{transform:translateX(4px)}
.cluster-header{display:flex;align-items:center;gap:14px;margin-bottom:12px;flex-wrap:wrap}
.cluster-badge{color:white;padding:4px 14px;border-radius:14px;font-size:.85em;font-weight:700}
.cluster-stats{font-size:.9em;color:#5C2E0B;font-weight:600}
.cluster-samples{display:flex;flex-wrap:wrap;gap:6px}
.cluster-sample-pill{background:#fdf6ec;color:#5C2E0B;padding:4px 12px;border-radius:12px;font-size:.85em;font-weight:600;font-family:'Consolas',monospace;border:1px solid #E8C39E}
.cluster-controls{display:flex;gap:14px;align-items:flex-end;flex-wrap:wrap;margin-bottom:20px}
.cluster-controls label{display:flex;flex-direction:column;gap:4px}
.cluster-controls label span{font-size:.85em;font-weight:600;color:#5C2E0B}
.cluster-controls select{padding:10px 14px;border:2px solid #E8C39E;border-radius:8px;background:white;font-size:.95em;min-width:220px}
.heatmap-wrapper{background:white;border-radius:12px;padding:20px;margin-top:24px;box-shadow:0 2px 8px rgba(0,0,0,.06);overflow-x:auto;width:100%}
.heatmap-wrapper h3{margin-bottom:8px;color:#5C2E0B}
.heatmap-grid{display:grid;gap:3px;margin-top:16px;padding:4px;width:100%;min-width:max-content}
.heatmap-cell{width:100%;height:100%;min-width:34px;min-height:34px;display:flex;align-items:center;justify-content:center;font-size:.75em;font-weight:700;border-radius:4px;cursor:default}
.heatmap-cell:hover{filter:brightness(1.15);z-index:2;box-shadow:0 2px 8px rgba(0,0,0,.25)}
.heatmap-col-header{background:transparent;color:#5C2E0B;font-size:.72em;writing-mode:vertical-rl;text-orientation:mixed;height:150px;min-width:44px;font-family:'Consolas',monospace;display:flex;align-items:flex-end;justify-content:center;padding-bottom:6px;white-space:nowrap;overflow:visible}
.heatmap-row-header{background:transparent;color:#5C2E0B;font-size:.78em;display:flex;align-items:center;justify-content:flex-start;padding-left:6px;font-family:'Consolas',monospace;white-space:nowrap;min-width:200px}
.heatmap-corner{background:transparent}
.heatmap-cell[data-tip]{position:relative;cursor:help}
.heatmap-cell[data-tip]:hover{filter:brightness(1.18);z-index:60}
.heatmap-cell[data-tip]::after{content:attr(data-tip);position:absolute;bottom:calc(100% + 8px);left:50%;transform:translateX(-50%);background:#1e293b;color:#fff;padding:6px 12px;border-radius:6px;font-size:12px;font-weight:500;white-space:nowrap;pointer-events:none;opacity:0;transition:opacity .15s ease;z-index:300;font-family:'Segoe UI',Tahoma,sans-serif;letter-spacing:.2px;box-shadow:0 4px 14px rgba(0,0,0,.28);writing-mode:horizontal-tb;text-orientation:mixed;}
.heatmap-cell[data-tip]::before{content:'';position:absolute;bottom:100%;left:50%;transform:translateX(-50%);border:6px solid transparent;border-top-color:#1e293b;pointer-events:none;opacity:0;transition:opacity .15s ease;z-index:300;}
.heatmap-cell[data-tip]:hover::after,.heatmap-cell[data-tip]:hover::before{opacity:1}
.heatmap-col-header[data-tip]::after,.heatmap-row-header[data-tip]::after{bottom:auto;top:50%;left:calc(100% + 10px);transform:translateY(-50%);}
.heatmap-col-header[data-tip]::before,.heatmap-row-header[data-tip]::before{display:none;}
.compare-panel{background:#fdf6ec;border-radius:12px;padding:20px;margin-bottom:24px}
.compare-select-row{display:flex;gap:12px;align-items:flex-end;flex-wrap:wrap;margin-bottom:16px}
.compare-select-row label{display:flex;flex-direction:column;gap:4px}
.compare-select-row label span{font-size:.85em;font-weight:600;color:#5C2E0B}
.compare-select-row select{padding:10px 14px;border:2px solid #E8C39E;border-radius:8px;background:white;font-size:.95em;min-width:220px}
.compare-select-row select:focus{outline:none;border-color:#8B4513;box-shadow:0 0 0 3px rgba(139,69,19,.15)}
.compare-verdict-box{padding:16px 20px;border-radius:10px;margin-bottom:16px;display:flex;gap:12px;align-items:flex-start;font-size:.95em}
.compare-verdict-box.danger{background:#fee2e2;border-left:5px solid #C62828;color:#7f1d1d}
.compare-verdict-box.warning{background:#fff3cd;border-left:5px solid #ffc107;color:#78350f}
.compare-verdict-box.info{background:#eff6ff;border-left:5px solid #1565C0;color:#1e3a8a}
.compare-warning{padding:20px;background:#fee2e2;color:#7f1d1d;border-radius:8px;border-left:4px solid #C62828}
.compare-section{background:white;border-radius:12px;padding:20px;margin-bottom:20px;box-shadow:0 2px 8px rgba(0,0,0,.06)}
.compare-section h3{margin-bottom:14px;color:#5C2E0B;font-size:1.15em;padding-bottom:8px;border-bottom:2px solid #E8C39E}
.compare-table{display:flex;flex-direction:column;gap:2px}
.compare-row{display:grid;grid-template-columns:180px 1fr 1fr;gap:12px;padding:10px 14px;border-radius:6px;align-items:center;font-size:.92em}
.compare-row.header{background:#f5e6d3;font-weight:700;border-bottom:2px solid #E8C39E}
.compare-row.diff{background:#fee2e2;border-left:3px solid #C62828}
.compare-row.match{background:#f0fdf4;border-left:3px solid #16a34a}
.compare-row .row-label{font-weight:600;color:#5C2E0B;font-size:.88em;text-transform:uppercase;letter-spacing:.3px}
.compare-row .row-val{color:#2c2c2c;word-break:break-word}
.compare-gene-block{margin-bottom:24px;padding-bottom:20px;border-bottom:1px dashed #E8C39E}
.compare-gene-block:last-child{border-bottom:none;margin-bottom:0;padding-bottom:0}
.compare-gene-block h4{margin-bottom:12px;color:#5C2E0B;font-size:1.02em;display:flex;align-items:center;gap:10px;flex-wrap:wrap}
.gene-count{font-size:.82em;font-weight:400;color:#5C2E0B;background:#fdf6ec;padding:3px 10px;border-radius:12px}
.gene-columns{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:14px}
.gene-col{background:#fdf6ec;border-radius:8px;padding:12px;border:1px solid #E8C39E}
.gene-col-head{font-size:.85em;font-weight:700;color:#5C2E0B;margin-bottom:8px;text-transform:uppercase;letter-spacing:.3px}
.compare-gene{display:inline-block;font-size:.82em;padding:3px 10px;border-radius:12px;margin:2px;font-family:'Consolas',monospace}
.compare-gene.shared{background:#dcfce7;color:#14532d;border:1px solid #86efac}
.compare-gene.only-a{background:#dbeafe;color:#1e3a8a;border:1px solid #93c5fd}
.compare-gene.only-b{background:#fce7f3;color:#831843;border:1px solid #f9a8d4}
.compare-full-details{background:#fdf6ec;border:1px solid #E8C39E;border-radius:10px;padding:12px 16px;margin-bottom:10px}
.compare-full-details[open]{background:white;border-color:#8B4513}
.compare-full-details summary{cursor:pointer;font-weight:700;color:#8B4513;list-style:none;display:flex;justify-content:space-between;align-items:center;gap:10px;flex-wrap:wrap;padding:4px 0}
.compare-full-details summary::-webkit-details-marker{display:none}
.compare-full-details summary::before{content:'▸';display:inline-block;margin-right:8px;color:#8B4513}
.compare-full-details[open] summary::before{content:'▾'}
.details-cat-name{font-size:1em}
.details-counts{font-size:.85em;font-weight:500;color:#5C2E0B;background:#f5e6d3;padding:3px 10px;border-radius:12px;margin-left:auto}
.compare-full-profile{display:grid;grid-template-columns:1fr 1fr;gap:16px;margin-top:14px;padding-top:14px;border-top:1px dashed #E8C39E}
.profile-col{background:#fdf6ec;border-radius:10px;padding:14px;border:1px solid #E8C39E}
.profile-head{display:flex;justify-content:space-between;align-items:center;gap:8px;margin-bottom:10px;padding-bottom:8px;border-bottom:2px solid #E8C39E;flex-wrap:wrap}
.profile-sample-label{font-weight:700;color:#5C2E0B;font-size:.92em;word-break:break-word}
.profile-sample-count{font-size:.78em;font-weight:600;color:#8B4513;background:#dcfce7;padding:2px 10px;border-radius:12px;white-space:nowrap}
.empty-note{font-size:.85em;color:#94a3b8;font-style:italic;padding:6px 0;display:inline-block}
@media (max-width:900px){.compare-full-profile{grid-template-columns:1fr}.compare-row{grid-template-columns:1fr}.compare-visuals-grid{grid-template-columns:1fr}}
.carb-section.cat-unclassified{border-left-color:#F57C00}
.carb-section.cat-ambiguous{border-left-color:#9E9E9E}
.cat-unclassified h3{color:#F57C00}
.cat-ambiguous h3{color:#9E9E9E}
</style>"""

    def _get_js(self, typing_json: str) -> str:
        return f"""
<script>
var sampleTyping = {typing_json};
var originalGenomeLists = {{}};
function escapeHtml(str) {{
    return String(str).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}}
function switchTab(tabName){{
    document.querySelectorAll('.tab-content').forEach(t=>t.classList.remove('active'));
    document.querySelectorAll('.tab-button').forEach(b=>b.classList.remove('active'));
    var content = document.getElementById(tabName+'-tab');
    var button = document.querySelector('.tab-button.'+tabName);
    if(content) content.classList.add('active');
    if(button) button.classList.add('active');
    if(event && event.currentTarget) event.currentTarget.classList.add('active');
    window.location.hash = tabName;
}}
function searchTable(tableId, searchId){{
    var input = document.getElementById(searchId);
    if(!input) return;
    var filter = input.value.toUpperCase();
    var table = document.getElementById(tableId);
    if(!table) return;
    var rows = table.tBodies[0].rows;
    for(var i=0;i<rows.length;i++){{
        var cells = rows[i].getElementsByTagName('td');
        var found = false;
        for(var j=0;j<cells.length;j++){{
            if(cells[j] && (cells[j].textContent||cells[j].innerText).toUpperCase().indexOf(filter)>-1){{ found = true; break; }}
        }}
        rows[i].style.display = found ? '' : 'none';
    }}
}}
function highlightGenome(tableId, searchId){{
    var el = document.getElementById(searchId);
    if(!el) return;
    var filter = el.value.toUpperCase().trim();
    var table = document.getElementById(tableId);
    if(!table) return;
    table.querySelectorAll('.genome-tag').forEach(function(t){{
        t.classList.remove('highlight');
        if(filter && t.textContent.toUpperCase().indexOf(filter)>-1) t.classList.add('highlight');
    }});
}}
function getTypingValue(genome, groupBy){{
    var info = sampleTyping[genome];
    if(!info) return "Unknown";
    function pick(k){{ return info[k] || 'ND'; }}
    switch(groupBy){{
        case 'MLST': return pick('MLST');
        case 'Otype': return pick('Otype');
        case 'T3SS': return pick('T3SS');
        case 'MLST-Otype': return pick('MLST')+' - '+pick('Otype');
        case 'MLST-T3SS': return pick('MLST')+' - '+pick('T3SS');
        case 'Otype-T3SS': return pick('Otype')+' - '+pick('T3SS');
        case 'MLST-Otype-T3SS': return pick('MLST')+' - '+pick('Otype')+' - '+pick('T3SS');
    }}
    return "Unknown";
}}
function groupRowGenomes(row, groupBy, originalList){{
    var genomesCell = null;
    for(var i=0;i<row.cells.length;i++){{
        if(row.cells[i].querySelector('.genome-list')){{ genomesCell = row.cells[i]; break; }}
    }}
    if(!genomesCell) return;
    var genomes = originalList.slice();
    if(genomes.length===0){{ genomesCell.innerHTML='<div class="genome-list">None</div>'; return; }}
    var groups = {{}};
    genomes.forEach(function(g){{
        var key = getTypingValue(g, groupBy);
        if(!groups[key]) groups[key]=[];
        groups[key].push(g);
    }});
    var html = '<div class="genome-list">';
    for(var key in groups){{
        var tags = groups[key].map(g=>'<span class="genome-tag">'+g+'</span>').join('');
        html += '<div class="genome-group"><div class="genome-group-header">'+key+'</div><div class="genome-group-tags">'+tags+'</div></div>';
    }}
    html += '</div>';
    genomesCell.innerHTML = html;
}}
function groupGenomesByTyping(tableId, groupBy){{
    var table = document.getElementById(tableId);
    if(!table) return;
    var tbody = table.tBodies[0];
    if(!tbody) return;
    var rows = tbody.rows;
    for(var i=0;i<rows.length;i++){{
        var row = rows[i];
        var nameCell = row.cells[0];
        if(!nameCell) continue;
        var name = nameCell.textContent.trim().replace(/⚠️/g,'').trim();
        if(!originalGenomeLists[name]){{
            var genomesCell = null;
            for(var j=0;j<row.cells.length;j++){{
                if(row.cells[j].querySelector('.genome-list')){{ genomesCell = row.cells[j]; break; }}
            }}
            if(genomesCell){{
                var tags = genomesCell.querySelectorAll('.genome-tag');
                originalGenomeLists[name] = Array.from(tags).map(t=>t.textContent.trim());
            }} else {{ originalGenomeLists[name] = []; }}
        }}
    }}
    if(!groupBy){{ resetGenomeList(tableId); return; }}
    for(var i=0;i<rows.length;i++){{
        var row = rows[i];
        var nameCell = row.cells[0];
        if(!nameCell) continue;
        var name = nameCell.textContent.trim().replace(/⚠️/g,'').trim();
        groupRowGenomes(row, groupBy, originalGenomeLists[name]||[]);
    }}
}}
function resetGenomeList(tableId){{
    var table = document.getElementById(tableId);
    if(!table) return;
    var tbody = table.tBodies[0];
    if(!tbody) return;
    var rows = tbody.rows;
    for(var i=0;i<rows.length;i++){{
        var row = rows[i];
        var nameCell = row.cells[0];
        if(!nameCell) continue;
        var name = nameCell.textContent.trim().replace(/⚠️/g,'').trim();
        var original = originalGenomeLists[name] || [];
        for(var j=0;j<row.cells.length;j++){{
            if(row.cells[j].querySelector('.genome-list')){{
                var tags = original.map(g=>'<span class="genome-tag">'+g+'</span>').join('');
                row.cells[j].innerHTML = '<div class="genome-list">'+tags+'</div>';
                break;
            }}
        }}
    }}
    var sel = document.querySelector('#'+tableId)?.closest('.tab-content')?.querySelector('.group-select');
    if(sel) sel.value = '';
}}
function sortTable(tableId, colIndex, type){{
    var table = document.getElementById(tableId);
    if(!table) return;
    var tbody = table.tBodies[0];
    var rows = Array.from(tbody.rows);
    var asc = table.getAttribute('data-sort-dir') !== 'asc';
    rows.sort(function(a,b){{
        var av = a.cells[colIndex].innerText.trim();
        var bv = b.cells[colIndex].innerText.trim();
        if(type==='number'){{
            av = parseFloat(av.replace(/,/g,''))||0;
            bv = parseFloat(bv.replace(/,/g,''))||0;
            return asc ? av-bv : bv-av;
        }}
        return asc ? av.localeCompare(bv) : bv.localeCompare(av);
    }});
    tbody.append.apply(tbody, rows);
    table.setAttribute('data-sort-dir', asc ? 'asc' : 'desc');
}}
function printSection(id){{
    var content = document.getElementById(id);
    if(!content) return;
    var w = window.open('', '_blank');
    var style = document.querySelector('style');
    w.document.write('<html><head><title>Print</title>');
    if(style) w.document.write('<style>'+style.textContent+'</style>');
    w.document.write('</head><body>'+content.innerHTML+'</body></html>');
    w.document.close(); w.print();
}}
function exportTableToCSV(tableId, filename){{
    var table = document.getElementById(tableId);
    if(!table) return;
    var rows = table.querySelectorAll('tr');
    var csv = [];
    for(var i=0;i<rows.length;i++){{
        var row=[], cols=rows[i].querySelectorAll('td, th');
        for(var j=0;j<cols.length;j++){{
            row.push('"'+(cols[j].innerText||'').replace(/"/g,'""')+'"');
        }}
        csv.push(row.join(','));
    }}
    var blob = new Blob([csv.join('\\n')],{{type:'text/csv'}});
    var a = document.createElement('a');
    a.download = filename; a.href = URL.createObjectURL(blob);
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
}}

// ---------- Compare tool ----------
function setCompareMode(mode){{
    document.querySelectorAll('.mode-btn').forEach(b => {{ b.classList.toggle('active', b.dataset.mode === mode); }});
    document.querySelectorAll('.compare-mode').forEach(m => m.classList.remove('active'));
    const el = document.getElementById('mode-' + mode);
    if (el) el.classList.add('active');
    if (mode === 'cluster') runCluster();
}}
function computeSimilarity(A, B) {{
    var typingKeys = Object.keys(A.typing);
    var matches = 0;
    typingKeys.forEach(function(k) {{
        if ((A.typing[k] || '') === (B.typing[k] || '')) matches++;
    }});
    var typingPct = typingKeys.length ? (matches / typingKeys.length) * 100 : 0;
    var cats = ['AMR','Virulence','BACMET','Plasmids','Mutations'];
    var catScores = [];
    cats.forEach(function(cat) {{
        var sa = new Set(A.genes[cat] || []);
        var sb = new Set(B.genes[cat] || []);
        var union = new Set([...sa, ...sb]);
        if (union.size === 0) return;
        var shared = 0;
        union.forEach(function(g) {{ if (sa.has(g) && sb.has(g)) shared++; }});
        catScores.push((shared / union.size) * 100);
    }});
    var genePct = catScores.length ? catScores.reduce((a,b)=>a+b,0)/catScores.length : 0;
    var overall = (typingPct * 0.5) + (genePct * 0.5);
    return {{
        typing: Math.round(typingPct), genes: Math.round(genePct), overall: Math.round(overall),
        typingMatches: matches, typingTotal: typingKeys.length
    }};
}}
function similarityColor(pct) {{
    if (pct >= 95) return '#059669';
    if (pct >= 90) return '#10b981';
    if (pct >= 80) return '#84cc16';
    if (pct >= 70) return '#eab308';
    if (pct >= 50) return '#f59e0b';
    if (pct >= 30) return '#f97316';
    return '#C62828';
}}
function drawGauge(pct) {{
    var color = similarityColor(pct);
    var radius = 60;
    var circumference = 2 * Math.PI * radius;
    var offset = circumference * (1 - pct / 100);
    return '<svg viewBox="0 0 160 160" class="similarity-gauge">' +
        '<circle cx="80" cy="80" r="' + radius + '" fill="none" stroke="#f5e6d3" stroke-width="14"/>' +
        '<circle cx="80" cy="80" r="' + radius + '" fill="none" stroke="' + color + '" stroke-width="14" ' +
        'stroke-dasharray="' + circumference + '" stroke-dashoffset="' + offset + '" stroke-linecap="round" ' +
        'transform="rotate(-90 80 80)"/>' +
        '<text x="80" y="76" text-anchor="middle" font-size="28" font-weight="700" fill="' + color + '">' + pct + '%</text>' +
        '<text x="80" y="98" text-anchor="middle" font-size="11" fill="#5C2E0B">overall</text></svg>';
}}
function typingValuePill(va, vb) {{
    var match = (va === vb);
    var color = match ? '#10b981' : '#C62828';
    var bg = match ? '#dcfce7' : '#fee2e2';
    var icon = match ? '✓' : '✗';
    return '<span class="typing-pill" style="background:' + bg + ';color:' + color + ';border-color:' + color + ';">' +
        icon + ' ' + escapeHtml(va) + '</span>';
}}
function runCompare() {{
    var data = window.PSEUDOSCOPE_COMPARE_DATA || {{}};
    var a = document.getElementById('compare-a').value;
    var b = document.getElementById('compare-b').value;
    if (!a || !b) return;
    if (a === b) {{
        document.getElementById('compare-output').innerHTML = '<div class="compare-warning">Pick two different samples.</div>';
        document.getElementById('compare-visuals').innerHTML = '';
        document.getElementById('compare-verdict').style.display = 'none';
        return;
    }}
    if (!data[a] || !data[b]) {{
        document.getElementById('compare-output').innerHTML = '<div class="compare-warning">Sample data not found.</div>';
        return;
    }}
    var A = data[a], B = data[b];
    var sim = computeSimilarity(A, B);

    var visualsHtml = '<div class="compare-visuals-grid">' +
        '<div class="visual-card gauge-card">' + drawGauge(sim.overall) + '</div>' +
        '<div class="visual-card metrics-card">' +
            '<div class="metric-block"><div class="metric-label">Typing match</div>' +
            '<div class="metric-bar-track"><div class="metric-bar-fill" style="width:' + sim.typing + '%;background:' + similarityColor(sim.typing) + ';"></div></div>' +
            '<div class="metric-value">' + sim.typingMatches + '/' + sim.typingTotal + ' fields (' + sim.typing + '%)</div></div>' +
            '<div class="metric-block"><div class="metric-label">Gene content match</div>' +
            '<div class="metric-bar-track"><div class="metric-bar-fill" style="width:' + sim.genes + '%;background:' + similarityColor(sim.genes) + ';"></div></div>' +
            '<div class="metric-value">' + sim.genes + '% Jaccard similarity</div></div>' +
        '</div></div>';
    document.getElementById('compare-visuals').innerHTML = visualsHtml;

    var typingKeys = Object.keys(A.typing);
    var typingRows = '';
    typingKeys.forEach(function(k) {{
        var va = A.typing[k] || 'ND';
        var vb = B.typing[k] || 'ND';
        var diff = va !== vb;
        typingRows += '<div class="compare-row ' + (diff ? 'diff' : 'match') + '">' +
            '<span class="row-label">' + k + '</span>' +
            '<span class="row-val">' + typingValuePill(va, vb) + '</span>' +
            '<span class="row-val">' + typingValuePill(vb, va) + '</span></div>';
    }});

    var geneCategories = ['AMR','Virulence','BACMET','Plasmids','Mutations'];
    var geneHtml = '';
    geneCategories.forEach(function(cat) {{
        var ga = new Set(A.genes[cat] || []);
        var gb = new Set(B.genes[cat] || []);
        var union = new Set([...ga, ...gb]);
        if (union.size === 0) return;
        var shared = [], onlyA = [], onlyB = [];
        union.forEach(function(g) {{
            if (ga.has(g) && gb.has(g)) shared.push(g);
            else if (ga.has(g)) onlyA.push(g);
            else onlyB.push(g);
        }});
        var total = union.size;
        var sharedPct = (shared.length/total*100);
        var aPct = (onlyA.length/total*100);
        var bPct = (onlyB.length/total*100);
        geneHtml += '<div class="compare-gene-block"><h4>' + cat + ' <span class="gene-count">' +
            shared.length + ' shared · ' + onlyA.length + ' only A · ' + onlyB.length + ' only B (' +
            sharedPct.toFixed(0) + '% match)</span></h4>' +
            '<div class="gene-proportion-bar">' +
                '<div class="gene-prop-segment shared" style="width:' + sharedPct + '%;" title="Shared"></div>' +
                '<div class="gene-prop-segment only-a" style="width:' + aPct + '%;" title="Only A"></div>' +
                '<div class="gene-prop-segment only-b" style="width:' + bPct + '%;" title="Only B"></div>' +
            '</div><div class="gene-columns">' +
            '<div class="gene-col"><div class="gene-col-head">Shared (' + shared.length + ')</div>' +
            (shared.map(g=>'<span class="compare-gene shared">' + escapeHtml(g) + '</span>').join('') || '<em>None</em>') + '</div>' +
            '<div class="gene-col"><div class="gene-col-head">Only in A (' + onlyA.length + ')</div>' +
            (onlyA.map(g=>'<span class="compare-gene only-a">' + escapeHtml(g) + '</span>').join('') || '<em>None</em>') + '</div>' +
            '<div class="gene-col"><div class="gene-col-head">Only in B (' + onlyB.length + ')</div>' +
            (onlyB.map(g=>'<span class="compare-gene only-b">' + escapeHtml(g) + '</span>').join('') || '<em>None</em>') + '</div>' +
            '</div></div>';
    }});

    var fullProfileHtml = '<div class="compare-section"><h3>📋 Full Gene Profiles</h3>' +
        '<p style="font-size:0.9em;color:#5C2E0B;margin-bottom:12px;">Click any category to expand both samples\\' complete gene lists side by side.</p>';
    geneCategories.forEach(function(cat) {{
        var ga = A.genes[cat] || [];
        var gb = B.genes[cat] || [];
        if (ga.length === 0 && gb.length === 0) return;
        var tagsA = ga.length ? ga.map(g=>'<span class="genome-tag">' + escapeHtml(g) + '</span>').join('') : '<em class="empty-note">None detected</em>';
        var tagsB = gb.length ? gb.map(g=>'<span class="genome-tag">' + escapeHtml(g) + '</span>').join('') : '<em class="empty-note">None detected</em>';
        fullProfileHtml += '<details class="compare-full-details"><summary>' +
            '<span class="details-cat-name">' + cat + '</span>' +
            '<span class="details-counts">Sample A: ' + ga.length + ' · Sample B: ' + gb.length + '</span></summary>' +
            '<div class="compare-full-profile">' +
                '<div class="profile-col"><div class="profile-head">' +
                    '<span class="profile-sample-label">' + escapeHtml(a) + '</span>' +
                    '<span class="profile-sample-count">' + ga.length + ' gene(s)</span>' +
                '</div><div class="genome-list">' + tagsA + '</div></div>' +
                '<div class="profile-col"><div class="profile-head">' +
                    '<span class="profile-sample-label">' + escapeHtml(b) + '</span>' +
                    '<span class="profile-sample-count">' + gb.length + ' gene(s)</span>' +
                '</div><div class="genome-list">' + tagsB + '</div></div>' +
            '</div></details>';
    }});
    fullProfileHtml += '</div>';

    var verdict = document.getElementById('compare-verdict');
    if (sim.overall >= 95) {{
        verdict.innerHTML = '<div class="compare-verdict-box danger"><i class="fas fa-exclamation-triangle"></i> ' +
            '<strong>Near-identical profile (' + sim.overall + '%)</strong> — likely the same strain. ' +
            'Investigate as possible <em>transmission pair</em> or duplicate sample.</div>';
    }} else if (sim.overall >= 85) {{
        verdict.innerHTML = '<div class="compare-verdict-box warning"><i class="fas fa-info-circle"></i> ' +
            '<strong>High similarity (' + sim.overall + '%)</strong> — very likely the same clonal lineage with small marker differences.</div>';
    }} else if (sim.overall >= 60) {{
        verdict.innerHTML = '<div class="compare-verdict-box warning"><i class="fas fa-info-circle"></i> ' +
            '<strong>Moderate similarity (' + sim.overall + '%)</strong> — related lineage, but not a recent transmission pair.</div>';
    }} else {{
        verdict.innerHTML = '<div class="compare-verdict-box info"><i class="fas fa-info-circle"></i> ' +
            '<strong>Low similarity (' + sim.overall + '%)</strong> — distinct lineages.</div>';
    }}
    verdict.style.display = 'block';

    var tableHtml = '<div class="compare-section"><h3>🧬 Typing Profile</h3>' +
        '<div class="compare-table"><div class="compare-row header">' +
        '<span class="row-label">Field</span>' +
        '<span class="row-val"><strong>' + escapeHtml(a) + '</strong></span>' +
        '<span class="row-val"><strong>' + escapeHtml(b) + '</strong></span>' +
        '</div>' + typingRows + '</div></div>';

    document.getElementById('compare-output').innerHTML =
        tableHtml +
        '<div class="compare-section"><h3>💊 Gene Content</h3>' +
        (geneHtml || '<p><em>No gene data available for these samples.</em></p>') +
        '</div>' + fullProfileHtml;
}}
function swapCompare() {{
    var a = document.getElementById('compare-a');
    var b = document.getElementById('compare-b');
    var tmp = a.value; a.value = b.value; b.value = tmp;
    runCompare();
}}
function resetCompare() {{
    var sel = document.getElementById('compare-a');
    if (!sel) return;
    var opts = sel.options;
    if (opts.length >= 2) {{
        document.getElementById('compare-a').value = opts[0].value;
        document.getElementById('compare-b').value = opts[1].value;
    }}
    document.getElementById('compare-output').innerHTML = '';
    document.getElementById('compare-visuals').innerHTML = '';
    document.getElementById('compare-verdict').style.display = 'none';
}}
function getPrecomputedSim(i, j) {{
    if (i === j) return 100;
    if (i > j) {{ var t = i; i = j; j = t; }}
    var sim = window.PSEUDOSCOPE_SIM;
    if (!sim || !sim.flat) {{
        var data = window.PSEUDOSCOPE_COMPARE_DATA || {{}};
        var sorted = Object.keys(data).sort();
        if (!sorted[i] || !sorted[j]) return 0;
        return computeSimilarity(data[sorted[i]], data[sorted[j]]).overall;
    }}
    var n = sim.n;
    return sim.flat[i * n - i * (i + 1) / 2 + (j - i - 1)];
}}

function renderHeatmapDOM(samples, n) {{
    var rowLabelWidth = 200;
    var cellSize = 34;
    var truncate = function(s, max) {{ return s.length > max ? s.substring(0, max - 1) + '…' : s; }};
    var html = '<div class="heatmap-wrapper"><h3>🔷 Pairwise Similarity Matrix</h3>' +
        '<p style="font-size:0.9em;color:#5C2E0B;margin-bottom:12px;">Hover any cell for details.</p>' +
        '<div class="heatmap-grid" style="grid-template-columns: ' +
        rowLabelWidth + 'px repeat(' + n + ', ' + cellSize + 'px);">';
    html += '<div class="heatmap-cell heatmap-corner"></div>';
    samples.forEach(function(s) {{
        html += '<div class="heatmap-cell heatmap-col-header" data-tip="' + escapeHtml(s) + '">' +
                escapeHtml(s) + '</div>';
    }});
    samples.forEach(function(s1, i) {{
        html += '<div class="heatmap-cell heatmap-row-header" data-tip="' + escapeHtml(s1) + '" ' +
                'style="justify-content:flex-start;padding-left:8px;width:' + rowLabelWidth + 'px;">' +
                escapeHtml(truncate(s1, 26)) + '</div>';
        samples.forEach(function(s2, j) {{
            var pct = (i === j) ? 100 : getPrecomputedSim(i, j);
            var color = (i === j) ? '#5C2E0B' : similarityColor(pct);
            var text = (i === j) ? '—' : Math.round(pct);
            var tip = (i === j) ? escapeHtml(s1) + ' (self)' :
                    escapeHtml(s1) + ' vs ' + escapeHtml(s2) + ' → ' + Math.round(pct) + '%';
            html += '<div class="heatmap-cell" data-tip="' + tip + '" style="background:' + color +
                    ';color:white;">' + text + '</div>';
        }});
    }});
    html += '</div></div>';
    document.getElementById('cluster-heatmap').innerHTML = html;
}}

function renderHeatmapCanvas(samples, n) {{
    var cellSize = Math.max(3, Math.min(12, Math.floor(700 / n)));
    var rowLabelW = 160;
    var headerH = 60;
    var canvasW = rowLabelW + n * cellSize + 20;
    var canvasH = headerH + n * cellSize + 20;
    var canvas = document.createElement('canvas');
    canvas.width = canvasW;
    canvas.height = canvasH;
    canvas.style.display = 'block';
    canvas.style.maxWidth = '100%';
    canvas.style.cursor = 'crosshair';
    var ctx = canvas.getContext('2d');
    ctx.font = Math.max(8, Math.min(11, cellSize - 1)) + 'px sans-serif';
    ctx.textBaseline = 'middle';
    for (var i = 0; i < n; i++) {{
        for (var j = 0; j < n; j++) {{
            var pct = (i === j) ? 100 : getPrecomputedSim(i, j);
            ctx.fillStyle = (i === j) ? '#5C2E0B' : similarityColor(pct);
            ctx.fillRect(rowLabelW + j * cellSize, headerH + i * cellSize, cellSize, cellSize);
        }}
    }}
    ctx.fillStyle = '#5C2E0B';
    ctx.textAlign = 'right';
    var rowStep = Math.max(1, Math.floor(n / 40));
    for (var i2 = 0; i2 < n; i2 += rowStep) {{
        var lb = samples[i2].length > 20 ? samples[i2].substring(0, 19) + '…' : samples[i2];
        ctx.fillText(lb, rowLabelW - 5, headerH + i2 * cellSize + cellSize / 2);
    }}
    ctx.save();
    ctx.translate(0, headerH - 5);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'left';
    var colStep = Math.max(1, Math.floor(n / 40));
    for (var j2 = 0; j2 < n; j2 += colStep) {{
        var lb2 = samples[j2].length > 15 ? samples[j2].substring(0, 14) + '…' : samples[j2];
        ctx.fillText(lb2, 5, rowLabelW + j2 * cellSize + cellSize / 2);
    }}
    ctx.restore();
    canvas.addEventListener('click', function(e) {{
        var rect = canvas.getBoundingClientRect();
        var x = (e.clientX - rect.left) * (canvas.width / rect.width) - rowLabelW;
        var y = (e.clientY - rect.top) * (canvas.height / rect.height) - headerH;
        var jj = Math.floor(x / cellSize);
        var ii = Math.floor(y / cellSize);
        if (ii >= 0 && ii < n && jj >= 0 && jj < n && ii !== jj) {{
            document.getElementById('compare-a').value = samples[ii];
            document.getElementById('compare-b').value = samples[jj];
            setCompareMode('pair');
            runCompare();
        }}
    }});
    var wrapper = document.getElementById('cluster-heatmap');
    wrapper.innerHTML =
        '<div class="heatmap-wrapper"><h3>🔷 Pairwise Similarity Matrix</h3>' +
        '<p style="font-size:0.9em;color:#5C2E0B;margin-bottom:12px;">' +
        '🎨 Canvas-rendered for ' + n + ' samples. ' +
        '<strong>Click any cell to open pair-compare.</strong></p>' +
        '<div style="overflow-x:auto;"></div></div>';
    wrapper.querySelector('div[style*="overflow"]').appendChild(canvas);
}}

function runCluster() {{
    var data = window.PSEUDOSCOPE_COMPARE_DATA || {{}};
    var samples = Object.keys(data).sort();
    if (samples.length < 2) {{
        document.getElementById('cluster-heatmap').innerHTML = '<div class="compare-warning">Need at least two samples.</div>';
        return;
    }}
    var thresholdEl = document.getElementById('cluster-threshold');
    var threshold = thresholdEl ? parseInt(thresholdEl.value) : 90;
    var matrix = {{}};
    samples.forEach(function(s) {{ matrix[s] = {{}}; }});
    samples.forEach(function(s1, i) {{
        samples.forEach(function(s2, j) {{
            if (i === j) {{ matrix[s1][s2] = 100; }}
            else {{ matrix[s1][s2] = getPrecomputedSim(i, j); }}
        }});
    }});
    var parent = {{}};
    samples.forEach(function(s) {{ parent[s] = s; }});
    function find(x) {{ while (parent[x] !== x) {{ parent[x] = parent[parent[x]]; x = parent[x]; }} return x; }}
    function union(a, b) {{ var ra = find(a), rb = find(b); if (ra !== rb) parent[ra] = rb; }}
    for (var i = 0; i < samples.length; i++) {{
        for (var j = i + 1; j < samples.length; j++) {{
            if (matrix[samples[i]][samples[j]] >= threshold) union(samples[i], samples[j]);
        }}
    }}
    var groups = {{}};
    samples.forEach(function(s) {{ var r = find(s); if (!groups[r]) groups[r] = []; groups[r].push(s); }});
    var clusters = Object.values(groups).filter(function(g) {{ return g.length >= 2; }});
    clusters.sort(function(a, b) {{ return b.length - a.length; }});

    var clusterHtml = '';
    if (clusters.length === 0) {{
        clusterHtml = '<div class="alert-box alert-info" style="border-left-color:#1565C0;">' +
            'No clusters detected at ≥' + threshold + '% similarity. Every sample is sufficiently distinct.</div>';
    }} else {{
        clusters.forEach(function(cluster, idx) {{
            var sum = 0, count = 0;
            for (var i = 0; i < cluster.length; i++) {{
                for (var j = i + 1; j < cluster.length; j++) {{ sum += matrix[cluster[i]][cluster[j]]; count++; }}
            }}
            var mean = count ? (sum / count).toFixed(1) : '100';
            var color = similarityColor(parseFloat(mean));
            clusterHtml += '<div class="cluster-card" style="border-left-color:' + color + ';">' +
                '<div class="cluster-header">' +
                    '<span class="cluster-badge" style="background:' + color + ';">Cluster ' + (idx + 1) + '</span>' +
                    '<span class="cluster-stats">' + cluster.length + ' samples · mean ' + mean + '% similar</span>' +
                '</div><div class="cluster-samples">' +
                    cluster.map(s => '<span class="cluster-sample-pill">' + escapeHtml(s) + '</span>').join('') +
                '</div></div>';
        }});
    }}

    document.getElementById('cluster-clusters').innerHTML = clusterHtml;

    if (samples.length <= 60) {{
        renderHeatmapDOM(samples, samples.length);
    }} else {{
        renderHeatmapCanvas(samples, samples.length);
    }}

    var summary = document.getElementById('cluster-summary');
    if (summary) summary.textContent = clusters.length + ' cluster(s) detected at ≥' + threshold + '%';
}}

document.addEventListener('DOMContentLoaded', function(){{
    var hash = window.location.hash.substring(1);
    var target = hash ? document.querySelector('.tab-button.'+hash) : document.querySelector('.tab-button');
    if(target) target.click();
    document.querySelectorAll('.data-table').forEach(function(table){{
        var headers = table.querySelectorAll('th');
        headers.forEach(function(h, idx){{
            var type = h.getAttribute('data-sort') || 'string';
            h.style.cursor = 'pointer';
            h.addEventListener('click', function(){{ sortTable(table.id, idx, type); }});
            var icon = document.createElement('span');
            icon.className='sort-icon'; icon.innerHTML='⇅';
            h.appendChild(icon);
        }});
    }});
    document.querySelectorAll('.accordion-header').forEach(function(header){{
        header.addEventListener('click', function(){{
            var content = this.nextElementSibling;
            content.style.display = content.style.display === 'block' ? 'none' : 'block';
        }});
    }});
    document.querySelectorAll('.collapsible-header').forEach(function(h){{
        h.addEventListener('click', function(){{
            h.parentElement.classList.toggle('open');
        }});
    }});
    document.querySelectorAll('.copy-btn').forEach(function(b){{
        b.addEventListener('click', function(){{
            var c = this.getAttribute('data-citation') || '';
            navigator.clipboard.writeText(c).then(()=>{{
                var t = this.innerHTML;
                this.innerHTML='✓ Copied!';
                setTimeout(()=>{{this.innerHTML=t;}},2000);
            }});
        }});
    }});
}});
</script>"""

    def _precompute_similarity(self, compare_data: Dict[str, Dict]) -> Dict[str, Any]:
        """Precompute pairwise similarity for the compare/cluster tab.

        Returns a flat upper-triangle array of length n*(n-1)/2 indexed by
        i*n - i*(i+1)/2 + (j-i-1), so the JS heatmap can read without
        recomputing similarity in the browser. Scales to thousands of samples.
        """
        samples = sorted(compare_data.keys())
        n = len(samples)
        if n < 2:
            return {'samples': samples, 'flat': [], 'n': n}
        cats = ('AMR', 'Virulence', 'BACMET', 'Plasmids', 'Mutations')
        gene_sets, typing_data = {}, {}
        for s in samples:
            d = compare_data[s]
            gene_sets[s] = {cat: set(d['genes'].get(cat, [])) for cat in cats}
            typing_data[s] = d['typing']
        typing_keys = list(typing_data[samples[0]].keys()) if samples else []
        n_keys = len(typing_keys)
        flat = []
        for i in range(n):
            s1 = samples[i]
            for j in range(i + 1, n):
                s2 = samples[j]
                matches = sum(
                    1 for k in typing_keys
                    if typing_data[s1].get(k) == typing_data[s2].get(k)
                )
                typing_pct = (matches / n_keys * 100) if n_keys else 0
                cat_scores = []
                for cat in cats:
                    sa, sb = gene_sets[s1][cat], gene_sets[s2][cat]
                    union = sa | sb
                    if not union:
                        continue
                    cat_scores.append(len(sa & sb) / len(union) * 100)
                gene_pct = (sum(cat_scores) / len(cat_scores)) if cat_scores else 0
                flat.append(round(typing_pct * 0.5 + gene_pct * 0.5, 1))
        return {'samples': samples, 'flat': flat, 'n': n}

    # ------------------------------------------------------------------
    # Main report assembly
    # ------------------------------------------------------------------
    def generate_main_report(self, integrated_data: Dict[str, Any], output_dir: Path) -> str:
        print("\n🎨 Generating PSEUDOSCOPE ULTIMATE HTML report...")
        samples_data = integrated_data.get('samples', {})
        patterns = integrated_data.get('patterns', {})
        gene_centric = integrated_data.get('gene_centric', {})
        metadata = integrated_data.get('metadata', {})
        self._current_samples_data = samples_data

        sample_typing_js = {}
        for sample, data in samples_data.items():
            t = data.get('typing', {})
            sample_typing_js[sample] = {
                "MLST":  t.get('MLST', 'ND'),
                "Otype": t.get('O_Type', 'ND'),
                "T3SS":  t.get('T3SS', 'ND'),
            }

        css = self._get_css()
        js = self._get_js(json.dumps(sample_typing_js))

        TABS = [
            ('summary', 'Summary', 'chart-pie'),
            ('samples', 'Sample Overview', 'list-alt'),
            ('qc', 'FASTA QC', 'chart-line'),
            ('mlst', 'MLST', 'code-branch'),
            ('serotype', 'O-Serotype', 'tag'),
            ('t3ss', 'ExoU/ExoS Pathotype', 'syringe'),
            ('amr', 'AMR Genes', 'biohazard'),
            ('virulence', 'Virulence', 'virus'),
            ('plasmid', 'Plasmid Replicons', 'dna'),
            ('bacmet', 'Bacmet2', 'flask'),
            ('mutation', 'Mutations', 'dna'),
            ('mge', 'Mobile Genetic Elements', 'mobile-alt'),
            ('patterns', 'Patterns', 'project-diagram'),
            ('compare', 'Compare Isolates', 'balance-scale'),
            ('citation', 'Citation', 'book'),
            ('guide', 'Guide', 'question-circle'),
            ('calltoaction', 'Call to Action', 'globe'),
            ('export', 'Export', 'download'),
        ]

        nav_html = ''
        for i, (tid, title, icon) in enumerate(TABS):
            active = ' active' if i == 0 else ''
            nav_html += (f'<button class="tab-button {tid}{active}" onclick="switchTab(\'{tid}\')">'
                         f'<i class="fas fa-{icon}"></i> {title}</button>')

        section_methods = {
            'summary': self._sec_summary,
            'samples': self._sec_sample_overview,
            'qc': self._sec_qc,
            'mlst': self._sec_mlst,
            'serotype': self._sec_serotype,
            't3ss': self._sec_t3ss,
            'amr': self._sec_amr,
            'virulence': self._sec_virulence,
            'plasmid': self._sec_plasmid,
            'bacmet': self._sec_bacmet,
            'mutation': self._sec_mutation,
            'mge': self._sec_mge,
            'patterns': self._sec_patterns,
            'compare': self._sec_compare,
            'citation': self._sec_citation,
            'guide': self._sec_guide,
            'calltoaction': self._sec_calltoaction,
            'export': self._sec_export,
        }

        tabs_html = ''
        for i, (tid, title, icon) in enumerate(TABS):
            active = ' active' if i == 0 else ''
            color = self.tab_colors.get(tid, '#8B4513')
            content = section_methods[tid](integrated_data)
            tabs_html += f'''
            <div id="{tid}-tab" class="tab-content{active}">
                <h2 class="section-header" style="border-color:{color};">
                    <span><i class="fas fa-{icon}" style="color:{color};"></i> {title}</span>
                    <button class="print-section-btn" onclick="printSection('{tid}-tab')">
                        <i class="fas fa-print"></i> Print</button>
                </h2>
                {content}
            </div>'''

        total_samples = len(samples_data)
        n_st = len(patterns.get('mlst_distribution', {}))
        n_o = len(patterns.get('o_type_distribution', {}))
        n_exou = len(patterns.get('t3ss_to_samples', {}).get('ExoU', []))
        n_amr = sum(len(g) for g in gene_centric.get('amr_databases', {}).values())
        n_vir = sum(len(g) for g in gene_centric.get('virulence_databases', {}).values())
        n_mge = len(integrated_data.get('mge_data', {}).get('per_sample', {}))
        n_mut = len(integrated_data.get('mutation_data', {}).get('mutations', []))

        dash_html = f'''
        <div class="dashboard-grid">
            <div class="dashboard-card" onclick="switchTab('summary')" style="border-left-color:#8B4513;">
                <i class="fas fa-vial fa-2x" style="color:#8B4513;"></i>
                <div class="card-number">{total_samples}</div>
                <div class="card-label">Total Samples</div></div>
            <div class="dashboard-card" onclick="switchTab('mlst')" style="border-left-color:#1565C0;">
                <i class="fas fa-code-branch fa-2x" style="color:#1565C0;"></i>
                <div class="card-number">{n_st}</div>
                <div class="card-label">Unique STs</div></div>
            <div class="dashboard-card" onclick="switchTab('serotype')" style="border-left-color:#6A1B9A;">
                <i class="fas fa-tag fa-2x" style="color:#6A1B9A;"></i>
                <div class="card-number">{n_o}</div>
                <div class="card-label">O-Types</div></div>
            <div class="dashboard-card" onclick="switchTab('t3ss')" style="border-left-color:#D32F2F;">
                <i class="fas fa-syringe fa-2x" style="color:#D32F2F;"></i>
                <div class="card-number">{n_exou}</div>
                <div class="card-label">ExoU Isolates</div></div>
            <div class="dashboard-card" onclick="switchTab('amr')" style="border-left-color:#C62828;">
                <i class="fas fa-biohazard fa-2x" style="color:#C62828;"></i>
                <div class="card-number">{n_amr}</div>
                <div class="card-label">AMR Genes</div></div>
            <div class="dashboard-card" onclick="switchTab('virulence')" style="border-left-color:#E65100;">
                <i class="fas fa-virus fa-2x" style="color:#E65100;"></i>
                <div class="card-number">{n_vir}</div>
                <div class="card-label">Virulence Genes</div></div>
            <div class="dashboard-card" onclick="switchTab('mge')" style="border-left-color:#16A085;">
                <i class="fas fa-mobile-alt fa-2x" style="color:#16A085;"></i>
                <div class="card-number">{n_mge}</div>
                <div class="card-label">MGE-Profiled</div></div>
            <div class="dashboard-card" onclick="switchTab('mutation')" style="border-left-color:#00BCD4;">
                <i class="fas fa-dna fa-2x" style="color:#00BCD4;"></i>
                <div class="card-number">{n_mut}</div>
                <div class="card-label">Mutations</div></div>
        </div>'''

        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PseudoScope Ultimate P. aeruginosa Report v3.0.0</title>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    {css}
    {js}
</head>
<body>
<div class="container">
    <div class="main-header">
        <h1><i class="fas fa-dna"></i> PseudoScope Ultimate P. aeruginosa Analysis Report</h1>
        <p>Gene-Centric Cross-Genome Analysis with Full Typing, T3SS Pathotyping, MGE, and Interaction</p>
        <div class="metadata-bar">
            <div class="metadata-item"><i class="fas fa-calendar"></i><span>Generated: {metadata.get('analysis_date', 'Unknown')}</span></div>
            <div class="metadata-item"><i class="fas fa-database"></i><span>Samples: {total_samples}</span></div>
            <div class="metadata-item"><i class="fas fa-code-branch"></i><span>Tool: PseudoScope Ultimate v3.0.0</span></div>
            <div class="metadata-item"><i class="fas fa-university"></i><span>University of Ghana Medical School</span></div>
        </div>
    </div>
    {dash_html}
    <div class="tab-navigation">{nav_html}</div>
    {tabs_html}
    <div class="footer">
<p><strong>PSEUDOSCOPE</strong> - Gene-Centric Cross-Genome Integration</p>
<p class="timestamp">Generated on {metadata.get('analysis_date', 'Unknown')}</p>
<div class="authorship">
<p><strong>Technical Support &amp; Inquiries:</strong></p>
<p>Maintainer: Brown Beckley</p>
<p>GitHub: <a href="https://github.com/bbeckley-hub" style="color:#E8C39E;">bbeckley-hub</a></p>
<p>Email: <a href="mailto:brownbeckley94@gmail.com" style="color:#E8C39E;">brownbeckley94@gmail.com</a></p>
<p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
<p style="margin-top:10px;font-size:11px;opacity:0.8;">
Powered by the PseudoScope gene-centric integration engine
</p>
</div>
</div>
</div>
</body>
</html>'''
        output_file = output_dir / "genius_pseudomonas_ultimate_gene_centric_report.html"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"    ✅ HTML report saved: {output_file}")
        return str(output_file)

    # =========================================================================
    # SECTIONS
    # =========================================================================

    def _sec_summary(self, d):
        samples_data = d['samples']
        patterns = d['patterns']
        gene_centric = d['gene_centric']

        total = len(samples_data)

        def clean_dist(c):
            return Counter({k: v for k, v in c.items()
                            if k not in ('ND', 'Not Assigned', '', 'nan', 'None')})

        st_dist = clean_dist(patterns.get('mlst_distribution', Counter()))
        o_dist = clean_dist(patterns.get('o_type_distribution', Counter()))
        t3ss_dist = patterns.get('t3ss_distribution', Counter())

        n_st = len(st_dist)
        n_o = len(o_dist)
        n_exou = t3ss_dist.get('ExoU', 0)
        n_exos = t3ss_dist.get('ExoS', 0)
        n_both = t3ss_dist.get('ExoU+ExoS', 0)
        n_neither = t3ss_dist.get('Neither', 0)

        top_st = st_dist.most_common(1)
        top_st = top_st[0] if top_st else ('—', 0)
        top_o = o_dist.most_common(1)
        top_o = top_o[0] if top_o else ('—', 0)

        all_amr, all_vir = [], []
        for genes in gene_centric.get('amr_databases', {}).values():
            all_amr.extend(genes)
        for genes in gene_centric.get('virulence_databases', {}).values():
            all_vir.extend(genes)

        n_amr_unique = len({g['gene'] for g in all_amr})
        n_vir_unique = len({g['gene'] for g in all_vir})
        n_mge = len(d.get('mge_data', {}).get('per_sample', {}))
        n_mut = len(d.get('mutation_data', {}).get('mutations', []))
        n_high_risk = len(patterns.get('high_risk_combinations', []))

        exou_pct = (n_exou / total * 100) if total else 0

        hero = f'''
        <div class="summary-hero">
            <h2><i class="fas fa-dna"></i> P. aeruginosa Cohort at a Glance</h2>
            <p>A gene-centric cross-genome analysis of <strong>{total}</strong> isolate{'s' if total != 1 else ''},
            profiling population structure, T3SS pathotype, resistome, mobilome, and co-occurrence patterns.
            Every gene-centric table below groups hits by typing context, so you can trace which clones
            carry which resistance and virulence determinants.</p>
            <div>
                <span class="hero-badge"><i class="fas fa-syringe"></i> {exou_pct:.0f}% ExoU+</span>
                <span class="hero-badge"><i class="fas fa-project-diagram"></i> {n_high_risk} high-risk combos</span>
                <span class="hero-badge"><i class="fas fa-mobile-alt"></i> {n_mge} MGE-profiled</span>
            </div>
        </div>'''

        def kpi(value, label, sub, color, icon, tab):
            return (f'<div class="kpi-card" '
                    f'style="background:linear-gradient(135deg,{color} 0%,{color}cc 100%);" '
                    f'onclick="switchTab(\'{tab}\')" title="Click to open">'
                    f'<i class="fas fa-{icon} kpi-icon"></i>'
                    f'<div class="kpi-value">{value}</div>'
                    f'<div class="kpi-label">{label}</div>'
                    f'<div class="kpi-sub">{sub}</div></div>')

        kpis = ''.join([
            kpi(total, 'Isolates', 'in the cohort', '#8B4513', 'vial', 'samples'),
            kpi(n_st, 'Sequence Types', f'dominant: {top_st[0]}', '#1565C0', 'code-branch', 'mlst'),
            kpi(n_o, 'O-Types', f'dominant: {top_o[0]}', '#6A1B9A', 'tag', 'serotype'),
            kpi(n_exou, 'ExoU+', f'cytotoxic pathotype ({exou_pct:.0f}%)', '#C62828', 'syringe', 't3ss'),
            kpi(n_amr_unique, 'AMR Genes', 'unique across DBs', '#E65100', 'biohazard', 'amr'),
            kpi(n_vir_unique, 'Virulence Genes', 'from VFDB', '#6A1B9A', 'virus', 'virulence'),
            kpi(n_mut, 'Mutations', 'point mutations', '#00BCD4', 'dna', 'mutation'),
        ])

        t3ss_segments = [
            ('ExoU (cytotoxic)',  n_exou,    '#C62828'),
            ('ExoS (invasive)',   n_exos,    '#1565C0'),
            ('ExoU+ExoS',         n_both,    '#6A1B9A'),
            ('Neither',           n_neither, '#9E9E9E'),
        ]
        donut = self._donut_chart(t3ss_segments, size=210, thickness=34,
                                  center_text=str(total), center_sub='isolates')
        legend = ''.join(
            f'<div class="legend-item"><span class="legend-swatch" style="background:{c};"></span>'
            f'<strong>{l}</strong> — {v} ({(v/total*100) if total else 0:.0f}%)</div>'
            for l, v, c in t3ss_segments
        )

        visual_grid = f'''
        <div class="visual-grid" style="grid-template-columns:minmax(280px,480px);justify-content:start;">
            <div class="visual-card-big">
                <h3><i class="fas fa-syringe" style="color:#D32F2F;"></i> T3SS Pathotype Split</h3>
                <div class="donut-wrapper">
                    {donut}
                    <div class="legend-list">{legend}</div>
                </div>
            </div>
        </div>'''

        top_amr = sorted(all_amr, key=lambda x: x['count'], reverse=True)[:6]
        top_vir = sorted(all_vir, key=lambda x: x['count'], reverse=True)[:6]
        amr_bars = self._horizontal_bars(
            [(g['gene'], g['count'], '#E65100') for g in top_amr], total)
        vir_bars = self._horizontal_bars(
            [(g['gene'], g['count'], '#6A1B9A') for g in top_vir], total)

        top_genes_grid = f'''
        <div class="visual-grid" style="grid-template-columns:1fr 1fr;">
            <div class="visual-card-big">
                <h3><i class="fas fa-biohazard" style="color:#E65100;"></i> Top 6 AMR Genes</h3>
                {amr_bars}
            </div>
            <div class="visual-card-big">
                <h3><i class="fas fa-virus" style="color:#6A1B9A;"></i> Top 6 Virulence Genes</h3>
                {vir_bars}
            </div>
        </div>'''

        tour_items = [
            ('samples',      'Sample Overview',         'All typing on one row per isolate',              '#2E7D32', 'list-alt'),
            ('qc',           'FASTA QC',                'Assembly metrics + fastANI species check',      '#00695C', 'chart-line'),
            ('mlst',         'MLST',                    'Pasteur scheme — ST distribution & combinations', '#1565C0', 'code-branch'),
            ('serotype',     'O-Serotype',              'PAST serotyping & vaccine-relevant types',      '#6A1B9A', 'tag'),
            ('t3ss',         'ExoU/ExoS Pathotype',     'The primary virulence axis in P. aeruginosa',   '#D32F2F', 'syringe'),
            ('amr',          'AMR Genes',               'Multi-database resistance profiling',           '#C62828', 'biohazard'),
            ('virulence',    'Virulence',               'VFDB gene-centric tables with filters',         '#E65100', 'virus'),
            ('plasmid',      'Plasmid Replicons',       'Horizontal transfer markers',                   '#3F51B5', 'dna'),
            ('bacmet',       'Bacmet2',                 'Biocide & heavy-metal resistance',              '#607D8B', 'flask'),
            ('mutation',     'Mutations',               'Point mutations from AMRFinderPlus',            '#00BCD4', 'dna'),
            ('mge',          'Mobile Genetic Elements', 'mobileOG-db profiling by typing',               '#16A085', 'mobile-alt'),
            ('patterns',     'Patterns',                'Multi-way typing combinations',                 '#AD1457', 'project-diagram'),
            ('compare',      'Compare Isolates',        'Pairwise similarity & cluster detection',       '#0d9488', 'balance-scale'),
            ('citation',     'Citation',                'Full acknowledgements & suggested text',        '#00838F', 'book'),
            ('guide',        'Guide',                   'AI assistant guide + usage tips',               '#4527A0', 'question-circle'),
            ('export',       'Export',                  'Download tables as CSV / JSON',                 '#37474F', 'download'),
        ]
        tour_cards = ''.join(
            f'''<div class="tour-card" style="border-left-color:{c};" onclick="switchTab('{tid}')">
                <div class="tour-icon" style="background:{c};"><i class="fas fa-{icon}"></i></div>
                <div class="tour-body">
                    <div class="tour-title">{title}</div>
                    <div class="tour-desc">{desc}</div>
                </div>
            </div>'''
            for tid, title, desc, c, icon in tour_items
        )

        tour_section = f'''
        <h3 style="margin-top:36px;color:#5C2E0B;">
            <i class="fas fa-compass"></i> Tour the Report
        </h3>
        <p style="color:#666;font-size:.92em;margin-bottom:8px;">
            Click any card to jump directly to that tab.
        </p>
        <div class="tour-grid">{tour_cards}</div>'''

        stats_rows = f'''
        <tr><td>Total Isolates</td><td><strong>{total}</strong></td><td>Complete genomic analysis</td></tr>
        <tr><td>Unique Sequence Types</td><td><strong>{n_st}</strong></td><td>MLST Pasteur scheme</td></tr>
        <tr><td>Unique O-Types</td><td><strong>{n_o}</strong></td><td>PAST serotyping</td></tr>
        <tr><td>ExoU+ (cytotoxic)</td><td><span class="badge badge-exou">{n_exou}</span></td><td>Acute infection pathotype — higher mortality</td></tr>
        <tr><td>ExoS+ (invasive)</td><td><span class="badge badge-exos">{n_exos}</span></td><td>Chronic infection pathotype — CF persistence</td></tr>
        <tr><td>AMR Genes (unique)</td><td><strong>{n_amr_unique}</strong></td><td>Across all AMR databases</td></tr>
        <tr><td>Virulence Genes (unique)</td><td><strong>{n_vir_unique}</strong></td><td>From VFDB</td></tr>
        <tr><td>Plasmid Replicons</td><td><strong>{sum(len(g) for g in gene_centric.get('plasmid_databases', {}).values())}</strong></td><td>PlasmidFinder replicon families</td></tr>
        <tr><td>Bacmet2 Genes</td><td><strong>{sum(len(g) for g in gene_centric.get('bacmet_databases', {}).values())}</strong></td><td>Biocide & heavy-metal resistance</td></tr>
        <tr><td>Unique Point Mutations</td><td><strong>{n_mut}</strong></td><td>AMRFinderPlus mutation calls</td></tr>
        <tr><td>MGE-Profiled Isolates</td><td><strong>{n_mge}</strong></td><td>mobileOG-db protein family hits</td></tr>
        <tr><td>High-Risk AMR + Virulence Combos</td><td><span class="badge badge-critical">{n_high_risk}</span></td><td>Isolates with both critical AMR and virulence genes</td></tr>
        '''
        stats_section = f'''
        <div class="collapsible-section">
            <div class="collapsible-header">
                <span><i class="fas fa-chart-bar"></i> Full Statistics Table</span>
                <i class="fas fa-chevron-down collapsible-toggle"></i>
            </div>
            <div class="collapsible-content">
                <div class="scrollable-table">
                    <table class="data-table">
                        <thead><tr><th>Metric</th><th>Count</th><th>Details</th></tr></thead>
                        <tbody>{stats_rows}</tbody>
                    </table>
                </div>
            </div>
        </div>'''

        how_to = '''
        <div class="collapsible-section">
            <div class="collapsible-header">
                <span><i class="fas fa-lightbulb"></i> How to Use the Grouping Feature</span>
                <i class="fas fa-chevron-down collapsible-toggle"></i>
            </div>
            <div class="collapsible-content">
                <p style="color:#444;font-size:.95em;line-height:1.6;">
                    Every gene-centric table (AMR, Virulence, Bacmet2, Plasmids, Mutations) has a
                    <strong>"Group genomes by"</strong> dropdown above it. Choose any combination of
                    <strong>MLST &middot; O-type &middot; T3SS pathotype</strong>
                    and the genome column reorganises instantly. Click <strong>Reset</strong> to return to a flat list.
                </p>
                <p style="color:#444;font-size:.95em;line-height:1.6;margin-top:10px;">
                    <strong>Why this matters:</strong> you can instantly ask
                    <em>"Is blaVIM restricted to ST111?"</em> or
                    <em>"Do ExoU+ isolates cluster in specific O-types?"</em>
                    — the answer is one dropdown away.
                </p>
            </div>
        </div>'''

        return f'''
        {hero}
        <div class="kpi-grid">{kpis}</div>
        {visual_grid}
        {top_genes_grid}
        {tour_section}
        {stats_section}
        {how_to}
        '''

    def _sec_sample_overview(self, d):
        samples_data = d['samples']
        rows = ''
        for sample, data in sorted(samples_data.items()):
            t = data.get('typing', {})
            mlst = t.get('MLST', 'ND')
            otype = t.get('O_Type', 'ND')
            t3ss = t.get('T3SS', 'ND')

            vir_genes = []
            for db in ('vfdb', 'ecoli_vf', 'pa_vf', 'pseudomonas_vf'):
                vir_genes.extend(data.get('abricate_databases', {}).get(db, []))
            vir_count = len(vir_genes)
            mge_hits = data.get('mge_hits', 0)

            t3ss_class = 't3ss-exou' if t3ss == 'ExoU' else \
                         't3ss-exos' if t3ss == 'ExoS' else \
                         't3ss-both' if t3ss == 'ExoU+ExoS' else 't3ss-neither'
            t3ss_badge = f'<span class="typing-badge {t3ss_class}">{esc(t3ss)}</span>'

            vir_tags = ''.join(f'<span class="genome-tag">{esc(g)}</span>' for g in vir_genes)
            vir_details = (
                f'<details class="vir-details"><summary>{vir_count} gene(s)</summary>'
                f'<div class="vir-gene-tags">{vir_tags or "<em>None</em>"}</div></details>'
            )

            rows += f'''<tr>
                <td><strong>{esc(sample)}</strong></td>
                <td>{esc(mlst)}</td>
                <td>{esc(otype)}</td>
                <td>{t3ss_badge}</td>
                <td>{vir_details}</td>
                <td>{mge_hits}</td>
            </tr>'''

        return f'''
        {self._alert('info', 'fa-info-circle',
            '<h3>🧬 Sample Overview – Population Snapshot</h3>'
            '<p>Master reference table with all typing layers on a single row per isolate. '
            'Click the virulence column to expand full gene lists inline. '
            'Carbapenemase classification is on its own dedicated tab for accuracy.</p>')}
        <input type="text" class="search-box" id="search-samples-table"
               onkeyup="searchTable('samples-table','search-samples-table')"
               placeholder="🔍 Search samples by ID, ST, O-type, T3SS...">
        <div class="action-buttons">
            <button class="action-btn btn-primary"
                onclick="exportTableToCSV('samples-table','sample_overview.csv')">
                <i class="fas fa-download"></i> Export CSV</button>
            <button class="action-btn btn-light"
                onclick="document.getElementById('search-samples-table').value='';searchTable('samples-table','search-samples-table')">
                <i class="fas fa-sync"></i> Clear</button>
        </div>
        <div class="scrollable-table">
            <table id="samples-table" class="data-table">
                <thead><tr>
                    <th data-sort="string">Sample</th>
                    <th data-sort="string">MLST</th>
                    <th data-sort="string">O-Type</th>
                    <th data-sort="string">T3SS Pathotype</th>
                    <th data-sort="number">Virulence</th>
                    <th data-sort="number">MobileOG Hits</th>
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    def _sec_qc(self, d):
        qc_data = d.get('qc_data', {})
        if not qc_data:
            return self._alert('warning', 'fa-exclamation-circle',
                               '<h3>No QC Data Available</h3>'
                               '<p>The FASTA_QC_summary.html file was not found.</p>')

        credit = self._credit_bar('#00695C', '📊', 'FASTA QC & Species Confirmation',
            'Assembly metrics computed with <strong>Biopython</strong>. '
            'Species confirmation performed with <strong>fastANI</strong> '
            '(Jain et al., <em>Nat Commun</em> 2018) against a curated P. aeruginosa + close-relative panel.')

        all_metrics = set()
        for m in qc_data.values():
            all_metrics.update(m.keys())
        metric_list = sorted(all_metrics)

        species_cols = [c for c in metric_list if any(k in c for k in
                        ('Species', 'ANI', 'Confirmed', 'Contamination'))]
        numeric_cols = [c for c in metric_list if c not in species_cols]

        header_cells = ''.join(f'<th data-sort="number">{esc(m)}</th>' for m in numeric_cols)
        species_header = ''.join(f'<th data-sort="string">{esc(m)}</th>' for m in species_cols)

        rows = ''
        for sample, metrics in sorted(qc_data.items()):
            cells = ''
            for m in numeric_cols:
                v = metrics.get(m, 'ND')
                if isinstance(v, float):
                    v = f"{v:,.0f}" if v > 1000 else f"{v:.2f}"
                cells += f'<td>{v}</td>'
            species_cells = ''
            for m in species_cols:
                v = metrics.get(m, 'ND')
                v_str = str(v)
                if 'Confirmed' in m:
                    if v_str == 'Yes':
                        species_cells += '<td><span style="color:#28a745;font-weight:bold;">✓ Yes</span></td>'
                    elif v_str == 'No':
                        species_cells += '<td><span style="color:#C62828;font-weight:bold;">✗ No</span></td>'
                    else:
                        species_cells += f'<td>{esc(v)}</td>'
                elif m == 'Contamination Suspected':
                    if v_str == 'Yes':
                        species_cells += '<td><span style="color:#e67e22;font-weight:bold;">⚠️ Yes</span></td>'
                    else:
                        species_cells += '<td>No</td>'
                else:
                    species_cells += f'<td>{esc(v)}</td>'
            rows += f'<tr><td><strong>{esc(sample)}</strong></td>{cells}{species_cells}</tr>'

        return f'''
        {credit}
        {self._alert('info', 'fa-chart-line',
            '<h3>📏 FASTA Quality Control + Species Confirmation</h3>'
            '<p><strong>Assembly metrics:</strong> Number of contigs (lower better, &lt;200 typical draft, &lt;100 excellent), N50 (higher better, &gt;50kb good), GC% (P. aeruginosa ~60–70%), total length (~5–8 Mb expected).</p>'
            '<p><strong>Species confirmation:</strong> fastANI computes ANI against a reference panel. A genome is confirmed when the best match is P. aeruginosa at ≥95% ANI. Contamination flagged when second-best exceeds 90% ANI, is within 5% of the best, and is not P. aeruginosa.</p>')}
        <input type="text" class="search-box" id="search-qc-table"
               onkeyup="searchTable('qc-table','search-qc-table')" placeholder="🔍 Search sample...">
        <div class="action-buttons">
            <button class="action-btn btn-primary"
                onclick="exportTableToCSV('qc-table','fasta_qc.csv')">
                <i class="fas fa-download"></i> Export QC Data</button>
        </div>
        <div class="master-scrollable-container">
            <table id="qc-table" class="data-table">
                <thead><tr>
                    <th data-sort="string">Sample</th>
                    {header_cells}{species_header}
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    def _sec_mlst(self, d):
        P = d['patterns']
        credit = self._credit_bar('#1565C0', '🧬', 'MLST Typing (Pasteur scheme)',
            'Powered by <a href="https://github.com/tseemann/mlst" target="_blank" '
            'style="color:#1565C0;font-weight:bold;">Prof. Torsten Seemann\'s mlst</a> '
            'and the <a href="https://pubmlst.org/" target="_blank" style="color:#1565C0;font-weight:bold;">'
            'PubMLST database</a> (Jolley et al., <em>Wellcome Open Res</em> 2018).<br>'
            '<span style="font-size:0.9em;color:#6c757d;"><i class="fas fa-info-circle"></i> '
            'Due to licensing on the PubMLST database, alleles added after December 31, 2024 '
            'may not be included.</span>')
        return f'''
        {credit}
        {self._alert('info', 'fa-code-branch',
            '<h3>🔬 MLST (Pasteur scheme) – Tracking High-Risk Clones</h3>'
            '<p>MLST uses seven housekeeping genes (<em>acsA, aroE, guaA, mutL, nuoD, ppsA, trpE</em>). '
            'Globally disseminated high-risk clones include ST111 (VIM-2, ExoU), ST175 (VIM-2, XDR), '
            'ST235 (KPC/NDM/VIM, ExoU), ST244, ST277 (SPM-1, Brazil), ST395 (KPC, ExoS).</p>'
            f'<p><strong>{len(P.get("mlst_distribution", {}))} unique STs</strong> identified.</p>')}
        <h3>📊 ST Distribution</h3>
        {self._distribution_table('mlst-dist-table', P.get('mlst_distribution', Counter()), 'ST')}
        {self._combo_table('mlst-otype-table', 'ST – O-Type Combination', P.get('mlst_otype', {}))}
        {self._combo_table('mlst-t3ss-table', 'ST – T3SS Pathotype Combination', P.get('mlst_t3ss', {}))}'''

    def _sec_serotype(self, d):
        P = d['patterns']
        o_dist = P.get('o_type_distribution', Counter())
        credit = self._credit_bar('#6A1B9A', '🏷️', 'O-Serotyping (PAST)',
            'Powered by <strong>pasty</strong> (Robert A. Petit III) — a modern, '
            'standalone wrapper for the original PAST tool '
            '(<a href="https://cge.food.dtu.dk/services/PAst/" target="_blank" '
            'style="color:#6A1B9A;font-weight:bold;">CGE, Technical University of Denmark</a>). '
            'The O-antigen is a major surface antigen and an investigational vaccine target.')
        return f'''
        {credit}
        {self._alert('info', 'fa-tag',
            '<h3>🏷️ O-Serotype Distribution</h3>'
            '<p>Common clinical serotypes include O1, O2, O5, O6, O11. O11 is strongly associated with epidemic multidrug-resistant clones (e.g., ST235, ST111).</p>'
            f'<p><strong>{len(o_dist)} unique O-types</strong> identified.</p>')}
        <h3>📊 O-Type Distribution</h3>
        {self._distribution_table('otype-dist-table', o_dist, 'O-Type')}
        {self._combo_table('otype-t3ss-table', 'O-Type – T3SS Pathotype', P.get('otype_t3ss', {}))}'''

    def _sec_t3ss(self, d):
        P = d['patterns']
        samples_data = d['samples']
        t3ss_dist = P.get('t3ss_distribution', Counter())
        total_samples = len(samples_data)

        effector_counts = Counter()
        apparatus_counts = Counter()
        regulator_counts = Counter()
        for sample, data in samples_data.items():
            genes = data.get('amrfinder', {}).get('all_genes', [])
            for db in ('vfdb', 'ecoli_vf', 'pa_vf', 'pseudomonas_vf'):
                genes.extend(data.get('abricate_databases', {}).get(db, []))
            genes_lower = [str(g).lower() for g in genes]
            for eff in ('exou', 'exos', 'exot', 'exoy', 'toxa'):
                if any(eff in g for g in genes_lower):
                    effector_counts[eff] += 1
            for struct in ('pscc', 'pscj', 'pcrv', 'pcrh', 'pcrg'):
                if any(struct in g for g in genes_lower):
                    apparatus_counts[struct] += 1
            for reg in ('exsa', 'exsc', 'exsd'):
                if any(reg in g for g in genes_lower):
                    regulator_counts[reg] += 1

        credit = self._credit_bar('#D32F2F', '💉',
            'Type III Secretion System (T3SS) – The Primary Virulence Axis',
            'The T3SS is <em>the</em> defining virulence determinant of <em>Pseudomonas aeruginosa</em>. '
            'Every clinical isolate is essentially one of two pathotypes: <strong>ExoU+</strong> (cytotoxic, acute) '
            'or <strong>ExoS+</strong> (invasive, chronic). The effectors are mutually exclusive and clinically actionable.')

        effector_rows = ''
        for eff, name, role in [
            ('exou', 'ExoU', 'Phospholipase A2 — cytotoxic, acute infections, higher mortality in bacteraemia and VAP'),
            ('exos', 'ExoS', 'ADP-ribosyltransferase — invasive, chronic infections, CF persistence'),
            ('exot', 'ExoT', 'ADP-ribosyltransferase — anti-phagocytic, nearly universal'),
            ('exoy', 'ExoY', 'Adenylate cyclase — disrupts host cytoskeleton, nearly universal'),
            ('toxa', 'Exotoxin A', 'EF-2 inhibitor — systemic toxicity, tissue damage, nearly universal'),
        ]:
            count = effector_counts.get(eff, 0)
            pct = (count / total_samples * 100) if total_samples else 0
            badge = ''
            if eff == 'exou':
                badge = '<span class="badge badge-exou">Cytotoxic</span>'
            elif eff == 'exos':
                badge = '<span class="badge badge-exos">Invasive</span>'
            effector_rows += (f'<tr><td><strong>{name}</strong></td><td>{badge}</td>'
                              f'<td>{count} / {total_samples} ({pct:.1f}%)</td>'
                              f'<td>{role}</td></tr>')

        apparatus_rows = ''
        for s, name in [('pscc', 'PscC'), ('pscj', 'PscJ'), ('pcrv', 'PcrV'),
                        ('pcrh', 'PcrH'), ('pcrg', 'PcrG')]:
            count = apparatus_counts.get(s, 0)
            pct = (count / total_samples * 100) if total_samples else 0
            apparatus_rows += (f'<tr><td><strong>{name}</strong></td>'
                               f'<td>{count} / {total_samples} ({pct:.1f}%)</td></tr>')

        regulator_rows = ''
        for r, name in [('exsa', 'ExsA'), ('exsc', 'ExsC'), ('exsd', 'ExsD')]:
            count = regulator_counts.get(r, 0)
            pct = (count / total_samples * 100) if total_samples else 0
            regulator_rows += (f'<tr><td><strong>{name}</strong></td>'
                               f'<td>{count} / {total_samples} ({pct:.1f}%)</td></tr>')

        clinical_box = '''
        <div class="alert-box" style="border-left-color:#D32F2F; background:#fee2e2; border-radius:8px; padding:20px; margin:20px 0;">
            <i class="fas fa-stethoscope fa-2x" style="color:#D32F2F;"></i>
            <div>
                <h3 style="color:#7f1d1d; margin-bottom:10px;">🩺 ExoU vs ExoS — Why It Matters Clinically</h3>
                <p style="margin:6px 0;"><strong>ExoU+ (Cytotoxic pathotype)</strong></p>
                <ul style="margin:6px 0 12px 20px;">
                    <li>Phospholipase A2 activity → rapid cell lysis, tissue destruction</li>
                    <li>Associated with <strong>acute infections</strong>: ventilator-associated pneumonia, bacteraemia, keratitis</li>
                    <li>Higher <strong>mortality</strong> in bloodstream infections and ICU settings</li>
                    <li>Often found in <strong>MDR/XDR</strong> clones (ST235, ST111)</li>
                </ul>
                <p style="margin:6px 0;"><strong>ExoS+ (Invasive pathotype)</strong></p>
                <ul style="margin:6px 0 12px 20px;">
                    <li>ADP-ribosyltransferase activity → host cell invasion, immune evasion</li>
                    <li>Associated with <strong>chronic infections</strong>: cystic fibrosis, chronic wounds</li>
                    <li>Better <strong>biofilm formation</strong> and persistence</li>
                    <li>Common in environmental isolates and early CF colonisation</li>
                </ul>
                <p style="margin:6px 0;"><strong>Rare: neither ExoU nor ExoS detected</strong> — could indicate an assembly gap in the T3SS locus, an uncommon T3SS variant, or a genuinely effectors-negative lineage (rare).</p>
            </div>
        </div>'''

        return f'''
        {credit}
        {self._alert('info', 'fa-syringe',
            '<h3>💉 ExoU/ExoS Pathotype – The T3SS Virulence Axis</h3>'
            '<p>The T3SS injects effector proteins directly into host cells. Two effectors are mutually exclusive '
            'and define the pathotype: <strong>ExoU</strong> (cytotoxic) and <strong>ExoS</strong> (invasive). '
            'This is the single most important virulence feature in P. aeruginosa epidemiology.</p>')}
        {clinical_box}
        <h3>📊 Pathotype Distribution</h3>
        {self._distribution_table('t3ss-dist-table', t3ss_dist, 'T3SS Pathotype')}
        <h3 style="margin-top:30px;">📈 T3SS Effectors – Population Prevalence</h3>
        <div class="scrollable-table">
            <table class="data-table">
                <thead><tr>
                    <th>Effector</th><th>Pathotype</th><th>Prevalence</th><th>Role</th>
                </tr></thead>
                <tbody>{effector_rows}</tbody>
            </table>
        </div>
        <h3 style="margin-top:30px;">🔧 T3SS Apparatus – Population Prevalence</h3>
        <div class="scrollable-table">
            <table class="data-table">
                <thead><tr><th>Component</th><th>Prevalence</th></tr></thead>
                <tbody>{apparatus_rows}</tbody>
            </table>
        </div>
        <h3 style="margin-top:30px;">🎛️ T3SS Regulators – Population Prevalence</h3>
        <div class="scrollable-table">
            <table class="data-table">
                <thead><tr><th>Regulator</th><th>Prevalence</th></tr></thead>
                <tbody>{regulator_rows}</tbody>
            </table>
        </div>
        <h3 style="margin-top:30px;">🔗 Cross-Typing Combinations</h3>
        {self._combo_table('st-t3ss-table', 'ST – T3SS Pathotype', P.get('mlst_t3ss', {}))}
        {self._combo_table('otype-t3ss-table2', 'O-Type – T3SS Pathotype', P.get('otype_t3ss', {}))}'''

    def _sec_amr(self, d):
        gene_centric = d['gene_centric']
        amr_dbs = gene_centric.get('amr_databases', {})
        total_samples = len(d['samples'])
        all_genes = []
        for genes in amr_dbs.values():
            all_genes.extend(genes)
        all_genes.sort(key=lambda x: x['count'], reverse=True)

        credit = self._credit_bar('#C62828', '🧬', 'AMR Detection Tools & Databases',
            '<strong>ABRicate</strong> by '
            '<a href="https://github.com/tseemann/abricate" target="_blank" style="color:#C62828;font-weight:bold;">Prof. Torsten Seemann</a> '
            '→ CARD, ARG-ANNOT, NCBI AMR, MEGARes, ResFinder.<br>'
            '<strong>AMRFinderPlus</strong> by '
            '<a href="https://github.com/ncbi/amr" target="_blank" style="color:#C62828;font-weight:bold;">NCBI</a> '
            '→ comprehensive resistance gene + point-mutation detection.')

        why_multi_db = '''
        <div class="alert-box" style="border-left-color:#00695c; background:#e8f5e9; border-radius:8px; padding:18px 22px; margin:20px 0;">
            <div style="display:flex; gap:15px; align-items:flex-start;">
                <i class="fas fa-info-circle fa-2x" style="color:#00695c; margin-top:3px;"></i>
                <div>
                    <h4 style="margin:0 0 10px 0; color:#00695c; font-size:1.1em;">🔬 Why Do We Use Multiple Databases?</h4>
                    <p style="margin:6px 0; font-size:0.95em; line-height:1.5;">PseudoScope screens every genome against <strong>five independent AMR databases</strong> plus <strong>AMRFinderPlus</strong> — because <strong>no single database is comprehensive</strong>.</p>
                    <ul style="margin:6px 0 10px 20px; font-size:0.93em; line-height:1.6;">
                        <li>Some researchers pick one database and only report hits from that one. This is <strong>fast but incomplete</strong>.</li>
                        <li>Conversely, single-database hits with weak support can be <strong>false positives</strong> that a second database would have flagged.</li>
                        <li>The result: <em>biased prevalence estimates</em> and <em>missed resistance signals</em>.</li>
                    </ul>
                    <p style="margin:10px 0 6px 0; font-size:0.95em; line-height:1.5;"><strong>✅ Our approach — report everything, provenance preserved:</strong></p>
                    <ul style="margin:6px 0 10px 20px; font-size:0.93em; line-height:1.6;">
                        <li>All hits from all databases are kept <strong>separate and unmerged</strong>.</li>
                        <li>Cross-database agreement becomes a <strong>confidence signal</strong>.</li>
                    </ul>
                </div>
            </div>
        </div>
        '''

        intrinsic_primer = '''
        <div class="alert-box" style="border-left-color:#6f42c1; background:#f3e8ff; border-radius:8px; padding:20px; margin:20px 0;">
            <i class="fas fa-shield-alt fa-2x" style="color:#6f42c1;"></i>
            <div>
                <h4 style="margin:0 0 12px 0; color:#6f42c1; font-size:1.1em;">🧬 P. aeruginosa's Triple Resistance Architecture</h4>
                <ul style="margin:10px 0 10px 20px; font-size:0.93em; line-height:1.6;">
                    <li><strong style="color:#7c3aed;">1. Intrinsic resistance</strong> — chromosomal ampC, OprD porin loss, 12 RND efflux systems, low outer-membrane permeability.</li>
                    <li><strong style="color:#e11d48;">2. Acquired resistance</strong> — carbapenemases, ESBLs, 16S rRNA methylases, mcr genes. See the <strong>Carbapenemase & OXA</strong> tab for accurate classification.</li>
                    <li><strong style="color:#e67e22;">3. Adaptive resistance</strong> — biofilm formation, persister cells, small-colony variants.</li>
                </ul>
                <p style="margin:8px 0 0 0; font-size:0.92em; background:#fff3cd; padding:8px 14px; border-radius:4px; border-left:3px solid #ffc107;">
                    <i class="fas fa-lightbulb" style="color:#856404;"></i>
                    <strong>Why this matters:</strong> A P. aeruginosa isolate can be MDR <em>without</em> any acquired resistance gene. Conversely, a susceptible-looking isolate can become resistant overnight by acquiring a single plasmid.
                </p>
            </div>
        </div>
        '''

        confidence_box = '''
        <div class="alert-box" style="border-left-color:#0891b2; background:#e0f2fe; border-radius:8px; padding:16px 20px; margin:20px 0;">
            <i class="fas fa-layer-group fa-2x" style="color:#0891b2;"></i>
            <div>
                <h4 style="margin:0 0 10px 0; color:#0891b2; font-size:1.05em;">🎯 Cross-Database Confidence Tiers</h4>
                <ul style="margin:6px 0 0 20px; font-size:0.93em; line-height:1.6;">
                    <li><span style="color:#16a34a; font-weight:bold;">🟢 High confidence</span> — gene found in <strong>3 or more databases</strong>.</li>
                    <li><span style="color:#f59e0b; font-weight:bold;">🟡 Moderate confidence</span> — gene found in <strong>2 databases</strong>.</li>
                    <li><span style="color:#C62828; font-weight:bold;">🔴 Low confidence / investigate</span> — gene found in <strong>only 1 database</strong>.</li>
                </ul>
            </div>
        </div>
        '''

        caveat = '''
        <div class="alert-box alert-warning" style="border-left-color:#ffc107;">
            <i class="fas fa-exclamation-triangle fa-2x"></i>
            <div>
                <h4 style="margin:0 0 8px 0; color:#856404;">⚠️ Gene Presence ≠ Phenotypic Resistance</h4>
                <ul style="margin:6px 0 12px 20px; font-size:0.92em; line-height:1.6;">
                    <li><strong>Silent / truncated genes</strong> — blaVIM can be present but non-functional.</li>
                    <li><strong>Expression regulation</strong> — MexXY is inducible; AmpC derepressed only when AmpR is mutated.</li>
                    <li><strong>Porin loss matters</strong> — carbapenem resistance often requires both a carbapenemase AND OprD loss.</li>
                    <li><strong>Naming ambiguity</strong> — ResFinder appends _1 to primary alleles.</li>
                </ul>
                <p style="background:#fff8e1; padding:8px 12px; border-radius:4px; border-left:3px solid #f59e0b; font-size:0.92em;">
                    <strong>Clinical bottom line:</strong> AST remains the gold standard for treatment decisions.
                </p>
            </div>
        </div>
        '''

        filter_buttons = self._filter_buttons('amr-table', [
            ('Carbapenemases (KPC)', 'KPC', 'btn-danger', 'fa-biohazard'),
            ('Carbapenemases (NDM)', 'NDM', 'btn-danger', 'fa-biohazard'),
            ('Carbapenemases (VIM)', 'VIM', 'btn-danger', 'fa-biohazard'),
            ('Carbapenemases (IMP)', 'IMP', 'btn-danger', 'fa-biohazard'),
            ('GES', 'GES', 'btn-warning', 'fa-biohazard'),
            ('SPM', 'SPM', 'btn-warning', 'fa-biohazard'),
            ('OXA-48-like', 'OXA-48', 'btn-danger', 'fa-biohazard'),
            ('OXA-23-like', 'OXA-23', 'btn-danger', 'fa-biohazard'),
            ('OXA-24/40-like', 'OXA-24', 'btn-danger', 'fa-biohazard'),
            ('OXA-58-like', 'OXA-58', 'btn-danger', 'fa-biohazard'),
            ('OXA-2-like', 'OXA-2', 'btn-warning', 'fa-biohazard'),
            ('OXA-10-like', 'OXA-10', 'btn-warning', 'fa-biohazard'),
            ('OXA-50 (intrinsic)', 'OXA-50', 'btn-secondary', 'fa-info-circle'),
            ('ESBL (PER)', 'PER', 'btn-warning', 'fa-biohazard'),
            ('ESBL (VEB)', 'VEB', 'btn-warning', 'fa-biohazard'),
            ('Colistin (mcr)', 'mcr', 'btn-danger', 'fa-biohazard'),
            ('16S methylase (armA)', 'armA', 'btn-warning', 'fa-biohazard'),
            ('16S methylase (rmt)', 'rmt', 'btn-warning', 'fa-biohazard'),
            ('Aminoglycosides (aac)', 'aac', 'btn-info', 'fa-syringe'),
            ('Efflux pumps (mex)', 'mex', 'btn-secondary', 'fa-filter'),
        ])

        info_block = self._gene_family_info('#C62828',
            'Role of each AMR gene family in <em>P. aeruginosa</em>:', [
                ('blaKPC/NDM/VIM/IMP', 'Carbapenemases — destroy last-resort carbapenems. Often plasmid-borne.'),
                ('blaGES/PER/VEB', 'ESBLs — extended-spectrum β-lactamases.'),
                ('blaOXA-48/23/24/58-like', 'Acquired OXA carbapenemases — clinically decisive.'),
                ('blaOXA-2-like / blaOXA-10-like', 'Narrow to extended-spectrum Class D. ESBL variants affect ceftazidime/cefepime.'),
                ('blaOXA-50-like', 'Intrinsic chromosomal. Not an acquired resistance marker.'),
                ('mcr-1 to mcr-10', 'Plasmid-mediated colistin resistance.'),
                ('armA / rmtA-F', '16S rRNA methyltransferases — high-level pan-aminoglycoside resistance.'),
                ('aac/aph/ant', 'Aminoglycoside-modifying enzymes.'),
                ('mexAB-oprM', 'Intrinsic RND efflux system.'),
                ('mexXY', 'Inducible RND efflux — aminoglycosides and zwitterionic β-lactams.'),
                ('gyrA/parC mutations', 'QRDR mutations — fluoroquinolone resistance (see Mutations tab).'),
            ])

        return f'''
        {credit}
        {self._alert('info', 'fa-biohazard',
            '<h3>🧬 AMR Genes – Gene-Centric View + Grouping</h3>'
            '<p>Each gene is shown with all genomes carrying it. Use the grouping dropdown to reveal clone associations.</p>')}
        {filter_buttons}
        {self._gene_table('amr-table', 'All AMR Genes Across Databases',
                          all_genes, total_samples, self.analyzer.critical_amr_genes)}
        <details class="edu-details" style="margin-top:30px;background:#fdf6ec;border-left:4px solid #00695C;padding:14px 18px;border-radius:8px;box-shadow:0 2px 8px rgba(139,69,19,.06);">
            <summary style="cursor:pointer;font-weight:700;color:#5C2E0B;font-size:1.05em;">
                <i class="fas fa-book-medical" style="color:#00695C;"></i> Interpretation Guide — click to expand
            </summary>
            <div style="padding-top:16px;">
                {why_multi_db}
                {intrinsic_primer}
                {confidence_box}
                {caveat}
                {info_block}
            </div>
        </details>
        <details class="edu-details" style="margin-top:16px;background:#fdf6ec;border-left:4px solid #C62828;padding:14px 18px;border-radius:8px;box-shadow:0 2px 8px rgba(139,69,19,.06);">
            <summary style="cursor:pointer;font-weight:700;color:#5C2E0B;font-size:1.05em;">
                <i class="fas fa-database" style="color:#C62828;"></i> Database Provenance — click to expand
            </summary>
            <div style="padding-top:16px;">
                {self._database_cards(amr_dbs, {'amrfinder': 'AMRfinder'})}
            </div>
        </details>'''

    def _sec_virulence(self, d):
        gene_centric = d['gene_centric']
        vir_dbs = gene_centric.get('virulence_databases', {})
        total_samples = len(d['samples'])
        all_genes = []
        for genes in vir_dbs.values():
            all_genes.extend(genes)
        all_genes.sort(key=lambda x: x['count'], reverse=True)

        credit = self._credit_bar('#E65100', '🦠', 'Virulence Detection',
            '<strong>ABRicate</strong> by Prof. Torsten Seemann → '
            '<a href="http://www.mgc.ac.cn/VFs/" target="_blank" style="color:#E65100;font-weight:bold;">'
            'VFDB (Virulence Factor Database)</a>.')

        filter_buttons = self._filter_buttons('vir-table', [
            ('T3SS effectors (ExoU)', 'exoU', 'btn-danger', 'fa-syringe'),
            ('T3SS effectors (ExoS)', 'exoS', 'btn-danger', 'fa-syringe'),
            ('T3SS effectors (ExoT)', 'exoT', 'btn-warning', 'fa-syringe'),
            ('T3SS effectors (ExoY)', 'exoY', 'btn-warning', 'fa-syringe'),
            ('Exotoxin A (toxA)', 'toxA', 'btn-danger', 'fa-biohazard'),
            ('Phospholipase C (plcH)', 'plcH', 'btn-info', 'fa-filter'),
            ('Elastase (lasB)', 'lasB', 'btn-warning', 'fa-filter'),
            ('Elastase (lasA)', 'lasA', 'btn-warning', 'fa-filter'),
            ('Protease (aprA)', 'aprA', 'btn-info', 'fa-filter'),
            ('Pyocyanin (phz)', 'phz', 'btn-info', 'fa-filter'),
            ('Alginate (alg)', 'alg', 'btn-info', 'fa-filter'),
            ('Pyoverdine (pvd)', 'pvd', 'btn-info', 'fa-filter'),
            ('Type IV pili (pil)', 'pil', 'btn-secondary', 'fa-filter'),
            ('Flagella (fli)', 'fli', 'btn-secondary', 'fa-filter'),
            ('T6SS (hcp)', 'hcp', 'btn-secondary', 'fa-filter'),
            ('T6SS (vgrG)', 'vgrG', 'btn-secondary', 'fa-filter'),
        ])

        info_block = self._gene_family_info('#E65100',
            'Role of each virulence family in <em>P. aeruginosa</em>:', [
                ('exoU, exoS, exoT, exoY', 'T3SS effectors — ExoU is cytotoxic, ExoS is invasive, ExoT/ExoY modulatory.'),
                ('toxA (exotoxin A)', 'EF-2 ADP-ribosylation — inhibits host protein synthesis.'),
                ('plcH, plcN', 'Phospholipase C — tissue destruction.'),
                ('lasA, lasB', 'Elastases — degrade elastin and collagen.'),
                ('aprA', 'Alkaline protease — immune evasion.'),
                ('phz (pyocyanin)', 'Blue phenazine pigment — ROS generation, tissue damage.'),
                ('alg (alginate)', 'Alginate biosynthesis — mucoid phenotype, chronic CF.'),
                ('pvd (pyoverdine)', 'Siderophore — iron acquisition, biofilm signalling.'),
                ('pil (Type IV pili)', 'Twitching motility, adhesion, DNA uptake.'),
                ('fli (flagella)', 'Motility, adhesion, immune stimulation.'),
                ('T6SS (hcp, vgrG)', 'Type VI Secretion System — interbacterial competition.'),
            ])

        return f'''
        {credit}
        {self._alert('info', 'fa-virus',
            '<h3>🦠 Virulence Factors – Gene-Centric View</h3>'
            '<p>Use the grouping dropdown to see which clones carry specific virulence factors.</p>')}
        {filter_buttons}
        {info_block}
        {self._gene_table('vir-table', 'All Virulence Genes',
                          all_genes, total_samples, self.analyzer.critical_virulence_genes)}
        {self._database_cards(vir_dbs)}'''

    def _sec_plasmid(self, d):
        gene_centric = d['gene_centric']
        plasmid_dbs = gene_centric.get('plasmid_databases', {})
        total_samples = len(d['samples'])
        all_genes = []
        for genes in plasmid_dbs.values():
            all_genes.extend(genes)
        all_genes.sort(key=lambda x: x['count'], reverse=True)

        credit = self._credit_bar('#3F51B5', '🔌', 'Plasmid Replicon Typing',
            '<strong>ABRicate</strong> → '
            '<a href="https://genepi.food.dtu.dk/PlasmidFinder/" target="_blank" style="color:#3F51B5;font-weight:bold;">'
            'PlasmidFinder</a> from the Center for Genomic Epidemiology (DTU).')

        filter_buttons = self._filter_buttons('plasmid-table', [
            ('IncP', 'IncP', 'btn-info', 'fa-plug'),
            ('IncQ', 'IncQ', 'btn-info', 'fa-plug'),
            ('IncU', 'IncU', 'btn-info', 'fa-plug'),
            ('IncHI', 'IncHI', 'btn-secondary', 'fa-plug'),
            ('Col', 'Col', 'btn-secondary', 'fa-plug'),
        ])

        info_block = self._gene_family_info('#3F51B5',
            'Plasmid replicon families in <em>P. aeruginosa</em>:', [
                ('IncP-2', 'Large, broad-host-range plasmids — often carry carbapenemases (blaIMP, blaVIM).'),
                ('IncQ', 'Broad-host-range plasmids — aminoglycoside resistance, mobQ family.'),
                ('IncU', 'Common in P. aeruginosa MDR plasmids — variable resistance cargo.'),
                ('IncHI', 'Large conjugative plasmids — heavy metal and antibiotic resistance.'),
                ('Col', 'Small colicin plasmids — sometimes carry toxin-antitoxin systems.'),
                ('pUM505', 'Well-known P. aeruginosa plasmid family — carries mer operon and resistance genes.'),
            ])

        return f'''
        {credit}
        {self._alert('info', 'fa-dna',
            '<h3>🔌 Plasmid Replicons – Horizontal Gene Transfer Markers</h3>'
            '<p>Replicon families reveal which clones carry which plasmids.</p>')}
        {filter_buttons}
        {info_block}
        {self._gene_table('plasmid-table', 'Plasmid Replicons', all_genes, total_samples)}
        {self._database_cards(plasmid_dbs)}'''

    def _sec_bacmet(self, d):
        gene_centric = d['gene_centric']
        bacmet_dbs = gene_centric.get('bacmet_databases', {})
        total_samples = len(d['samples'])
        all_genes = []
        for genes in bacmet_dbs.values():
            all_genes.extend(genes)
        all_genes.sort(key=lambda x: x['count'], reverse=True)

        credit = self._credit_bar('#607D8B', '🧪', 'BacMet2 – Biocide & Heavy Metal Resistance',
            '<strong>ABRicate</strong> → '
            '<a href="https://bacmet.biomedicine.gu.se/" target="_blank" style="color:#607D8B;font-weight:bold;">'
            'BacMet2 database</a>.')

        filter_buttons = self._filter_buttons('bacmet-table', [
            ('Biocides (qac)', 'qac', 'btn-info', 'fa-filter'),
            ('Chlorhexidine (cep)', 'cep', 'btn-info', 'fa-filter'),
            ('Formaldehyde (form)', 'form', 'btn-info', 'fa-filter'),
            ('Mercury (mer)', 'mer', 'btn-warning', 'fa-filter'),
            ('Arsenic (ars)', 'ars', 'btn-warning', 'fa-filter'),
            ('Silver (sil)', 'sil', 'btn-warning', 'fa-filter'),
            ('Copper (cop)', 'cop', 'btn-warning', 'fa-filter'),
            ('Chromium (chr)', 'chr', 'btn-warning', 'fa-filter'),
            ('Cadmium (cad)', 'cad', 'btn-warning', 'fa-filter'),
            ('Zinc (znt)', 'znt', 'btn-warning', 'fa-filter'),
            ('Cobalt-Zinc-Cadmium (czc)', 'czc', 'btn-warning', 'fa-filter'),
            ('Nickel (nik)', 'nik', 'btn-warning', 'fa-filter'),
            ('Lead (pbr)', 'pbr', 'btn-warning', 'fa-filter'),
            ('Oxidative stress (soxR)', 'soxR', 'btn-secondary', 'fa-filter'),
            ('Envelope stress (cpxR)', 'cpxR', 'btn-secondary', 'fa-filter'),
            ('Efflux (mex)', 'mex', 'btn-secondary', 'fa-filter'),
        ])

        info_block = self._gene_family_info('#607D8B',
            'Environmental co-selection markers in <em>P. aeruginosa</em>:', [
                ('qac family', 'Quaternary ammonium compounds — hospital disinfectants.'),
                ('cep', 'Chlorhexidine resistance — antiseptic tolerance.'),
                ('form', 'Formaldehyde resistance.'),
                ('mer', 'Mercury resistance — often co-located with antibiotic resistance.'),
                ('ars', 'Arsenic resistance.'),
                ('sil / cop', 'Silver and copper resistance.'),
                ('czc / cad / znt', 'Cobalt-zinc-cadmium efflux — heavy metal co-selection.'),
                ('soxR / cpxR', 'Stress regulators that upregulate multidrug efflux pumps.'),
            ])

        return f'''
        {credit}
        {self._alert('info', 'fa-flask',
            '<h3>🧪 BacMet2 – Biocide & Heavy Metal Resistance</h3>'
            '<p>Environmental co-selection markers. Critical for infection control.</p>')}
        {filter_buttons}
        {info_block}
        {self._gene_table('bacmet-table', 'All BacMet2 Genes', all_genes, total_samples)}
        {self._database_cards(bacmet_dbs)}'''

    def _sec_mutation(self, d):
        mut = d.get('mutation_data', {})
        mutations = mut.get('mutations', [])
        total_samples = len(d['samples'])
        if not mutations:
            return self._alert('warning', 'fa-dna',
                '<h3>No Mutation Data Available</h3>'
                '<p>The mutation_summary.html file was not found.</p>')

        credit = self._credit_bar('#00BCD4', '🧬', 'AMRFinderPlus – Point Mutations',
            '<strong>AMRFinderPlus</strong> by NCBI — comprehensive database of resistance genes and point mutations.')

        filter_buttons = self._filter_buttons('mutation-table', [
            ('Quinolone (gyrA)', 'gyrA', 'btn-warning', 'fa-biohazard'),
            ('Quinolone (parC)', 'parC', 'btn-warning', 'fa-biohazard'),
            ('Rifampin (rpoB)', 'rpoB', 'btn-warning', 'fa-biohazard'),
            ('23S rRNA', '23S', 'btn-info', 'fa-filter'),
            ('mprF', 'mprF', 'btn-info', 'fa-filter'),
            ('rplC', 'rplC', 'btn-info', 'fa-filter'),
            ('rplD', 'rplD', 'btn-info', 'fa-filter'),
            ('oprD', 'oprD', 'btn-info', 'fa-filter'),
        ])

        info_block = self._gene_family_info('#00BCD4',
            'Clinical relevance of key mutations in <em>P. aeruginosa</em>:', [
                ('gyrA / parC (quinolones)', 'QRDR mutations — reduce fluoroquinolone susceptibility.'),
                ('rpoB (rifampin)', 'Rifampin resistance.'),
                ('23S rRNA', 'Macrolide and linezolid resistance.'),
                ('mprF', 'Daptomycin non-susceptibility.'),
                ('oprD', 'Porin loss — carbapenem resistance independent of carbapenemases.'),
                ('mexR / nalC / nalD', 'Efflux pump regulators — MexAB-OprM overexpression.'),
                ('ampR / ampD', 'AmpC derepression.'),
            ])

        rows = ''
        for m in mutations:
            tags = ''.join(f'<span class="genome-tag">{esc(g)}</span>' for g in m['genomes'])
            rows += (f'<tr><td><strong>{esc(m["gene"])}</strong></td>'
                     f'<td>{esc(m["class"])}</td>'
                     f'<td><strong>{len(m["genomes"])}</strong> '
                     f'({len(m["genomes"]) / total_samples * 100:.1f}%)</td>'
                     f'<td><div class="genome-list">{tags}</div></td></tr>')

        return f'''
        {credit}
        {self._alert('info', 'fa-dna',
            '<h3>🧬 Point Mutations – Gene-Centric View</h3>'
            '<p>Mutations in gyrA, parC, rpoB, 23S rRNA, and efflux regulators confer resistance independently of acquired genes.</p>')}
        {filter_buttons}
        {info_block}
        {self._grouping_dropdown('mutation-table')}
        <input type="text" class="search-box" id="search-mutation-table"
               onkeyup="searchTable('mutation-table','search-mutation-table')"
               placeholder="🔍 Search mutation by gene or mutation name...">
        <input type="text" class="search-box" id="highlight-mutation-table"
               onkeyup="highlightGenome('mutation-table','highlight-mutation-table')"
               placeholder="🔍 Highlight genomes containing specific text">
        <div class="action-buttons">
            <button class="action-btn btn-primary"
                onclick="exportTableToCSV('mutation-table','mutations.csv')">
                <i class="fas fa-download"></i> Export Mutations</button>
        </div>
        <style>
          #mutation-table {{ table-layout: fixed; width: 100%; }}
          #mutation-table th:nth-child(1),
          #mutation-table td:nth-child(1) {{ width: 22%; word-break: break-word; overflow-wrap: anywhere; }}
          #mutation-table th:nth-child(2),
          #mutation-table td:nth-child(2) {{ width: 18%; word-break: break-word; overflow-wrap: anywhere; }}
          #mutation-table th:nth-child(3),
          #mutation-table td:nth-child(3) {{ width: 12%; text-align: center; white-space: nowrap; }}
          #mutation-table th:nth-child(4),
          #mutation-table td:nth-child(4) {{ width: 48%; }}
          #mutation-table td:nth-child(4) .genome-list {{ max-height: 240px; }}
        </style>
        <div class="master-scrollable-container">
            <table id="mutation-table" class="data-table">
                <thead><tr>
                    <th data-sort="string">Gene</th>
                    <th data-sort="string">Class</th>
                    <th data-sort="number">Count</th>
                    <th data-sort="string">Genomes (scrollable, groupable)</th>
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    def _sec_mge(self, d):
        mge = d.get('mge_data', {})
        stats = mge.get('stats', {})
        per_sample = mge.get('per_sample', {})
        analysis = mge.get('analysis', {})
        note = mge.get('note', '')

        credit = self._credit_bar('#16A085', '🧬', 'Mobile Genetic Elements – mobileOG-db',
            '<strong>mobileOG-db</strong> (Brown et al., <em>Appl Environ Microbiol</em> 2022) — '
            'a manually curated database of protein families mediating the life cycle of bacterial MGEs.')

        why_mge = '''
        <div class="alert-box" style="border-left-color:#16A085; background:#e0f2f1; border-radius:8px; padding:20px; margin:20px 0;">
            <i class="fas fa-dna fa-2x" style="color:#16A085;"></i>
            <div>
                <h3 style="color:#00695c; margin-bottom:10px;">🧬 Why Mobile Genetic Elements Matter</h3>
                <ul style="margin:10px 0 10px 20px;">
                    <li><strong>ICEs</strong> — PAPI-1 and PAPI-2 carry pathogenicity islands.</li>
                    <li><strong>Plasmids</strong> — carry blaVIM, blaIMP, blaNDM on mobile cassettes.</li>
                    <li><strong>Integrons</strong> — cassette-capture systems that accumulate resistance genes.</li>
                    <li><strong>IS elements</strong> — move resistance genes into new genomic contexts.</li>
                    <li><strong>Prophages</strong> — sometimes carry virulence factors.</li>
                </ul>
            </div>
        </div>'''

        category_cards = self._mge_category_cards()

        if not per_sample:
            return f'''
            {credit}
            {why_mge}
            {category_cards}
            {self._alert('warning', 'fa-mobile-alt',
                '<h3>No MGE Data Available</h3>'
                '<p>The pseudoscope_mge_summary.html file was not found or could not be parsed.</p>')}'''

        cards = ''
        cards += self._stat_card(stats.get('Total Samples', len(per_sample)),
                                 'Total Samples', '#8B4513', 'fa-vial')
        cards += self._stat_card(stats.get('With mobileOG hits', '—'),
                                 'With mobileOG hits', '#16A085', 'fa-check')
        cards += self._stat_card(stats.get('Without mobileOG hits', '—'),
                                 'Without mobileOG hits', '#6c757d', 'fa-times')
        cards += self._stat_card(stats.get('Total mobileOG hits', analysis.get('grand_total', '—')),
                                 'Total mobileOG hits', '#E65100', 'fa-database')
        cards += self._stat_card(stats.get('Key MGE-associated genes', '—'),
                                 'Key MGE-associated genes', '#6A1B9A', 'fa-star')
        cards += self._stat_card(stats.get('Runtime', '—'),
                                 'Runtime', '#C62828', 'fa-clock')

        cat_rows = ''
        for c in analysis.get('category_totals', []):
            cat_rows += (f'<tr><td><strong>{esc(c["category"])}</strong></td>'
                         f'<td>{c["total"]:,}</td>'
                         f'<td>{c["pct_of_all"]}%</td>'
                         f'<td>{c["mean"]}</td>'
                         f'<td>{c["max"]}</td>'
                         f'<td>{esc(c["max_sample"])}</td></tr>')

        all_cols = ['Sample']
        seen = set()
        for s in per_sample.values():
            for k in s.keys():
                if k not in seen:
                    seen.add(k)
                    all_cols.append(k)
        sample_header = ''.join(
            f'<th data-sort="{"string" if c == "Sample" else "number"}">{esc(c)}</th>'
            for c in all_cols)
        sample_rows = ''
        for sample, vals in sorted(per_sample.items()):
            cells = f'<td><strong>{esc(sample)}</strong></td>'
            for c in all_cols[1:]:
                cells += f'<td>{vals.get(c, 0)}</td>'
            sample_rows += f'<tr>{cells}</tr>'

        agg_html = self._build_mge_aggregate(per_sample)

        top_html = ''
        for col, label in [('Key MGE-assoc.', 'Key MGE-associated'),
                           ('IS-assoc.', 'IS-associated'),
                           ('Plasmid-assoc.', 'Plasmid-associated')]:
            top = sorted(per_sample.items(),
                         key=lambda x: x[1].get(col, 0), reverse=True)[:5]
            rows_ = ''.join(f'<tr><td><strong>{esc(s)}</strong></td>'
                            f'<td>{v.get(col, 0)}</td></tr>' for s, v in top)
            top_html += f'''<h4>Top 5 — {label}</h4>
            <table class="data-table">
                <thead><tr><th>Sample</th><th>{label} hits</th></tr></thead>
                <tbody>{rows_}</tbody></table>'''

        note_html = ''
        if note:
            note_html = (f'<div class="alert-box alert-warning" style="font-size:.9em;">'
                         f'<strong>Interpretation note:</strong> {note}</div>')

        return f'''
        {credit}
        {self._alert('info', 'fa-mobile-alt',
            '<h3>🧬 Mobile Genetic Elements – mobileOG Profile</h3>'
            '<p>Per-sample counts of mobileOG-db protein family hits across functional categories. '
            'Individual matches represent MGE-associated protein signatures.</p>')}
        <div class="stats-grid">{cards}</div>
        {note_html}
        {why_mge}
        <h3 style="margin-top:30px;">📚 Understanding the Categories</h3>
        {category_cards}
        <h3 style="margin-top:30px;">📊 Category Totals (aggregated across all samples)</h3>
        <div class="scrollable-table">
            <table class="data-table">
                <thead><tr>
                    <th data-sort="string">Category</th>
                    <th data-sort="number">Total hits</th>
                    <th data-sort="number">% of all hits</th>
                    <th data-sort="number">Mean/sample</th>
                    <th data-sort="number">Max/sample</th>
                    <th data-sort="string">Sample with max</th>
                </tr></thead>
                <tbody>{cat_rows}</tbody>
            </table>
        </div>
        <h3 style="margin-top:30px;">🧬 Per-sample mobileOG Profile</h3>
        <input type="text" class="search-box" id="search-mge-table"
               onkeyup="searchTable('mge-table','search-mge-table')"
               placeholder="🔍 Search sample...">
        <div class="action-buttons">
            <button class="action-btn btn-primary"
                onclick="exportTableToCSV('mge-table', 'mge_profile.csv')">
                <i class="fas fa-download"></i> Export MGE CSV</button>
        </div>
        <div class="master-scrollable-container">
            <table id="mge-table" class="data-table">
                <thead><tr>{sample_header}</tr></thead>
                <tbody>{sample_rows}</tbody>
            </table>
        </div>
        <h3 style="margin-top:30px;">📈 Aggregate MGE Profile by Typing</h3>
        {agg_html}
        <h3 style="margin-top:30px;">🏆 Top-N Views</h3>
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:20px;">
            {top_html}
        </div>'''

    def _mge_category_cards(self) -> str:
        cards = [
            ('Integr.', 'Integrase', '#ef4444', 'fa-cut',
             'Enzymes catalysing integration and excision of mobile elements.',
             'High integrase counts suggest active chromosomal integration of foreign DNA.'),
            ('Transfer', 'Transfer / Conjugation', '#f97316', 'fa-share-alt',
             'Type IV secretion systems (T4SS), relaxases, and coupling proteins.',
             'The molecular machinery of horizontal gene transfer.'),
            ('Stab.', 'Stability / Maintenance', '#eab308', 'fa-anchor',
             'Partitioning systems (parAB), toxin-antitoxin cassettes.',
             'A plasmid with strong stability systems persists even without selection.'),
            ('Phage', 'Bacteriophage', '#22c55e', 'fa-virus',
             'Complete phage-related proteins including capsids, tails, lysogeny modules.',
             'Phages can carry virulence genes and mobilise chromosomal DNA.'),
            ('Repl.', 'Replication', '#14b8a6', 'fa-sync',
             'Replication initiators, primases, helicases.',
             'Only elements with a functional replication module persist extrachromosomally.'),
            ('IS-assoc.', 'Insertion Sequences', '#0891b2', 'fa-exchange-alt',
             'Transposases and IS elements — the smallest autonomous MGEs.',
             'IS elements frequently disrupt genes or carry resistance genes into new contexts.'),
            ('ICE-assoc.', 'Integrative Conjugative Elements', '#3b82f6', 'fa-project-diagram',
             'Self-transmissible elements that integrate and conjugate.',
             'ICEs move large multi-gene cassettes — key vehicles for pathogenicity islands.'),
            ('Plasmid-assoc.', 'Plasmid', '#8b5cf6', 'fa-plug',
             'Plasmid backbone, replication, partitioning, and conjugal transfer proteins.',
             'Plasmids are the single most important vehicle for acquired resistance in P. aeruginosa.'),
            ('Phage-assoc.', 'Phage-associated (lysogeny)', '#d946ef', 'fa-dna',
             'Prophage and lysogeny-associated proteins.',
             'Prophage carriage often correlates with virulence gene acquisition.'),
            ('Key MGE-assoc.', 'Key MGE signatures', '#C62828', 'fa-star',
             'Curated subset of high-confidence mobileOG-db entries.',
             'The cleanest single indicator of MGE activity in a genome.'),
        ]
        html = ('<div style="display:grid;'
                'grid-template-columns:repeat(auto-fit,minmax(340px,1fr));'
                'gap:15px;margin:20px 0;">')
        for key, label, color, icon, definition, importance in cards:
            html += f'''
            <div style="background:white;padding:16px;border-radius:10px;
                        border-left:4px solid {color};
                        box-shadow:0 2px 8px rgba(0,0,0,.06);">
                <div style="display:flex;align-items:center;gap:10px;margin-bottom:8px;">
                    <i class="fas {icon}" style="color:{color};font-size:1.4em;"></i>
                    <div>
                        <strong style="color:{color};font-size:1.05em;">{label}</strong>
                        <div style="font-size:0.8em;color:#888;font-family:monospace;">{key}</div>
                    </div>
                </div>
                <p style="font-size:0.88em;color:#333;margin:0 0 8px 0;">{definition}</p>
                <p style="font-size:0.85em;color:#555;margin:0;padding-top:8px;
                          border-top:1px dashed #e0e0e0;">
                    <i class="fas fa-bullseye" style="color:{color};"></i>
                    <em>{importance}</em>
                </p>
            </div>'''
        html += '</div>'
        return html

    def _build_mge_aggregate(self, per_sample):
        samples_data = getattr(self, '_current_samples_data', {})
        if not samples_data or not per_sample:
            return ''
        fields = [
            ('MLST', 'MLST (ST)'),
            ('O_Type', 'O-type'),
            ('T3SS', 'T3SS Pathotype'),
        ]
        sample_cols = list(next(iter(per_sample.values())).keys())
        blocks = ''
        for field_key, field_label in fields:
            groups = defaultdict(list)
            for sample, vals in per_sample.items():
                sd = samples_data.get(sample, {})
                typ = sd.get('typing', {}).get(field_key, 'ND')
                if typ in ('Not Assigned', 'ND', ''):
                    continue
                groups[typ].append(vals)
            if not groups:
                continue
            rows = ''
            for grp, records in sorted(groups.items(), key=lambda x: -len(x[1])):
                means = {}
                for c in sample_cols:
                    vals_ = [r.get(c, 0) for r in records
                             if isinstance(r.get(c, 0), (int, float))]
                    means[c] = round(sum(vals_) / len(vals_), 1) if vals_ else 0
                cells = ''.join(f'<td>{means[c]}</td>' for c in sample_cols)
                rows += (f'<tr><td><strong>{esc(grp)}</strong></td>'
                         f'<td>{len(records)}</td>{cells}</tr>')
            header = ''.join(f'<th>{esc(c)}</th>' for c in sample_cols)
            blocks += f'''
            <h4 style="margin-top:20px;">Mean MGE profile by {field_label}</h4>
            <div class="master-scrollable-container">
                <table class="data-table">
                    <thead><tr>
                        <th data-sort="string">{field_label}</th>
                        <th data-sort="number">N samples</th>
                        {header}
                    </tr></thead>
                    <tbody>{rows}</tbody>
                </table>
            </div>'''
        return blocks

    def _sec_patterns(self, d):
        P = d['patterns']
        html = self._alert('info', 'fa-project-diagram',
            '<h3>🔍 Cross-Genome Pattern Discovery</h3>'
            '<p>Multi-way typing combinations, high-risk combinations, and gene co-occurrence.</p>')
        html += self._combo_table('triple-table', 'Triple Typing (ST – O-type – T3SS)',
                                  P.get('mlst_otype_t3ss', {}))
        high_risk = P.get('high_risk_combinations', [])
        if high_risk:
            rows = ''
            for c in high_risk:
                rows += (f'<tr><td><strong>{esc(c["sample"])}</strong></td>'
                         f'<td>{esc(c["mlst"])}</td><td>{esc(c["o_type"])}</td>'
                         f'<td>{esc(c["t3ss"])}</td>'
                         f'<td>{", ".join(esc(x) for x in c["critical_amr_genes"])}</td>'
                         f'<td>{", ".join(esc(x) for x in c["critical_virulence_genes"])}</td></tr>')
            html += f'''<h3>⚠️ High-Risk Combinations</h3>
            <div class="master-scrollable-container"><table class="data-table">
                <thead><tr><th>Sample</th><th>MLST</th><th>O-type</th><th>T3SS</th>
                <th>Critical AMR</th><th>Critical Virulence</th></tr></thead>
                <tbody>{rows}</tbody></table></div>'''
        cooc = P.get('gene_cooccurrence', {})
        if cooc:
            pairs = []
            for g1, partners in cooc.items():
                for g2, cnt in partners.items():
                    pairs.append((g1, g2, cnt))
            pairs.sort(key=lambda x: x[2], reverse=True)
            rows = ''.join(f'<tr><td>{esc(g1)}</td><td>{esc(g2)}</td><td>{cnt}</td></tr>'
                           for g1, g2, cnt in pairs[:500])
            html += f'''<h3>📈 Gene Co-occurrence (Top 500)</h3>
            <div class="master-scrollable-container"><table class="data-table">
                <thead><tr><th>Gene 1</th><th>Gene 2</th><th>Co-occurrence</th></tr></thead>
                <tbody>{rows}</tbody></table></div>'''
        return html

    def _sec_compare(self, d):
        samples_data = d.get('samples', {})
        mutation_data = d.get('mutation_data', {})
        mutation_list = mutation_data.get('mutations', [])

        if len(samples_data) < 2:
            return self._alert('warning', 'fa-balance-scale',
                '<h3>Not Enough Samples</h3>'
                f'<p>Comparison requires at least two samples. This dataset has '
                f'{len(samples_data)}.</p>')

        sample_mutations: Dict[str, List[str]] = defaultdict(list)
        for m in mutation_list:
            label = f"{m.get('gene', '')} {m.get('mutation', '')}".strip()
            for genome in m.get('genomes', []):
                sample_mutations[genome].append(label)

        compare_data: Dict[str, Dict[str, Any]] = {}

        def db_genes(sample: str, db: str) -> List[str]:
            return list(samples_data.get(sample, {})
                        .get('abricate_databases', {}).get(db, []))

        for sample, sd in samples_data.items():
            t = sd.get('typing', {})
            amr_genes = sorted(set(
                sd.get('amrfinder', {}).get('all_genes', [])
                + db_genes(sample, 'card')
                + db_genes(sample, 'resfinder')
                + db_genes(sample, 'argannot')
                + db_genes(sample, 'megares')
                + db_genes(sample, 'ncbi')
            ))
            virulence_genes = sorted(set(
                db_genes(sample, 'vfdb')
                + db_genes(sample, 'ecoli_vf')
                + db_genes(sample, 'pa_vf')
                + db_genes(sample, 'pseudomonas_vf')
            ))
            bacmet_genes = sorted(set(db_genes(sample, 'bacmet2')))
            plasmid_genes = sorted(set(db_genes(sample, 'plasmidfinder')))
            mutations = sorted(set(sample_mutations.get(sample, [])))

            compare_data[sample] = {
                'typing': {
                    'MLST':   t.get('MLST', 'ND'),
                    'O-type': t.get('O_Type', 'ND'),
                    'T3SS Pathotype': t.get('T3SS', 'ND'),
                },
                'genes': {
                    'AMR':       amr_genes,
                    'Virulence': virulence_genes,
                    'BACMET':    bacmet_genes,
                    'Plasmids':  plasmid_genes,
                    'Mutations': mutations,
                },
            }

        sorted_samples = sorted(compare_data.keys())
        options_a = ''.join(
            f'<option value="{esc(s)}"{" selected" if i == 0 else ""}>{esc(s)}</option>'
            for i, s in enumerate(sorted_samples)
        )
        options_b = ''.join(
            f'<option value="{esc(s)}"{" selected" if i == 1 else ""}>{esc(s)}</option>'
            for i, s in enumerate(sorted_samples)
        )
        blob_json = json.dumps(compare_data, default=str, ensure_ascii=False)
        precomp = self._precompute_similarity(compare_data)
        precomp_json = json.dumps(precomp, default=str, ensure_ascii=False)

        return f'''
        <div class="alert-box alert-info">
            <i class="fas fa-balance-scale fa-2x"></i>
            <div>
                <h3>⚖️ Compare Isolates</h3>
                <p>Two modes: <strong>Pair Compare</strong> for a deep side-by-side of two isolates, and
                <strong>Cluster Detection</strong> to find groups of highly similar samples.</p>
            </div>
        </div>
        <div class="mode-toggle">
            <button class="mode-btn active" data-mode="pair" onclick="setCompareMode('pair')">
                <i class="fas fa-balance-scale"></i> Pair Compare
            </button>
            <button class="mode-btn" data-mode="cluster" onclick="setCompareMode('cluster')">
                <i class="fas fa-project-diagram"></i> Cluster Detection
            </button>
        </div>
        <div id="mode-pair" class="compare-mode active">
            <div class="compare-panel">
                <div class="compare-select-row">
                    <label for="compare-a">
                        <span>Sample A</span>
                        <select id="compare-a">{options_a}</select>
                    </label>
                    <button class="action-btn btn-primary" onclick="runCompare()">
                        <i class="fas fa-balance-scale"></i> Compare</button>
                    <button class="action-btn btn-light" onclick="swapCompare()">
                        <i class="fas fa-exchange-alt"></i> Swap</button>
                    <button class="action-btn btn-light" onclick="resetCompare()">
                        <i class="fas fa-sync"></i> Reset</button>
                    <label for="compare-b">
                        <span>Sample B</span>
                        <select id="compare-b">{options_b}</select>
                    </label>
                </div>
                <div id="compare-verdict" style="display:none;"></div>
                <div id="compare-visuals"></div>
                <div id="compare-output"></div>
            </div>
        </div>
        <div id="mode-cluster" class="compare-mode">
            <div class="compare-panel">
                <div class="cluster-controls">
                    <label>
                        <span>Similarity threshold:</span>
                        <select id="cluster-threshold" onchange="runCluster()">
                            <option value="95">Tight clusters (≥95%)</option>
                            <option value="90" selected>Moderate clusters (≥90%)</option>
                            <option value="85">Loose clusters (≥85%)</option>
                            <option value="80">Broad clusters (≥80%)</option>
                        </select>
                    </label>
                    <button class="action-btn btn-primary" onclick="runCluster()">
                        <i class="fas fa-project-diagram"></i> Detect Clusters</button>
                    <span id="cluster-summary" style="font-size:.9em;color:#5C2E0B;font-weight:600;"></span>
                </div>
                <div id="cluster-clusters"></div>
                <div id="cluster-heatmap"></div>
            </div>
        </div>
        <script>
        window.PSEUDOSCOPE_COMPARE_DATA = {blob_json};
        window.PSEUDOSCOPE_SIM = {precomp_json};
        </script>'''

    def _sec_citation(self, d):
        # Category-based colors — one for the main tool, one for databases, one for tools
        DATABASES = {
            'PubMLST', 'CARD', 'ResFinder', 'VFDB', 'PlasmidFinder',
            'BacMet2', 'MEGARes', 'mobileOG-db',
        }

        def fmt(name, text, link_type, link_url):
            if name == 'PseudoScope':
                color, border = '#FFF8E1', '#B8860B'   # amber — main tool
            elif name in DATABASES:
                color, border = '#E0F2F1', '#00695C'   # teal — databases
            else:
                color, border = '#E3F2FD', '#1565C0'   # blue — tools
            if link_type == 'doi':
                link_html = f'<a href="https://doi.org/{link_url}" target="_blank" style="color:#0D6EFD;font-weight:600;">🔗 {link_url}</a>'
            else:
                link_html = f'<a href="{link_url}" target="_blank" style="color:{border};font-weight:600;">🔗 {link_url}</a>'
            return f'''
            <div class="citation-card" style="background:{color};border-left:5px solid {border};padding:15px;border-radius:8px;margin-bottom:15px;">
                <strong>{name}</strong> – {text} {link_html}
                <button class="copy-btn" data-citation="{text} {link_url}">📋 Copy</button>
            </div>'''

        html = self._alert('info', 'fa-quote-right',
            '<h3>📚 How to Cite PseudoScope and Its Dependencies</h3>'
            '<p>If you use this tool in your research, please cite the main tool and the relevant third‑party tools and databases.</p>')

        html += fmt("PseudoScope",
            "Beckley B, Amarh V. PseudoScope: a species‑specific bioinformatics suite for rapid and accessible Pseudomonas aeruginosa genomic analysis. GitHub 2026.",
            "url", "https://github.com/bbeckley-hub/pseudoscope")
        html += fmt("MLST (Pasteur scheme)",
            "Seemann T. MLST: Scan contig files against PubMLST typing schemes. GitHub.",
            "url", "https://github.com/tseemann/mlst")
        html += fmt("PubMLST",
            "Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. Wellcome Open Res. 2018;3:124.",
            "doi", "10.12688/wellcomeopenres.14826.1")
        html += fmt("PAST (Original)",
            "Thomsen MCF, et al. A web-based tool for serotyping of Pseudomonas aeruginosa. CGE, Technical University of Denmark.",
            "url", "https://cge.food.dtu.dk/services/PAst/")
        html += fmt("pasty",
            "Petit RA III. pasty: A tool for in silico serogrouping of Pseudomonas aeruginosa isolates. GitHub.",
            "url", "https://github.com/rpetit3/pasty")
        html += fmt("camlhmp",
            "Petit RA III. camlhmp: Classification through yAML Heuristic Mapping Protocol. GitHub.",
            "url", "https://github.com/rpetit3/camlhmp")
        html += fmt("AMRFinderPlus",
            "Feldgarden M, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021;11(1):12728.",
            "doi", "10.1038/s41598-021-91456-0")
        html += fmt("ABRicate",
            "Seemann T. ABRicate: mass screening of contigs for antibiotic resistance genes. GitHub.",
            "url", "https://github.com/tseemann/abricate")
        html += fmt("CARD",
            "McArthur AG, et al. The comprehensive antibiotic resistance database. Antimicrob Agents Chemother. 2013;57(7):3348-57.",
            "doi", "10.1128/AAC.00419-13")
        html += fmt("ResFinder",
            "Florensa AF, et al. ResFinder – an open online resource for identification of antimicrobial resistance genes. Microb Genom. 2022;8(1):000748.",
            "doi", "10.1099/mgen.0.000748")
        html += fmt("VFDB",
            "Chen L, et al. VFDB 2012 update: toward the genetic diversity and molecular evolution of bacterial virulence factors. Nucleic Acids Res. 2012;40(Database issue):D641-5.",
            "doi", "10.1093/nar/gkr989")
        html += fmt("PlasmidFinder",
            "Carattoli A, et al. In silico detection and typing of plasmids using PlasmidFinder. Antimicrob Agents Chemother. 2014;58(7):3895-903.",
            "doi", "10.1128/AAC.02412-14")
        html += fmt("BacMet2",
            "Pal C, et al. BacMet: antibacterial biocide and metal resistance genes database. Nucleic Acids Res. 2014;42(Database issue):D737-43.",
            "doi", "10.1093/nar/gkt1252")
        html += fmt("MEGARes",
            "Doster E, et al. MEGARes 2.0: a database for classification of antimicrobial drug, biocide and metal resistance determinants. Nucleic Acids Res. 2020;48(D1):D561-D569.",
            "doi", "10.1093/nar/gkz1010")
        html += fmt("mobileOG-db",
            "Brown CL, et al. mobileOG-db: a Manually Curated Database of Protein Families Mediating the Life Cycle of Bacterial Mobile Genetic Elements. Appl Environ Microbiol. 2022;88(18):e00991-22.",
            "doi", "10.1128/aem.00991-22")
        html += fmt("Prodigal",
            "Hyatt D, et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics. 2010;11:119.",
            "doi", "10.1186/1471-2105-11-119")
        html += fmt("DIAMOND",
            "Buchfink B, Xie C, Huson DH. Fast and sensitive protein alignment using DIAMOND. Nat Methods. 2015;12:59-60.",
            "doi", "10.1038/nmeth.3176")
        html += fmt("fastANI",
            "Jain C, et al. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. Nat Commun. 2018;9(1):5114.",
            "doi", "10.1038/s41467-018-07641-9")
        html += fmt("Biopython",
            "Cock PJ, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009;25(11):1422-3.",
            "doi", "10.1093/bioinformatics/btp163")

        html += '''
        <div class="alert-box alert-success" style="margin-top:20px;">
            <i class="fas fa-hand-peace"></i>
            <div><strong>Suggested acknowledgement:</strong><br>
            "Genomic analysis was performed using PseudoScope [Beckley &amp; Amarh, 2026], which integrates MLST [Seemann, 2018] using the PubMLST database [Jolley et al., 2018], ABRicate [Seemann, 2018], AMRFinderPlus [Feldgarden et al., 2021], PAST serotyping [Thomsen et al., DTU] via the pasty wrapper [Petit et al., 2023], and fastANI [Jain et al., 2018] for comprehensive P. aeruginosa characterization. Antimicrobial resistance genes were identified using the CARD [McArthur et al., 2013] and ResFinder [Florensa et al., 2022] databases. For biocide and heavy metal resistance genes, BacMet [Pal et al., 2014] was used. Virulence and plasmid screening were performed with ABRicate using the VFDB [Chen et al., 2012] and PlasmidFinder [Carattoli et al., 2014] databases. Mobile genetic element profiling used mobileOG-db [Brown et al., 2022], Prodigal [Hyatt et al., 2010], and DIAMOND [Buchfink et al., 2015]. FASTA QC was performed using Biopython [Cock et al., 2009]."
            </div>
        </div>'''
        return html

    def _sec_guide(self, d):
        return '''
        <div class="alert-box alert-info">
            <i class="fas fa-robot fa-2x"></i>
            <div>
                <h3>🤖 AI Assistant Guide – Unleash the Power of AI for P. aeruginosa Genomic Epidemiology</h3>
                <p>Use large language models (LLMs) like <strong>ChatGPT, Claude, or Gemini</strong> to interact with your <em>P. aeruginosa</em> dataset.</p>
            </div>
        </div>
        <div style="margin:20px 0;">
            <div class="database-section">
                <h4><i class="fas fa-upload"></i> How to Feed This Report to AI</h4>
                <ol>
                    <li><strong>Upload the JSON file</strong> — <code>genius_pseudomonas_ultimate_gene_centric_report.json</code>.</li>
                    <li><strong>Upload the HTML report</strong> — modern AI tools parse HTML tables directly.</li>
                    <li><strong>Copy-paste specific tables</strong> — for quick insights, paste a table into chat.</li>
                </ol>
                <p style="margin-top:10px;background:#e8f5e9;padding:10px;border-radius:5px;">
                    <i class="fas fa-info-circle"></i> <strong>Pro tip:</strong> tell the AI: <em>"You are a bioinformatician analysing P. aeruginosa genomes. The attached data contains typing, AMR, virulence, Bacmet2, mutation, and MGE information."</em>
                </p>
            </div>
            <div class="database-section">
                <h4><i class="fas fa-chart-line"></i> Scientifically Relevant Questions to Ask</h4>
                <div style="display:grid;grid-template-columns:1fr 1fr;gap:10px;">
                    <div style="background:#fdf6ec;padding:10px;border-radius:8px;">
                        <strong>🧬 Epidemiology &amp; Clonality</strong>
                        <ul style="margin-top:5px;font-size:.9em;">
                            <li>What are the most common MLST sequence types?</li>
                            <li>Which STs are associated with carbapenemase production?</li>
                            <li>What is the O-type distribution, and does it correlate with ST?</li>
                        </ul>
                    </div>
                    <div style="background:#fdf6ec;padding:10px;border-radius:8px;">
                        <strong>💉 T3SS Pathotyping</strong>
                        <ul style="margin-top:5px;font-size:.9em;">
                            <li>How many isolates are ExoU+ vs ExoS+?</li>
                            <li>Do ExoU+ isolates cluster in specific STs?</li>
                            <li>Are carbapenemase-positive isolates more likely to be ExoU+?</li>
                        </ul>
                    </div>
                    <div style="background:#fdf6ec;padding:10px;border-radius:8px;">
                        <strong>💊 Carbapenemase &amp; OXA</strong>
                        <ul style="margin-top:5px;font-size:.9em;">
                            <li>How many isolates carry true acquired carbapenemases?</li>
                            <li>Are OXA-2-like or OXA-10-like ESBL variants present?</li>
                            <li>Any <em>mcr</em>-positive isolates?</li>
                        </ul>
                    </div>
                    <div style="background:#fdf6ec;padding:10px;border-radius:8px;">
                        <strong>🦠 Virulence &amp; MGE</strong>
                        <ul style="margin-top:5px;font-size:.9em;">
                            <li>Which isolates carry <em>exoU</em>? Are they MDR?</li>
                            <li>Do high-MGE-hit isolates carry more AMR genes?</li>
                            <li>Is plasmid replicon content correlated with ST?</li>
                        </ul>
                    </div>
                </div>
            </div>
            <div class="database-section">
                <h4><i class="fas fa-balance-scale"></i> Scientific Rigour &amp; Ethical AI Use</h4>
                <ul>
                    <li><strong>AI is your co-pilot, not the pilot.</strong></li>
                    <li><strong>Verify, verify, verify.</strong> Cross-check critical calls with primary literature.</li>
                    <li><strong>No patient-identifiable data.</strong> Only upload aggregated, de-identified genomic data.</li>
                    <li><strong>Transparency in publications.</strong> Mention AI-assisted pattern discovery in methods.</li>
                    <li><strong>AI hallucination is real.</strong> Treat every AI statement as a hypothesis.</li>
                </ul>
            </div>
            <div class="database-section" style="background:#fff3cd;border-left:6px solid #ffc107;">
                <h4><i class="fas fa-smile-wink"></i> A (Mostly Serious) AI Survival Guide</h4>
                <ul>
                    <li><strong>If the AI says "I don't know"</strong> — trust it.</li>
                    <li><strong>If the AI says "It is widely known"</strong> — ask for a reference.</li>
                    <li><strong>If the AI offers a completely novel evolutionary theory</strong> — check your coffee.</li>
                    <li><strong>Remember:</strong> AI won't take your job — but a microbiologist who knows how to use AI might! 😉</li>
                </ul>
            </div>
        </div>'''


    def _sec_calltoaction(self, d):
            return '''
            <div class="alert-box alert-info">
            <i class="fas fa-globe fa-2x"></i>
            <div>
                <h3>The Global Burden of AMR and Our Call to Action</h3>
                <p>Antimicrobial resistance (AMR) is one of the top global public health threats, with an estimated <strong>1.27 million direct deaths annually</strong>. <em>Pseudomonas aeruginosa</em> is a WHO critical-priority pathogen. Genomic surveillance is our best tool. Tracking AMR and virulence determinants is essential to inform treatment guidelines and infection control.</p>
                <p>We developed <strong>PSEUDOSCOPE</strong> to empower researchers and clinicians – especially in low-resource settings – to analyse their own sequencing data without extensive bioinformatics expertise.</p>
            </div>
            </div>

            <div style="background:#e8f5e9;padding:20px;border-radius:12px;margin:20px 0;">
            <h3><i class="fas fa-bacterium"></i> ESCAPE AMR – Our Ongoing Project (ESKAPE Pathogens)</h3>
            <p><strong>PSEUDOSCOPE</strong> is one of the first modules of a larger initiative called <strong>ESCAPE AMR</strong> (formerly ESKAPE). We target the notorious <strong>ESKAPE pathogens</strong>:</p>
            <ul>
                <li><strong>E</strong>nterococcus faecium</li>
                <li><strong>S</strong>taphylococcus aureus</li>
                <li><strong>K</strong>lebsiella pneumoniae</li>
                <li><strong>A</strong>cinetobacter baumannii</li>
                <li><strong>P</strong>seudomonas aeruginosa</li>
                <li><strong>E</strong>nterobacter species (including E. coli)</li>
            </ul>
            <p>These bacteria "escape" the effects of antibiotics – hence the name. But we believe the name is also a global call to action:</p>
            <div style="background:#fff3e0;padding:15px;border-radius:8px;margin:15px 0;">
                <p><strong>🔹 E</strong>veryone must join forces – researchers, clinicians, policymakers, and citizens.<br>
                <strong>🔹 S</strong>mart surveillance is our first line of defence. No more guessing – we need genomic data.<br>
                <strong>🔹 K</strong>nowledge must be shared openly. No paywalls, no closed silos.<br>
                <strong>🔹 A</strong>frica bears a heavy AMR burden, but African solutions are already emerging.<br>
                <strong>🔹 P</strong>revention is cheaper than cure. Let's stop resistant infections before they spread.<br>
                <strong>🔹 E</strong>very day we delay, more lives are at stake. The time to act is now, not tomorrow.</p>
            </div>
            <p><i class="fas fa-laugh-squint"></i> <strong>"We didn't choose the name ESKAPE because it sounds cool (though it does). We chose it because it reminds us every single day: we must ESCAPE the AMR crisis – together, urgently, and with the best science we have."</strong><br>
            — Brown Beckley, lead developer (who secretly hopes this pun makes you smile, not roll your eyes 😉)</p>
            </div>

            <div style="text-align:center;margin:40px 0;">
            <i class="fas fa-star" style="font-size:3em;color:#ffc107;"></i>
            <h3>🤝 We Invite You to Contribute!</h3>
            <p><strong>If you find this tool useful, please:</strong></p>
            <div class="action-buttons" style="justify-content:center;">
                <a href="https://github.com/bbeckley-hub/pseudoscope" target="_blank" class="action-btn btn-primary" style="text-decoration:none;"><i class="fab fa-github"></i> ⭐ Star us on GitHub</a>
                <a href="mailto:brownbeckley94@gmail.com" class="action-btn btn-success" style="text-decoration:none;"><i class="fas fa-envelope"></i> Contact the team</a>
                <a href="https://github.com/bbeckley-hub/pseudoscope/issues" target="_blank" class="action-btn btn-warning" style="text-decoration:none;"><i class="fas fa-bug"></i> Report issues</a>
            </div>
            <p style="margin-top:20px;"><i class="fas fa-chalkboard-user"></i> <strong>We welcome collaborations</strong> to adapt this tool for other pathogens and to improve AMR surveillance in Africa and beyond.</p>
            <p><i class="fas fa-hand-holding-heart"></i> If you are a funder or organisation interested in supporting the <strong>ESCAPE AMR</strong> project, please reach out. Together we can build a free, open-source ecosystem for genomic surveillance of all ESKAPE pathogens.</p>
            </div>

            <div style="background:#f8f9fa;padding:15px;border-radius:8px;">
            <i class="fas fa-quote-left"></i> "AMR is a silent pandemic, but we have the tools to fight it – if we share them, if we teach each other, and if we act with urgency. Let's escape the era of untreatable infections."<br>
            <strong>— The ESCAPE AMR Team, University of Ghana Medical School</strong>
            </div>'''

    def _sec_export(self, d):
        exports = [
            ('samples-table', 'Sample Overview CSV', 'fa-table', 'sample_overview.csv'),
            ('amr-table', 'AMR Genes CSV', 'fa-biohazard', 'amr_genes.csv'),
            ('vir-table', 'Virulence Genes CSV', 'fa-virus', 'virulence_genes.csv'),
            ('bacmet-table', 'Bacmet2 Genes CSV', 'fa-flask', 'bacmet_genes.csv'),
            ('plasmid-table', 'Plasmid Replicons CSV', 'fa-plug', 'plasmid_replicons.csv'),
            ('mutation-table', 'Mutations CSV', 'fa-dna', 'mutations.csv'),
            ('qc-table', 'FASTA QC CSV', 'fa-chart-line', 'fasta_qc.csv'),
            ('mge-table', 'MGE Profile CSV', 'fa-mobile-alt', 'mge_profile.csv'),
        ]
        cards = ''
        for tid, label, icon, fname in exports:
            cards += f'''<div class="dashboard-card" style="border-left-color:#37474F;"
                onclick="exportTableToCSV('{tid}', '{fname}')">
                <i class="fas {icon} fa-2x" style="color:#37474F;"></i>
                <div class="card-label">{label}</div></div>'''
        cards += '''<div class="dashboard-card" style="border-left-color:#37474F;"
            onclick="location.href='genius_pseudomonas_ultimate_gene_centric_report.json'">
            <i class="fas fa-file-code fa-2x" style="color:#37474F;"></i>
            <div class="card-label">Complete JSON Data</div></div>'''
        return f'''
        {self._alert('info', 'fa-download',
            '<h3>📥 Export Data</h3>'
            '<p>Download any table as CSV, or the full dataset as JSON.</p>')}
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:20px;margin-top:20px;">
            {cards}
        </div>'''


# =============================================================================
# ORCHESTRATOR
# =============================================================================
class PseudomonasUltimateReporter:
    """Coordinates file discovery, parsing, analysis, and report generation."""

    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "PSEUDOSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = PseudomonasHTMLParser()
        self.analyzer = PseudomonasDataAnalyzer()
        self.generator = PseudomonasHTMLGenerator(self.analyzer)
        self.metadata = {
            "tool_name": "PseudoScope Ultimate P. aeruginosa Reporter",
            "version": "3.0.0",
            "author": "Brown Beckley <brownbeckley94@gmail.com>",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir),
        }

    def find_files(self) -> Dict[str, Any]:
        print("🔍 Searching for PseudoScope reports...")
        found = {
            'mlst_html': None, 'qc_html': None, 'amrfinder_html': None,
            'abricate_htmls': [], 'past_html': None,
            'mutation_html': None, 'mge_html': None,
        }
        for f in self.input_dir.glob('**/*.html'):
            name = f.name.lower()
            if 'mlst_summary' in name and found['mlst_html'] is None:
                found['mlst_html'] = f
                print(f"    ✓ MLST: {f.name}")
            elif 'fasta_qc_summary' in name and found['qc_html'] is None:
                found['qc_html'] = f
                print(f"    ✓ FASTA QC: {f.name}")
            elif 'amrfinder_summary' in name and found['amrfinder_html'] is None:
                found['amrfinder_html'] = f
                print(f"    ✓ AMRfinder: {f.name}")
            elif 'past_summary' in name and found['past_html'] is None:
                found['past_html'] = f
                print(f"    ✓ PAST: {f.name}")
            elif 'mutation_summary' in name and found['mutation_html'] is None:
                found['mutation_html'] = f
                print(f"    ✓ Mutation: {f.name}")
            elif ('mge_summary' in name or 'mge_ultimate' in name) and found['mge_html'] is None:
                found['mge_html'] = f
                print(f"    ✓ MGE: {f.name}")
            elif any(f'{db}_summary_report' in name or f'{db}_summary' in name
                     for db in self.parser.abricate_databases):
                if f not in found['abricate_htmls']:
                    found['abricate_htmls'].append(f)
                    print(f"    ✓ ABRicate: {f.name}")
        return found

    def integrate(self, files: Dict) -> Dict[str, Any]:
        print("\n🔗 Integrating data...")
        data = {
            'metadata': self.metadata,
            'samples': {},
            'patterns': {},
            'gene_centric': {},
            'qc_data': {},
            'mutation_data': {},
            'mge_data': {},
        }

        mlst_data = self.parser.parse_mlst_report(files['mlst_html']) if files['mlst_html'] else {}
        past_data = self.parser.parse_past_report(files['past_html']) if files['past_html'] else {}
        if files['qc_html']:
            data['qc_data'] = self.parser.parse_qc_report(files['qc_html'])
        if files['mutation_html']:
            data['mutation_data'] = self.parser.parse_mutation_summary_html(files['mutation_html'])

        amr_by_sample, amr_freq = {}, {}
        if files['amrfinder_html']:
            amr_by_sample, amr_freq = self.parser.parse_amrfinder_report(files['amrfinder_html'])

        abricate_by_sample = defaultdict(dict)
        abricate_freq = {}
        for f in files['abricate_htmls']:
            db_name = 'unknown'
            for db in self.parser.abricate_databases:
                if db in f.name.lower():
                    db_name = db
                    break
            if db_name == 'unknown':
                continue
            genes_by_sample, gene_freq = self.parser.parse_abricate_report(f)
            for s, gs in genes_by_sample.items():
                abricate_by_sample[s][db_name] = gs
            abricate_freq[db_name] = gene_freq

        if files['mge_html']:
            mge_raw = self.parser.parse_mge_summary_html(files['mge_html'])
            mge_raw['analysis'] = self.analyzer.compute_mge_stats(mge_raw)
            data['mge_data'] = mge_raw

        all_samples = set(mlst_data) | set(past_data) | set(amr_by_sample) \
                      | set(abricate_by_sample) | set(data['qc_data']) \
                      | set(data['mge_data'].get('per_sample', {}))
        all_samples = sorted(all_samples)
        if not all_samples:
            print("❌ No samples found in any report")
            return {}

        print(f"📊 Found {len(all_samples)} unique samples")

        for sample in all_samples:
            amr_genes = amr_by_sample.get(sample, [])
            vir_genes = []
            for vir_db in ('vfdb', 'ecoli_vf', 'pa_vf', 'pseudomonas_vf'):
                vir_genes.extend(abricate_by_sample.get(sample, {}).get(vir_db, []))
            all_genes_for_derivation = list(set(amr_genes + vir_genes))

            t3ss = self.analyzer.derive_t3ss_pathotype(all_genes_for_derivation)

            mge_rec = data['mge_data'].get('per_sample', {}).get(sample, {})
            mge_hits = mge_rec.get('mobileOG Hits', 0) if isinstance(mge_rec, dict) else 0
            if not isinstance(mge_hits, (int, float)):
                mge_hits = 0

            data['samples'][sample] = {
                'typing': {
                    'MLST':   mlst_data.get(sample, {}).get('ST', 'ND'),
                    'O_Type': past_data.get(sample, {}).get('O_Type', 'ND'),
                    'T3SS':   t3ss,
                },
                'amrfinder': {'all_genes': amr_genes,
                              'critical_genes': [], 'high_risk_genes': []},
                'abricate_databases': dict(abricate_by_sample.get(sample, {})),
                'mge_hits': mge_hits,
            }

        data['gene_frequencies'] = {'amrfinder': amr_freq, 'abricate': abricate_freq}

        print("\n🧠 Building gene-centric + pattern tables...")
        data['gene_centric'] = self.analyzer.create_gene_centric_tables(data)
        data['patterns'] = self.analyzer.create_cross_genome_patterns(data)
        self.generator._current_samples_data = data['samples']
        return data

    def write_json(self, data: Dict[str, Any]) -> Path:
        print("\n📝 Writing JSON report...")
        out = self.output_dir / "genius_pseudomonas_ultimate_gene_centric_report.json"

        def serial(obj):
            if obj is None or isinstance(obj, (str, int, float, bool)):
                return obj
            if isinstance(obj, (list, tuple, set)):
                return [serial(x) for x in obj]
            if isinstance(obj, dict):
                return {str(k): serial(v) for k, v in obj.items()}
            if isinstance(obj, (Counter, defaultdict)):
                return {str(k): serial(v) for k, v in obj.items()}
            if isinstance(obj, Path):
                return str(obj)
            if hasattr(obj, 'isoformat'):
                return obj.isoformat()
            return str(obj)

        with open(out, 'w', encoding='utf-8') as f:
            json.dump(serial(data), f, indent=2, ensure_ascii=False)
        print(f"    ✅ JSON saved: {out}")
        return out

    def write_csvs(self, data: Dict[str, Any]):
        print("\n📊 Writing CSV reports...")
        rows = []
        for sample, s in data['samples'].items():
            t = s['typing']
            rows.append({
                'Sample': sample,
                'MLST': t['MLST'],
                'O_Type': t['O_Type'],
                'T3SS_Pathotype': t['T3SS'],
                'Virulence_Gene_Count': len(
                    [g for db in ('vfdb', 'ecoli_vf', 'pa_vf', 'pseudomonas_vf')
                     for g in s.get('abricate_databases', {}).get(db, [])]),
                'MobileOG_Hits': s.get('mge_hits', 0),
            })
        pd.DataFrame(rows).to_csv(self.output_dir / "sample_overview.csv", index=False)

        gc = data.get('gene_centric', {})
        total = len(data['samples']) or 1
        for cat, fname in [('amr_databases', 'amr_genes.csv'),
                           ('virulence_databases', 'virulence_genes.csv'),
                           ('bacmet_databases', 'bacmet_genes.csv'),
                           ('plasmid_databases', 'plasmid_replicons.csv')]:
            out = []
            for db, genes in gc.get(cat, {}).items():
                for g in genes:
                    out.append({
                        'Gene': g['gene'],
                        'Database': g['database'],
                        'Count': g['count'],
                        'Percentage': f"{(g['count'] / total) * 100:.1f}%",
                        'Genomes': ';'.join(g.get('genomes', [])),
                    })
            if out:
                pd.DataFrame(out).to_csv(self.output_dir / fname, index=False)

        muts = data.get('mutation_data', {}).get('mutations', [])
        if muts:
            pd.DataFrame([{
                'Gene': m['gene'], 'Mutation': m['mutation'],
                'Class': m['class'],
                'Count': m['count'], 'Genomes': ';'.join(m['genomes'])
            } for m in muts]).to_csv(self.output_dir / "mutations.csv", index=False)

        if data.get('qc_data'):
            pd.DataFrame([{'Sample': s, **m} for s, m in data['qc_data'].items()]) \
                .to_csv(self.output_dir / "fasta_qc.csv", index=False)

        mge_ps = data.get('mge_data', {}).get('per_sample', {})
        if mge_ps:
            pd.DataFrame([{'Sample': s, **v} for s, v in mge_ps.items()]) \
                .to_csv(self.output_dir / "mge_profile.csv", index=False)

        print("    ✅ CSVs written")

    def run(self) -> bool:
        print("=" * 80)
        print("🧬 PSEUDOSCOPE ULTIMATE P. AERUGINOSA REPORTER v3.0.0")
        print("=" * 80)
        print(f"📁 Input:  {self.input_dir}")
        print(f"📁 Output: {self.output_dir}")

        files = self.find_files()
        if not any([files['mlst_html'], files['qc_html'], files['amrfinder_html'],
                    files['abricate_htmls'], files['past_html'],
                    files['mutation_html'], files['mge_html']]):
            print("❌ No recognized input files found")
            return False

        data = self.integrate(files)
        if not data:
            return False

        self.write_json(data)
        self.write_csvs(data)
        self.generator.generate_main_report(data, self.output_dir)

        n = len(data['samples'])
        n_exou = len(data['patterns'].get('t3ss_to_samples', {}).get('ExoU', []))
        mge_n = len(data['mge_data'].get('per_sample', {}))
        print("\n" + "=" * 80)
        print("✅ REPORT COMPLETE")
        print("=" * 80)
        print(f"   Samples:                {n}")
        print(f"   ExoU+ isolates:         {n_exou}")
        print(f"   MGE-profiled:           {mge_n}")
        print(f"   Output dir:             {self.output_dir}")
        print("=" * 80)
        return True


def main():
    parser = argparse.ArgumentParser(
        description='PseudoScope Ultimate P. aeruginosa Reporter v3.0.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python pseudoscope_ultimate_reporter.py -i /path/to/reports

Author: Brown Beckley <brownbeckley94@gmail.com>""")
    parser.add_argument('-i', '--input-dir', required=True,
                        help='Directory containing PseudoScope reports')
    parser.add_argument('-o', '--output-dir', help='Custom output directory')
    args = parser.parse_args()

    inp = Path(args.input_dir)
    if not inp.exists():
        print(f"❌ Input directory not found: {inp}")
        sys.exit(1)

    reporter = PseudomonasUltimateReporter(inp)
    if args.output_dir:
        reporter.output_dir = Path(args.output_dir)
        reporter.output_dir.mkdir(parents=True, exist_ok=True)

    if not reporter.run():
        sys.exit(1)


if __name__ == "__main__":
    main()