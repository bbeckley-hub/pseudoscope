#!/usr/bin/env python3
"""
PSEUDOSCOPE ULTIMATE SAMPLE-CENTRIC REPORTER (v3.0.0)
================================================================================
Hybrid gene-centric + sample-centric reporter for P. aeruginosa genomic data.

- Sample Overview tab shows all typing on one row per isolate.
- Sample-centric isolate boxes (lazy-loaded) for AMR, Virulence, Plasmids,
  Bacmet2, Mutations, and MGE. Each box loads its tables on click; browser
  stays responsive on 1000+ samples.
- Dashboard shows 5 KPI cards: Total Samples, AMR Genes, Virulence Genes,
  Mutations, MGE Hits.
- Supports PseudoScope brown palette; GitHub star CTA in footer.

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School, Department of Medical Biochemistry
Version: 3.0.0
License: MIT
================================================================================
"""

import os
import sys
import json
import re
import io
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any, Optional
from datetime import datetime
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

from bs4 import BeautifulSoup


def esc(v) -> str:
    if v is None:
        return ""
    return (str(v).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


# =============================================================================
# PARSER
# =============================================================================
class PseudomonasHTMLParser:
    """Parses PseudoScope MLST/PAST/QC HTML and pseudo_* TSV summaries."""

    def __init__(self):
        self.abricate_databases = [
            'ncbi', 'card', 'resfinder', 'vfdb', 'argannot',
            'plasmidfinder', 'megares', 'ecoli_vf', 'bacmet2'
        ]
        self.abricate_tsv_files = {
            'ncbi': 'pseudo_ncbi_abricate_summary.tsv',
            'card': 'pseudo_card_abricate_summary.tsv',
            'resfinder': 'pseudo_resfinder_abricate_summary.tsv',
            'vfdb': 'pseudo_vfdb_abricate_summary.tsv',
            'argannot': 'pseudo_argannot_abricate_summary.tsv',
            'plasmidfinder': 'pseudo_plasmidfinder_abricate_summary.tsv',
            'megares': 'pseudo_megares_abricate_summary.tsv',
            'ecoli_vf': 'pseudo_ecoli_vf_abricate_summary.tsv',
            'bacmet2': 'pseudo_bacmet2_abricate_summary.tsv',
        }

    def normalize_sample_id(self, sample_id: str) -> str:
        sample = str(sample_id)
        for ext in ('.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff',
                    '.txt', '.tsv', '.csv'):
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        return sample.strip()

    def clean_st(self, st_str: str) -> str:
        if not st_str or st_str == 'ND':
            return 'ND'
        st_str = str(st_str).strip()
        if not st_str.upper().startswith('ST'):
            st_str = 'ST' + st_str
        return st_str

    # ------------------------------------------------------------------
    # TSV loaders
    # ------------------------------------------------------------------
    def load_amrfinder_from_tsv(self, input_dir: Path) -> Tuple[Dict[str, List[Dict]], Dict[str, int]]:
        path = input_dir / 'pseudo_amrfinder_summary.tsv'
        if not path.exists():
            print(f"  ⚠️ AMRfinder TSV not found: {path.name}")
            return {}, {}
        try:
            df = pd.read_csv(path, sep='\t')
        except Exception as e:
            print(f"  ❌ Failed to read AMRfinder TSV: {e}")
            return {}, {}
        sample_col = 'Genome' if 'Genome' in df.columns else df.columns[0]
        details = defaultdict(list)
        gene_counts = Counter()
        for _, row in df.iterrows():
            sample = self.normalize_sample_id(row[sample_col])
            g = {k: (None if pd.isna(v) else v) for k, v in row.to_dict().items()}
            for key in ('Gene_Symbol', 'Element_symbol', 'Element symbol', 'gene'):
                if key in g:
                    g['gene'] = g.pop(key)
                    break
            else:
                for col in df.columns:
                    if 'gene' in col.lower():
                        g['gene'] = g.pop(col)
                        break
            details[sample].append(g)
            gene_counts[g.get('gene', 'unknown')] += 1
        print(f"  ✅ Loaded AMRfinder TSV: {len(details)} samples, {sum(gene_counts.values())} hits")
        return dict(details), dict(gene_counts)

    def load_abricate_from_tsv(self, input_dir: Path) -> Tuple[Dict[str, Dict[str, List[Dict]]], Dict[str, Dict[str, int]]]:
        details = defaultdict(lambda: defaultdict(list))
        counts = defaultdict(lambda: defaultdict(int))
        for db, fname in self.abricate_tsv_files.items():
            path = input_dir / fname
            if not path.exists():
                continue
            try:
                df = pd.read_csv(path, sep='\t')
            except Exception:
                continue
            sample_col = ('genome' if 'genome' in df.columns
                          else 'file' if 'file' in df.columns
                          else df.columns[0])
            for _, row in df.iterrows():
                sample = self.normalize_sample_id(row[sample_col])
                g = {k: (None if pd.isna(v) else v) for k, v in row.to_dict().items()}
                details[sample][db].append(g)
                gene = row.get('gene', 'unknown')
                counts[db][gene] += 1
        print(f"  ✅ Loaded ABRicate TSVs: {len(details)} samples across {len(counts)} databases")
        return dict(details), dict(counts)

    def load_mutations_from_tsv(self, input_dir: Path) -> Dict[str, List[Dict]]:
        path = input_dir / 'mutation_summary.tsv'
        if not path.exists():
            print(f"  ⚠️ Mutation TSV not found: {path.name}")
            return {}
        try:
            df = pd.read_csv(path, sep='\t')
        except Exception as e:
            print(f"  ❌ Failed to read mutation TSV: {e}")
            return {}
        by_sample = defaultdict(list)
        for _, row in df.iterrows():
            sample = self.normalize_sample_id(row.get('genome', ''))
            if not sample:
                continue

            def clean(v):
                return '' if pd.isna(v) else str(v)

            by_sample[sample].append({
                'gene':      clean(row.get('gene_symbol', '')),
                'mutation':  clean(row.get('element_name', '')),
                'class':     clean(row.get('class', '')),
                'contig':    clean(row.get('contig_id', '')),
                'start':     clean(row.get('start', '')),
                'stop':      clean(row.get('stop', '')),
                'strand':    clean(row.get('strand', '')),
                'coverage':  clean(row.get('coverage', '')),
                'identity':  clean(row.get('identity', '')),
                'accession': clean(row.get('accession', '')),
            })
        print(f"  ✅ Loaded mutations TSV: {len(by_sample)} samples")
        return dict(by_sample)

    def load_mge_from_tsv(self, input_dir: Path) -> Dict[str, List[Dict]]:
        """Load per-sample mobileOG-db hits from pseudoscope_mge_hits.tsv."""
        candidates = ['pseudoscope_mge_hits.tsv', 'pseudo_mge_hits.tsv',
                      'pseudoscope_mge_all_hits.tsv']
        path = None
        for cand in candidates:
            p = input_dir / cand
            if p.exists():
                path = p
                break
        if path is None:
            print(f"  ⚠️ MGE TSV not found (tried {', '.join(candidates)})")
            return {}
        try:
            df = pd.read_csv(path, sep='\t')
        except Exception as e:
            print(f"  ❌ Failed to read MGE TSV: {e}")
            return {}
        per_sample = defaultdict(list)
        for _, row in df.iterrows():
            sample = self.normalize_sample_id(row.get('Sample', ''))
            if not sample:
                continue

            def clean(v):
                return '' if pd.isna(v) else str(v)

            per_sample[sample].append({
                'rank':             clean(row.get('Rank', '')),
                'mobileog_id':      clean(row.get('mobileOG_ID', '')),
                'gene':             clean(row.get('Gene_Name', row.get('Gene', ''))),
                'contig':           clean(row.get('Contig', '')),
                'start':            clean(row.get('Start', '')),
                'stop':             clean(row.get('Stop', '')),
                'strand':           clean(row.get('Strand', '')),
                'major_category':   clean(row.get('Major_Category', '')),
                'minor_category':   clean(row.get('Minor_Category', '')),
                'percent_identity': clean(row.get('%_Identity', row.get('Percent_Identity', ''))),
                'align_over_qlen':  clean(row.get('Alignment_Length', row.get('Align_Over_Qlen', ''))),
                'coverage_percent': clean(row.get('%_Coverage', row.get('Coverage_Percent', ''))),
            })
        if per_sample:
            total = sum(len(v) for v in per_sample.values())
            print(f"  ✅ Loaded MGE hits: {len(per_sample)} samples, {total} total hits")
        return dict(per_sample)

    # ------------------------------------------------------------------
    # HTML parsing (MLST, PAST, QC)
    # ------------------------------------------------------------------
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
                if cols:
                    row_data = [col.get_text().strip() for col in cols]
                    if len(row_data) == len(headers):
                        data.append(row_data)
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
                raw = row[sample_col]
                if not raw:
                    continue
                sample = self.normalize_sample_id(raw)
                qc = {}
                for col in df.columns:
                    if col == sample_col:
                        continue
                    val = row[col]
                    if pd.isna(val) or val == '' or val == 'ND':
                        qc[col] = 'ND'
                    else:
                        cleaned = str(val).replace('%', '').replace(',', '').strip()
                        try:
                            qc[col] = float(cleaned)
                        except Exception:
                            qc[col] = str(val)
                results[sample] = qc
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
            headers = [c.get_text().strip() for c in rows[0].find_all(['th', 'td'])]
            data = []
            for row in rows[1:]:
                cols = row.find_all(['td', 'th'])
                if cols:
                    row_data = [c.get_text().strip() for c in cols]
                    if len(row_data) >= 2:
                        data.append(row_data)
            if not data:
                return {}
            df = pd.DataFrame(data)
            if len(df.columns) > len(headers):
                df = df.iloc[:, :len(headers)]
            df.columns = [c.strip() for c in headers[:len(df.columns)]]
            if 'Sample' in df.columns:
                df['norm'] = df['Sample'].apply(self.normalize_sample_id)
            else:
                df['norm'] = df.iloc[:, 0].apply(self.normalize_sample_id)
            results = {}
            for _, row in df.iterrows():
                sample = row['norm']
                st = 'ND'
                if 'ST' in df.columns:
                    st = self.clean_st(str(row['ST']) if pd.notna(row['ST']) else 'ND')
                allele_profile = 'ND'
                if 'Allele Profile' in df.columns and pd.notna(row['Allele Profile']):
                    allele_profile = str(row['Allele Profile'])
                results[sample] = {'ST': st, 'Allele_Profile': allele_profile}
            print(f"    ✓ Found {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing MLST: {e}")
            return {}

    def parse_past_report(self, file_path: Path) -> Dict[str, Dict]:
        print(f"  🧬 Parsing PAST serotype: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            df = self.parse_html_table(html_content, 0)
            if df.empty:
                return {}
            df.columns = [c.strip() for c in df.columns]
            sample_col = None
            for col in df.columns:
                if 'sample' in col.lower() or 'id' in col.lower():
                    sample_col = col
                    break
            if not sample_col:
                sample_col = df.columns[0]
            df['norm'] = df[sample_col].apply(self.normalize_sample_id)
            results = {}
            for _, row in df.iterrows():
                sample = row['norm']
                o = 'ND'
                for col in ('Type', 'O-type', 'O_Type'):
                    if col in df.columns and pd.notna(row.get(col)):
                        o = str(row[col])
                        break
                results[sample] = {'O_Type': o}
            print(f"    ✓ Found {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing PAST: {e}")
            return {}


# =============================================================================
# ANALYZER
# =============================================================================
class PseudomonasDataAnalyzer:
    """Patterns, gene-centric tables, T3SS/carbapenemase derivation."""

    def __init__(self):
        self.carbapenemase_families = {
            'KPC': ['kpc'], 'NDM': ['ndm'], 'VIM': ['vim'], 'IMP': ['imp'],
            'GES': ['ges'], 'SPM': ['spm'],
            'OXA-23-like': ['oxa-23'],
            'OXA-24/40-like': ['oxa-24', 'oxa-40'],
            'OXA-48-like': ['oxa-48', 'oxa-181', 'oxa-232', 'oxa-204', 'oxa-244'],
            'OXA-58-like': ['oxa-58'],
        }
        self.oxa_class_d_esbl_families = {
            'OXA-2-like': ['oxa-2', 'oxa-3', 'oxa-15', 'oxa-20',
                           'oxa-21', 'oxa-32', 'oxa-34', 'oxa-45'],
            'OXA-10-like': ['oxa-10', 'oxa-11', 'oxa-13', 'oxa-14',
                            'oxa-16', 'oxa-17', 'oxa-19', 'oxa-28',
                            'oxa-35', 'oxa-74', 'oxa-101', 'oxa-128'],
        }
        self.intrinsic_oxa_families = {
            'OXA-50-like': ['oxa-50', 'oxa-395', 'oxa-396'],
        }
        self.other_esbl_families = {
            'PER': ['per'], 'VEB': ['veb'], 'BEL': ['bel'],
            'CTX-M': ['ctx-m'], 'TEM': ['tem'], 'SHV': ['shv'],
        }
        self.colistin_resistance = {f'mcr-{i}' for i in range(1, 11)}
        self.methyltransferases = {
            'arma', 'rmta', 'rmtb', 'rmtc', 'rmtd', 'rmte', 'rmtf', 'rmtg', 'rmth',
        }
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
        }
        self.critical_amr_genes = set()
        for d in (self.carbapenemase_families, self.oxa_class_d_esbl_families,
                  self.intrinsic_oxa_families, self.other_esbl_families):
            for sub in d.values():
                self.critical_amr_genes.update(sub)
        self.critical_amr_genes.update(self.colistin_resistance)
        self.critical_amr_genes.update(self.methyltransferases)

        self.critical_virulence_genes = {
            'exou', 'exos', 'exot', 'exoy', 'exoa', 'toxa',
            'plch', 'plcn', 'lasb', 'lasa', 'apra',
        }

    def normalize_gene_name(self, raw: str) -> str:
        name = str(raw).lower().strip()
        name = name.replace("beta-lactamase ", "")
        if name.startswith("bla"):
            name = name[3:]
        tokens = re.split(r'[_/]', name)
        return tokens[0] if tokens else name

    def match_gene(self, input_gene: str) -> str:
        clean = self.normalize_gene_name(input_gene)
        if clean in self.critical_amr_genes:
            return clean
        for base in self.critical_amr_genes:
            if clean.startswith(base + "-") or clean == base:
                return base
        return ""

    def classify_carbapenemase_family(self, raw_gene: str) -> Tuple[str, str]:
        clean = self.normalize_gene_name(raw_gene)
        if not clean:
            return '', ''
        for fam, pats in self.carbapenemase_families.items():
            if any(clean == p or clean.startswith(p + '-') or clean.startswith(p) for p in pats):
                return fam, 'acquired_carbapenemase'
        for fam, pats in self.oxa_class_d_esbl_families.items():
            if any(clean == p or clean.startswith(p + '-') or clean.startswith(p) for p in pats):
                return fam, 'oxa_esbl'
        for fam, pats in self.intrinsic_oxa_families.items():
            if any(clean == p or clean.startswith(p + '-') or clean.startswith(p) for p in pats):
                return fam, 'intrinsic_oxa'
        for fam, pats in self.other_esbl_families.items():
            if any(clean == p or clean.startswith(p + '-') or clean.startswith(p) for p in pats):
                return fam, 'other_esbl'
        for pat in self.colistin_resistance:
            if clean == pat or clean.startswith(pat):
                return pat, 'colistin'
        for pat in self.methyltransferases:
            if clean == pat or clean.startswith(pat):
                return pat, 'methyltransferase'
        for prefix in self.unclassified_prefixes:
            if clean == prefix:
                return f'{prefix.upper()} (class-level)', 'ambiguous_class_level'
        m = re.match(r'^([a-z]+(?:-[a-z]+)?)-(\d+)$', clean)
        if m and m.group(1) in self.unclassified_prefixes:
            return f'{self.unclassified_prefixes[m.group(1)]} — {clean}', 'unclassified'
        return '', ''

    @staticmethod
    def derive_t3ss_pathotype(genes: List[str]) -> str:
        gl = [str(g).lower() for g in genes]
        exou = any('exou' in g for g in gl)
        exos = any('exos' in g for g in gl)
        if exou and exos:
            return 'ExoU+ExoS'
        if exou:
            return 'ExoU'
        if exos:
            return 'ExoS'
        return 'Neither'

    def create_gene_centric_tables(self, integrated_data: Dict[str, Any]) -> Dict[str, Any]:
        gc = {
            'amr_databases': {}, 'virulence_databases': {},
            'plasmid_databases': {}, 'bacmet_databases': {},
            'combined_gene_frequencies': [],
        }
        if 'amrfinder' in integrated_data.get('gene_frequencies', {}):
            amr = integrated_data['gene_frequencies']['amrfinder']
            lst = [{'gene': g, 'database': 'AMRfinder',
                    'frequency': str(d.get('count', 0)),
                    'count': int(d.get('count', 0)),
                    'genomes': d.get('genomes', [])}
                   for g, d in amr.items()]
            if lst:
                gc['amr_databases']['amrfinder'] = sorted(
                    lst, key=lambda x: x['count'], reverse=True)

        if 'abricate' in integrated_data.get('gene_frequencies', {}):
            for db, db_genes in integrated_data['gene_frequencies']['abricate'].items():
                lst = [{'gene': g, 'database': db.upper(),
                        'frequency': str(d.get('count', 0)),
                        'count': int(d.get('count', 0)),
                        'genomes': d.get('genomes', [])}
                       for g, d in db_genes.items()]
                if not lst:
                    continue
                lst.sort(key=lambda x: x['count'], reverse=True)
                if db in ('vfdb', 'ecoli_vf', 'pa_vf', 'pseudomonas_vf'):
                    gc['virulence_databases'][db] = lst
                elif db == 'plasmidfinder':
                    gc['plasmid_databases'][db] = lst
                elif db == 'bacmet2':
                    gc['bacmet_databases'][db] = lst
                else:
                    gc['amr_databases'][db] = lst

        all_genes = []
        for t in ('amr_databases', 'virulence_databases',
                  'plasmid_databases', 'bacmet_databases'):
            for gs in gc.get(t, {}).values():
                all_genes.extend(gs)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        gc['combined_gene_frequencies'] = all_genes
        return gc

    def create_cross_genome_patterns(self, integrated_data: Dict[str, Any]) -> Dict[str, Any]:
        samples = integrated_data.get('samples', {})
        P = {
            'mlst_distribution': Counter(),
            'o_type_distribution': Counter(),
            't3ss_distribution': Counter(),
            'st_to_samples': defaultdict(list),
            'o_type_to_samples': defaultdict(list),
            't3ss_to_samples': defaultdict(list),
            'gene_cooccurrence': defaultdict(Counter),
            'high_risk_combinations': [],
            'carbapenemase_by_sample': defaultdict(
                lambda: defaultdict(lambda: defaultdict(dict))),
            'carbapenemase_family_distribution': Counter(),
            'carbapenemase_category_distribution': Counter(),
        }
        COMBO_FIELDS = {
            'mlst_otype':      ['MLST', 'O_Type'],
            'mlst_t3ss':       ['MLST', 'T3SS'],
            'otype_t3ss':      ['O_Type', 'T3SS'],
            'mlst_otype_t3ss': ['MLST', 'O_Type', 'T3SS'],
        }
        combos = {k: defaultdict(list) for k in COMBO_FIELDS}

        sample_genes = defaultdict(list)
        gc = integrated_data.get('gene_centric', {})
        for t in ('amr_databases', 'virulence_databases',
                  'plasmid_databases', 'bacmet_databases'):
            for gs in gc.get(t, {}).values():
                for g in gs:
                    for genome in g['genomes']:
                        if g['gene'] not in sample_genes[genome]:
                            sample_genes[genome].append(g['gene'])

        for sample, data in samples.items():
            seen = set()
            sources = [(raw, 'AMRfinder')
                       for raw in data.get('amrfinder', {}).get('all_genes', [])]
            for db, gs in data.get('abricate_databases', {}).items():
                for raw in gs:
                    sources.append((raw, db.upper()))
            for raw, db_name in sources:
                fam, cat = self.classify_carbapenemase_family(raw)
                if not fam or not cat:
                    continue
                key = (fam, cat, raw)
                if key in seen:
                    continue
                seen.add(key)
                bucket = P['carbapenemase_by_sample'][sample][cat][fam]
                if raw not in bucket:
                    bucket[raw] = []
                if db_name not in bucket[raw]:
                    bucket[raw].append(db_name)
                if cat == 'acquired_carbapenemase':
                    P['carbapenemase_family_distribution'][fam] += 1
                P['carbapenemase_category_distribution'][cat] += 1

        def ok(v):
            return v and v not in ('Not Assigned', 'ND', '', 'nan')

        for sample, data in samples.items():
            t = data.get('typing', {})
            values = {
                'MLST':   t.get('MLST', 'ND'),
                'O_Type': t.get('O_Type', 'ND'),
                'T3SS':   t.get('T3SS', 'ND'),
            }
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
                    combos[key][label].append(sample)

            amr_genes = data.get('amrfinder', {}).get('all_genes', [])
            vir_genes = []
            for vdb in ('vfdb', 'ecoli_vf'):
                vir_genes.extend(data.get('abricate_databases', {}).get(vdb, []))
            crit_amr = [g for g in amr_genes if self.match_gene(g)]
            crit_vir = [g for g in vir_genes
                        if any(c in str(g).lower() for c in self.critical_virulence_genes)]
            if crit_amr and crit_vir:
                P['high_risk_combinations'].append({
                    'sample': sample,
                    'mlst': values['MLST'],
                    'o_type': values['O_Type'],
                    't3ss': values['T3SS'],
                    'critical_amr_genes': crit_amr,
                    'critical_virulence_genes': crit_vir,
                })

            genes = sample_genes.get(sample, [])
            for i, g1 in enumerate(genes):
                for g2 in genes[i + 1:]:
                    P['gene_cooccurrence'][g1][g2] += 1

        for k in combos:
            P[k] = dict(combos[k])
        return P


# =============================================================================
# HTML GENERATOR
# =============================================================================
class PseudomonasHTMLGenerator:
    """Builds the interactive multi-tab HTML report with lazy-loaded boxes."""

    def __init__(self, analyzer: PseudomonasDataAnalyzer):
        self.analyzer = analyzer
        self.tab_colors = {
            'summary': '#8B4513', 'samples': '#2E7D32',
            'amr': '#C62828', 'virulence': '#E65100',
            'plasmid': '#3F51B5', 'bacmet': '#607D8B',
            'mutation': '#00BCD4', 'mge': '#16A085',
            'citation': '#00838F', 'guide': '#4527A0',
            'calltoaction': '#2E7D32', 'export': '#37474F',
        }

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    def _credit_bar(self, color: str, icon: str, title: str, body: str) -> str:
        return f'''
        <div class="scientific-note" style="background:linear-gradient(135deg,#fdf6ec 0%,#f5e6d3 100%);border-left:6px solid {color};margin-bottom:20px;padding:15px;border-radius:8px;">
            <div style="display:flex;align-items:center;gap:12px;flex-wrap:wrap;">
                <span style="font-size:1.4em;">{icon}</span>
                <div>
                    <strong style="font-size:1.1em;color:{color};">{title}</strong><br>
                    <span style="font-size:.95em;color:#333;">{body}</span>
                </div>
            </div>
        </div>'''

    def _alert(self, kind: str, icon: str, html_body: str) -> str:
        return f'<div class="alert-box alert-{kind}"><i class="fas {icon} fa-2x"></i><div>{html_body}</div></div>'

    def _stat_card(self, value, label: str, color: str = '#8B4513', icon: str = '') -> str:
        icon_html = (f'<i class="fas {icon} fa-2x" style="opacity:.9;margin-bottom:8px;"></i>'
                     if icon else '')
        return f'''
        <div class="stat-card" style="background:linear-gradient(135deg,{color} 0%,{color}dd 100%);">
            {icon_html}
            <div class="stat-value">{value}</div>
            <div class="stat-label">{label}</div>
        </div>'''

    def _gene_family_info(self, border: str, title: str, items: list) -> str:
        html = (f'<div style="margin:10px 0 20px 0;background:#fdf6ec;padding:15px;'
                f'border-radius:8px;font-size:.9em;border-left:4px solid {border};">'
                f'<strong><i class="fas fa-info-circle"></i> {title}</strong><br>')
        for n, d in items:
            html += f'• <strong>{n}</strong> – {d}<br>'
        return html + '</div>'

    # ------------------------------------------------------------------
    # CSS
    # ------------------------------------------------------------------
    def _get_css(self) -> str:
        return """
<style>
:root{
  --brown:#8B4513;--brown-light:#A0522D;--brown-dark:#5C2E0B;
  --summary-color:#8B4513;--samples-color:#2E7D32;
  --amr-color:#C62828;--virulence-color:#E65100;--plasmid-color:#3F51B5;
  --bacmet-color:#607D8B;--mutation-color:#00BCD4;--mge-color:#16A085;
  --citation-color:#00838F;--guide-color:#4527A0;
  --calltoaction-color:#2E7D32;--export-color:#37474F;
}
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:'Segoe UI',Tahoma,Geneva,Verdana,sans-serif;line-height:1.6;color:#333;background:#f5f5f5;min-width:1200px}
.container{max-width:none;margin:0 auto;padding:20px;width:100%;overflow-x:auto}
.main-header{background:linear-gradient(135deg,#2A1408 0%,#5C2E0B 40%,#8B4513 100%);color:white;padding:30px;border-radius:15px;box-shadow:0 10px 30px rgba(0,0,0,.2);margin-bottom:30px;text-align:center}
.main-header h1{font-size:2.8em;margin-bottom:10px;color:white}
.metadata-bar{background:rgba(255,255,255,.1);padding:15px;border-radius:10px;margin:20px 0;display:flex;justify-content:space-around;flex-wrap:wrap;gap:15px}
.metadata-item{display:flex;align-items:center;gap:8px;font-size:.95em}
.dashboard-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:20px;margin-bottom:30px}
.dashboard-card{background:white;padding:22px;border-radius:12px;box-shadow:0 5px 20px rgba(0,0,0,.1);text-align:center;transition:all .3s ease;cursor:pointer;border-left:5px solid var(--brown);position:relative;overflow:hidden}
.dashboard-card:hover{transform:translateY(-10px);box-shadow:0 15px 30px rgba(0,0,0,.2)}
.card-number{font-size:2.6em;font-weight:bold;margin:10px 0;background:linear-gradient(90deg,#5C2E0B,#8B4513);-webkit-background-clip:text;-webkit-text-fill-color:transparent}
.card-label{font-size:.9em;color:#555;font-weight:600}
.tab-navigation{display:flex;gap:5px;margin-bottom:20px;flex-wrap:wrap;background:white;padding:15px;border-radius:12px;box-shadow:0 5px 20px rgba(0,0,0,.1);position:sticky;top:10px;z-index:100}
.tab-button{padding:11px 18px;background:#f5f5f5;border:none;border-radius:8px;cursor:pointer;font-weight:600;color:#666;transition:all .3s ease;display:flex;align-items:center;gap:8px;font-size:.88em}
.tab-button.active{color:white}
.tab-button.summary.active{background:var(--summary-color)}
.tab-button.samples.active{background:var(--samples-color)}
.tab-button.amr.active{background:var(--amr-color)}
.tab-button.virulence.active{background:var(--virulence-color)}
.tab-button.plasmid.active{background:var(--plasmid-color)}
.tab-button.bacmet.active{background:var(--bacmet-color)}
.tab-button.mutation.active{background:var(--mutation-color)}
.tab-button.mge.active{background:var(--mge-color)}
.tab-button.citation.active{background:var(--citation-color)}
.tab-button.guide.active{background:var(--guide-color)}
.tab-button.calltoaction.active{background:var(--calltoaction-color)}
.tab-button.export.active{background:var(--export-color)}
.tab-content{display:none;background:white;padding:30px;border-radius:15px;box-shadow:0 10px 30px rgba(0,0,0,.1);margin-bottom:30px;animation:fadeIn .5s ease;width:100%;overflow-x:auto}
.tab-content.active{display:block}
@keyframes fadeIn{from{opacity:0;transform:translateY(20px)}to{opacity:1;transform:translateY(0)}}
.section-header{color:#2c3e50;margin-bottom:25px;padding-bottom:15px;border-bottom:3px solid var(--brown);font-size:1.8em;display:flex;align-items:center;justify-content:space-between}
.data-table{width:100%;border-collapse:collapse;margin:20px 0;font-size:.95em;box-shadow:0 2px 10px rgba(0,0,0,.1);border-radius:8px;overflow:hidden}
.data-table th{background:#5C2E0B;color:white;padding:14px;text-align:left;font-weight:600;position:sticky;top:0;white-space:nowrap;cursor:pointer}
.data-table th:hover{background:#3D1E08}
.data-table td{padding:12px;border-bottom:1px solid #e0e0e0;vertical-align:top}
.data-table tr:hover td{background:#fdf6ec}
.scrollable-table{max-height:600px;overflow-y:auto;border:1px solid #e0e0e0;border-radius:8px;margin:20px 0;width:100%}
.master-scrollable-container{width:100%;overflow-x:auto;border:1px solid #e0e0e0;border-radius:8px;margin:20px 0}
.genome-list{display:flex;flex-wrap:wrap;gap:5px;max-height:220px;overflow-y:auto;padding:6px;background:#fdf6ec;border-radius:5px}
.genome-tag{display:inline-block;background:#e6ffe6;color:#5C2E0B;padding:3px 10px;border-radius:12px;font-size:.85em;border:1px solid #E8C39E;white-space:nowrap;margin:2px}
.genome-tag.highlight{background-color:#ffff99 !important;color:#000 !important;border:1px solid #ffc107}
.search-box{width:100%;padding:12px;margin-bottom:20px;border:2px solid #e0e0e0;border-radius:8px;font-size:1em;transition:all .3s ease}
.search-box:focus{outline:none;border-color:#8B4513;box-shadow:0 0 0 3px rgba(139,69,19,.1)}
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
.badge{display:inline-block;padding:5px 15px;border-radius:20px;font-size:.85em;font-weight:600;margin:2px}
.badge-critical{background:#8B0000;color:white}
.badge-exou{background:#C62828;color:white}
.badge-exos{background:#1565C0;color:white}
.alert-box{padding:20px;border-radius:10px;margin:20px 0;display:flex;align-items:flex-start;gap:20px;border-left:5px solid}
.alert-success{background:#d4edda;color:#155724;border-left-color:#28a745}
.alert-warning{background:#fff3cd;color:#856404;border-left-color:#ffc107}
.alert-danger{background:#f5e6d3;color:#3D1E08;border-left-color:#8B4513}
.alert-info{background:#d1ecf1;color:#0c5460;border-left-color:#17a2b8}
.scientific-note{border-radius:8px}
.database-section{margin:30px 0;padding:25px;border-radius:12px;background:#fdf6ec;box-shadow:0 3px 15px rgba(0,0,0,.08)}
.print-section-btn{background:#8B4513;color:white;border:none;border-radius:5px;padding:8px 15px;cursor:pointer;display:flex;align-items:center;gap:5px;font-size:.9em}
.print-section-btn:hover{background:#A0522D}
.footer{text-align:center;padding:30px;color:white;margin-top:40px;border-radius:15px;background:linear-gradient(135deg,#2A1408 0%,#5C2E0B 100%)}
.footer a{color:#E8C39E;text-decoration:none}
.footer a:hover{text-decoration:underline}
.footer .star-cta{display:inline-block;background:linear-gradient(135deg,#ffc107,#ff9800);color:#2A1408;font-weight:800;padding:10px 22px;border-radius:24px;margin-top:12px;text-decoration:none;box-shadow:0 3px 10px rgba(0,0,0,.25);transition:transform .2s}
.footer .star-cta:hover{transform:translateY(-2px);text-decoration:none}
.sort-icon{margin-left:5px;font-size:.8em;opacity:.6}
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
.copy-btn{background:#8B4513;color:white;border:none;padding:4px 14px;border-radius:16px;cursor:pointer;font-size:.82em;font-weight:600;transition:background .2s}
.copy-btn:hover{background:#A0522D}

.isolate-box{border:1px solid #E8C39E;border-radius:12px;margin-bottom:18px;padding:18px 20px;background:#fdf6ec;box-shadow:0 2px 8px rgba(139,69,19,.08);transition:box-shadow .2s}
.isolate-box:hover{box-shadow:0 4px 14px rgba(139,69,19,.15)}
.isolate-box .sample-header{display:flex;align-items:center;gap:12px;flex-wrap:wrap}
.isolate-box .sample-header h3{font-size:1.25em;margin:0;color:#5C2E0B}
.isolate-box .total-badge{background:#8B4513;color:white;padding:4px 14px;border-radius:20px;font-weight:bold;font-size:.85em}
.isolate-box .typing-info{display:flex;flex-wrap:wrap;gap:6px}
.isolate-box .toggle-btn{background:#8B4513;color:white;border:none;padding:8px 16px;border-radius:8px;cursor:pointer;font-weight:600;font-size:.85em;display:inline-flex;align-items:center;gap:6px;margin-left:auto;transition:background .2s}
.isolate-box .toggle-btn:hover{background:#A0522D}
.isolate-box .box-details{margin-top:14px}
.database-table-wrapper{margin:14px 0;overflow-x:auto;border:1px solid #E8C39E;border-radius:8px;background:white}
.database-table-wrapper table{width:100%;border-collapse:collapse;font-size:.85em;min-width:800px}
.database-table-wrapper table th{background:#5C2E0B;color:white;padding:8px 12px;text-align:left;white-space:nowrap}
.database-table-wrapper table td{padding:8px 12px;border-bottom:1px solid #f0f0f0;white-space:nowrap}
.database-table-wrapper table tr:hover{background:#fdf6ec}
.db-title{font-weight:bold;color:#8B4513;margin:10px 0 5px 0;font-size:1.05em;border-left:4px solid #8B4513;padding-left:10px}
.filter-controls{display:flex;flex-wrap:wrap;gap:10px;align-items:center;background:#fdf6ec;padding:15px;border-radius:8px;margin-bottom:20px}
.filter-controls select{padding:10px;border-radius:8px;border:2px solid #E8C39E;background:white;min-width:150px}
.results-counter{font-size:.9em;color:#5C2E0B;font-weight:600;padding:8px 12px;background:#f5e6d3;border-radius:6px;white-space:nowrap}

@media print{body *{visibility:hidden}.tab-content.active,.tab-content.active *{visibility:visible}.tab-content.active{position:absolute;left:0;top:0;width:100%;padding:20px;box-shadow:none;border-radius:0}.print-section-btn,.tab-navigation,.dashboard-grid,.search-box,.action-buttons,.filter-controls,.toggle-btn{display:none !important}.isolate-box .box-details{display:block !important}}
</style>"""

    # ------------------------------------------------------------------
    # JS
    # ------------------------------------------------------------------
    def _get_js(self, typing_json: str) -> str:
        return f"""
<script>
var sampleTyping = {typing_json};
window.PSEUDOSCOPE_BOX_DATA = window.PSEUDOSCOPE_BOX_DATA || {{}};

function escapeHtml(str) {{
    return String(str)
        .replace(/&/g,'&amp;').replace(/</g,'&lt;')
        .replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}}

function switchTab(tabName) {{
    document.querySelectorAll('.tab-content').forEach(t => t.classList.remove('active'));
    document.querySelectorAll('.tab-button').forEach(b => b.classList.remove('active'));
    var content = document.getElementById(tabName + '-tab');
    var button = document.querySelector('.tab-button.' + tabName);
    if (content) content.classList.add('active');
    if (button) button.classList.add('active');
    if (event && event.currentTarget) event.currentTarget.classList.add('active');
    window.location.hash = tabName;
}}

function searchTable(tableId, searchId) {{
    var input = document.getElementById(searchId);
    if (!input) return;
    var filter = input.value.toUpperCase();
    var table = document.getElementById(tableId);
    if (!table || !table.tBodies[0]) return;
    var rows = table.tBodies[0].rows;
    for (var i = 0; i < rows.length; i++) {{
        var found = false;
        var cells = rows[i].getElementsByTagName('td');
        for (var j = 0; j < cells.length; j++) {{
            if (cells[j] && (cells[j].textContent || cells[j].innerText).toUpperCase().indexOf(filter) > -1) {{
                found = true; break;
            }}
        }}
        rows[i].style.display = found ? '' : 'none';
    }}
}}

function highlightGenome(tableId, searchId) {{
    var el = document.getElementById(searchId);
    if (!el) return;
    var filter = el.value.toUpperCase().trim();
    var table = document.getElementById(tableId);
    if (!table) return;
    table.querySelectorAll('.genome-tag').forEach(function(t) {{
        t.classList.remove('highlight');
        if (filter && t.textContent.toUpperCase().indexOf(filter) > -1) t.classList.add('highlight');
    }});
}}

// ========== LAZY BOX RENDERING ==========
function toggleBoxDetails(tabId, sample) {{
    var box = document.querySelector('#' + tabId + '-tab .isolate-box[data-sample="' + CSS.escape(sample) + '"]');
    if (!box) return;
    var details = box.querySelector('.box-details');
    var btn = box.querySelector('.toggle-btn');
    var isOpen = details.style.display === 'block';
    if (isOpen) {{
        details.style.display = 'none';
        if (btn) btn.innerHTML = '<i class="fas fa-chevron-down"></i> Show Details';
        return;
    }}
    if (!details.innerHTML.trim()) {{
        if (tabId === 'mutation') {{
            details.innerHTML = buildMutationTable(sample);
        }} else if (tabId === 'mge') {{
            details.innerHTML = buildMgeTable(sample);
        }} else {{
            details.innerHTML = buildBoxContent(tabId, sample);
        }}
    }}
    details.style.display = 'block';
    if (btn) btn.innerHTML = '<i class="fas fa-chevron-up"></i> Hide Details';
}}

function buildBoxContent(tabId, sample) {{
    var root = window.PSEUDOSCOPE_BOX_DATA || {{}};
    var byTab = root[tabId] || {{}};
    var entry = byTab[sample];
    if (!entry) return '<p>No data available</p>';
    var html = '';
    if (entry.amrfinder && entry.amrfinder.length) {{
        html += buildDBTable('AMRfinder', entry.amrfinder, 'amrfinder');
    }}
    if (entry.abricate) {{
        for (var db in entry.abricate) {{
            var genes = entry.abricate[db];
            if (genes && genes.length) html += buildDBTable(db.toUpperCase(), genes, db);
        }}
    }}
    return html || '<p style="color:#888;font-style:italic;">No gene details</p>';
}}

function buildDBTable(dbName, genes, dbKey) {{
    if (!genes || !genes.length) return '';
    var keys = Object.keys(genes[0]);
    var priority = ['gene', 'product', 'coverage_percent', 'identity_percent',
                    'accession', 'contig', 'start', 'stop', 'class', 'subclass',
                    'scope', 'resistance'];
    var ordered = priority.filter(function(k) {{ return keys.indexOf(k) !== -1; }});
    keys.forEach(function(k) {{ if (ordered.indexOf(k) === -1) ordered.push(k); }});
    var html = '<div class="database-table-wrapper" data-db="' + dbKey + '">';
    html += '<div class="db-title">' + escapeHtml(dbName) + '</div>';
    html += '<table><thead><tr>';
    ordered.forEach(function(c) {{
        var display = c.replace(/_/g, ' ').replace(/\\b\\w/g, function(x) {{ return x.toUpperCase(); }});
        html += '<th>' + escapeHtml(display) + '</th>';
    }});
    html += '</tr></thead><tbody>';
    genes.forEach(function(g) {{
        html += '<tr>';
        ordered.forEach(function(c) {{
            var v = (g[c] == null) ? '' : g[c];
            html += '<td>' + escapeHtml(v) + '</td>';
        }});
        html += '</tr>';
    }});
    html += '</tbody></table></div>';
    return html;
}}

function buildMutationTable(sample) {{
    var root = window.PSEUDOSCOPE_BOX_DATA || {{}};
    var byTab = root['mutation'] || {{}};
    var entry = byTab[sample];
    if (!entry || !entry.mutations || !entry.mutations.length) return '<p>No mutations</p>';
    var cols = ['gene','mutation','class','contig','start','stop','strand','coverage','identity','accession'];
    var html = '<div class="database-table-wrapper" data-db="mutations"><div class="db-title">Mutations</div><table><thead><tr>';
    cols.forEach(function(c) {{
        var display = c.replace(/_/g, ' ').replace(/\\b\\w/g, function(x) {{ return x.toUpperCase(); }});
        html += '<th>' + display + '</th>';
    }});
    html += '</tr></thead><tbody>';
    entry.mutations.forEach(function(m) {{
        html += '<tr>';
        cols.forEach(function(c) {{ html += '<td>' + escapeHtml((m[c] == null) ? '' : m[c]) + '</td>'; }});
        html += '</tr>';
    }});
    html += '</tbody></table></div>';
    return html;
}}

function buildMgeTable(sample) {{
    var root = window.PSEUDOSCOPE_BOX_DATA || {{}};
    var byTab = root['mge'] || {{}};
    var entry = byTab[sample];
    if (!entry || !entry.mge || !entry.mge.length) return '<p>No MGE hits</p>';
    var cols = ['rank', 'mobileog_id', 'gene', 'contig', 'start', 'stop', 'strand',
                'major_category', 'minor_category', 'percent_identity',
                'align_over_qlen', 'coverage_percent'];
    var headers = ['Rank', 'mobileOG ID', 'Gene', 'Contig', 'Start', 'Stop', 'Strand',
                   'Major Category', 'Minor Category', '% Identity',
                   'Align / Qlen', 'Coverage %'];
    var html = '<div class="database-table-wrapper" data-db="mge">';
    html += '<div class="db-title">MGE (mobileOG-db)</div><table><thead><tr>';
    for (var h = 0; h < headers.length; h++) {{
        html += '<th>' + escapeHtml(headers[h]) + '</th>';
    }}
    html += '</tr></thead><tbody>';
    for (var i = 0; i < entry.mge.length; i++) {{
        html += '<tr>';
        for (var j = 0; j < cols.length; j++) {{
            var v = (entry.mge[i][cols[j]] == null) ? '' : entry.mge[i][cols[j]];
            html += '<td>' + escapeHtml(v) + '</td>';
        }}
        html += '</tr>';
    }}
    html += '</tbody></table></div>';
    return html;
}}

function expandAllBoxes(tabId) {{
    var boxes = document.querySelectorAll('#' + tabId + '-tab .isolate-box');
    var visible = [];
    boxes.forEach(function(b) {{ if (b.style.display !== 'none') visible.push(b); }});
    if (visible.length > 100) {{
        if (!confirm('You have ' + visible.length + ' visible samples. Expanding all may take a moment. Continue?')) return;
    }}
    var i = 0, batch = 20;
    function step() {{
        var end = Math.min(i + batch, visible.length);
        for (; i < end; i++) {{
            var box = visible[i];
            var sample = box.getAttribute('data-sample');
            var details = box.querySelector('.box-details');
            var btn = box.querySelector('.toggle-btn');
            if (details && details.style.display !== 'block') {{
                if (!details.innerHTML.trim()) {{
                    if (tabId === 'mutation') {{
                        details.innerHTML = buildMutationTable(sample);
                    }} else if (tabId === 'mge') {{
                        details.innerHTML = buildMgeTable(sample);
                    }} else {{
                        details.innerHTML = buildBoxContent(tabId, sample);
                    }}
                }}
                details.style.display = 'block';
                if (btn) btn.innerHTML = '<i class="fas fa-chevron-up"></i> Hide Details';
            }}
        }}
        if (i < visible.length) setTimeout(step, 30);
    }}
    step();
}}

function collapseAllBoxes(tabId) {{
    document.querySelectorAll('#' + tabId + '-tab .isolate-box').forEach(function(box) {{
        var details = box.querySelector('.box-details');
        if (details) details.style.display = 'none';
        var btn = box.querySelector('.toggle-btn');
        if (btn) btn.innerHTML = '<i class="fas fa-chevron-down"></i> Show Details';
    }});
}}

function filterBoxes(tabId) {{
    var searchEl = document.getElementById('search-' + tabId);
    var dbEl = document.getElementById('dbFilter-' + tabId);
    var search = (searchEl ? searchEl.value : '').toUpperCase();
    var dbFilter = dbEl ? dbEl.value : 'all';
    var boxes = document.querySelectorAll('#' + tabId + '-tab .isolate-box');
    var visible = 0;
    boxes.forEach(function(box) {{
        var sample = box.getAttribute('data-sample') || '';
        var show = (!search || sample.toUpperCase().indexOf(search) !== -1);
        box.style.display = show ? '' : 'none';
        if (show) visible++;
        if (show && dbEl) {{
            var wrappers = box.querySelectorAll('.database-table-wrapper');
            wrappers.forEach(function(w) {{
                var dbName = w.getAttribute('data-db') || '';
                var keep = (dbFilter === 'all' || dbName === dbFilter);
                w.style.display = keep ? '' : 'none';
                var title = w.previousElementSibling;
                if (title && title.classList.contains('db-title')) title.style.display = keep ? '' : 'none';
            }});
        }}
    }});
    var counter = document.getElementById('counter-' + tabId);
    if (counter) counter.textContent = visible + ' shown';
}}

function resetBoxFilters(tabId) {{
    var s = document.getElementById('search-' + tabId);
    var d = document.getElementById('dbFilter-' + tabId);
    if (s) s.value = '';
    if (d) d.value = 'all';
    filterBoxes(tabId);
}}

function sortTable(tableId, colIndex, type) {{
    var table = document.getElementById(tableId);
    if (!table || !table.tBodies[0]) return;
    var tbody = table.tBodies[0];
    var rows = Array.from(tbody.rows);
    var asc = table.getAttribute('data-sort-dir') !== 'asc';
    rows.sort(function(a, b) {{
        var av = a.cells[colIndex].innerText.trim();
        var bv = b.cells[colIndex].innerText.trim();
        if (type === 'number') {{
            av = parseFloat(av.replace(/,/g, '')) || 0;
            bv = parseFloat(bv.replace(/,/g, '')) || 0;
            return asc ? av - bv : bv - av;
        }}
        return asc ? av.localeCompare(bv) : bv.localeCompare(av);
    }});
    tbody.append.apply(tbody, rows);
    table.setAttribute('data-sort-dir', asc ? 'asc' : 'desc');
}}

function printSection(sectionId) {{
    var content = document.getElementById(sectionId);
    if (!content) return;
    content.querySelectorAll('.isolate-box').forEach(function(box) {{
        var tabId = box.getAttribute('data-tab');
        var sample = box.getAttribute('data-sample');
        var details = box.querySelector('.box-details');
        if (details && !details.innerHTML.trim()) {{
            if (tabId === 'mutation') {{
                details.innerHTML = buildMutationTable(sample);
            }} else if (tabId === 'mge') {{
                details.innerHTML = buildMgeTable(sample);
            }} else {{
                details.innerHTML = buildBoxContent(tabId, sample);
            }}
        }}
    }});
    var w = window.open('', '_blank');
    var style = document.querySelector('style');
    w.document.write('<html><head><title>Print</title>');
    if (style) w.document.write('<style>' + style.textContent + '</style>');
    w.document.write('</head><body>' + content.innerHTML + '</body></html>');
    w.document.close();
    w.print();
}}

function exportTableToCSV(tableId, filename) {{
    var table = document.getElementById(tableId);
    if (!table) return;
    var rows = table.querySelectorAll('tr');
    var csv = [];
    for (var i = 0; i < rows.length; i++) {{
        var row = [], cols = rows[i].querySelectorAll('td, th');
        for (var j = 0; j < cols.length; j++) {{
            row.push('"' + (cols[j].innerText || '').replace(/"/g, '""') + '"');
        }}
        csv.push(row.join(','));
    }}
    var blob = new Blob([csv.join('\\n')], {{ type: 'text/csv' }});
    var a = document.createElement('a');
    a.download = filename; a.href = URL.createObjectURL(blob);
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
}}

document.addEventListener('DOMContentLoaded', function() {{
    var hash = window.location.hash.substring(1);
    var target = hash ? document.querySelector('.tab-button.' + hash) : document.querySelector('.tab-button');
    if (target) target.click();
    document.querySelectorAll('.data-table').forEach(function(table) {{
        var headers = table.querySelectorAll('th');
        headers.forEach(function(h, idx) {{
            var type = h.getAttribute('data-sort') || 'string';
            h.style.cursor = 'pointer';
            h.addEventListener('click', function() {{ sortTable(table.id, idx, type); }});
            var icon = document.createElement('span');
            icon.className = 'sort-icon'; icon.innerHTML = '⇅';
            h.appendChild(icon);
        }});
    }});
    document.querySelectorAll('.copy-btn').forEach(function(b) {{
        b.addEventListener('click', function() {{
            var c = this.getAttribute('data-citation') || '';
            navigator.clipboard.writeText(c).then(() => {{
                var t = this.innerHTML;
                this.innerHTML = '✓ Copied!';
                setTimeout(() => {{ this.innerHTML = t; }}, 2000);
            }});
        }});
    }});
}});
</script>"""

    # ------------------------------------------------------------------
    # Main assembly
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
            ('amr', 'AMR', 'biohazard'),
            ('virulence', 'Virulence', 'virus'),
            ('plasmid', 'Plasmids', 'plug'),
            ('bacmet', 'Bacmet2', 'flask'),
            ('mutation', 'Mutations', 'dna'),
            ('mge', 'MGE', 'mobile-alt'),
            ('citation', 'Citation', 'book'),
            ('guide', 'Guide', 'question-circle'),
            ('calltoaction', 'Call to Action', 'globe'),
            ('export', 'Export', 'download'),
        ]

        nav = ''
        for i, (tid, title, icon) in enumerate(TABS):
            active = ' active' if i == 0 else ''
            nav += (f'<button class="tab-button {tid}{active}" onclick="switchTab(\'{tid}\')">'
                    f'<i class="fas fa-{icon}"></i> {title}</button>')

        methods = {
            'summary': self._sec_summary,
            'samples': self._sec_sample_overview,
            'mge': self._sec_mge,
            'citation': self._sec_citation,
            'guide': self._sec_guide,
            'calltoaction': self._sec_calltoaction,
            'export': self._sec_export,
        }

        tabs_html = ''
        for i, (tid, title, icon) in enumerate(TABS):
            active = ' active' if i == 0 else ''
            color = self.tab_colors.get(tid, '#8B4513')
            if tid == 'amr':
                content = self._gen_sample_boxes(
                    integrated_data, 'amr', 'AMR',
                    ['amrfinder', 'resfinder', 'card', 'argannot', 'megares', 'ncbi'])
            elif tid == 'virulence':
                content = self._gen_sample_boxes(
                    integrated_data, 'virulence', 'Virulence',
                    ['vfdb', 'ecoli_vf'])
            elif tid == 'plasmid':
                content = self._gen_sample_boxes(
                    integrated_data, 'plasmid', 'Plasmids', ['plasmidfinder'])
            elif tid == 'bacmet':
                content = self._gen_sample_boxes(
                    integrated_data, 'bacmet', 'Bacmet2', ['bacmet2'])
            elif tid == 'mutation':
                content = self._gen_mutation_boxes(integrated_data)
            else:
                content = methods[tid](integrated_data)

            tabs_html += f'''
            <div id="{tid}-tab" class="tab-content{active}">
                <h2 class="section-header" style="border-color:{color};">
                    <span><i class="fas fa-{icon}" style="color:{color};"></i> {title}</span>
                    <button class="print-section-btn" onclick="printSection('{tid}-tab')">
                        <i class="fas fa-print"></i> Print</button>
                </h2>
                {content}
            </div>'''

        total = len(samples_data)
        n_amr = sum(len(g) for g in gene_centric.get('amr_databases', {}).values())
        n_vir = sum(len(g) for g in gene_centric.get('virulence_databases', {}).values())
        n_mut = sum(len(v) for v in integrated_data.get('mutation_details', {}).values())
        mge_hits_total = sum(len(h) for h in integrated_data.get('mge_details', {}).values())

        dash = f'''
        <div class="dashboard-grid">
            <div class="dashboard-card" onclick="switchTab('samples')"><div class="card-number">{total}</div><div class="card-label">Total Samples</div></div>
            <div class="dashboard-card" onclick="switchTab('amr')"><div class="card-number">{n_amr}</div><div class="card-label">AMR Genes</div></div>
            <div class="dashboard-card" onclick="switchTab('virulence')"><div class="card-number">{n_vir}</div><div class="card-label">Virulence Genes</div></div>
            <div class="dashboard-card" onclick="switchTab('mutation')"><div class="card-number">{n_mut}</div><div class="card-label">Mutations</div></div>
            <div class="dashboard-card" onclick="switchTab('mge')"><div class="card-number">{mge_hits_total}</div><div class="card-label">MGE Hits</div></div>
        </div>'''

        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PseudoScope Ultimate P. aeruginosa Sample-Centric Report v3.0.0</title>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    {css}
    {js}
</head>
<body>
<div class="container">
    <div class="main-header">
        <h1><i class="fas fa-dna"></i> PseudoScope Ultimate P. aeruginosa Analysis Report</h1>
        <p>Sample-Centric Isolate Boxes with Lazy Loading</p>
        <div class="metadata-bar">
            <div class="metadata-item"><i class="fas fa-calendar"></i><span>Generated: {metadata.get('analysis_date','Unknown')}</span></div>
            <div class="metadata-item"><i class="fas fa-database"></i><span>Samples: {total}</span></div>
            <div class="metadata-item"><i class="fas fa-code-branch"></i><span>Tool: PseudoScope Ultimate v3.0.0</span></div>
            <div class="metadata-item"><i class="fas fa-university"></i><span>University of Ghana Medical School</span></div>
        </div>
    </div>
    {dash}
    <div class="tab-navigation">{nav}</div>
    {tabs_html}
    <div class="footer">
<p><strong>PSEUDOSCOPE</strong> - Sample-Centric Interactive Reports</p>
<p class="timestamp">Generated on {metadata.get('analysis_date','Unknown')}</p>
<div class="authorship">
<p><strong>Technical Support &amp; Inquiries:</strong></p>
<p>Maintainer: Brown Beckley</p>
<p>GitHub: <a href="https://github.com/bbeckley-hub" style="color:#E8C39E;">bbeckley-hub</a></p>
<p>Email: <a href="mailto:brownbeckley94@gmail.com" style="color:#E8C39E;">brownbeckley94@gmail.com</a></p>
<p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
<p style="margin-top:10px;font-size:11px;opacity:0.8;">
Powered by the PseudoScope sample-centric integration engine
</p>
</div>
</div>
</div>
</body>
</html>'''
        output_file = output_dir / "pseudoscope_ultimate_sample_centric_report.html"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"    ✅ HTML report saved: {output_file}")
        return str(output_file)

    # =========================================================================
    # SECTIONS
    # =========================================================================
    def _sec_summary(self, d):
        samples = d['samples']
        P = d['patterns']
        gc = d['gene_centric']
        total = len(samples)
        total_amr = sum(len(g) for g in gc.get('amr_databases', {}).values())
        total_vir = sum(len(g) for g in gc.get('virulence_databases', {}).values())
        total_plasmid = sum(len(g) for g in gc.get('plasmid_databases', {}).values())
        total_bacmet = sum(len(g) for g in gc.get('bacmet_databases', {}).values())
        mutation_details = d.get('mutation_details', {})
        n_mut = sum(len(v) for v in mutation_details.values())
        mge_details = d.get('mge_details', {})
        n_mge_samples = len(mge_details)
        mge_hits_total = sum(len(v) for v in mge_details.values())

        return f'''
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>📊 Analysis Overview</h3>
                <p>This report analyses <strong>{total}</strong> <em>Pseudomonas aeruginosa</em> genomes.</p>
                <ul>
                    <li><strong>Sample Overview</strong> — every isolate on one row with typing summary.</li>
                    <li><strong>Sample-centric, lazy-loaded</strong> gene tabs (AMR, Virulence, Plasmids, Bacmet2, Mutations, MGE) — each isolate in its own collapsed box; tables built on click. <strong>Browser-safe on 1000+ samples.</strong></li>
                </ul>
            </div>
        </div>
        <h3>Key Statistics</h3>
        <div class="scrollable-table">
            <table class="data-table">
                <thead><tr><th>Metric</th><th>Count</th><th>Notes</th></tr></thead>
                <tbody>
                    <tr><td>Total Samples</td><td><strong>{total}</strong></td><td>Complete genomic analysis</td></tr>
                    <tr><td>Unique STs</td><td><strong>{len(P.get('mlst_distribution', {}))}</strong></td><td>MLST Pasteur scheme</td></tr>
                    <tr><td>Unique O-Types</td><td><strong>{len(P.get('o_type_distribution', {}))}</strong></td><td>PAST serotyping</td></tr>
                    <tr><td>AMR genes (unique)</td><td><strong>{total_amr}</strong></td><td>All AMR databases</td></tr>
                    <tr><td>Virulence genes (unique)</td><td><strong>{total_vir}</strong></td><td>VFDB</td></tr>
                    <tr><td>Plasmid replicons</td><td><strong>{total_plasmid}</strong></td><td>PlasmidFinder</td></tr>
                    <tr><td>Bacmet2 genes</td><td><strong>{total_bacmet}</strong></td><td>Biocide &amp; metal resistance</td></tr>
                    <tr><td>Total mutation events</td><td><strong>{n_mut}</strong></td><td>AMRFinderPlus point mutations</td></tr>
                    <tr><td>MGE-profiled samples</td><td><strong>{n_mge_samples}</strong></td><td>mobileOG-db hits: {mge_hits_total}</td></tr>
                </tbody>
            </table>
        </div>
        <h3 style="margin-top:24px;">📘 How to Use the Isolate Boxes</h3>
        <ol style="margin-left:20px;">
            <li>Open any gene-centric tab (AMR, Virulence, Plasmids, Bacmet2, Mutations, MGE).</li>
            <li>Each isolate shows a header (sample name, count, typing badges).</li>
            <li>Click <strong>Show Details</strong> inside a box to load its tables on demand.</li>
            <li>Filter by sample name (search box) or by database (dropdown).</li>
            <li>Use <strong>Expand All (visible)</strong> for batched, non-blocking expansion.</li>
        </ol>'''

    def _sec_sample_overview(self, d):
        samples = d['samples']
        rows = ''
        for sample, data in sorted(samples.items()):
            t = data.get('typing', {})
            mlst = t.get('MLST', 'ND')
            otype = t.get('O_Type', 'ND')
            t3ss = t.get('T3SS', 'ND')
            vir_count = 0
            for db in ('vfdb', 'ecoli_vf'):
                vir_count += len(data.get('abricate_databases', {}).get(db, []))
            t3ss_class = ('t3ss-exou' if t3ss == 'ExoU'
                          else 't3ss-exos' if t3ss == 'ExoS'
                          else 't3ss-both' if t3ss == 'ExoU+ExoS'
                          else 't3ss-neither')
            rows += f'''<tr>
                <td><strong>{esc(sample)}</strong></td>
                <td>{esc(mlst)}</td>
                <td>{esc(otype)}</td>
                <td><span class="typing-badge {t3ss_class}">{esc(t3ss)}</span></td>
                <td>{vir_count}</td>
            </tr>'''

        return f'''
        {self._alert('info', 'fa-list-alt',
            '<h3>🧬 Sample Overview – Population Snapshot</h3>'
            '<p>Master reference table showing MLST, O-type, T3SS pathotype, and virulence gene count for every isolate.</p>')}
        <input type="text" class="search-box" id="search-samples-table"
               onkeyup="searchTable('samples-table','search-samples-table')"
               placeholder="🔍 Search samples...">
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('samples-table','sample_overview.csv')">
                <i class="fas fa-download"></i> Export CSV</button>
        </div>
        <div class="scrollable-table">
            <table id="samples-table" class="data-table">
                <thead><tr>
                    <th data-sort="string">Sample</th>
                    <th data-sort="string">MLST</th>
                    <th data-sort="string">O-Type</th>
                    <th data-sort="string">T3SS Pathotype</th>
                    <th data-sort="number">Virulence Count</th>
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    def _gen_sample_boxes(self, d, tab_id, title, db_list):
        samples = d.get('samples', {})
        amr_details = d.get('amrfinder_details', {})
        abricate_details = d.get('abricate_details', {})

        relevant = []
        for sample in samples:
            has = False
            if 'amrfinder' in db_list and amr_details.get(sample):
                has = True
            for db in db_list:
                if db != 'amrfinder' and abricate_details.get(sample, {}).get(db):
                    has = True
            if has:
                relevant.append(sample)
        relevant.sort()

        if not relevant:
            return self._alert('warning', 'fa-exclamation-circle',
                f'<h3>No {title} Data Available</h3><p>No samples with {title} genes were found.</p>')

        box_data = {}
        for sample in relevant:
            entry = {'amrfinder': [], 'abricate': {}}
            if 'amrfinder' in db_list:
                entry['amrfinder'] = amr_details.get(sample, [])
            for db in db_list:
                if db != 'amrfinder':
                    entry['abricate'][db] = abricate_details.get(sample, {}).get(db, [])
            box_data[sample] = entry
        box_json = json.dumps(box_data, default=str, ensure_ascii=False)

        db_options = '<option value="all">All Databases</option>'
        for db in db_list:
            disp = 'AMRfinder' if db == 'amrfinder' else db.upper()
            db_options += f'<option value="{db}">{disp}</option>'

        html = f'''
        {self._alert('info', 'fa-info-circle',
            f'<h3>🧬 {title} – Interactive Isolate Boxes (Lazy-Loaded)</h3>'
            '<p>Each box shows one isolate. Click <strong>Show Details</strong> to load its tables on demand — the report opens instantly even with thousands of samples.</p>')}
        <div class="filter-controls">
            <input type="text" class="search-box" id="search-{tab_id}"
                   onkeyup="filterBoxes('{tab_id}')"
                   placeholder="🔍 Search sample..." style="max-width:320px;margin-bottom:0;">
            <select id="dbFilter-{tab_id}" onchange="filterBoxes('{tab_id}')">{db_options}</select>
            <button class="action-btn btn-info" onclick="expandAllBoxes('{tab_id}')">
                <i class="fas fa-expand-alt"></i> Expand All (visible)</button>
            <button class="action-btn btn-light" onclick="collapseAllBoxes('{tab_id}')">
                <i class="fas fa-compress-alt"></i> Collapse All</button>
            <button class="action-btn btn-success" onclick="resetBoxFilters('{tab_id}')">
                <i class="fas fa-sync"></i> Clear</button>
            <span class="results-counter" id="counter-{tab_id}">{len(relevant)} shown</span>
        </div>
        <div id="box-container-{tab_id}">
        '''

        for sample in relevant:
            t = samples.get(sample, {}).get('typing', {})
            mlst = t.get('MLST', 'ND')
            otype = t.get('O_Type', 'ND')
            t3ss = t.get('T3SS', 'ND')
            total_genes = len(amr_details.get(sample, [])) if 'amrfinder' in db_list else 0
            for db in db_list:
                if db != 'amrfinder':
                    total_genes += len(abricate_details.get(sample, {}).get(db, []))
            t3ss_cls = ('t3ss-exou' if t3ss == 'ExoU'
                        else 't3ss-exos' if t3ss == 'ExoS'
                        else 't3ss-both' if t3ss == 'ExoU+ExoS'
                        else 't3ss-neither')
            html += f'''
            <div class="isolate-box" data-sample="{esc(sample)}" data-tab="{tab_id}">
                <div class="sample-header">
                    <h3><i class="fas fa-microbe"></i> {esc(sample)}</h3>
                    <span class="total-badge">Total Genes: {total_genes}</span>
                    <div class="typing-info">
                        <span class="typing-badge">ST: {esc(mlst)}</span>
                        <span class="typing-badge">O-Type: {esc(otype)}</span>
                        <span class="typing-badge {t3ss_cls}">{esc(t3ss)}</span>
                    </div>
                    <button class="toggle-btn" onclick="toggleBoxDetails('{tab_id}', '{esc(sample)}')">
                        <i class="fas fa-chevron-down"></i> Show Details
                    </button>
                </div>
                <div class="box-details" style="display:none;"></div>
            </div>'''

        html += '</div>'
        html += f'''
        <script>
        window.PSEUDOSCOPE_BOX_DATA = window.PSEUDOSCOPE_BOX_DATA || {{}};
        window.PSEUDOSCOPE_BOX_DATA["{tab_id}"] = {box_json};
        </script>'''
        return html

    def _gen_mutation_boxes(self, d):
        samples = d.get('samples', {})
        mutation_details = d.get('mutation_details', {})
        relevant = [s for s in samples if mutation_details.get(s)]
        relevant.sort()

        if not relevant:
            return self._alert('warning', 'fa-dna',
                '<h3>No Mutation Data Available</h3><p>No point mutations were found.</p>')

        box_data = {s: {'mutations': mutation_details.get(s, [])} for s in relevant}
        box_json = json.dumps(box_data, default=str, ensure_ascii=False)

        credit = self._credit_bar('#00BCD4', '🧬', 'AMRFinderPlus – Point Mutations',
            'AMRFinderPlus by NCBI — comprehensive database for resistance genes and point mutations.')

        info = self._gene_family_info('#00BCD4',
            'Clinical relevance of key mutations in P. aeruginosa:', [
                ('gyrA / parC', 'QRDR mutations — fluoroquinolone resistance.'),
                ('rpoB', 'Rifampin resistance.'),
                ('23S rRNA', 'Macrolide and linezolid resistance.'),
                ('oprD', 'Porin loss — carbapenem resistance mechanism independent of carbapenemases.'),
                ('mexR / nalC / nalD', 'Efflux pump regulators — MexAB-OprM overexpression.'),
                ('ampR / ampD', 'AmpC derepression.'),
            ])

        html = f'''
        {credit}
        {self._alert('info', 'fa-dna',
            '<h3>🧬 Point Mutations – Sample-Centric View (Lazy-Loaded)</h3>'
            '<p>Each box shows one isolate. Click <strong>Show Details</strong> to load its mutation table.</p>')}
        {info}
        <div class="filter-controls">
            <input type="text" class="search-box" id="search-mutation"
                   onkeyup="filterBoxes('mutation')"
                   placeholder="🔍 Search sample..." style="max-width:320px;margin-bottom:0;">
            <button class="action-btn btn-info" onclick="expandAllBoxes('mutation')">
                <i class="fas fa-expand-alt"></i> Expand All (visible)</button>
            <button class="action-btn btn-light" onclick="collapseAllBoxes('mutation')">
                <i class="fas fa-compress-alt"></i> Collapse All</button>
            <button class="action-btn btn-success" onclick="resetBoxFilters('mutation')">
                <i class="fas fa-sync"></i> Clear</button>
            <span class="results-counter" id="counter-mutation">{len(relevant)} shown</span>
        </div>
        <div id="box-container-mutation">
        '''
        for sample in relevant:
            t = samples.get(sample, {}).get('typing', {})
            mlst = t.get('MLST', 'ND')
            otype = t.get('O_Type', 'ND')
            t3ss = t.get('T3SS', 'ND')
            n_mut = len(mutation_details.get(sample, []))
            t3ss_cls = ('t3ss-exou' if t3ss == 'ExoU'
                        else 't3ss-exos' if t3ss == 'ExoS'
                        else 't3ss-both' if t3ss == 'ExoU+ExoS'
                        else 't3ss-neither')
            html += f'''
            <div class="isolate-box" data-sample="{esc(sample)}" data-tab="mutation">
                <div class="sample-header">
                    <h3><i class="fas fa-microbe"></i> {esc(sample)}</h3>
                    <span class="total-badge">Total Mutations: {n_mut}</span>
                    <div class="typing-info">
                        <span class="typing-badge">ST: {esc(mlst)}</span>
                        <span class="typing-badge">O-Type: {esc(otype)}</span>
                        <span class="typing-badge {t3ss_cls}">{esc(t3ss)}</span>
                    </div>
                    <button class="toggle-btn" onclick="toggleBoxDetails('mutation', '{esc(sample)}')">
                        <i class="fas fa-chevron-down"></i> Show Details
                    </button>
                </div>
                <div class="box-details" style="display:none;"></div>
            </div>'''
        html += '</div>'
        html += f'''
        <script>
        window.PSEUDOSCOPE_BOX_DATA = window.PSEUDOSCOPE_BOX_DATA || {{}};
        window.PSEUDOSCOPE_BOX_DATA["mutation"] = {box_json};
        </script>'''
        return html

    def _sec_mge(self, d):
        """MGE tab — sample-centric lazy-loaded isolate boxes."""
        samples_data = d.get('samples', {})
        mge_details = d.get('mge_details', {})
        relevant = [s for s in samples_data if mge_details.get(s)]
        relevant.sort()

        if not relevant:
            return self._alert('warning', 'fa-mobile-alt',
                '<h3>No MGE Data Available</h3>'
                '<p>pseudoscope_mge_hits.tsv was not found or contains no rows.</p>')

        blob = {s: {'mge': mge_details.get(s, [])} for s in relevant}
        blob_json = json.dumps(blob, default=str, ensure_ascii=False)

        credit = self._credit_bar('#16A085', '🧬', 'Mobile Genetic Elements – mobileOG-db',
            '<strong>mobileOG-db</strong> by '
            '<a href="https://github.com/clb21565/mobileOG-db" target="_blank" '
            'style="color:#16A085;font-weight:bold;">Brown, Mullet and colleagues</a>.<br>'
            'Please cite: Brown CL, et al. mobileOG-db: a Manually Curated Database of '
            'Protein Families Mediating the Life Cycle of Bacterial Mobile Genetic Elements. '
            '<em>Appl Environ Microbiol</em>. 2022;88(18):e00991-22. '
            '<a href="https://doi.org/10.1128/aem.00991-22" target="_blank" '
            'style="color:#16A085;font-weight:bold;">doi:10.1128/aem.00991-22</a>')

        info = self._gene_family_info('#16A085',
            'MGE categories and why they matter:', [
                ('replication/recombination/repair', 'Core machinery for element replication and host DNA integration.'),
                ('transfer', 'Conjugal transfer systems (T4SS) — primary route for horizontal gene transfer.'),
                ('integration/excision', 'Integrases and transposases moving genetic material between replicons.'),
                ('stability/defense', 'Partitioning and toxin-antitoxin systems ensuring persistence.'),
                ('phage', 'Bacteriophage structural and lysogeny proteins.'),
                ('plasmid', 'Plasmid backbone proteins — important vehicle for acquired resistance.'),
                ('IS', 'Insertion sequence transposases.'),
                ('ICE', 'Integrative conjugative elements moving large gene cassettes.'),
            ])

        html = f'''
        {credit}
        {self._alert('info', 'fa-mobile-alt',
            '<h3>🧬 Mobile Genetic Elements – Sample-Centric View (Lazy-Loaded)</h3>'
            '<p>mobileOG-db protein family hits for each isolate. Click <strong>Show Details</strong> '
            'on any box to load its MGE table on demand.</p>')}
        <div class="filter-controls">
            <input type="text" class="search-box" id="search-mge"
                   onkeyup="filterBoxes('mge')"
                   placeholder="🔍 Search sample by name..." style="max-width:320px;margin-bottom:0;">
            <button class="action-btn btn-info" onclick="expandAllBoxes('mge')">
                <i class="fas fa-expand-alt"></i> Expand All (visible)</button>
            <button class="action-btn btn-light" onclick="collapseAllBoxes('mge')">
                <i class="fas fa-compress-alt"></i> Collapse All</button>
            <button class="action-btn btn-success" onclick="resetBoxFilters('mge')">
                <i class="fas fa-sync"></i> Clear</button>
            <span class="results-counter" id="counter-mge">{len(relevant)} shown</span>
        </div>
        <div id="box-container-mge">
        '''

        for sample in relevant:
            sd = samples_data.get(sample, {})
            t = sd.get('typing', {})
            mlst = t.get('MLST', 'ND')
            otype = t.get('O_Type', 'ND')
            t3ss = t.get('T3SS', 'ND')
            n_hits = len(mge_details.get(sample, []))
            t3ss_cls = ('t3ss-exou' if t3ss == 'ExoU'
                        else 't3ss-exos' if t3ss == 'ExoS'
                        else 't3ss-both' if t3ss == 'ExoU+ExoS'
                        else 't3ss-neither')
            html += f'''
            <div class="isolate-box" data-sample="{esc(sample)}" data-tab="mge">
                <div class="sample-header">
                    <h3><i class="fas fa-microbe"></i> {esc(sample)}</h3>
                    <span class="total-badge">Total MGE hits: {n_hits}</span>
                    <div class="typing-info">
                        <span class="typing-badge">ST: {esc(mlst)}</span>
                        <span class="typing-badge">O-Type: {esc(otype)}</span>
                        <span class="typing-badge {t3ss_cls}">{esc(t3ss)}</span>
                    </div>
                    <button class="toggle-btn" onclick="toggleBoxDetails('mge', '{esc(sample)}')">
                        <i class="fas fa-chevron-down"></i> Show Details
                    </button>
                </div>
                <div class="box-details" style="display:none;"></div>
            </div>'''

        html += '</div>'
        html += f'''
        <script>
        window.PSEUDOSCOPE_BOX_DATA = window.PSEUDOSCOPE_BOX_DATA || {{}};
        window.PSEUDOSCOPE_BOX_DATA["mge"] = {blob_json};
        </script>'''
        html += '<h3 style="margin-top:30px;"><i class="fas fa-book-medical"></i> MGE Category Guide</h3>'
        html += '<p style="color:#666;font-size:.9em;margin-bottom:10px;">Expand below for MGE category descriptions.</p>'
        html += info
        return html

    def _sec_citation(self, d):
        colors = ['#E3F2FD','#E8F5E9','#FFF3E0','#FCE4EC','#F3E5F5',
                  '#E0F7FA','#FFF8E1','#FBE9E7','#EDE7F6','#E0F2F1',
                  '#F1F8E9','#FFFDE7','#FCE4EC','#E3F2FD','#E8EAF6']
        i = [0]

        def fmt(name, text, kind, url):
            c = colors[i[0] % len(colors)]
            i[0] += 1
            if kind == 'doi':
                link = f'<a href="https://doi.org/{url}" target="_blank" style="color:#0D6EFD;font-weight:600;">🔗 {url}</a>'
            else:
                link = f'<a href="{url}" target="_blank" style="color:#8B4513;font-weight:600;">🔗 {url}</a>'
            return (f'<div style="background:{c};border-left:5px solid #00838F;padding:14px;'
                    f'border-radius:8px;margin-bottom:12px;">'
                    f'<strong>{name}</strong> – {text} {link} '
                    f'<button class="copy-btn" data-citation="{text} {url}">📋 Copy</button></div>')

        html = self._alert('info', 'fa-quote-right',
            '<h3>📚 How to Cite PseudoScope and Its Dependencies</h3>'
            '<p>If you use this tool in your research, please cite the main tool and the relevant third-party tools and databases.</p>')

        html += fmt("PseudoScope",
            "Beckley B, Amarh V. PseudoScope: a species-specific bioinformatics suite for rapid and accessible Pseudomonas aeruginosa genomic analysis. GitHub 2026.",
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
        html += fmt("Biopython",
            "Cock PJ, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009;25(11):1422-3.",
            "doi", "10.1093/bioinformatics/btp163")

        html += '''
        <div class="alert-box alert-success" style="margin-top:20px;">
            <i class="fas fa-hand-peace"></i>
            <div><strong>Suggested acknowledgement:</strong><br>
            "Genomic analysis was performed using PseudoScope [Beckley &amp; Amarh, 2026], which integrates MLST [Seemann, 2018] using the PubMLST database [Jolley et al., 2018], ABRicate [Seemann, 2018], AMRFinderPlus [Feldgarden et al., 2021], PAST serotyping [Thomsen et al., DTU] via the pasty wrapper [Petit et al., 2023], and fastANI [Jain et al., 2018] for comprehensive P. aeruginosa characterization. Antimicrobial resistance genes were identified using the CARD [McArthur et al., 2013] and ResFinder [Florensa et al., 2022] databases. For biocide and heavy metal resistance genes, BacMet [Pal et al., 2014] was used. Virulence and plasmid screening were performed with ABRicate using the VFDB [Chen et al., 2012] and PlasmidFinder [Carattoli et al., 2014] databases. Mobile genetic element profiling used mobileOG-db [Brown et al., 2022], Prodigal [Hyatt et al., 2010], and DIAMOND [Buchfink et al., 2015]."
            </div>
        </div>'''
        return html

    def _sec_guide(self, d):
        return '''
        <div class="alert-box alert-info">
            <i class="fas fa-question-circle fa-2x"></i>
            <div><h3>📘 PseudoScope Sample-Centric Report – Guide</h3>
            <p>This guide explains how to get the most out of each tab and the lazy-loading boxes.</p></div>
        </div>
        <div class="database-section">
            <h4><i class="fas fa-cogs"></i> How to Use the Report</h4>
            <ul>
                <li><strong>Summary</strong> — overview + key statistics + how-to.</li>
                <li><strong>Sample Overview</strong> — one row per isolate with typing summary.</li>
                <li><strong>AMR</strong> — resistance genes as lazy-loaded isolate boxes.</li>
                <li><strong>Virulence</strong> — VFDB virulence genes as lazy-loaded isolate boxes.</li>
                <li><strong>Plasmids</strong> — plasmid replicons as lazy-loaded isolate boxes.</li>
                <li><strong>Bacmet2</strong> — biocide &amp; heavy-metal resistance as lazy-loaded isolate boxes.</li>
                <li><strong>Mutations</strong> — point mutations as lazy-loaded isolate boxes.</li>
                <li><strong>MGE</strong> — mobileOG-db hits as lazy-loaded isolate boxes.</li>
                <li><strong>Citation, Call to Action, Export</strong> — citations, ESCAPE AMR project, downloads.</li>
            </ul>
        </div>
        <div class="database-section" style="background:#e8f5e9;border-left:5px solid #28a745;">
            <h4><i class="fas fa-lightbulb"></i> Pro Tips</h4>
            <ul>
                <li><strong>Lazy boxes</strong> — click "Show Details" on any isolate box to load its gene/mutation/MGE tables on demand. This keeps the report responsive even with 1000+ samples.</li>
                <li><strong>Expand All (visible)</strong> — batch-opens the visible boxes at 20 per tick. Browser stays responsive.</li>
                <li><strong>Filter</strong> — search by sample name or restrict to one database via the dropdown.</li>
                <li><strong>Export</strong> — every typing table is exportable as CSV. Use the JSON file for AI-assisted analysis.</li>
                <li><strong>Print</strong> — printing a tab auto-expands all its isolate boxes first.</li>
            </ul>
        </div>'''

    def _sec_calltoaction(self, d):
        return '''
        <div class="alert-box alert-info">
            <i class="fas fa-globe fa-2x"></i>
            <div>
                <h3>The Global Burden of AMR and Our Call to Action</h3>
                <p>Antimicrobial resistance kills an estimated <strong>1.27 million people directly each year</strong>. <em>Pseudomonas aeruginosa</em> is a WHO critical-priority pathogen, especially carbapenem-resistant. Genomic surveillance is our best tool.</p>
                <p>We developed <strong>PseudoScope</strong> as part of the <strong>ESCAPE AMR</strong> project.</p>
            </div>
        </div>
        <div style="background:#e8f5e9;padding:20px;border-radius:12px;margin:20px 0;">
            <h3><i class="fas fa-bacterium"></i> ESCAPE AMR – Our Vision</h3>
            <ul>
                <li><strong>E</strong>veryone must join forces.</li>
                <li><strong>S</strong>mart surveillance is our first line of defence.</li>
                <li><strong>K</strong>nowledge must be shared openly.</li>
                <li><strong>A</strong>frica bears a heavy AMR burden — African solutions are emerging.</li>
                <li><strong>P</strong>revention is cheaper than cure.</li>
                <li><strong>E</strong>very day we delay, more lives are at stake.</li>
            </ul>
        </div>
        <div style="text-align:center;margin:40px 0;">
            <i class="fas fa-star" style="font-size:3em;color:#ffc107;"></i>
            <h3>🤝 We Invite You to Contribute!</h3>
            <div class="action-buttons" style="justify-content:center;">
                <a href="https://github.com/bbeckley-hub/pseudoscope" target="_blank" class="action-btn btn-primary" style="text-decoration:none;"><i class="fab fa-github"></i> ⭐ Star us on GitHub</a>
                <a href="mailto:brownbeckley94@gmail.com" class="action-btn btn-success" style="text-decoration:none;"><i class="fas fa-envelope"></i> Contact us</a>
                <a href="https://github.com/bbeckley-hub/pseudoscope/issues" target="_blank" class="action-btn btn-warning" style="text-decoration:none;"><i class="fas fa-bug"></i> Report issues</a>
            </div>
            <p style="margin-top:20px;"><i class="fas fa-chalkboard-user"></i> <strong>We welcome collaborations</strong> to adapt this tool for other pathogens.</p>
        </div>'''

    def _sec_export(self, d):
        cards = ''
        for tid, label, icon, fname in [
            ('samples-table', 'Sample Overview CSV', 'fa-table', 'sample_overview.csv'),
        ]:
            cards += f'''
            <div class="dashboard-card" onclick="exportTableToCSV('{tid}', '{fname}')">
                <div style="font-size:2em;color:#37474F;"><i class="fas {icon}"></i></div>
                <div class="card-label">{label}</div>
            </div>'''
        cards += '''
        <div class="dashboard-card" onclick="location.href='pseudoscope_ultimate_sample_centric_report.json'">
            <div style="font-size:2em;color:#37474F;"><i class="fas fa-file-code"></i></div>
            <div class="card-label">Complete JSON Data</div>
        </div>'''
        return f'''
        {self._alert('info', 'fa-download',
            '<h3>📥 Export Data</h3>'
            '<p>Typing tables export to CSV directly. Lazy-loaded gene/mutation/MGE tables export via the JSON file, which also powers AI-assisted analysis.</p>')}
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:20px;margin-top:20px;">
            {cards}
        </div>'''


# =============================================================================
# ORCHESTRATOR
# =============================================================================
class PseudomonasUltimateReporter:
    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "PSEUDOSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS"
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

    def find_files(self):
        print("🔍 Searching for PseudoScope reports...")
        found = {'mlst_html': None, 'past_html': None, 'qc_html': None}
        for f in self.input_dir.glob('**/*.html'):
            n = f.name.lower()
            if 'mlst_summary' in n and not found['mlst_html']:
                found['mlst_html'] = f
                print(f"    ✓ MLST: {f.name}")
            elif 'past_summary' in n and not found['past_html']:
                found['past_html'] = f
                print(f"    ✓ PAST: {f.name}")
            elif ('fasta_qc_summary' in n or 'fastaqc_summary' in n) and not found['qc_html']:
                found['qc_html'] = f
                print(f"    ✓ FASTA QC: {f.name}")
        return found

    def integrate(self, files):
        print("\n🔗 Integrating data...")
        data = {
            'metadata': self.metadata,
            'samples': {},
            'patterns': {},
            'gene_centric': {},
            'qc_data': {},
            'amrfinder_details': {},
            'abricate_details': {},
            'mutation_details': {},
            'mge_details': {},
        }
        mlst_data = self.parser.parse_mlst_report(files['mlst_html']) if files['mlst_html'] else {}
        past_data = self.parser.parse_past_report(files['past_html']) if files['past_html'] else {}
        if files['qc_html']:
            data['qc_data'] = self.parser.parse_qc_report(files['qc_html'])

        amr_details, amr_freq = self.parser.load_amrfinder_from_tsv(self.input_dir)
        abricate_details, abricate_freq = self.parser.load_abricate_from_tsv(self.input_dir)
        mutation_details = self.parser.load_mutations_from_tsv(self.input_dir)
        mge_details = self.parser.load_mge_from_tsv(self.input_dir)

        data['amrfinder_details'] = amr_details
        data['abricate_details'] = abricate_details
        data['mutation_details'] = mutation_details
        data['mge_details'] = mge_details

        all_samples = (set(mlst_data) | set(past_data) | set(amr_details)
                       | set(abricate_details) | set(data['qc_data']) | set(mutation_details)
                       | set(mge_details))
        all_samples = sorted(all_samples)
        if not all_samples:
            print("❌ No samples found in any report")
            return {}
        print(f"📊 Found {len(all_samples)} unique samples")

        for sample in all_samples:
            amr_genes = [g.get('gene', '') for g in amr_details.get(sample, [])
                         if g.get('gene')]
            vir_genes = []
            for vir_db in ('vfdb', 'ecoli_vf', 'pa_vf', 'pseudomonas_vf'):
                for g in abricate_details.get(sample, {}).get(vir_db, []):
                    if g.get('gene'):
                        vir_genes.append(g['gene'])
            all_genes_for_t3ss = list(set(amr_genes + vir_genes))
            t3ss = self.analyzer.derive_t3ss_pathotype(all_genes_for_t3ss)

            data['samples'][sample] = {
                'typing': {
                    'MLST':   mlst_data.get(sample, {}).get('ST', 'ND'),
                    'O_Type': past_data.get(sample, {}).get('O_Type', 'ND'),
                    'T3SS':   t3ss,
                },
                'amrfinder': {'all_genes': amr_genes},
                'abricate_databases': {
                    db: [g.get('gene', '') for g in gs if g.get('gene')]
                    for db, gs in abricate_details.get(sample, {}).items()
                },
            }

        amr_gene_freq = {}
        for sample, genes in amr_details.items():
            for g in genes:
                gene = g.get('gene', '')
                if not gene:
                    continue
                if gene not in amr_gene_freq:
                    amr_gene_freq[gene] = {'count': 0, 'genomes': []}
                if sample not in amr_gene_freq[gene]['genomes']:
                    amr_gene_freq[gene]['genomes'].append(sample)
                    amr_gene_freq[gene]['count'] += 1

        abricate_gene_freq = {}
        for sample, dbs in abricate_details.items():
            for db, gs in dbs.items():
                if db not in abricate_gene_freq:
                    abricate_gene_freq[db] = {}
                for g in gs:
                    gene = g.get('gene', '')
                    if not gene:
                        continue
                    if gene not in abricate_gene_freq[db]:
                        abricate_gene_freq[db][gene] = {'count': 0, 'genomes': []}
                    if sample not in abricate_gene_freq[db][gene]['genomes']:
                        abricate_gene_freq[db][gene]['genomes'].append(sample)
                        abricate_gene_freq[db][gene]['count'] += 1

        data['gene_frequencies'] = {'amrfinder': amr_gene_freq, 'abricate': abricate_gene_freq}
        data['gene_centric'] = self.analyzer.create_gene_centric_tables(data)
        data['patterns'] = self.analyzer.create_cross_genome_patterns(data)
        return data

    def write_json(self, data):
        print("\n📝 Writing JSON report...")
        out = self.output_dir / "pseudoscope_ultimate_sample_centric_report.json"

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

    def write_csvs(self, data):
        print("\n📊 Writing CSV reports...")
        rows = []
        for sample, s in data['samples'].items():
            t = s['typing']
            vir_count = 0
            for db in ('vfdb', 'ecoli_vf'):
                vir_count += len(s.get('abricate_databases', {}).get(db, []))
            rows.append({
                'Sample': sample,
                'MLST': t['MLST'],
                'O_Type': t['O_Type'],
                'T3SS_Pathotype': t['T3SS'],
                'Virulence_Count': vir_count,
            })
        pd.DataFrame(rows).to_csv(self.output_dir / "sample_overview.csv", index=False)

        gc = data.get('gene_centric', {})
        total = len(data['samples']) or 1
        for cat, fname in [('amr_databases', 'amr_genes.csv'),
                           ('virulence_databases', 'virulence_genes.csv'),
                           ('bacmet_databases', 'bacmet_genes.csv'),
                           ('plasmid_databases', 'plasmid_replicons.csv')]:
            out = []
            for db, gs in gc.get(cat, {}).items():
                for g in gs:
                    out.append({'Gene': g['gene'], 'Database': g['database'],
                                'Count': g['count'],
                                'Percentage': f"{(g['count'] / total) * 100:.1f}%",
                                'Genomes': ';'.join(g.get('genomes', []))})
            if out:
                pd.DataFrame(out).to_csv(self.output_dir / fname, index=False)

        muts = data.get('mutation_details', {})
        if muts:
            rows = []
            for sample, mlist in muts.items():
                for m in mlist:
                    rows.append({
                        'Sample': sample,
                        'Gene': m.get('gene', ''),
                        'Mutation': m.get('mutation', ''),
                        'Class': m.get('class', ''),
                        'Contig': m.get('contig', ''),
                        'Start': m.get('start', ''),
                        'Stop': m.get('stop', ''),
                        'Strand': m.get('strand', ''),
                        'Coverage': m.get('coverage', ''),
                        'Identity': m.get('identity', ''),
                        'Accession': m.get('accession', ''),
                    })
            pd.DataFrame(rows).to_csv(self.output_dir / "mutations.csv", index=False)

        mge = data.get('mge_details', {})
        if mge:
            rows = []
            for sample, hits in mge.items():
                for h in hits:
                    rows.append({'Sample': sample, **h})
            if rows:
                pd.DataFrame(rows).to_csv(self.output_dir / "mge_hits.csv", index=False)

        carb_rows = []
        for sample, cats in data['patterns'].get('carbapenemase_by_sample', {}).items():
            for cat, fams in cats.items():
                for fam, raw_map in fams.items():
                    for raw, dbs in raw_map.items():
                        carb_rows.append({
                            'Sample': sample, 'Category': cat, 'Family': fam,
                            'Raw_Gene': raw, 'Databases': ';'.join(dbs),
                            'Database_Count': len(dbs),
                        })
        if carb_rows:
            pd.DataFrame(carb_rows).to_csv(
                self.output_dir / "carbapenemase_classification.csv", index=False)

        print("    ✅ CSV reports written")

    def run(self):
        print("=" * 80)
        print("🧬 PSEUDOSCOPE ULTIMATE P. AERUGINOSA SAMPLE-CENTRIC REPORTER v3.0.0")
        print("=" * 80)
        print(f"📁 Input:  {self.input_dir}")
        print(f"📁 Output: {self.output_dir}")
        files = self.find_files()
        if not any(files.values()):
            print("⚠️ No MLST/PAST/QC HTML found — will rely on TSVs only")
        data = self.integrate(files)
        if not data:
            return False
        self.write_json(data)
        self.write_csvs(data)
        self.generator.generate_main_report(data, self.output_dir)
        n = len(data['samples'])
        n_mut = len(data.get('mutation_details', {}))
        n_mge = len(data.get('mge_details', {}))
        mge_hits = sum(len(v) for v in data.get('mge_details', {}).values())
        print("\n" + "=" * 80)
        print("✅ REPORT COMPLETE")
        print("=" * 80)
        print(f"   Samples:               {n}")
        print(f"   Samples with mutations: {n_mut}")
        print(f"   Samples with MGE hits: {n_mge} ({mge_hits} total hits)")
        print(f"   Output dir:            {self.output_dir}")
        print("=" * 80)
        return True


def main():
    p = argparse.ArgumentParser(
        description='PseudoScope Ultimate P. aeruginosa Sample-Centric Reporter v3.0.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python pseudoscope_sample_centric.py -i /path/to/reports

Author: Brown Beckley <brownbeckley94@gmail.com>""")
    p.add_argument('-i', '--input-dir', required=True)
    p.add_argument('-o', '--output-dir')
    args = p.parse_args()
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