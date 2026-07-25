#!/usr/bin/env python3
"""
GENIUS PSEUDOMONAS AERUGINOSA SAMPLE-CENTRIC REPORTER (v1.0.0)
================================================================================
Hybrid reporter for P. aeruginosa genomic data.
- MLST and serotype parsing from HTML (only)
- AMR, Virulence, Plasmid, Bacmet, and Mutation data loaded strictly from TSV summaries
  with prefix "pseudo_"
- Gene‑centric tables with dynamic grouping by ST, O‑type, or ST‑O
- Sample‑centric interactive boxes for each isolate (AMR, Virulence, BACMET, Plasmids, Mutations)
- Educational content, filter buttons, and export capabilities
- Corrected citations for PAST serotyping and pasty wrapper

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School, Department of Medical Biochemistry
Version: 1.0.0
Date: 2026-07-24
License: MIT
================================================================================
"""

import os
import sys
import json
import re
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


class PseudomonasHTMLParser:
    """
    Parses HTML reports from PseudoScope pipeline for MLST and serotype.
    Loads TSV summaries for AMRfinder, ABRicate, and mutations.
    TSV files are expected to have the prefix "pseudo_".
    """

    def __init__(self):
        self.abricate_databases = [
            'ncbi', 'card', 'resfinder', 'vfdb', 'argannot',
            'plasmidfinder', 'megares', 'ecoli_vf', 'bacmet2'
        ]

    def normalize_sample_id(self, sample_id: str) -> str:
        sample = str(sample_id)
        for ext in ['.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv']:
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

    def load_amrfinder_from_tsv(self, input_dir: Path) -> Tuple[Dict[str, List[Dict]], Dict[str, int]]:
        tsv_path = input_dir / 'pseudo_amrfinder_summary.tsv'
        if not tsv_path.exists():
            print(f"  ⚠️ AMRfinder TSV not found: {tsv_path}")
            return {}, {}
        df = pd.read_csv(tsv_path, sep='\t')
        amr_details = defaultdict(list)
        gene_counts = Counter()
        sample_col = 'Genome' if 'Genome' in df.columns else 'sample'
        if sample_col not in df.columns:
            sample_col = df.columns[0]
        for _, row in df.iterrows():
            sample = str(row[sample_col])
            sample = self.normalize_sample_id(sample)
            gene_dict = row.to_dict()
            gene_dict = {k: (None if pd.isna(v) else v) for k, v in gene_dict.items()}
            if 'Gene_Symbol' in gene_dict:
                gene_dict['gene'] = gene_dict.pop('Gene_Symbol')
            elif 'Element_symbol' in gene_dict:
                gene_dict['gene'] = gene_dict.pop('Element_symbol')
            elif 'Element symbol' in gene_dict:
                gene_dict['gene'] = gene_dict.pop('Element symbol')
            elif 'gene' not in gene_dict:
                for col in df.columns:
                    if 'gene' in col.lower():
                        gene_dict['gene'] = gene_dict.pop(col)
                        break
            amr_details[sample].append(gene_dict)
            gene = gene_dict.get('gene', 'unknown')
            gene_counts[gene] += 1
        return dict(amr_details), dict(gene_counts)

    def load_abricate_from_tsv(self, input_dir: Path) -> Tuple[Dict[str, Dict[str, List[Dict]]], Dict[str, Dict[str, int]]]:
        abricate_details = defaultdict(lambda: defaultdict(list))
        abricate_gene_counts = defaultdict(lambda: defaultdict(int))
        for db in self.abricate_databases:
            fname = f"pseudo_{db}_abricate_summary.tsv"
            path = input_dir / fname
            if not path.exists():
                continue
            df = pd.read_csv(path, sep='\t')
            sample_col = 'genome' if 'genome' in df.columns else 'file'
            if sample_col not in df.columns:
                sample_col = df.columns[0]
            for _, row in df.iterrows():
                sample = str(row[sample_col])
                sample = self.normalize_sample_id(sample)
                gene_dict = row.to_dict()
                gene_dict = {k: (None if pd.isna(v) else v) for k, v in gene_dict.items()}
                abricate_details[sample][db].append(gene_dict)
                gene = row['gene']
                abricate_gene_counts[db][gene] += 1
        return dict(abricate_details), dict(abricate_gene_counts)

    def load_mutations_from_tsv(self, input_dir: Path) -> Dict[str, List[Dict]]:
        tsv_path = input_dir / 'mutation_summary.tsv'  
        if not tsv_path.exists():
            print(f"  ⚠️ Mutation TSV not found: {tsv_path}")
            return {}
        df = pd.read_csv(tsv_path, sep='\t')
        mutations_by_sample = defaultdict(list)
        for _, row in df.iterrows():
            sample = str(row['genome'])
            sample = self.normalize_sample_id(sample)
            def clean_val(v):
                if pd.isna(v):
                    return ''
                return str(v)
            mut_dict = {
                'gene': clean_val(row.get('gene_symbol', '')),
                'mutation': clean_val(row.get('element_name', '')),
                'class': clean_val(row.get('class', '')),
                'subclass': clean_val(row.get('subclass', '')),
                'contig': clean_val(row.get('contig_id', '')),
                'start': clean_val(row.get('start', '')),
                'stop': clean_val(row.get('stop', '')),
                'strand': clean_val(row.get('strand', '')),
                'coverage': clean_val(row.get('coverage', '')),
                'identity': clean_val(row.get('identity', '')),
                'accession': clean_val(row.get('accession', ''))
            }
            mutations_by_sample[sample].append(mut_dict)
        return dict(mutations_by_sample)

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
            if not data:
                return pd.DataFrame()
            return pd.DataFrame(data, columns=headers)
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
                        except:
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
                    st = self.clean_st(st_raw)
                elif 'st' in [col.lower() for col in df.columns]:
                    st_col = [col for col in df.columns if col.lower() == 'st'][0]
                    st_raw = str(row[st_col]) if pd.notna(row[st_col]) else 'ND'
                    st = self.clean_st(st_raw)
                allele_profile = 'ND'
                if 'Allele Profile' in df.columns:
                    allele_profile = str(row['Allele Profile']) if pd.notna(row['Allele Profile']) else 'ND'
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


class PseudomonasDataAnalyzer:
    def __init__(self):
        self.critical_amr_genes = {
            'blaKPC', 'blaNDM', 'blaVIM', 'blaIMP', 'blaOXA', 'blaGES',
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6',
            'armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD'
        }
        self.critical_virulence_genes = {
            'exoU', 'exoS', 'exoT', 'exoY', 'exoA', 'toxA',
            'pld', 'plcH', 'plcN', 'lasB', 'lasA', 'aprA'
        }

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
            gene_list = []
            for gene, data in amr_data.items():
                gene_list.append({
                    'gene': gene,
                    'database': 'AMRfinder',
                    'frequency': str(data.get('count', 0)),
                    'count': data.get('count', 0),
                    'genomes': data.get('genomes', [])
                })
            if gene_list:
                gene_centric['amr_databases']['amrfinder'] = sorted(gene_list, key=lambda x: x['count'], reverse=True)

        if 'abricate' in integrated_data.get('gene_frequencies', {}):
            for db_name, db_genes in integrated_data['gene_frequencies']['abricate'].items():
                gene_list = []
                for gene, data in db_genes.items():
                    gene_list.append({
                        'gene': gene,
                        'database': db_name.upper(),
                        'frequency': str(data.get('count', 0)),
                        'count': data.get('count', 0),
                        'genomes': data.get('genomes', [])
                    })
                if gene_list:
                    gene_list.sort(key=lambda x: x['count'], reverse=True)
                    if db_name in ['vfdb', 'ecoli_vf', 'pa_vf', 'pseudomonas_vf']:
                        gene_centric['virulence_databases'][db_name] = gene_list
                    elif db_name == 'plasmidfinder':
                        gene_centric['plasmid_databases'][db_name] = gene_list
                    elif db_name == 'bacmet2':
                        gene_centric['bacmet_databases'][db_name] = gene_list
                    else:
                        gene_centric['amr_databases'][db_name] = gene_list

        all_genes = []
        for db_type in ['amr_databases', 'virulence_databases', 'plasmid_databases', 'bacmet_databases']:
            for genes in gene_centric.get(db_type, {}).values():
                all_genes.extend(genes)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        gene_centric['combined_gene_frequencies'] = all_genes

        return gene_centric

    def create_cross_genome_patterns(self, integrated_data: Dict[str, Any]) -> Dict[str, Any]:
        patterns = {
            'st_distribution': Counter(),
            'o_type_distribution': Counter(),
            'st_o_combinations': defaultdict(list),
            'st_to_samples': defaultdict(list),
            'gene_cooccurrence': defaultdict(Counter)
        }
        samples_data = integrated_data.get('samples', {})

        sample_genes = defaultdict(list)
        gene_centric = integrated_data.get('gene_centric', {})
        for db_type in ['amr_databases', 'virulence_databases', 'plasmid_databases', 'bacmet_databases']:
            for genes in gene_centric.get(db_type, {}).values():
                for gene_data in genes:
                    for genome in gene_data['genomes']:
                        if gene_data['gene'] not in sample_genes[genome]:
                            sample_genes[genome].append(gene_data['gene'])

        for sample, data in samples_data.items():
            st = data.get('mlst', {}).get('ST', 'ND')
            o_type = data.get('serotype', {}).get('O_Type', 'ND')
            if st != 'ND':
                patterns['st_distribution'][st] += 1
                patterns['st_to_samples'][st].append(sample)
            if o_type != 'ND':
                patterns['o_type_distribution'][o_type] += 1
            if st != 'ND' and o_type != 'ND':
                patterns['st_o_combinations'][f"{st} - {o_type}"].append(sample)

            genes = sample_genes.get(sample, [])
            for i, g1 in enumerate(genes):
                for g2 in genes[i+1:]:
                    patterns['gene_cooccurrence'][g1][g2] += 1

        return patterns


class PseudomonasHTMLGenerator:
    def __init__(self, data_analyzer: PseudomonasDataAnalyzer):
        self.data_analyzer = data_analyzer

    def generate_main_report(self, integrated_data: Dict[str, Any], output_dir: Path) -> str:
        print("\n🎨 Generating P. aeruginosa SAMPLE-CENTRIC HTML report...")
        html = self._create_sample_centric_html(
            metadata=integrated_data.get('metadata', {}),
            samples_data=integrated_data.get('samples', {}),
            patterns=integrated_data.get('patterns', {}),
            gene_centric=integrated_data.get('gene_centric', {}),
            qc_data=integrated_data.get('qc_data', {}),
            integrated_data=integrated_data
        )
        output_file = output_dir / "genius_pseudomonas_sample_centric_report.html"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"    ✅ HTML report saved: {output_file}")
        return str(output_file)

    def _create_sample_centric_html(self, **kwargs) -> str:
        css = """
        <style>
        :root {
            --brown: #8B4513;
            --brown-light: #A0522D;
            --brown-dark: #5C2E0B;
            --summary-color: var(--brown);
            --samples-color: #2E7D32;
            --qc-color: #00695C;
            --mlst-color: #1565C0;
            --serotype-color: #6A1B9A;
            --amr-color: #C62828;
            --virulence-color: #E65100;
            --plasmid-color: #3F51B5;
            --bacmet-color: #607D8B;
            --mutation-color: #00BCD4;
            --patterns-color: #AD1457;
            --citation-color: #00838F;
            --guide-color: #4527A0;
            --calltoaction-color: #2E7D32;
            --export-color: #37474F;
        }
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; line-height: 1.6; color: #333; background: #f5f5f5; min-width: 1200px; }
        .container { max-width: none; margin: 0 auto; padding: 20px; width: 100%; overflow-x: auto; }
        .main-header { background: linear-gradient(135deg, var(--brown-dark) 0%, var(--brown) 100%); color: white; padding: 30px; border-radius: 15px; margin-bottom: 30px; text-align: center; }
        .main-header h1 { font-size: 2.8em; margin-bottom: 10px; color: white; }
        .metadata-bar { background: rgba(255,255,255,0.1); padding: 15px; border-radius: 10px; margin: 20px 0; display: flex; justify-content: space-around; flex-wrap: wrap; gap: 15px; }
        .dashboard-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; margin-bottom: 30px; }
        .dashboard-card { background: white; padding: 25px; border-radius: 12px; box-shadow: 0 5px 20px rgba(0,0,0,0.1); text-align: center; transition: all 0.3s ease; cursor: pointer; border-left: 5px solid; }
        .dashboard-card:hover { transform: translateY(-10px); }
        .card-number { font-size: 3em; font-weight: bold; margin: 15px 0; background: linear-gradient(90deg, var(--brown-dark), var(--brown)); -webkit-background-clip: text; -webkit-text-fill-color: transparent; }
        .tab-navigation { display: flex; gap: 5px; margin-bottom: 20px; flex-wrap: wrap; background: white; padding: 15px; border-radius: 12px; position: sticky; top: 10px; z-index: 100; box-shadow: 0 5px 20px rgba(0,0,0,0.1); }
        .tab-button { padding: 12px 25px; background: #f5f5f5; border: none; border-radius: 8px; cursor: pointer; font-weight: 600; color: #666; transition: all 0.3s ease; display: flex; align-items: center; gap: 8px; }
        .tab-button.active { color: white; }
        .tab-button.summary.active { background: var(--summary-color); }
        .tab-button.samples.active { background: var(--samples-color); }
        .tab-button.qc.active { background: var(--qc-color); }
        .tab-button.mlst.active { background: var(--mlst-color); }
        .tab-button.serotype.active { background: var(--serotype-color); }
        .tab-button.amr.active { background: var(--amr-color); }
        .tab-button.virulence.active { background: var(--virulence-color); }
        .tab-button.plasmid.active { background: var(--plasmid-color); }
        .tab-button.bacmet.active { background: var(--bacmet-color); }
        .tab-button.mutation.active { background: var(--mutation-color); }
        .tab-button.patterns.active { background: var(--patterns-color); }
        .tab-button.citation.active { background: var(--citation-color); }
        .tab-button.guide.active { background: var(--guide-color); }
        .tab-button.calltoaction.active { background: var(--calltoaction-color); }
        .tab-button.export.active { background: var(--export-color); }
        .tab-content { display: none; background: white; padding: 30px; border-radius: 15px; margin-bottom: 30px; width: 100%; overflow-x: auto; animation: fadeIn 0.5s ease; box-shadow: 0 10px 30px rgba(0,0,0,0.1); }
        .tab-content.active { display: block; }
        @keyframes fadeIn { from { opacity: 0; transform: translateY(20px); } to { opacity: 1; transform: translateY(0); } }
        .section-header { color: #2c3e50; margin-bottom: 25px; padding-bottom: 15px; border-bottom: 3px solid var(--brown); font-size: 1.8em; display: flex; justify-content: space-between; }
        .data-table { width: 100%; border-collapse: collapse; margin: 20px 0; font-size: 0.95em; box-shadow: 0 2px 10px rgba(0,0,0,0.1); border-radius: 8px; overflow: hidden; }
        .data-table th { background: var(--brown-dark); color: white; padding: 15px; text-align: left; cursor: pointer; white-space: nowrap; }
        .data-table th:hover { background: var(--brown); }
        .data-table td { padding: 12px; border-bottom: 1px solid #e0e0e0; vertical-align: top; word-wrap: break-word; }
        .data-table tr:hover { background: #f8f9fa; }
        .sort-icon { margin-left: 5px; font-size: 0.8em; opacity: 0.6; }
        .search-box { width: 100%; padding: 12px; margin-bottom: 20px; border: 2px solid #e0e0e0; border-radius: 8px; font-size: 1em; transition: 0.3s; }
        .search-box:focus { outline: none; border-color: var(--brown); box-shadow: 0 0 0 3px rgba(139,69,19,0.1); }
        .action-buttons { display: flex; gap: 10px; margin: 20px 0; flex-wrap: wrap; clear: both; }
        .action-btn { padding: 10px 20px; border: none; border-radius: 8px; cursor: pointer; font-weight: 600; display: inline-flex; align-items: center; gap: 8px; transition: 0.3s; }
        .action-btn:hover { transform: translateY(-2px); box-shadow: 0 5px 15px rgba(0,0,0,0.2); }
        .btn-primary { background: var(--brown); color: white; }
        .btn-danger { background: #dc3545; color: white; }
        .btn-warning { background: #ffc107; color: black; }
        .btn-info { background: #17a2b8; color: white; }
        .btn-secondary { background: #6c757d; color: white; }
        .btn-light { background: #f8f9fa; color: #212529; border: 1px solid #dee2e6; }
        .badge { display: inline-block; padding: 5px 15px; border-radius: 20px; font-size: 0.85em; font-weight: 600; margin: 2px; }
        .badge-critical { background: #8B0000; color: white; }
        .alert-box { padding: 20px; border-radius: 10px; margin: 20px 0; display: flex; align-items: center; gap: 20px; border-left: 5px solid; }
        .alert-info { background: #d1ecf1; color: #0c5460; border-left-color: #17a2b8; }
        .alert-danger { background: #f8d7da; color: #721c24; border-left-color: #dc3545; }
        .alert-success { background: #d4edda; color: #155724; border-left-color: #28a745; }
        .alert-warning { background: #fff3cd; color: #856404; border-left-color: #ffc107; }
        .genome-list { display: flex; flex-wrap: wrap; gap: 5px; max-height: 200px; overflow-y: auto; padding: 5px; background: #f8f9fa; border-radius: 5px; }
        .genome-tag { background: #e3f2fd; color: #1976d2; padding: 3px 10px; border-radius: 12px; font-size: 0.85em; border: 1px solid #bbdefb; white-space: nowrap; margin: 2px; }
        .genome-tag.highlight { background-color: #ffff99 !important; color: #000 !important; border: 1px solid #ffc107; }
        .genome-group { margin-bottom: 10px; width: 100%; }
        .genome-group-header { font-weight: bold; background: #e0e0e0; padding: 4px 8px; border-radius: 4px; margin: 5px 0; font-size: 0.85em; display: inline-block; }
        .genome-group-tags { display: flex; flex-wrap: wrap; gap: 5px; margin-left: 10px; }
        .grouping-controls { background: #f0f7f0; padding: 12px; border-radius: 8px; margin: 15px 0; display: flex; flex-wrap: wrap; gap: 10px; align-items: center; border-left: 4px solid var(--brown); }
        .grouping-controls label { font-weight: bold; margin-right: 5px; }
        .group-btn { background: white; border: 1px solid var(--brown); color: var(--brown); padding: 6px 12px; border-radius: 20px; cursor: pointer; font-size: 0.85em; transition: all 0.2s; }
        .group-btn:hover { background: var(--brown); color: white; }
        .group-btn.active { background: var(--brown); color: white; }
        .scrollable-table { overflow-x: auto; width: 100%; border: 1px solid #e0e0e0; border-radius: 8px; }
        .scrollable-table table { margin: 0; }
        .citation-card { padding: 15px; border-radius: 8px; margin-bottom: 15px; border-left: 5px solid; transition: transform 0.2s; }
        .citation-card:hover { transform: translateX(5px); }
        .citation-card .copy-btn { color: white; border: none; padding: 5px 15px; border-radius: 20px; cursor: pointer; font-size: 0.8em; margin-left: 10px; }
        .citation-card a { text-decoration: none; font-weight: 600; }
        .citation-card a.doi { color: #0D6EFD; }
        .citation-card a.url { color: #2E7D32; }
        .citation-card a:hover { text-decoration: underline; }
        .footer { text-align: center; padding: 30px; background: linear-gradient(135deg, var(--brown-dark), var(--brown)); color: white; border-radius: 15px; margin-top: 40px; }
        .footer a { color: #ffd700; text-decoration: none; }
        .footer a:hover { text-decoration: underline; }
        .scientific-note { background: #f9f9e0; padding: 15px; border-radius: 10px; margin-bottom: 20px; border-left: 5px solid var(--brown); font-size: 0.95em; }
        .features-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 20px; margin: 20px 0; }
        .feature-item { background: #f8f9fa; padding: 20px; border-radius: 12px; border-left: 5px solid var(--brown); transition: all 0.3s; box-shadow: 0 2px 10px rgba(0,0,0,0.05); }
        .feature-item:hover { transform: translateY(-5px); box-shadow: 0 8px 25px rgba(0,0,0,0.15); }
        .feature-item h4 { color: var(--brown-dark); margin-bottom: 8px; }
        .feature-item i { color: var(--brown); margin-right: 8px; }
        .sample-box { border: 1px solid #ddd; border-radius: 12px; margin-bottom: 30px; padding: 20px; background: #fafafa; box-shadow: 0 2px 8px rgba(0,0,0,0.06); }
        .sample-box .sample-header { display: flex; align-items: center; gap: 15px; flex-wrap: wrap; margin-bottom: 15px; border-bottom: 2px solid #e0e0e0; padding-bottom: 10px; }
        .sample-box .sample-header h3 { font-size: 1.4em; margin: 0; }
        .sample-box .sample-header .total-badge { background: var(--brown); color: white; padding: 4px 16px; border-radius: 20px; font-weight: bold; font-size: 0.9em; }
        .sample-box .sample-header .typing-info { display: flex; flex-wrap: wrap; gap: 8px; margin-left: auto; }
        .sample-box .sample-header .typing-badge { display: inline-block; padding: 3px 10px; border-radius: 12px; font-size: 0.8em; font-weight: 600; background: #e0e0e0; color: #333; border: 1px solid #ccc; }
        .sample-box .database-table-wrapper { margin: 15px 0; overflow-x: auto; }
        .sample-box .database-table-wrapper table { width: 100%; border-collapse: collapse; font-size: 0.85em; min-width: 800px; }
        .sample-box .database-table-wrapper table th { background: var(--brown-dark); color: white; padding: 8px 12px; text-align: left; white-space: nowrap; }
        .sample-box .database-table-wrapper table td { padding: 8px 12px; border-bottom: 1px solid #e0e0e0; white-space: nowrap; }
        .sample-box .database-table-wrapper table tr:hover { background: #f1f1f1; }
        .sample-box .db-title { font-weight: bold; color: var(--brown); margin: 10px 0 5px 0; font-size: 1.1em; border-left: 4px solid var(--brown); padding-left: 10px; }
        .filter-controls { display: flex; flex-wrap: wrap; gap: 10px; align-items: center; background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; }
        .filter-controls select { padding: 10px; border-radius: 8px; border: 2px solid #ddd; background: white; min-width: 150px; }
        @media print { .tab-navigation, .dashboard-grid, .search-box, .action-buttons, .print-section-btn { display: none; } .tab-content.active { display: block; } }
        </style>
        """

        js = """
        <script>
        var sampleTyping = {};
        var originalGenomeLists = {};

        function switchTab(tabName) {
            document.querySelectorAll('.tab-content').forEach(t => t.classList.remove('active'));
            document.querySelectorAll('.tab-button').forEach(b => b.classList.remove('active'));
            document.getElementById(tabName + '-tab').classList.add('active');
            event.currentTarget.classList.add('active');
            window.location.hash = tabName;
        }

        function searchTable(tableId, searchId) {
            const filter = document.getElementById(searchId).value.toUpperCase();
            const rows = document.getElementById(tableId).getElementsByTagName('tr');
            for (let i = 1; i < rows.length; i++) {
                let found = false;
                for (let cell of rows[i].getElementsByTagName('td')) {
                    if (cell.textContent.toUpperCase().indexOf(filter) > -1) { found = true; break; }
                }
                rows[i].style.display = found ? '' : 'none';
            }
        }

        function highlightGenome(tableId, searchId) {
            const filter = document.getElementById(searchId).value.toUpperCase().trim();
            const table = document.getElementById(tableId);
            const allTags = table.querySelectorAll('.genome-tag');
            allTags.forEach(tag => tag.classList.remove('highlight'));
            if (filter === '') return;
            allTags.forEach(tag => {
                if (tag.textContent.toUpperCase().indexOf(filter) > -1) {
                    tag.classList.add('highlight');
                }
            });
        }

        function getTypingValue(genome, groupBy) {
            var info = sampleTyping[genome];
            if (!info) return "Unknown";
            if (groupBy === "ST") return info.ST;
            if (groupBy === "Otype") return info.Otype;
            if (groupBy === "ST-O") return info.ST + " - " + info.Otype;
            return "Unknown";
        }

        function groupRowGenomes(row, groupBy, originalList) {
            let genomesCell = null;
            for (let i = 0; i < row.cells.length; i++) {
                if (row.cells[i].querySelector('.genome-list')) {
                    genomesCell = row.cells[i];
                    break;
                }
            }
            if (!genomesCell) {
                console.warn("Could not find genomes cell in row");
                return;
            }
            var genomes = originalList.slice();
            if (genomes.length === 0) {
                genomesCell.innerHTML = '<div class="genome-list">None</div>';
                return;
            }
            var groups = {};
            genomes.forEach(function(genome) {
                var key = getTypingValue(genome, groupBy);
                if (!groups[key]) groups[key] = [];
                groups[key].push(genome);
            });
            var html = '<div class="genome-list">';
            for (var key in groups) {
                var tags = groups[key].map(g => `<span class="genome-tag">${g}</span>`).join('');
                html += `<div class="genome-group"><div class="genome-group-header">${key}</div><div class="genome-group-tags">${tags}</div></div>`;
            }
            html += '</div>';
            genomesCell.innerHTML = html;
        }

        function groupGenomesByTyping(tableId, groupBy) {
            var table = document.getElementById(tableId);
            if (!table) {
                console.error("Table not found:", tableId);
                return;
            }
            var tbody = table.tBodies[0];
            if (!tbody) {
                console.error("No tbody found in table", tableId);
                return;
            }
            var rows = tbody.rows;
            for (var i = 0; i < rows.length; i++) {
                var row = rows[i];
                var geneNameCell = row.cells[0];
                if (!geneNameCell) continue;
                var geneName = geneNameCell.textContent.trim().replace(/⚠️/g, '').trim();
                if (!originalGenomeLists[geneName]) {
                    var genomesCell = null;
                    for (var j = 0; j < row.cells.length; j++) {
                        if (row.cells[j].querySelector('.genome-list')) {
                            genomesCell = row.cells[j];
                            break;
                        }
                    }
                    if (genomesCell) {
                        var tags = genomesCell.querySelectorAll('.genome-tag');
                        var genomes = Array.from(tags).map(tag => tag.textContent.trim());
                        originalGenomeLists[geneName] = genomes;
                    } else {
                        originalGenomeLists[geneName] = [];
                    }
                }
            }
            for (var i = 0; i < rows.length; i++) {
                var row = rows[i];
                var geneNameCell = row.cells[0];
                if (!geneNameCell) continue;
                var geneName = geneNameCell.textContent.trim().replace(/⚠️/g, '').trim();
                var original = originalGenomeLists[geneName] || [];
                groupRowGenomes(row, groupBy, original);
            }
            var container = table.closest('.tab-content');
            if (container) {
                var btns = container.querySelectorAll('.group-btn');
                btns.forEach(btn => btn.classList.remove('active'));
                var activeBtn = container.querySelector(`.group-btn[data-group="${groupBy}"]`);
                if (activeBtn) activeBtn.classList.add('active');
            }
        }

        function resetGenomeList(tableId) {
            var table = document.getElementById(tableId);
            if (!table) return;
            var tbody = table.tBodies[0];
            if (!tbody) return;
            var rows = tbody.rows;
            for (var i = 0; i < rows.length; i++) {
                var row = rows[i];
                var geneNameCell = row.cells[0];
                if (!geneNameCell) continue;
                var geneName = geneNameCell.textContent.trim().replace(/⚠️/g, '').trim();
                var original = originalGenomeLists[geneName] || [];
                var genomesCell = null;
                for (var j = 0; j < row.cells.length; j++) {
                    if (row.cells[j].querySelector('.genome-list')) {
                        genomesCell = row.cells[j];
                        break;
                    }
                }
                if (genomesCell) {
                    var tags = original.map(g => `<span class="genome-tag">${g}</span>`).join('');
                    genomesCell.innerHTML = `<div class="genome-list">${tags}</div>`;
                }
            }
            var container = table.closest('.tab-content');
            if (container) {
                var btns = container.querySelectorAll('.group-btn');
                btns.forEach(btn => btn.classList.remove('active'));
            }
        }

        function sortTable(tableId, colIndex, type) {
            const table = document.getElementById(tableId);
            const tbody = table.tBodies[0];
            const rows = Array.from(tbody.rows);
            const isAscending = table.getAttribute('data-sort-dir') !== 'asc';
            rows.sort((a,b) => {
                let aVal = a.cells[colIndex].innerText.trim();
                let bVal = b.cells[colIndex].innerText.trim();
                if (type === 'number') {
                    aVal = parseFloat(aVal.replace(/,/g,'')) || 0;
                    bVal = parseFloat(bVal.replace(/,/g,'')) || 0;
                    return isAscending ? aVal - bVal : bVal - aVal;
                } else {
                    return isAscending ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
                }
            });
            tbody.append(...rows);
            table.setAttribute('data-sort-dir', isAscending ? 'asc' : 'desc');
            const headers = table.querySelectorAll('th');
            headers.forEach((th,idx) => { let icon = th.querySelector('.sort-icon'); if(icon) icon.innerHTML = '⇅'; });
            const currentHeader = headers[colIndex];
            const icon = currentHeader.querySelector('.sort-icon');
            if(icon) icon.innerHTML = isAscending ? '↑' : '↓';
        }

        function exportTableToCSV(tableId, filename) {
            const rows = document.getElementById(tableId).querySelectorAll('tr');
            const csv = [];
            for (let row of rows) {
                const cols = row.querySelectorAll('td, th');
                const rowData = Array.from(cols).map(c => '"' + c.innerText.replace(/"/g, '""') + '"');
                csv.push(rowData.join(','));
            }
            const blob = new Blob([csv.join('\\n')], {type: 'text/csv'});
            const a = document.createElement('a');
            a.href = URL.createObjectURL(blob);
            a.download = filename;
            a.click();
            URL.revokeObjectURL(a.href);
        }

        function printSection(sectionId) {
            const content = document.getElementById(sectionId);
            const win = window.open('', '_blank');
            win.document.write('<html><head><title>Print</title><style>' + document.querySelector('style').textContent + '</style></head><body>');
            win.document.write(content.innerHTML);
            win.document.write('</body></html>');
            win.document.close();
            win.print();
        }

        function filterBoxes(tabId) {
            var search = document.getElementById('search-' + tabId).value.toUpperCase();
            var dbFilter = document.getElementById('dbFilter-' + tabId).value;
            var boxes = document.querySelectorAll('#' + tabId + '-tab .sample-box');
            boxes.forEach(function(box) {
                var sample = box.getAttribute('data-sample') || '';
                var show = true;
                if (search && sample.toUpperCase().indexOf(search) === -1) show = false;
                var dbWrappers = box.querySelectorAll('.database-table-wrapper');
                if (dbFilter !== 'all') {
                    dbWrappers.forEach(function(wrapper) {
                        var dbName = wrapper.getAttribute('data-db') || '';
                        if (dbName === dbFilter) {
                            wrapper.style.display = '';
                        } else {
                            wrapper.style.display = 'none';
                        }
                    });
                    var dbTitles = box.querySelectorAll('.db-title');
                    dbTitles.forEach(function(title) {
                        var wrapper = title.nextElementSibling;
                        if (wrapper && wrapper.style.display === 'none') {
                            title.style.display = 'none';
                        } else {
                            title.style.display = '';
                        }
                    });
                } else {
                    dbWrappers.forEach(function(wrapper) { wrapper.style.display = ''; });
                    box.querySelectorAll('.db-title').forEach(function(title) { title.style.display = ''; });
                }
                box.style.display = show ? '' : 'none';
            });
        }

        function resetBoxFilters(tabId) {
            document.getElementById('search-' + tabId).value = '';
            document.getElementById('dbFilter-' + tabId).value = 'all';
            filterBoxes(tabId);
        }

        document.addEventListener('DOMContentLoaded', function() {
            const hash = window.location.hash.substring(1);
            if(hash) {
                let btn = document.querySelector(`.tab-button.${hash}`);
                if(btn) btn.click();
            } else document.querySelector('.tab-button').click();

            document.querySelectorAll('.data-table').forEach(table => {
                const headers = table.querySelectorAll('th');
                headers.forEach((th, idx) => {
                    const type = th.getAttribute('data-sort') || 'string';
                    th.style.cursor = 'pointer';
                    th.addEventListener('click', () => sortTable(table.id, idx, type));
                    let icon = document.createElement('span');
                    icon.className = 'sort-icon';
                    icon.innerHTML = '⇅';
                    th.appendChild(icon);
                });
            });

            document.querySelectorAll('.copy-btn').forEach(btn => {
                btn.addEventListener('click', function() {
                    const citation = this.getAttribute('data-citation');
                    navigator.clipboard.writeText(citation).then(() => {
                        const originalText = this.innerHTML;
                        this.innerHTML = '✓ Copied!';
                        setTimeout(() => { this.innerHTML = originalText; }, 2000);
                    });
                });
            });
        });
        </script>
        """

        samples_data = kwargs['samples_data']
        sample_typing_js = {}
        for sample, data in samples_data.items():
            st = data.get('mlst', {}).get('ST', 'ND')
            otype = data.get('serotype', {}).get('O_Type', 'ND')
            sample_typing_js[sample] = {'ST': st, 'Otype': otype}
        sample_typing_json = json.dumps(sample_typing_js)
        js = js.replace('var sampleTyping = {};', f'var sampleTyping = {sample_typing_json};')

        metadata = kwargs['metadata']
        samples = kwargs['samples_data']
        patterns = kwargs['patterns']
        gene_centric = kwargs['gene_centric']
        qc_data = kwargs.get('qc_data', {})
        total_amr = sum(len(genes) for genes in gene_centric.get('amr_databases', {}).values())
        total_vir = sum(len(genes) for genes in gene_centric.get('virulence_databases', {}).values())
        total_plasmid = sum(len(genes) for genes in gene_centric.get('plasmid_databases', {}).values())
        total_bacmet = sum(len(genes) for genes in gene_centric.get('bacmet_databases', {}).values())
        mutation_data = kwargs.get('integrated_data', {}).get('mutation_data', {})
        total_mutations = len(mutation_data.get('mutations', []))

        html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>GENIUS P. aeruginosa Sample-Centric Report</title>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
{css}{js}</head>
<body><div class="container">
<div class="main-header"><h1><i class="fas fa-dna"></i> GENIUS P. aeruginosa Sample-Centric Analysis Report</h1>
<p>Interactive Isolate Boxes for AMR, Virulence, Plasmids, Bacmet, and Mutations with Typing Badges</p>
<div class="metadata-bar"><div class="metadata-item"><i class="fas fa-calendar"></i> {metadata.get('analysis_date','Unknown')}</div>
<div class="metadata-item"><i class="fas fa-database"></i> Samples: {len(samples)}</div>
<div class="metadata-item"><i class="fas fa-user-md"></i>Pathogen: P. aeruginosa</div>
<div class="metadata-item"><i class="fas fa-university"></i> University of Ghana Medical School</div></div></div>
<div class="dashboard-grid">
    <div class="dashboard-card card-summary" onclick="switchTab('summary')"><div class="card-number">{len(samples)}</div><div class="card-label">Total Samples</div><i class="fas fa-vial fa-2x"></i></div>
    <div class="dashboard-card card-mlst" onclick="switchTab('mlst')"><div class="card-number">{len(patterns.get('st_distribution',{}))}</div><div class="card-label">Unique STs</div><i class="fas fa-code-branch fa-2x"></i></div>
    <div class="dashboard-card card-serotype" onclick="switchTab('serotype')"><div class="card-number">{len(patterns.get('o_type_distribution',{}))}</div><div class="card-label">O‑types</div><i class="fas fa-tag fa-2x"></i></div>
    <div class="dashboard-card card-amr" onclick="switchTab('amr')"><div class="card-number">{total_amr}</div><div class="card-label">AMR Genes</div><i class="fas fa-biohazard fa-2x"></i></div>
    <div class="dashboard-card card-virulence" onclick="switchTab('virulence')"><div class="card-number">{total_vir}</div><div class="card-label">Virulence Genes</div><i class="fas fa-virus fa-2x"></i></div>
    <div class="dashboard-card card-plasmid" onclick="switchTab('plasmid')"><div class="card-number">{total_plasmid}</div><div class="card-label">Plasmid Replicons</div><i class="fas fa-dna fa-2x"></i></div>
    <div class="dashboard-card card-bacmet" onclick="switchTab('bacmet')"><div class="card-number">{total_bacmet}</div><div class="card-label">Bacmet2 Genes</div><i class="fas fa-flask fa-2x"></i></div>
    <div class="dashboard-card card-mutation" onclick="switchTab('mutation')"><div class="card-number">{total_mutations}</div><div class="card-label">Mutations</div><i class="fas fa-dna fa-2x"></i></div>
</div>
<div class="tab-navigation">
    <button class="tab-button summary active" onclick="switchTab('summary')"><i class="fas fa-chart-pie"></i> Summary</button>
    <button class="tab-button samples" onclick="switchTab('samples')"><i class="fas fa-list-alt"></i> Sample Overview</button>
    <button class="tab-button qc" onclick="switchTab('qc')"><i class="fas fa-chart-line"></i> FASTA QC</button>
    <button class="tab-button mlst" onclick="switchTab('mlst')"><i class="fas fa-code-branch"></i> MLST</button>
    <button class="tab-button serotype" onclick="switchTab('serotype')"><i class="fas fa-tag"></i> O‑Serotype</button>
    <button class="tab-button amr" onclick="switchTab('amr')"><i class="fas fa-biohazard"></i> AMR</button>
    <button class="tab-button virulence" onclick="switchTab('virulence')"><i class="fas fa-virus"></i> Virulence</button>
    <button class="tab-button plasmid" onclick="switchTab('plasmid')"><i class="fas fa-dna"></i> Plasmids</button>
    <button class="tab-button bacmet" onclick="switchTab('bacmet')"><i class="fas fa-flask"></i> Bacmet2 </button>
    <button class="tab-button mutation" onclick="switchTab('mutation')"><i class="fas fa-dna"></i> Mutations</button>
    <button class="tab-button patterns" onclick="switchTab('patterns')"><i class="fas fa-project-diagram"></i> Patterns</button>
    <button class="tab-button citation" onclick="switchTab('citation')"><i class="fas fa-book"></i> Citation</button>
    <button class="tab-button guide" onclick="switchTab('guide')"><i class="fas fa-question-circle"></i> Guide</button>
    <button class="tab-button calltoaction" onclick="switchTab('calltoaction')"><i class="fas fa-globe"></i> Call to Action</button>
    <button class="tab-button export" onclick="switchTab('export')"><i class="fas fa-download"></i> Export</button>
</div>
<div id="summary-tab" class="tab-content active">{self._summary_section(kwargs)}</div>
<div id="samples-tab" class="tab-content">{self._samples_section(kwargs)}</div>
<div id="qc-tab" class="tab-content">{self._qc_section(kwargs)}</div>
<div id="mlst-tab" class="tab-content">{self._mlst_section(kwargs)}</div>
<div id="serotype-tab" class="tab-content">{self._serotype_section(kwargs)}</div>
<div id="amr-tab" class="tab-content">{self._sample_centric_section(kwargs, 'amr', 'AMR', ['amrfinder', 'card', 'resfinder', 'ncbi', 'megares', 'argannot'])}</div>
<div id="virulence-tab" class="tab-content">{self._sample_centric_section(kwargs, 'virulence', 'Virulence', ['vfdb'])}</div>
<div id="plasmid-tab" class="tab-content">{self._sample_centric_section(kwargs, 'plasmid', 'Plasmids', ['plasmidfinder'])}</div>
<div id="bacmet-tab" class="tab-content">{self._sample_centric_section(kwargs, 'bacmet', 'Bacmet2', ['bacmet2'])}</div>
<div id="mutation-tab" class="tab-content">{self._mutation_boxes_section(kwargs)}</div>
<div id="patterns-tab" class="tab-content">{self._patterns_section(kwargs)}</div>
<div id="citation-tab" class="tab-content">{self._citation_section()}</div>
<div id="guide-tab" class="tab-content">{self._guide_section()}</div>
<div id="calltoaction-tab" class="tab-content">{self._calltoaction_section()}</div>
<div id="export-tab" class="tab-content">{self._export_section()}</div>
<div class="footer">
    <h3>GENIUS P. aeruginosa Sample-Centric Reporter v1.0.0</h3>
    <p>University of Ghana Medical School | Department of Medical Biochemistry</p>
    <p>Author: Brown Beckley &lt;brownbeckley94@gmail.com&gt;</p>
    <p>GitHub: <a href="https://github.com/bbeckley-hub/pseudoscope" target="_blank">https://github.com/bbeckley-hub/pseudoscope</a></p>
    <p>Generated on {metadata.get('analysis_date','Unknown')}</p>
</div>
</div></body></html>"""
        return html

    def _summary_section(self, kwargs):
        samples = kwargs['samples_data']
        patterns = kwargs['patterns']
        total = len(samples)
        st_unique = len(patterns.get('st_distribution', {}))
        o_unique = len(patterns.get('o_type_distribution', {}))
        gene_centric = kwargs['gene_centric']
        amr_total = sum(len(genes) for genes in gene_centric.get('amr_databases', {}).values())
        vir_total = sum(len(genes) for genes in gene_centric.get('virulence_databases', {}).values())
        plasmid_total = sum(len(genes) for genes in gene_centric.get('plasmid_databases', {}).values())
        bacmet_total = sum(len(genes) for genes in gene_centric.get('bacmet_databases', {}).values())
        mutation_data = kwargs.get('integrated_data', {}).get('mutation_data', {})
        mutation_count = len(mutation_data.get('mutations', []))

        features = [
            ("🧬 Sample-Centric Boxes", "Each isolate is displayed in its own box with typing badges (ST, O‑type) and detailed tables for AMR, Virulence, Plasmids, Bacmet, and Mutations."),
            ("🔗 Dynamic Grouping", "Group genomes by MLST, O‑type, or ST‑O combinations to instantly see which clones harbour specific genes or mutations."),
            ("📚 Comprehensive Database Integration", "Combines AMRfinder, ABRicate (CARD, ResFinder, VFDB, PlasmidFinder, BacMet2, etc.) for maximum sensitivity and consensus."),
            ("🔍 Interactive Filters", "Filter sample boxes by isolate name or by database to focus on specific datasets."),
            ("🧩 ST‑O Combination Analysis", "Identify epidemic clones and their associated resistance/virulence profiles."),
            ("📊 FASTA QC Metrics", "Assembly quality (N50, contig count) ensures reliable gene calls."),
            ("🧬 Mutation Tracking", "Point mutations in gyrA, parC, rpoB, 23S rRNA, etc. – detect resistance even without acquired genes."),
            ("📤 Export & AI‑Ready", "Export any table as CSV or download the complete JSON for upload to ChatGPT, Claude, or Gemini for further analysis."),
        ]

        feature_html = ""
        colours = ['#e3f2fd', '#e8f5e9', '#fff3e0', '#fce4ec', '#f3e5f5', '#e0f7fa', '#fff8e1', '#fbe9e7']
        for i, (title, desc) in enumerate(features):
            colour = colours[i % len(colours)]
            feature_html += f'''
            <div class="feature-item" style="background: {colour}; border-left-color: var(--brown);">
                <h4><i class="fas fa-{title.split()[0].lower()}"></i> {title}</h4>
                <p>{desc}</p>
            </div>
            '''

        return f"""
        <div class="scientific-note"><i class="fas fa-flask"></i> <strong>Why this matters – Surveillance of <em>P. aeruginosa</em></strong><br>
        <em>Pseudomonas aeruginosa</em> is a leading cause of healthcare‑associated infections, especially in immunocompromised patients, cystic fibrosis, and burn wounds. 
        It exhibits both intrinsic and acquired resistance to many antibiotics, and its virulence arsenal (T3SS, exotoxin A, biofilms) makes infections difficult to treat.
        <br><br>
        This report provides a sample‑centric view – each isolate is presented in its own box with typing badges and detailed gene/mutation tables. This allows rapid assessment of individual isolates and their resistance/virulence profiles.
        </div>
        <div class="alert-box alert-info"><i class="fas fa-info-circle fa-2x"></i><div><h3>Analysis Overview</h3><p>Sample‑centric report for <strong>{total}</strong> P. aeruginosa genomes. Each isolate has its own box with typing badges and detailed tables. Use the sample box tabs to explore AMR, virulence, plasmids, bacmet, and mutations per isolate.</p></div></div>
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="switchTab('amr')"><i class="fas fa-biohazard"></i> View AMR Boxes</button>
            <button class="action-btn btn-success" onclick="switchTab('virulence')"><i class="fas fa-virus"></i> View Virulence Boxes</button>
            <button class="action-btn btn-info" onclick="switchTab('plasmid')"><i class="fas fa-dna"></i> View Plasmid Boxes</button>
            <button class="action-btn btn-secondary" onclick="switchTab('bacmet')"><i class="fas fa-flask"></i> View Bacmet Boxes</button>
            <button class="action-btn btn-warning" onclick="switchTab('mutation')"><i class="fas fa-dna"></i> View Mutation Boxes</button>
        </div>
        <h3>Key Statistics</h3>
        <div class="scrollable-table"><table id="summary-stats" class="data-table"><thead><tr><th data-sort="string">Metric</th><th data-sort="string">Count</th><th>Details</th></tr></thead><tbody>
        <tr><td>Total Samples</td><td><strong>{total}</strong></td><td>Complete genomic analysis</td></tr>
        <tr><td>Unique STs</td><td><strong>{st_unique}</strong></td><td>MLST</td></tr>
        <tr><td>Unique O‑types</td><td><strong>{o_unique}</strong></td><td>PAST serotyping</td></tr>
        <tr><td>AMR Genes</td><td><strong>{amr_total}</strong></td><td>All databases (AMRfinder + ABRicate)</td></tr>
        <tr><td>Virulence Genes</td><td><strong>{vir_total}</strong></td><td>VFDB</td></tr>
        <tr><td>Plasmid Replicons</td><td><strong>{plasmid_total}</strong></td><td>PlasmidFinder database</td></tr>
        <tr><td>Bacmet2 Genes</td><td><strong>{bacmet_total}</strong></td><td>Biocide & heavy metal resistance</td></tr>
        <tr><td>Unique Mutations</td><td><strong>{mutation_count}</strong></td><td>AMRfinderPlus point mutations</td></tr>
        </tbody></table></div>
        <h3>✨ Key Features of This Report</h3>
        <div class="features-grid">{feature_html}</div>
        """

    def _samples_section(self, kwargs):
        samples = kwargs['samples_data']
        qc_data = kwargs.get('qc_data', {})
        html = """
        <div class="scientific-note"><i class="fas fa-database"></i> <strong>Population Structure Overview</strong><br>
        Understanding the genetic background of your <em>P. aeruginosa</em> isolates is critical for outbreak tracking and risk assessment.
        <ul>
            <li><strong>MLST (Sequence Type)</strong>: Multilocus sequence typing defines STs. High‑risk clones (e.g., ST111, ST235, ST244, ST277) are often associated with carbapenem resistance and severe infections.</li>
            <li><strong>O‑Serotype</strong>: The lipopolysaccharide O‑antigen is a vaccine target. Prevalent types include O1, O2, O5, O6, O11. Serotype O11 is linked to epidemic clones.</li>
            <li><strong>N50</strong>: A key assembly quality metric (the shortest contig length such that half the genome is in contigs of that size or larger). Higher N50 = better assembly.</li>
            <li><strong>Virulence Count</strong>: Number of detected virulence genes from VFDB. Higher counts may indicate greater pathogenic potential.</li>
        </ul>
        </div>
        <input type="text" class="search-box" id="search-samples" onkeyup="searchTable('samples-table','search-samples')" placeholder="🔍 Search samples...">
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="exportTableToCSV('samples-table','sample_overview.csv')"><i class="fas fa-download"></i> Export CSV</button></div>
        <div class="scrollable-table"><table id="samples-table" class="data-table"><thead><tr><th data-sort="string">Sample</th><th data-sort="string">ST</th><th data-sort="string">O‑Type</th><th data-sort="number">N50</th><th data-sort="number">Virulence Count</th></tr></thead><tbody>
        """
        for sample, data in samples.items():
            st = data.get('mlst',{}).get('ST','ND')
            o_type = data.get('serotype',{}).get('O_Type','ND')
            n50 = qc_data.get(sample, {}).get('N50', 'ND')
            if isinstance(n50, float):
                n50 = f"{n50:,.0f}"
            vir_cnt = len(data.get('virulence_genes',[]))
            html += f'<tr><td><strong>{sample}</strong></td><td>{st}</td><td>{o_type}</td><td>{n50}</td><td>{vir_cnt}</td></tr>'
        html += '</tbody></table></div>'
        return html

    def _qc_section(self, kwargs):
        qc_data = kwargs.get('qc_data', {})
        if not qc_data:
            return '<div class="alert-box alert-warning"><i class="fas fa-exclamation-circle"></i><div>No FASTA QC data found.</div></div>'

        html = """
        <div class="scientific-note" style="background: linear-gradient(135deg, #f8f9fa 0%, #d1ecf1 100%); border-left: 6px solid #17a2b8; margin-bottom: 20px;">
            <div style="display: flex; align-items: center; gap: 12px; flex-wrap: wrap;">
                <span style="font-size: 1.4em;">🐍</span>
                <div>
                    <strong style="font-size: 1.2em; color: #0c5460;">FASTA QC – Powered by Biopython</strong><br>
                    <span style="font-size: 0.95em; color: #333;">
                        <strong>Biopython</strong> by <a href="https://biopython.org/" target="_blank" style="color: #17a2b8; font-weight: bold;">The Biopython Consortium</a> 
                        <i class="fas fa-arrow-right"></i> Used for parsing FASTA files and computing quality metrics (GC%, N50, contig counts, etc.).
                        <br>
                        <span style="color: #6c757d;">Biopython is an open‑source toolkit for computational biology, essential for robust and reproducible assembly QC.</span>
                    </span>
                </div>
            </div>
        </div>
        """

        html += """
        <div class="scientific-note"><i class="fas fa-chart-line"></i> <strong>FASTA Quality Control – Why It Matters</strong><br>
        Assembly quality directly affects gene detection. Poor assemblies (low N50, high contig count) may miss genes or create false positives.
        <ul>
            <li><strong>N50</strong>: Ideally >20 kb for accurate gene annotation.</li>
            <li><strong>Contig count</strong>: Fewer contigs (e.g., <200) indicate a more complete assembly.</li>
            <li><strong>GC%</strong>: <em>P. aeruginosa</em> typically 66–67%. Deviations may indicate contamination.</li>
        </ul>
        </div>
        <input type="text" class="search-box" id="search-qc" onkeyup="searchTable('qc-table','search-qc')" placeholder="🔍 Search sample...">
        <div class="scrollable-table"><table id="qc-table" class="data-table"><thead><tr><th data-sort="string">Sample</th>"""
        metrics = set()
        for d in qc_data.values():
            metrics.update(d.keys())
        metrics = sorted(metrics)

        for m in metrics:
            html += f'<th data-sort="number">{m}</th>'
        html += '</tr></thead><tbody>'
        for sample, vals in sorted(qc_data.items()):
            html += f'<tr><td><strong>{sample}</strong></td>'
            for m in metrics:
                v = vals.get(m, 'ND')
                if isinstance(v, float):
                    v = f"{v:,.0f}" if v > 1000 else f"{v:.2f}"
                html += f'<td>{v}</td>'
            html += '</tr>'
        html += '</tbody></table></div>'
        return html

    def _mlst_section(self, kwargs):
        patterns = kwargs['patterns']
        st_dist = patterns.get('st_distribution', Counter())
        st_to_samples = patterns.get('st_to_samples', defaultdict(list))
        total = sum(st_dist.values())

        html = """
        <div class="scientific-note" style="background: linear-gradient(135deg, #f8f9fa 0%, #cce5ff 100%); border-left: 6px solid #007bff; margin-bottom: 20px;">
            <div style="display: flex; align-items: center; gap: 12px; flex-wrap: wrap;">
                <span style="font-size: 1.4em;">🧬</span>
                <div>
                    <strong style="font-size: 1.2em; color: #004085;">MLST – Powered by PubMLST & Seemann's MLST Scripts</strong><br>
                    <span style="font-size: 0.95em; color: #333;">
                        <strong>PubMLST</strong> database by <a href="https://pubmlst.org/" target="_blank" style="color: #007bff; font-weight: bold;">Keith Jolley & the PubMLST team</a> 
                        <i class="fas fa-arrow-right"></i> Typing schemes integrated via 
                        <a href="https://github.com/tseemann/mlst" target="_blank" style="color: #dc3545; font-weight: bold;">Seemann's MLST scripts</a>.
                        <br>
                        <span style="color: #6c757d;">
                            <i class="fas fa-exclamation-triangle" style="color: #ffc107;"></i> 
                            <strong>Important:</strong> Due to licensing restrictions on the PubMLST database, we cannot include alleles added after <strong>December 31, 2024</strong>. 
                            To use the latest alleles, you can update your local database by following the 
                            <a href="https://github.com/tseemann/mlst#updating-the-databases" target="_blank" style="color: #007bff;">instructions from Seemann's MLST</a>. 
                        </span>
                    </span>
                </div>
            </div>
        </div>
        """

        html += """
        <div class="scientific-note"><i class="fas fa-code-branch"></i> <strong>MLST – Tracking High‑Risk Clones</strong><br>
        Multi‑locus sequence typing uses seven housekeeping genes (<em>acsA, aroE, guaA, mutL, nuoD, ppsA, trpE</em>). 
        Globally disseminated high‑risk clones include:
        <ul>
            <li><strong>ST111</strong>: Often carries VIM‑2 carbapenemase and ExoU; associated with ventilator‑associated pneumonia.</li>
            <li><strong>ST235</strong>: Major carbapenemase producer (KPC, NDM, VIM) with ExoU; pandemic clone.</li>
            <li><strong>ST244</strong>: Common in environmental and some clinical settings; often less virulent.</li>
            <li><strong>ST277</strong>: Associated with multidrug resistance and cystic fibrosis.</li>
        </ul>
        </div>
        <h3>ST Distribution (Counts)</h3>
        <input type="text" class="search-box" id="search-mlst-counts" onkeyup="searchTable('mlst-counts-table','search-mlst-counts')" placeholder="🔍 Search ST...">
        <div class="scrollable-table"><table id="mlst-counts-table" class="data-table"><thead><tr><th data-sort="string">ST</th><th data-sort="number">Count</th><th data-sort="number">Percentage</th><th data-sort="string">Associated O‑types</th></tr></thead><tbody>
        """
        for st, cnt in st_dist.most_common():
            pct = cnt/total*100
            o_types = set()
            for combo in patterns.get('st_o_combinations',{}):
                if f"{st} - " in combo:
                    o_types.add(combo.split(" - ")[1])
            o_str = ', '.join(o_types) if o_types else 'ND'
            html += f'<tr><td><strong>{st}</strong></td><td>{cnt}</td><td>{pct:.1f}%</td><td>{o_str}</td></tr>'
        html += '</tbody></table></div>'

        html += """
        <h3>ST-Sample List</h3>
        <p>Each ST is shown with all samples that belong to it. Use the search box to highlight specific samples across the table.</p>
        <input type="text" class="search-box" id="search-mlst-tags" onkeyup="highlightGenome('mlst-tags-table','search-mlst-tags')" placeholder="🔍 Highlight sample...">
        <div class="scrollable-table"><table id="mlst-tags-table" class="data-table"><thead>
        <tr><th data-sort="string">ST</th><th data-sort="number">Count</th><th data-sort="string">Samples</th>
        </thead><tbody>
        """
        for st, samples_list in sorted(st_to_samples.items(), key=lambda x: len(x[1]), reverse=True):
            genome_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples_list)
            html += f'<tr><td><strong>{st}</strong></td><td>{len(samples_list)}</td><td><div class="genome-list">{genome_tags}</div></td></tr>'
        html += '</tbody></table></div>'
        return html

    def _serotype_section(self, kwargs):
        patterns = kwargs['patterns']
        samples = kwargs['samples_data']
        o_dist = patterns.get('o_type_distribution', Counter())
        total = sum(o_dist.values())
        html = f"""
        <div class="scientific-note"><i class="fas fa-tag"></i> <strong>O‑Serotype (PAST) – A Vaccine Target</strong><br>
        The O‑antigen is a major surface antigen and a component of the investigational vaccine for <em>P. aeruginosa</em>. 
        Common serotypes in clinical isolates include O1, O2, O5, O6, and O11. Serotype O11 is often associated with epidemic multidrug‑resistant clones.
        <br><br>
        Knowing the serotype distribution can guide vaccine development and infection control.
        <br><br>
        <span style="font-size: 0.9em; color: #666;"><i class="fas fa-code-branch"></i> PAST tool originally developed by the Technical University of Denmark (CGE). This implementation uses <strong>pasty</strong> by <a href="https://github.com/rpetit3/pasty" target="_blank" style="color: #dc3545;">Dr. Robert A. Petit III</a>, which makes the original PAst tool more user‑friendly and accessible for genomic surveillance. <a href="https://cge.food.dtu.dk/services/PAst/" target="_blank" style="color: #dc3545;">Original PAst‑1.0</a></span>
        </div>
        <div class="alert-box alert-info"><i class="fas fa-info-circle"></i><div><h3>O‑Type Distribution with Genome Tags</h3><p>Each O‑type is shown with all samples as genome tags. Use the search box to highlight specific isolates across the table.</p></div></div>
        <input type="text" class="search-box" id="search-serotype" onkeyup="highlightGenome('serotype-table','search-serotype')" placeholder="🔍 Highlight genome tags in O‑type table...">
        <div class="scrollable-table"><table id="serotype-table" class="data-table"><thead><tr><th data-sort="string">O‑Type</th><th data-sort="number">Count</th><th data-sort="number">Percentage</th><th data-sort="string">Sample</th></tr></thead><tbody>
        """
        for o, cnt in o_dist.most_common():
            if o == 'ND':
                continue
            pct = cnt/total*100
            sample_list = [s for s, d in samples.items() if d.get('serotype', {}).get('O_Type') == o]
            genome_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in sample_list)
            html += f'<tr><td><strong>{o}</strong></td><td>{cnt}</td><td>{pct:.1f}%</td><td><div class="genome-list">{genome_tags}</div></td></tr>'
        html += '</tbody></table></div>'
        return html

    def _grouping_controls(self, table_id: str) -> str:
        return f"""
        <div class="grouping-controls">
            <strong><i class="fas fa-layer-group"></i> Group genomes by:</strong>
            <button class="group-btn" data-group="ST" onclick="groupGenomesByTyping('{table_id}', 'ST')">MLST (ST)</button>
            <button class="group-btn" data-group="Otype" onclick="groupGenomesByTyping('{table_id}', 'Otype')">O‑type</button>
            <button class="group-btn" data-group="ST-O" onclick="groupGenomesByTyping('{table_id}', 'ST-O')">ST‑O</button>
            <button class="group-btn" onclick="resetGenomeList('{table_id}')">Reset (flat list)</button>
        </div>
        """

    def _sample_centric_section(self, kwargs, tab_id, title, db_list):
        samples_data = kwargs.get('samples_data', {})
        amr_details = kwargs.get('integrated_data', {}).get('amr_details', {})
        abricate_details = kwargs.get('integrated_data', {}).get('abricate_details', {})

        relevant_samples = []
        for sample in samples_data:
            has_data = False
            if 'amrfinder' in db_list and sample in amr_details and amr_details[sample]:
                has_data = True
            for db in db_list:
                if db != 'amrfinder' and sample in abricate_details and db in abricate_details[sample] and abricate_details[sample][db]:
                    has_data = True
            if has_data:
                relevant_samples.append(sample)
        relevant_samples.sort()

        if not relevant_samples:
            return f"""
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div><h3>No {title} Data Available</h3><p>No samples with {title} genes were found. Please check that the corresponding TSV files exist and contain data.</p></div>
            </div>
            """

        html = f"""
        <div class="scientific-note" style="background: linear-gradient(135deg, #f8f9fa 0%, #fff3cd 100%); border-left: 6px solid #17a2b8; margin-bottom: 20px;">
            <div style="display: flex; align-items: center; gap: 12px; flex-wrap: wrap;">
                <span style="font-size: 1.4em;">🧬</span>
                <div>
                    <strong style="font-size: 1.2em; color: #0c5460;">{title} Detection – Powered by ABRicate & AMRFinderPlus</strong><br>
                    <span style="font-size: 0.95em; color: #333;">
                        <strong>ABRicate</strong> by <a href="https://github.com/tseemann/abricate" target="_blank" style="color: #dc3545; font-weight: bold;">Prof. Torsten Seemann</a> 
                        <i class="fas fa-arrow-right"></i> Integrated with multiple databases.
                        <br>
                        <strong>AMRFinderPlus</strong> by <a href="https://github.com/ncbi/amr" target="_blank" style="color: #dc3545; font-weight: bold;">NCBI</a> 
                        <i class="fas fa-arrow-right"></i> Comprehensive resistance gene and point mutation detection.
                        <br>
                        <span style="color: #6c757d;">Isolate-centric view of your data</span>
                    </span>
                </div>
            </div>
        </div>
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>🧬 {title} – Interactive Sample Boxes</h3>
                <p>Each box represents one isolate. Inside, you will find separate tables for each database with full gene details – horizontally scrollable.</p>
                <p>Use the filters below to search by sample name or to show only a specific database.</p>
            </div>
        </div>
        <div class="filter-controls">
            <input type="text" class="search-box" id="search-{tab_id}" onkeyup="filterBoxes('{tab_id}')" placeholder="🔍 Search sample...">
            <select id="dbFilter-{tab_id}" onchange="filterBoxes('{tab_id}')">
                <option value="all">All Databases</option>
        """
        for db in db_list:
            display = db.upper()
            if db == 'amrfinder':
                display = 'AMRfinder'
            html += f'<option value="{db}">{display}</option>'
        html += f"""
            </select>
            <button class="action-btn btn-success" onclick="resetBoxFilters('{tab_id}')"><i class="fas fa-sync"></i> Clear Filters</button>
        </div>
        <div id="box-container-{tab_id}">
        """

        for sample in relevant_samples:
            typing = samples_data.get(sample, {}).get('typing', {})
            mlst = samples_data.get(sample, {}).get('mlst', {}).get('ST', 'ND')
            otype = samples_data.get(sample, {}).get('serotype', {}).get('O_Type', 'ND')

            html += f'<div class="sample-box" data-sample="{sample}">'
            total_genes = 0
            if 'amrfinder' in db_list and sample in amr_details:
                total_genes += len(amr_details[sample])
            for db in db_list:
                if db != 'amrfinder' and sample in abricate_details and db in abricate_details[sample]:
                    total_genes += len(abricate_details[sample][db])
            html += f"""
                <div class="sample-header">
                    <h3><i class="fas fa-microbe"></i> {sample}</h3>
                    <span class="total-badge">Total Genes: {total_genes}</span>
                    <div class="typing-info">
                        <span class="typing-badge">ST: {mlst}</span>
                        <span class="typing-badge">O-Type: {otype}</span>
                    </div>
                </div>
            """

            if 'amrfinder' in db_list and sample in amr_details and amr_details[sample]:
                html += self._make_database_table('AMRfinder', amr_details[sample], 'amrfinder')

            for db in db_list:
                if db == 'amrfinder':
                    continue
                if sample in abricate_details and db in abricate_details[sample] and abricate_details[sample][db]:
                    html += self._make_database_table(db.upper(), abricate_details[sample][db], db)

            html += '</div>'

        html += """
        </div>
        """
        return html

    def _make_database_table(self, db_name: str, genes: List[Dict], db_key: str) -> str:
        if not genes:
            return ''
        keys = list(genes[0].keys()) if genes else []
        priority = ['gene', 'product', 'coverage_percent', 'identity_percent', 'accession', 'contig', 'start', 'stop', 'class', 'subclass', 'scope', 'resistance']
        ordered = []
        for p in priority:
            if p in keys:
                ordered.append(p)
        for k in keys:
            if k not in ordered:
                ordered.append(k)

        html = f"""
                <div class="database-table-wrapper" data-db="{db_key}">
                    <div class="db-title">{db_name}</div>
                    <table>
                        <thead>
                            <tr>
        """
        for col in ordered:
            display = col.replace('_', ' ').title()
            html += f'<th>{display}</th>'
        html += '</tr></thead><tbody>'
        for gene_dict in genes:
            html += '<tr>'
            for col in ordered:
                val = gene_dict.get(col, '')
                if val is None:
                    val = ''
                html += f'<td>{val}</td>'
            html += '</tr>'
        html += """
                        </tbody>
                    </table>
                </div>
        """
        return html

    def _mutation_boxes_section(self, kwargs):
        samples_data = kwargs.get('samples_data', {})
        integrated_data = kwargs.get('integrated_data', {})
        mutation_details = integrated_data.get('mutation_details', {})
        relevant_samples = [s for s in samples_data if s in mutation_details and mutation_details[s]]
        relevant_samples.sort()

        if not relevant_samples:
            return """
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div><h3>No Mutation Data Available</h3><p>No mutations found for any sample. Please check that pseudo_mutation_summary.tsv exists and contains data.</p></div>
            </div>
            """

        html = f"""
        <div class="scientific-note" style="background: linear-gradient(135deg, #f8f9fa 0%, #d1ecf1 100%); border-left: 6px solid #17a2b8; margin-bottom: 20px;">
            <div style="display: flex; align-items: center; gap: 12px; flex-wrap: wrap;">
                <span style="font-size: 1.4em;">🧬</span>
                <div>
                    <strong style="font-size: 1.2em; color: #0c5460;">Point Mutations – Powered by AMRFinderPlus</strong><br>
                    <span style="font-size: 0.95em; color: #333;">
                        <strong>AMRFinderPlus</strong> by <a href="https://github.com/ncbi/amr" target="_blank" style="color: #dc3545; font-weight: bold;">NCBI</a> 
                        <i class="fas fa-arrow-right"></i> Detects point mutations in key resistance genes (gyrA, parC, rpoB, 23S rRNA, etc.).
                        <br>
                        <span style="color: #6c757d;">Isolate-centric view of your data.</span>
                    </span>
                </div>
            </div>
        </div>
        <div class="alert-box alert-info">
            <i class="fas fa-dna fa-2x"></i>
            <div>
                <h3>🧬 Point Mutations – Sample‑Centric Boxes</h3>
                <p>Each box represents one isolate. Inside, you will find a table of all detected point mutations with full details – horizontally scrollable.</p>
                <p>Use the search bar below to filter boxes by sample name.</p>
            </div>
        </div>
        <div class="filter-controls">
            <input type="text" class="search-box" id="search-mutation" onkeyup="filterBoxes('mutation')" placeholder="🔍 Search sample...">
            <button class="action-btn btn-success" onclick="resetBoxFilters('mutation')"><i class="fas fa-sync"></i> Clear Filters</button>
        </div>
        <div id="box-container-mutation">
        """

        for sample in relevant_samples:
            typing = samples_data.get(sample, {}).get('typing', {})
            mlst = samples_data.get(sample, {}).get('mlst', {}).get('ST', 'ND')
            otype = samples_data.get(sample, {}).get('serotype', {}).get('O_Type', 'ND')
            mutations = mutation_details.get(sample, [])
            total_mutations = len(mutations)

            html += f'<div class="sample-box" data-sample="{sample}">'
            html += f"""
                <div class="sample-header">
                    <h3><i class="fas fa-microbe"></i> {sample}</h3>
                    <span class="total-badge">Total Mutations: {total_mutations}</span>
                    <div class="typing-info">
                        <span class="typing-badge">ST: {mlst}</span>
                        <span class="typing-badge">O-Type: {otype}</span>
                    </div>
                </div>
            """

            if mutations:
                html += self._make_mutation_table(mutations)

            html += '</div>'

        html += """
        </div>
        """
        return html

    def _make_mutation_table(self, mutations: List[Dict]) -> str:
        if not mutations:
            return ''
        cols = ['gene', 'mutation', 'class', 'subclass', 'contig', 'start', 'stop', 'strand', 'coverage', 'identity', 'accession']
        html = """
        <div class="database-table-wrapper" data-db="mutations">
            <div class="db-title">Mutations</div>
            <table>
                <thead>
                    <tr>
        """
        for col in cols:
            display = col.replace('_', ' ').title()
            html += f'<th>{display}</th>'
        html += '</tr></thead><tbody>'
        for mut in mutations:
            html += '<tr>'
            for col in cols:
                val = mut.get(col, '')
                if val is None or (isinstance(val, float) and val != val):
                    val = ''
                html += f'<td>{val}</td>'
            html += '</tr>'
        html += """
                </tbody>
            </table>
        </div>
        """
        return html

    def _patterns_section(self, kwargs):
        patterns = kwargs['patterns']
        html = """
        <div class="scientific-note"><i class="fas fa-project-diagram"></i> <strong>Cross‑Genome Patterns – Tracking Clones and Co‑occurrence</strong><br>
        This section shows:
        <ul>
            <li><strong>ST‑O combination table</strong>: Which Sequence Types are associated with which O‑serotypes? Certain ST‑O pairs (e.g., ST235‑O11) define epidemic clones.</li>
            <li><strong>Gene co‑occurrence</strong>: Which pairs of genes (AMR, virulence, plasmid, bacmet) are frequently found together in the same genome. High co‑occurrence may indicate co‑selection or physical linkage (e.g., on the same plasmid).</li>
        </ul>
        Use these tables to identify outbreak clones and potential co‑selection risks.
        </div>
        <div class="alert-box alert-info"><i class="fas fa-table"></i><div><h3>ST‑O Combination Table</h3><p>Each combination of Sequence Type and O‑serotype with the list of samples as tags. Use the search box to highlight specific isolates.</p></div></div>
        <input type="text" class="search-box" id="search-sto" onkeyup="highlightGenome('sto-table','search-sto')" placeholder="🔍 Highlight sample in ST-O combinations...">
        <div class="scrollable-table"><table id="sto-table" class="data-table"><thead><tr><th data-sort="string">ST - O</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>
        """
        sto = patterns.get('st_o_combinations', {})
        for combo, samples_list in sorted(sto.items(), key=lambda x: len(x[1]), reverse=True):
            genome_tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples_list)
            html += f'<tr><td><strong>{combo}</strong></td><td>{len(samples_list)}</td><td><div class="genome-list">{genome_tags}</div></td></tr>'
        html += '</tbody></table></div>'

        cooc = patterns.get('gene_cooccurrence', {})
        cooc_list = []
        for g1, partners in cooc.items():
            for g2, cnt in partners.items():
                if g1 < g2:
                    cooc_list.append((g1, g2, cnt))
        cooc_list.sort(key=lambda x: x[2], reverse=True)
        if cooc_list:
            html += '<h3>Top 100 Gene Co‑occurrences</h3>'
            html += '<input type="text" class="search-box" id="search-cooc" onkeyup="searchTable(\'cooc-table\',\'search-cooc\')" placeholder="🔍 Search co‑occurrence pairs...">'
            html += '<div class="scrollable-table"><table id="cooc-table" class="data-table"><thead><tr><th data-sort="string">Gene 1</th><th data-sort="string">Gene 2</th><th data-sort="number">Co‑occurrence Count</th></tr></thead><tbody>'
            for g1, g2, cnt in cooc_list[:100]:
                html += f'<tr><td>{g1}</td><td>{g2}</td><td>{cnt}</td></tr>'
            html += '</tbody></table></div>'
        else:
            html += '<p>No gene co‑occurrence data available.</p>'
        return html

    def _citation_section(self):
        colours = [
            '#E3F2FD', '#E8F5E9', '#FFF3E0', '#FCE4EC', '#F3E5F5',
            '#E0F7FA', '#FFF8E1', '#FBE9E7', '#EDE7F6', '#E0F2F1',
            '#F1F8E9', '#FFFDE7', '#FCE4EC', '#E3F2FD', '#E8EAF6'
        ]
        colour_index = 0

        def format_citation(name, citation_text, link_type, link_url):
            nonlocal colour_index
            colour = colours[colour_index % len(colours)]
            colour_index += 1
            if link_type == 'doi':
                link_html = f'<a href="https://doi.org/{link_url}" target="_blank" class="doi"><i class="fas fa-external-link-alt"></i> {link_url}</a>'
            else:
                link_html = f'<a href="{link_url}" target="_blank" class="url"><i class="fas fa-external-link-alt"></i> {link_url}</a>'
            return f"""
            <div class="citation-card" style="background: {colour}; border-left-color: var(--citation-color);">
                <strong>{name}</strong> – {citation_text} {link_html}
                <button class="copy-btn" data-citation='{citation_text} {link_url}' style="background: #006064;">📋 Copy</button>
            </div>
            """
        html = """
        <div class="alert-box alert-info"><i class="fas fa-quote-right fa-2x"></i><div><h3>📚 How to Cite PseudoScope and Its Dependencies</h3><p>If you use this tool in your research, please cite the main tool and the relevant third‑party tools and databases.</p></div></div>
        """
        html += format_citation(
            "PseudoScope",
            "Beckley B, Amarh V. PseudoScope: a species‑specific bioinformatics suite for rapid and accessible Pseudomonas aeruginosa genomic analysis. Github 2026",
            "url",
            "https://github.com/brown-beckley/pseudoscope"
        )
        html += format_citation(
            "MLST",
            "Seemann T. MLST: Scan contig files against PubMLST typing schemes. GitHub.",
            "url",
            "https://github.com/tseemann/mlst"
        )
        html += format_citation(
            "PAST (P. aeruginosa serotyping) – Original Tool",
            "Thomsen MCF, et al. A web‑based tool for serotyping of Pseudomonas aeruginosa. CGE, Technical University of Denmark.",
            "url",
            "https://cge.food.dtu.dk/services/PAst/"
        )
        html += format_citation(
            "pasty",
            "Petit RA III, pasty: A tool easily taken advantage of for in silico serogrouping of Pseudomonas aeruginosa isolates. GitHub.",
            "url",
            "https://github.com/rpetit3/pasty"
        )
        html += format_citation(
                    "camlhmp",
                    "Petit RA III, camlhmp: Classification through yAML Heuristic Mapping Protocol (GitHub).",
                    "url",
                    "https://github.com/rpetit3/camlhmp"
        )
        html += format_citation(
            "AMRFinderPlus",
            "Feldgarden M, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021;11(1):12728.",
            "doi",
            "10.1038/s41598-021-91456-0"
        )
        html += format_citation(
            "ABRicate",
            "Seemann T. ABRicate: mass screening of contigs for antibiotic resistance genes. GitHub.",
            "url",
            "https://github.com/tseemann/abricate"
        )
        html += format_citation(
            "CARD",
            "McArthur AG, et al. The comprehensive antibiotic resistance database. Antimicrob Agents Chemother. 2013;57(7):3348-57.",
            "doi",
            "10.1128/AAC.00419-13"
        )
        html += format_citation(
            "ResFinder",
            "Florensa AF, et al. ResFinder – an open online resource for identification of antimicrobial resistance genes in next‑generation sequencing data and prediction of phenotypes from genotypes. Microb Genom. 2022;8(1):000748.",
            "doi",
            "10.1099/mgen.0.000748"
        )
        html += format_citation(
            "VFDB",
            "Chen L, et al. VFDB 2012 update: toward the genetic diversity and molecular evolution of bacterial virulence factors. Nucleic Acids Res. 2012;40(Database issue):D641-5.",
            "doi",
            "10.1093/nar/gkr989"
        )
        html += format_citation(
            "PlasmidFinder",
            "Carattoli A, et al. In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. Antimicrob Agents Chemother. 2014;58(7):3895-903.",
            "doi",
            "10.1128/AAC.02412-14"
        )
        html += format_citation(
            "BacMet",
            "Pal C, et al. BacMet: antibacterial biocide and metal resistance genes database. Nucleic Acids Res. 2014;42(Database issue):D737-43.",
            "doi",
            "10.1093/nar/gkt1252"
        )
        html += format_citation(
            "MEGARes",
            "Doster E, et al. MEGARes 2.0: a database for classification of antimicrobial drug, biocide and metal resistance determinants in metagenomic sequence data. Nucleic Acids Res. 2020;48(D1):D561-D569.",
            "doi",
            "10.1093/nar/gkz1010"
        )
        html += format_citation(
            "Biopython",
            "Cock PJ, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009;25(11):1422-3.",
            "doi",
            "10.1093/bioinformatics/btp163"
        )
        html += """
        <div class="alert-box alert-success"><i class="fas fa-hand-peace"></i><div><strong>Suggested acknowledgement:</strong><br>
        "Genomic analysis was performed using PseudoScope [Beckley & Amarh, 2026] which integrates MLST [Seemann, 2018] using the PubMLST database [Jolley et al., 2018], ABRicate [Seemann, 2018], AMRFinderPlus [Feldgarden et al., 2021], and PAST serotyping [Thomsen et al., DTU] via the pasty wrapper [Petit et al., 2023] for comprehensive P. aeruginosa characterization. Antimicrobial resistance genes were identified using the CARD [McArthur et al., 2013] and ResFinder [Florensa et al., 2022] databases. For biocide and heavy metal resistance genes, BacMet [Pal et al., 2014] was used. Virulence and plasmid screening were performed with ABRicate using the VFDB [Chen et al., 2012] and PlasmidFinder [Carattoli et al., 2014] databases."
        </div></div>
        """
        return html

    def _guide_section(self):
        return """
        <div class="alert-box alert-info"><i class="fas fa-question-circle fa-2x"></i><div><h3>📘 Welcome to GENIUS P. aeruginosa Sample-Centric Report – Your Guide</h3><p>This guide explains how to get the most out of each tab, why multiple databases are used, and some practical tips.</p></div></div>
        <div style="background: #f0f7f0; padding:20px; border-radius:12px; margin:20px 0;">
        <h4><i class="fas fa-cogs"></i> How to Use the Report</h4>
        <ul>
            <li><strong>Summary</strong> – Overview of your dataset and key statistics.</li>
            <li><strong>Sample Overview</strong> – List of all samples with ST, O‑type, and virulence counts.</li>
            <li><strong>FASTA QC</strong> – Assembly quality metrics. Good assemblies = reliable gene calls.</li>
            <li><strong>MLST</strong> – ST distribution and sample lists. High‑risk clones are highlighted in the notes.</li>
            <li><strong>O‑Serotype</strong> – Serotype distribution with genome tags. Use highlight to track specific isolates.</li>
            <li><strong>AMR</strong> – Each isolate is displayed in a box with typing badges and a detailed table of AMR genes. Filter by sample name or database. Horizontally scrollable tables.</li>
            <li><strong>Virulence</strong> – Same as AMR but for virulence factors. Identify hypervirulent clones.</li>
            <li><strong>Plasmids</strong> – Plasmid replicons per isolate.</li>
            <li><strong>Bacmet2 </strong> – Biocide and heavy metal resistance genes per isolate – important for hospital hygiene.</li>
            <li><strong>Mutations</strong> – Point mutations per isolate.</li>
            <li><strong>Patterns</strong> – ST‑O combinations and gene co‑occurrence. Find outbreak clones and co‑selection links.</li>
            <li><strong>Citation</strong> – All required citations, copy‑able with one click, with distinct coloured cards.</li>
            <li><strong>Call to Action</strong> – Learn about the global AMR burden and how you can contribute.</li>
            <li><strong>Export</strong> – Download tables as CSV and full JSON for AI or downstream analysis.</li>
        </ul>
        </div>
        <div style="background: #fff0f0; padding:20px; border-radius:12px; margin:20px 0; border-left: 5px solid #dc3545;">
        <h4><i class="fas fa-database"></i> Why Multiple Databases?</h4>
        <p>Some users ask: "Why use AMRfinder <strong>and</strong> ABRicate with CARD, ResFinder, and NCBI?"</p>
        <p><strong>Because no single database is perfect.</strong> Each has strengths and biases. For example:</p>
        <ul>
            <li><strong>AMRfinder</strong> – Curated by NCBI, excellent for clinically relevant genes, but may miss some environmental or novel variants.</li>
            <li><strong>CARD</strong> – Comprehensive with detailed resistance mechanisms, but often lag behind newly discovered genes.</li>
            <li><strong>ResFinder</strong> – Frequently updated and widely used in clinical labs, but focuses on acquired resistance.</li>
            <li><strong>NCBI</strong> – Large but sometimes includes non‑specific hits.</li>
        </ul>
        <p>By combining them, we <strong>maximise sensitivity</strong> and provide a consensus view. If a gene appears in multiple databases, it's a strong signal. If only one, it's worth double‑checking.</p>
        <p>Also, some databases perform better for certain species or gene families. For <em>P. aeruginosa</em>, CARD and ResFinder are generally excellent, but AMRfinder catches many carbapenemase variants. We give you all the data, you decide what to trust.</p>
        </div>
        <div style="background: #e8f5e9; padding:20px; border-radius:12px; margin:20px 0; border-left: 5px solid #28a745;">
        <h4><i class="fas fa-lightbulb"></i> Pro Tips</h4>
        <ul>
            <li><strong>Filtering sample boxes</strong> – In AMR, Virulence, Plasmid, Bacmet, and Mutation tabs, use the search bar to filter by sample name, and use the dropdown to show only specific databases.</li>
            <li><strong>Highlighting</strong> – In gene‑centric tables (MLST, Serotype, Patterns), use the "Highlight sample" search boxes to find a particular isolate – it will turn yellow!</li>
            <li><strong>Export</strong> – You can export any table as CSV. The JSON file is ideal for AI tools – upload it to ChatGPT/Claude and ask questions about your data.</li>
            <li><strong>Print</strong> – Each section has a print button to save as PDF for reports.</li>
        </ul>
        </div>
        """

    def _calltoaction_section(self):
        return """
        <div class="alert-box alert-info"><i class="fas fa-globe fa-2x"></i><div>
        <h3>The Global Burden of AMR and Our Call to Action</h3>
        <p>Antimicrobial resistance (AMR) kills an estimated <strong>1.27 million people directly each year</strong>. <em>Pseudomonas aeruginosa</em> is a WHO critical‑priority pathogen, especially when carbapenem‑resistant. Genomic surveillance is our best tool to track resistant clones, understand transmission, and guide infection control.</p>
        <p>We developed <strong>PseudoScope</strong> as part of the <strong>ESCAPE AMR</strong> project – an open‑source initiative targeting all ESKAPE pathogens. Our goal is to empower researchers, especially in low‑resource settings, with user‑friendly genomic analysis tools.</p>
        </div></div>

        <div style="background:#e8f5e9; padding:20px; border-radius:12px; margin:20px 0;">
        <h3><i class="fas fa-bacterium"></i> ESCAPE AMR – Our Vision</h3>
        <p>We believe the name “ESCAPE” is also a call to action:</p>
        <ul>
            <li><strong>E</strong>veryone must join forces – researchers, clinicians, policymakers.</li>
            <li><strong>S</strong>mart surveillance is our first line of defence – sequence, don’t guess.</li>
            <li><strong>K</strong>nowledge must be shared openly – no paywalls, no silos.</li>
            <li><strong>A</strong>frica bears a heavy AMR burden, but African solutions are emerging.</li>
            <li><strong>P</strong>revention is cheaper than cure – let’s stop resistant infections before they spread.</li>
            <li><strong>E</strong>very day we delay, more lives are at stake – act now.</li>
        </ul>
        </div>

        <div style="text-align:center; margin:40px 0;">
            <i class="fas fa-star" style="font-size:3em; color:#ffc107;"></i>
            <h3>🤝 We Invite You to Contribute!</h3>
            <p><strong>If you find this tool useful, please:</strong></p>
            <div class="action-buttons" style="justify-content:center;">
                <a href="https://github.com/bbeckley-hub/pseudoscope" target="_blank" class="action-btn btn-primary" style="text-decoration:none;"><i class="fab fa-github"></i> ⭐ Star us on GitHub</a>
                <a href="mailto:brownbeckley94@gmail.com" class="action-btn btn-success" style="text-decoration:none;"><i class="fas fa-envelope"></i> Contact the team</a>
                <a href="https://github.com/bbeckley-hub/pseudoscope/issues" target="_blank" class="action-btn btn-warning" style="text-decoration:none;"><i class="fas fa-bug"></i> Report issues</a>
            </div>
            <p style="margin-top:20px;"><i class="fas fa-chalkboard-user"></i> <strong>We welcome collaborations</strong> to adapt this tool for other pathogens and to improve AMR surveillance in Africa and beyond.</p>
            <p><i class="fas fa-hand-holding-heart"></i> If you are a funder or organisation interested in supporting the <strong>ESCAPE AMR</strong> project, please reach out. Together we can build a free, open‑source ecosystem for genomic surveillance of all ESKAPE pathogens.</p>
        </div>
        """

    def _export_section(self):
        return """
        <div class="alert-box alert-info"><i class="fas fa-download"></i><div>Export data tables as CSV or download complete JSON.</div></div>
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('samples-table','sample_overview.csv')">Sample Overview CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('summary-stats','summary_stats.csv')">Summary Stats CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('mlst-counts-table','mlst_counts.csv')">MLST Counts CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('serotype-table','serotype_counts.csv')">Serotype Counts CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('qc-table','fasta_qc.csv')">FASTA QC CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('sto-table','st_o_combinations.csv')">ST‑O Combinations CSV</button>
            <button class="action-btn btn-primary" onclick="exportTableToCSV('cooc-table','gene_cooccurrence.csv')">Co‑occurrence CSV</button>
            <button class="action-btn btn-success" onclick="location.href='genius_pseudomonas_sample_centric_report.json'">Download JSON</button>
        </div>
        """


class GeniusPseudomonasSampleCentricReporter:
    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "GENIUS_PSEUDOMONAS_SAMPLE_CENTRIC_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = PseudomonasHTMLParser()
        self.analyzer = PseudomonasDataAnalyzer()
        self.html_generator = PseudomonasHTMLGenerator(self.analyzer)
        self.metadata = {
            "tool_name": "GENIUS P. aeruginosa Sample-Centric Reporter",
            "version": "1.0.0",
            "author": "Brown Beckley <brownbeckley94@gmail.com>",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir)
        }

    def find_html_files(self) -> Dict[str, List[Path]]:
        print("🔍 Searching for HTML reports...")
        html_files = {
            'mlst': [], 'past': [], 'qc': []
        }
        # Only need HTML for MLST, PAST, and QC. AMR/Abricate/Mutation are TSV-only.
        for html_file in self.input_dir.glob("**/*.html"):
            name = html_file.name.lower()
            if 'mlst' in name:
                html_files['mlst'].append(html_file)
            elif 'past' in name or 'serotype' in name:
                html_files['past'].append(html_file)
            elif 'fasta_qc' in name or 'qc_summary' in name:
                html_files['qc'].append(html_file)
        for k, v in html_files.items():
            if v:
                print(f"  📁 {k.upper()}: {len(v)} files found")
        return html_files

    def integrate_all_data(self, html_files: Dict[str, List[Path]]) -> Dict[str, Any]:
        print("\n🔗 Integrating data (MLST/serotype from HTML; AMR/Virulence/Plasmid/Bacmet/Mutations from TSV strictly)...")
        integrated = {
            'metadata': self.metadata,
            'samples': {},
            'patterns': {},
            'gene_centric': {},
            'qc_data': {},
            'mutation_data': {},
            'amr_details': {},
            'abricate_details': {},
            'mutation_details': {}
        }

        if html_files['qc']:
            integrated['qc_data'] = self.parser.parse_qc_report(html_files['qc'][0])

        mlst_data = {}
        if html_files['mlst']:
            mlst_data = self.parser.parse_mlst_report(html_files['mlst'][0])
        else:
            print("  ⚠️ No MLST HTML found; MLST data will be empty.")

        past_data = {}
        if html_files['past']:
            past_data = self.parser.parse_past_report(html_files['past'][0])
        else:
            print("  ⚠️ No PAST HTML found; serotype data will be empty.")

        # ---- Strict TSV loading for AMRfinder ----
        amr_details, _ = self.parser.load_amrfinder_from_tsv(self.input_dir)
        if amr_details:
            print(f"  ✅ Loaded AMRfinder data for {len(amr_details)} samples from TSV")
        else:
            print("  ⚠️ AMRfinder TSV not found or empty; AMR tab will be empty.")

        # ---- Strict TSV loading for ABRicate ----
        abricate_details, _ = self.parser.load_abricate_from_tsv(self.input_dir)
        if abricate_details:
            print(f"  ✅ Loaded ABRicate data for {len(abricate_details)} samples from TSV")
        else:
            print("  ⚠️ ABRicate TSV not found or empty; virulence/plasmid/bacmet tabs will be empty.")

        # ---- Strict TSV loading for Mutations ----
        mutation_details = self.parser.load_mutations_from_tsv(self.input_dir)
        mutation_data = {}
        if mutation_details:
            print(f"  ✅ Loaded mutation data for {len(mutation_details)} samples from TSV")
            # Convert to gene-centric mutation list for the mutation tab (if needed)
            gene_mutation_dict = {}
            genome_counts = defaultdict(int)
            for sample, muts in mutation_details.items():
                for m in muts:
                    gene = m.get('gene', '')
                    mutation = m.get('mutation', '')
                    if gene and mutation:
                        key = (gene, mutation)
                        if key not in gene_mutation_dict:
                            gene_mutation_dict[key] = {
                                'gene': gene,
                                'mutation': mutation,
                                'class': m.get('class', ''),
                                'subclass': m.get('subclass', ''),
                                'genomes': []
                            }
                        gene_mutation_dict[key]['genomes'].append(sample)
                        genome_counts[sample] += 1
            mutations_list = []
            for (gene, mutation), data in gene_mutation_dict.items():
                data['count'] = len(data['genomes'])
                mutations_list.append(data)
            mutations_list.sort(key=lambda x: x['count'], reverse=True)
            mutation_data = {
                'mutations': mutations_list,
                'genome_mutation_counts': dict(genome_counts)
            }
        else:
            print("  ⚠️ Mutation TSV not found or empty; mutation tab will be empty.")

        # Collect all samples from all sources
        all_samples = set()
        all_samples.update(mlst_data.keys())
        all_samples.update(past_data.keys())
        all_samples.update(amr_details.keys())
        all_samples.update(abricate_details.keys())
        all_samples.update(integrated['qc_data'].keys())

        if not all_samples:
            print("❌ No samples found in any report!")
            return {}

        print(f"📊 Found {len(all_samples)} unique samples")

        # Build sample entries
        for sample in all_samples:
            virulence_genes = []
            for db in ['vfdb', 'ecoli_vf', 'pa_vf', 'pseudomonas_vf']:
                if db in abricate_details.get(sample, {}):
                    for gene_dict in abricate_details[sample][db]:
                        virulence_genes.append(gene_dict.get('gene', 'unknown'))

            integrated['samples'][sample] = {
                'mlst': mlst_data.get(sample, {'ST': 'ND', 'Allele_Profile': 'ND'}),
                'serotype': past_data.get(sample, {'O_Type': 'ND'}),
                'amr_genes': [d.get('gene', 'unknown') for d in amr_details.get(sample, [])],
                'virulence_genes': list(set(virulence_genes))
            }

        integrated['amr_details'] = amr_details
        integrated['abricate_details'] = abricate_details
        integrated['mutation_details'] = mutation_details

        # Build gene_frequencies from the details (counts and genome lists)
        amr_gene_freq = {}
        for sample, genes in amr_details.items():
            for gene_dict in genes:
                gene = gene_dict.get('gene', 'unknown')
                if gene not in amr_gene_freq:
                    amr_gene_freq[gene] = {'count': 0, 'genomes': []}
                if sample not in amr_gene_freq[gene]['genomes']:
                    amr_gene_freq[gene]['genomes'].append(sample)
                    amr_gene_freq[gene]['count'] += 1

        abricate_gene_freq = {}
        for sample, dbs in abricate_details.items():
            for db, genes in dbs.items():
                if db not in abricate_gene_freq:
                    abricate_gene_freq[db] = {}
                for gene_dict in genes:
                    gene = gene_dict.get('gene', 'unknown')
                    if gene not in abricate_gene_freq[db]:
                        abricate_gene_freq[db][gene] = {'count': 0, 'genomes': []}
                    if sample not in abricate_gene_freq[db][gene]['genomes']:
                        abricate_gene_freq[db][gene]['genomes'].append(sample)
                        abricate_gene_freq[db][gene]['count'] += 1

        integrated['gene_frequencies'] = {
            'amrfinder': amr_gene_freq,
            'abricate': abricate_gene_freq
        }

        print("\n🧠 Processing gene-centric analysis...")
        integrated['gene_centric'] = self.analyzer.create_gene_centric_tables(integrated)
        integrated['patterns'] = self.analyzer.create_cross_genome_patterns(integrated)

        integrated['mutation_data'] = mutation_data

        return integrated

    def generate_json_report(self, data: Dict[str, Any]) -> Path:
        print("\n📝 Generating JSON report...")
        out = self.output_dir / "genius_pseudomonas_sample_centric_report.json"
        with open(out, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, default=str)
        print(f"    ✅ JSON saved: {out}")
        return out

    def generate_csv_reports(self, data: Dict[str, Any]):
        print("\n📊 Generating CSV reports...")
        samples_df = pd.DataFrame([{
            'Sample': s,
            'ST': d['mlst']['ST'],
            'O_Type': d['serotype']['O_Type'],
            'Virulence_Count': len(d['virulence_genes'])
        } for s, d in data['samples'].items()])
        samples_df.to_csv(self.output_dir / "sample_overview.csv", index=False)

        amr_rows = []
        for db, genes in data['gene_centric'].get('amr_databases', {}).items():
            for g in genes:
                amr_rows.append({
                    'Gene': g['gene'],
                    'Database': g['database'],
                    'Count': g['count'],
                    'Genomes': ';'.join(g['genomes'])
                })
        if amr_rows:
            pd.DataFrame(amr_rows).to_csv(self.output_dir / "amr_genes.csv", index=False)

        vir_rows = []
        for db, genes in data['gene_centric'].get('virulence_databases', {}).items():
            for g in genes:
                vir_rows.append({
                    'Gene': g['gene'],
                    'Database': g['database'],
                    'Count': g['count'],
                    'Genomes': ';'.join(g['genomes'])
                })
        if vir_rows:
            pd.DataFrame(vir_rows).to_csv(self.output_dir / "virulence_genes.csv", index=False)

        plasmid_rows = []
        for db, genes in data['gene_centric'].get('plasmid_databases', {}).items():
            for g in genes:
                plasmid_rows.append({
                    'Replicon': g['gene'],
                    'Database': g['database'],
                    'Count': g['count'],
                    'Genomes': ';'.join(g['genomes'])
                })
        if plasmid_rows:
            pd.DataFrame(plasmid_rows).to_csv(self.output_dir / "plasmid_replicons.csv", index=False)

        bacmet_rows = []
        for db, genes in data['gene_centric'].get('bacmet_databases', {}).items():
            for g in genes:
                bacmet_rows.append({
                    'Gene': g['gene'],
                    'Database': g['database'],
                    'Count': g['count'],
                    'Genomes': ';'.join(g['genomes'])
                })
        if bacmet_rows:
            pd.DataFrame(bacmet_rows).to_csv(self.output_dir / "bacmet_genes.csv", index=False)

        mutation_data = data.get('mutation_data', {})
        mutations = mutation_data.get('mutations', [])
        if mutations:
            mut_rows = []
            for m in mutations:
                mut_rows.append({
                    'Gene': m['gene'],
                    'Mutation': m['mutation'],
                    'Class': m['class'],
                    'Subclass': m['subclass'],
                    'Count': m['count'],
                    'Genomes': ';'.join(m['genomes'])
                })
            pd.DataFrame(mut_rows).to_csv(self.output_dir / "mutations.csv", index=False)

        cooc = data['patterns'].get('gene_cooccurrence', {})
        cooc_list = []
        for g1, partners in cooc.items():
            for g2, cnt in partners.items():
                if g1 < g2:
                    cooc_list.append((g1, g2, cnt))
        cooc_list.sort(key=lambda x: x[2], reverse=True)
        if cooc_list:
            pd.DataFrame(cooc_list[:100], columns=['Gene1', 'Gene2', 'Count']).to_csv(self.output_dir / "gene_cooccurrence.csv", index=False)

        print("    ✅ CSV reports generated.")

    def run(self):
        print("=" * 80)
        print("🧠 GENIUS P. AERUGINOSA SAMPLE-CENTRIC REPORTER v1.0.0")
        print("=" * 80)
        print(f"📁 Input directory: {self.input_dir}")

        html_files = self.find_html_files()
        data = self.integrate_all_data(html_files)
        if not data:
            print("❌ No data integrated!")
            return False

        print("\n" + "=" * 80)
        print("📊 GENERATING SAMPLE-CENTRIC REPORTS")
        print("=" * 80)
        self.generate_json_report(data)
        self.generate_csv_reports(data)
        self.html_generator.generate_main_report(data, self.output_dir)

        total_samples = len(data['samples'])
        patterns = data['patterns']
        gene_centric = data['gene_centric']
        total_amr = sum(len(genes) for genes in gene_centric.get('amr_databases', {}).values())
        total_vir = sum(len(genes) for genes in gene_centric.get('virulence_databases', {}).values())
        total_plasmid = sum(len(genes) for genes in gene_centric.get('plasmid_databases', {}).values())
        total_bacmet = sum(len(genes) for genes in gene_centric.get('bacmet_databases', {}).values())
        mutation_count = len(data.get('mutation_data', {}).get('mutations', []))

        print("\n" + "=" * 80)
        print("✅ SAMPLE-CENTRIC ANALYSIS COMPLETE!")
        print("=" * 80)
        print(f"📁 Output directory: {self.output_dir}")
        print(f"📄 Files generated:")
        print(f"   • genius_pseudomonas_sample_centric_report.html (Interactive sample boxes)")
        print(f"   • genius_pseudomonas_sample_centric_report.json (Complete data)")
        print(f"   • sample_overview.csv")
        print(f"   • amr_genes.csv")
        print(f"   • virulence_genes.csv")
        print(f"   • plasmid_replicons.csv")
        print(f"   • bacmet_genes.csv")
        if mutation_count > 0:
            print(f"   • mutations.csv")
        print(f"   • gene_cooccurrence.csv")
        print(f"\n📈 ANALYSIS SUMMARY:")
        print(f"   • {total_samples} total samples analyzed")
        print(f"   • {total_amr} AMR genes")
        print(f"   • {total_vir} virulence genes")
        print(f"   • {total_plasmid} plasmid replicons")
        print(f"   • {total_bacmet} Bacmet2 genes")
        print(f"   • {mutation_count} unique point mutations")
        print(f"   • {len(patterns.get('st_o_combinations', {}))} unique ST‑O combinations")
        print("\n🎯 Next steps:")
        print("   1. Open genius_pseudomonas_sample_centric_report.html in your browser")
        print("   2. Explore the sample boxes in AMR, Virulence, Plasmid, Bacmet, and Mutation tabs")
        print("   3. Use the filter controls to focus on specific isolates or databases")
        print("   4. Check ST‑O combinations and gene co‑occurrence in Patterns tab")
        print("   5. Export data using the Export tab or individual CSV buttons")
        print("   6. Use the AI Guide to ask ChatGPT/Claude about your JSON data")
        print("\n" + "=" * 80)
        return True


def main():
    parser = argparse.ArgumentParser(
        description='GENIUS P. aeruginosa Sample-Centric Reporter (v1.0.0) - Hybrid TSV/HTML with Sample Boxes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python p_sample_centric.py -i /path/to/html/reports

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School, Department of Medical Biochemistry
        """
    )
    parser.add_argument('-i', '--input-dir', required=True,
                        help='Directory containing HTML report files and TSV summaries')
    parser.add_argument('-o', '--output-dir',
                        help='Custom output directory (default: input_dir/GENIUS_PSEUDOMONAS_SAMPLE_CENTRIC_REPORTS)')

    args = parser.parse_args()
    input_dir = Path(args.input_dir)

    if not input_dir.exists():
        print(f"❌ Input directory not found: {input_dir}")
        sys.exit(1)

    reporter = GeniusPseudomonasSampleCentricReporter(input_dir)

    if args.output_dir:
        reporter.output_dir = Path(args.output_dir)
        reporter.output_dir.mkdir(parents=True, exist_ok=True)

    success = reporter.run()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()