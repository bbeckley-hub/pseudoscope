#!/usr/bin/env python3
"""
PseudoScope - MLST Module for Pseudomonas aeruginosa (single scheme)
Author: Brown Beckley <brownbeckley94@gmail.com>
GitHub: bbeckley-hub
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2025-12-29
"""

import os
import sys
import json
import glob
import argparse
import subprocess
import random
import shutil
import pandas as pd  
from pathlib import Path
from typing import Dict, List, Any
from datetime import datetime

class PseudoMLSTAnalyzer:
    """Single‑scheme MLST analyzer for P. aeruginosa (paeruginosa)."""

    def __init__(self, base_dir: Path):
        self.base_dir = base_dir
        self.mlst_bin = self._find_mlst_binary()
        self.scheme = "paeruginosa"
        self.scheme_display = "OXFORD"
        self.db_path = base_dir / "db" / "pubmlst" / self.scheme
        self._check_db()

        # ASCII art for reports
        self.ascii_art = r"""
██████╗ ███████╗███████╗██╗   ██╗██████╗  ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██╔════╝██║   ██║██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
██████╔╝███████╗█████╗  ██║   ██║██║  ██║██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔═══╝ ╚════██║██╔══╝  ██║   ██║██║  ██║██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║     ███████║███████╗╚██████╔╝██████╔╝╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝     ╚══════╝╚══════╝ ╚═════╝ ╚═════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
"""

        # Science quotes
        self.science_quotes = [
            {"text": "Pseudomonas aeruginosa: the master of adaptation.", "author": "Unknown"},
            {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"},
            {"text": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"},
            {"text": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"},
            {"text": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie"},
            {"text": "Research is what I'm doing when I don't know what I'm doing.", "author": "Wernher von Braun"},
            {"text": "The science of today is the technology of tomorrow.", "author": "Edward Teller"},
        ]

    def _find_mlst_binary(self) -> Path:
        """Locate mlst binary in bin/ directory or system PATH."""
        local_bin = self.base_dir / "bin" / "mlst"
        if local_bin.exists() and local_bin.is_file():
            return local_bin
        system_bin = shutil.which("mlst")
        if system_bin:
            return Path(system_bin)
        raise FileNotFoundError("mlst binary not found in bin/ or in PATH")

    def _check_db(self):
        """Verify that the scheme database exists."""
        if not self.db_path.exists():
            raise FileNotFoundError(f"MLST scheme database not found: {self.db_path}")
        tfa_files = list(self.db_path.glob("*.tfa"))
        if not tfa_files:
            raise FileNotFoundError(f"No allele files (*.tfa) found in {self.db_path}")
        print(f"✅ Found MLST scheme '{self.scheme}' with {len(tfa_files)} loci")

    def get_random_quote(self):
        return random.choice(self.science_quotes)

    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find FASTA files from a file, directory, or glob pattern."""
        if os.path.isfile(input_path):
            return [Path(input_path)]

        if os.path.isdir(input_path):
            patterns = ["**/*.fna", "**/*.fasta", "**/*.fa", "**/*.fn", "**/*.gb", "**/*.gbk", "**/*.gbff"]
            files = []
            for pat in patterns:
                files.extend(glob.glob(os.path.join(input_path, pat), recursive=True))
        else:
            files = glob.glob(input_path, recursive=True)

        fasta_files = []
        for f in files:
            path = Path(f)
            if path.is_file() and self._is_fasta(path):
                fasta_files.append(path)
        return sorted(set(fasta_files))

    def _is_fasta(self, path: Path) -> bool:
        """Quick check for FASTA format."""
        try:
            with open(path, 'r') as f:
                first = f.readline().strip()
                return first.startswith('>')
        except:
            return False

    def run_mlst(self, input_file: Path, out_dir: Path) -> Dict[str, Any]:
        """Run mlst on a single FASTA file and return parsed results."""
        print(f"🔬 Processing: {input_file.name}")

        sample_dir = out_dir / input_file.stem
        sample_dir.mkdir(parents=True, exist_ok=True)

        raw_out = sample_dir / "mlst_raw_output.txt"

        cmd = [
            str(self.mlst_bin),
            str(input_file),
            "--scheme", self.scheme,
            "--csv",
            "--nopath"
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            with open(raw_out, 'w') as f:
                f.write(result.stdout)

            mlst_data = self._parse_mlst_output(result.stdout, input_file.name)
            self._generate_reports(mlst_data, sample_dir)

            st_display = mlst_data.get('st', 'UNKNOWN')
            print(f"✅ {input_file.name} -> ST{st_display}")
            return mlst_data

        except subprocess.CalledProcessError as e:
            print(f"❌ MLST failed for {input_file.name}")
            with open(raw_out, 'w') as f:
                f.write(f"ERROR: {e.stderr}\n{e.stdout}")
            return self._empty_result(input_file.name)

    def _parse_mlst_output(self, stdout: str, filename: str) -> Dict[str, Any]:
        """Parse mlst --csv output (single line after header)."""
        lines = [l.strip() for l in stdout.strip().split('\n') if l.strip()]
        if not lines:
            return self._empty_result(filename)

        for line in reversed(lines):
            if ',' in line and not line.startswith('file'):
                parts = [p.strip() for p in line.split(',')]
                break
        else:
            return self._empty_result(filename)

        if len(parts) < 3:
            return self._empty_result(filename)

        st_raw = parts[2]
        if st_raw == '-' or st_raw == '' or st_raw == '0':
            st = "UNKNOWN"
        else:
            st = st_raw.lstrip('ST')

        alleles = {}
        allele_parts = []
        for item in parts[3:]:
            if '(' in item and ')' in item:
                gene, allele = item.split('(')[0], item.split('(')[1].rstrip(')')
                alleles[gene] = allele
                allele_parts.append(f"{gene}({allele})")
            else:
                alleles[item] = "?"
                allele_parts.append(f"{item}(?)")

        return {
            "sample": Path(filename).stem,
            "original_filename": filename,
            "st": st,
            "scheme": self.scheme,
            "scheme_display": self.scheme_display,
            "alleles": alleles,
            "allele_profile": '-'.join(allele_parts) if allele_parts else "",
            "confidence": "HIGH" if st != "UNKNOWN" else "LOW",
            "mlst_assigned": st != "UNKNOWN",
            "pubmlst_link": f"https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_isolates&page=query&scheme_id=1&ST={st}" if st.isdigit() else "https://pubmlst.org/organisms/pseudomonas-aeruginosa"
        }

    def _empty_result(self, filename: str) -> Dict[str, Any]:
        return {
            "sample": Path(filename).stem,
            "original_filename": filename,
            "st": "UNKNOWN",
            "scheme": self.scheme,
            "scheme_display": self.scheme_display,
            "alleles": {},
            "allele_profile": "",
            "confidence": "LOW",
            "mlst_assigned": False,
            "pubmlst_link": "https://pubmlst.org/organisms/pseudomonas-aeruginosa"
        }

    def _generate_reports(self, data: Dict[str, Any], out_dir: Path):
        """Create HTML, text, TSV, and JSON reports."""
        self._write_text_report(data, out_dir)
        self._write_tsv_report(data, out_dir)
        self._write_json_report(data, out_dir)
        self._write_html_report(data, out_dir)

    def _write_text_report(self, data: Dict, out_dir: Path):
        lines = [
            "PSEUDOSCOPE - MLST Analysis Report",
            "======================================",
            f"Sample: {data['sample']}",
            f"Original File: {data['original_filename']}",
            f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"MLST Scheme: {data['scheme_display']}",
            "",
            "MLST TYPING RESULTS:",
            "-------------------",
            f"Sequence Type (ST): {data['st']}",
            f"Confidence: {data['confidence']}",
            f"MLST Status: {'Assigned' if data['mlst_assigned'] else 'Not Assigned'}",
            "",
            f"Allele Profile ({data['scheme_display']} scheme):",
            data['allele_profile'],
            "",
            "Detailed Alleles:",
        ]
        if data['alleles']:
            for g, a in data['alleles'].items():
                lines.append(f"  {g}: {a}")
        else:
            lines.append("  No alleles detected")
        lines.append("")
        lines.append(f"PubMLST Link: {data['pubmlst_link']}")

        with open(out_dir / "mlst_report.txt", 'w') as f:
            f.write('\n'.join(lines))

    def _write_tsv_report(self, data: Dict, out_dir: Path):
        tsv_line = f"{data['sample']}\t{data['original_filename']}\t{data['scheme_display'].lower()}\t{data['scheme']}\t{data['st']}\t{'Assigned' if data['mlst_assigned'] else 'Not Assigned'}\t{data['confidence']}\t{data['allele_profile']}\n"
        with open(out_dir / "mlst_report.tsv", 'w') as f:
            f.write("Sample\tOriginal_File\tScheme\tMLST_Database\tST\tMLST_Status\tConfidence\tAllele_Profile\n")
            f.write(tsv_line)

    def _write_json_report(self, data: Dict, out_dir: Path):
        with open(out_dir / "mlst_report.json", 'w') as f:
            json.dump({
                "metadata": {
                    "sample": data['sample'],
                    "original_filename": data['original_filename'],
                    "analysis_date": datetime.now().isoformat(),
                    "scheme": data['scheme'],
                    "scheme_display": data['scheme_display'],
                    "version": "1.0.0",
                    "tool": "PseudoScope"
                },
                "mlst_results": {
                    "sequence_type": data['st'],
                    "confidence": data['confidence'],
                    "mlst_assigned": data['mlst_assigned'],
                    "allele_profile": data['allele_profile'],
                    "alleles": data['alleles']
                }
            }, f, indent=2)

    def _write_html_report(self, data: Dict, out_dir: Path):
        """Individual HTML report with ASCII art, quote, and footer."""
        quote = self.get_random_quote()
        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PseudoScope MLST Report - {data['sample']}</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #2c0b0e 0%, #4a1c20 50%, #6b2b2f 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #fff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
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
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 5px; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; }}
        .card {{
            background: rgba(255,255,255,0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
        }}
        h2 {{ color: #a71d2a; border-bottom: 3px solid #dc3545; padding-bottom: 10px; margin-bottom: 20px; }}
        h3 {{ color: #a71d2a; margin: 20px 0 10px; }}
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
        }}
        .metric {{
            background: linear-gradient(135deg, #dc3545 0%, #a71d2a 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}
        .metric .label {{ font-size: 14px; opacity: 0.9; }}
        .metric .value {{ font-size: 28px; font-weight: bold; word-wrap: break-word; }}
        .profile-box {{
            background: #f8fafc;
            padding: 15px;
            border-radius: 8px;
            border-left: 4px solid #dc3545;
            font-family: monospace;
            white-space: pre-wrap;
            word-wrap: break-word;
        }}
        .allele-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
            gap: 10px;
            margin-top: 15px;
        }}
        .allele-item {{
            background: linear-gradient(135deg, #dc3545 0%, #a71d2a 100%);
            color: white;
            padding: 12px;
            border-radius: 6px;
            text-align: center;
            word-wrap: break-word;
        }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0,0,0,0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        .footer a {{ color: #dc3545; text-decoration: none; }}
        .footer a:hover {{ text-decoration: underline; }}
        .timestamp {{ color: #fbbf24; }}
    </style>
    <script>
        const quotes = {json.dumps(self.science_quotes)};
        let idx = 0;
        function rotateQuote() {{
            const q = document.getElementById('quote-text');
            const a = document.getElementById('quote-author');
            if (q && a) {{
                idx = (idx + 1) % quotes.length;
                q.innerText = `"${{quotes[idx].text}}"`;
                a.innerText = `— ${{quotes[idx].author}}`;
            }}
        }}
        setInterval(rotateQuote, 10000);
    </script>
</head>
<body>
<div class="container">
    <div class="ascii-container"><div class="ascii-art">{self.ascii_art}</div></div>
    <div class="quote-container">
        <div id="quote-text" class="quote-text">"{quote['text']}"</div>
        <div id="quote-author" class="quote-author">— {quote['author']}</div>
    </div>
    <div class="card">
        <h2>📊 Sample Information</h2>
        <div class="grid">
            <div class="metric"><span class="label">Sample</span><div class="value">{data['sample']}</div></div>
            <div class="metric"><span class="label">Analysis Date</span><div class="value">{datetime.now().strftime('%Y-%m-%d')}</div></div>
            <div class="metric"><span class="label">MLST Scheme</span><div class="value">{data['scheme_display']}</div></div>
        </div>
    </div>
    <div class="card">
        <h2>🧬 MLST Results</h2>
        <div class="grid">
            <div class="metric"><span class="label">Sequence Type (ST)</span><div class="value">ST{data['st']}</div></div>
            <div class="metric"><span class="label">Confidence</span><div class="value" style="color: {'#28a745' if data['confidence']=='HIGH' else '#dc3545'};">{data['confidence']}</div></div>
        </div>
        <h3>Allele Profile</h3>
        <div class="profile-box">{data['allele_profile']}</div>
        <h3>Alleles</h3>
        <div class="allele-grid">
'''
        if data['alleles']:
            for g, a in data['alleles'].items():
                html += f'            <div class="allele-item"><div style="font-size:14px;">{g}</div><div style="font-size:20px; font-weight:bold;">{a}</div></div>\n'
        else:
            html += '            <div class="allele-item" style="grid-column:1/-1;">No alleles detected</div>\n'
        html += f'''        </div>
        <p style="text-align:center; margin-top:20px;"><a href="{data['pubmlst_link']}" target="_blank" style="background:#dc3545; color:white; padding:8px 16px; border-radius:4px; text-decoration:none;">🔗 View on PubMLST</a></p>
    </div>
    <div class="footer">
        <p><strong>PSEUDOSCOPE</strong> — P. aeruginosa MLST Module</p>
        <p>Author: <a href="mailto:brownbeckley94@gmail.com">Brown Beckley</a> | GitHub: <a href="https://github.com/bbeckley-hub" target="_blank">bbeckley-hub</a></p>
        <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
        <p class="timestamp">Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
</div>
</body>
</html>'''
        with open(out_dir / "mlst_report.html", 'w', encoding='utf-8') as f:
            f.write(html)

    def process(self, input_path: str, out_dir: Path, batch: bool = False):
        fasta_files = self.find_fasta_files(input_path)
        if not fasta_files:
            print("❌ No FASTA files found.")
            return

        if not batch and len(fasta_files) > 1:
            print(f"⚠️ Found {len(fasta_files)} files. Use --batch to process all, or specify a single file.")
            return

        out_dir.mkdir(parents=True, exist_ok=True)

        all_results = []
        for f in fasta_files:
            res = self.run_mlst(f, out_dir)
            all_results.append(res)

        # Always create summary files (even for a single sample)
        if all_results:
            self._create_summary(all_results, out_dir)

        print(f"\n✅ Analysis complete. Results in: {out_dir}")

    def _create_summary(self, results: List[Dict], out_dir: Path):
        """Generate summary CSV/JSON/HTML for any number of samples (1 or more)."""
        df = pd.DataFrame([
            {
                "Sample": r['sample'],
                "ST": r['st'],
                "Confidence": r['confidence'],
                "Allele_Profile": r['allele_profile']
            } for r in results
        ])

        csv_path = out_dir / "mlst_summary.csv"
        df.to_csv(csv_path, index=False)
        print(f"✅ Summary CSV: {csv_path}")

        json_path = out_dir / "mlst_summary.json"
        with open(json_path, 'w') as f:
            json.dump({
                "metadata": {
                    "analysis_date": datetime.now().isoformat(),
                    "scheme": self.scheme,
                    "samples_analyzed": len(results)
                },
                "results": results
            }, f, indent=2)
        print(f"✅ Summary JSON: {json_path}")

        self._write_batch_html(results, out_dir)

    def _write_batch_html(self, results: List[Dict], out_dir: Path):
        """Batch/single summary HTML with ASCII art, quotes, search, scrollable table."""
        quote = self.get_random_quote()
        rows = ''
        for i, r in enumerate(results, 1):
            rows += f'''<tr>
                <td>{i}</td>
                <td>{r['sample']}</td>
                <td>ST{r['st']}</td>
                <td class="{'assigned' if r['st']!='UNKNOWN' else 'unknown'}">{r['st']}</td>
                <td class="conf-{r['confidence'].lower()}">{r['confidence']}</td>
                <td><code class="profile">{r['allele_profile']}</code></td>
            </tr>\n'''

        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PseudoScope MLST Summary</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #2c0b0e 0%, #4a1c20 50%, #6b2b2f 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #fff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1600px; margin: 0 auto; }}
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
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 5px; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; }}
        .stats {{
            display: flex;
            justify-content: space-around;
            flex-wrap: wrap;
            gap: 15px;
            margin-bottom: 20px;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #dc3545 0%, #a71d2a 100%);
            padding: 20px;
            border-radius: 10px;
            text-align: center;
            flex: 1 1 200px;
        }}
        .stat-card .value {{ font-size: 32px; font-weight: bold; }}
        .stat-card .label {{ font-size: 14px; opacity: 0.9; }}
        .controls {{
            display: flex;
            gap: 10px;
            margin: 20px 0;
            flex-wrap: wrap;
        }}
        #searchInput {{
            flex: 1;
            padding: 10px;
            border: none;
            border-radius: 5px;
            font-size: 16px;
        }}
        button {{
            padding: 10px 20px;
            background: #dc3545;
            color: white;
            border: none;
            border-radius: 5px;
            cursor: pointer;
            font-size: 16px;
        }}
        button:hover {{ background: #c82333; }}
        .table-container {{
            max-height: 500px;
            overflow-y: auto;
            border-radius: 10px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.3);
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
            color: #1f2937;
            font-size: 14px;
        }}
        th {{
            background: #dc3545;
            color: white;
            padding: 12px;
            position: sticky;
            top: 0;
            z-index: 10;
        }}
        td, th {{
            padding: 10px 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        tr:hover {{ background-color: #f8f9fa; }}
        .assigned {{ color: #28a745; font-weight: bold; }}
        .unknown {{ color: #dc3545; font-weight: bold; }}
        .conf-high {{ color: #28a745; }}
        .conf-low {{ color: #dc3545; }}
        .profile {{ white-space: pre-wrap; word-wrap: break-word; max-width: 300px; }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0,0,0,0.3);
            border-radius: 10px;
        }}
        .footer a {{ color: #dc3545; text-decoration: none; }}
        .footer a:hover {{ text-decoration: underline; }}
        .timestamp {{ color: #fbbf24; }}
    </style>
    <script>
        function filterTable() {{
            const input = document.getElementById('searchInput');
            const filter = input.value.toLowerCase();
            const rows = document.querySelectorAll('#mlstTable tbody tr');
            rows.forEach(row => {{
                const text = row.textContent.toLowerCase();
                row.style.display = text.includes(filter) ? '' : 'none';
            }});
        }}
        const quotes = {json.dumps(self.science_quotes)};
        let idx = 0;
        function rotateQuote() {{
            const q = document.getElementById('quote-text');
            const a = document.getElementById('quote-author');
            if (q && a) {{
                idx = (idx + 1) % quotes.length;
                q.innerText = `"${{quotes[idx].text}}"`;
                a.innerText = `— ${{quotes[idx].author}}`;
            }}
        }}
        setInterval(rotateQuote, 10000);
    </script>
</head>
<body>
<div class="container">
    <div class="ascii-container"><div class="ascii-art">{self.ascii_art}</div></div>
    <div class="quote-container">
        <div id="quote-text" class="quote-text">"{quote['text']}"</div>
        <div id="quote-author" class="quote-author">— {quote['author']}</div>
    </div>
    <div class="stats">
        <div class="stat-card"><div class="value">{len(results)}</div><div class="label">Samples</div></div>
        <div class="stat-card"><div class="value">{sum(1 for r in results if r['st']!='UNKNOWN')}</div><div class="label">ST Assigned</div></div>
        <div class="stat-card"><div class="value">{sum(1 for r in results if r['st']=='UNKNOWN')}</div><div class="label">Unknown ST</div></div>
    </div>
    <div class="controls">
        <input type="text" id="searchInput" placeholder="🔍 Search by sample, ST, allele profile..." onkeyup="filterTable()">
        <button onclick="document.getElementById('searchInput').value=''; filterTable();">Clear</button>
    </div>
    <div class="table-container">
        <table id="mlstTable">
            <thead>
                <tr><th>#</th><th>Sample</th><th>ST</th><th>Assigned</th><th>Confidence</th><th>Allele Profile</th></tr>
            </thead>
            <tbody>
                {rows}
            </tbody>
        </table>
    </div>
    <div class="footer">
        <p><strong>PSEUDOSCOPE</strong> — P. aeruginosa MLST Module</p>
        <p>Author: <a href="mailto:brownbeckley94@gmail.com">Brown Beckley</a> | GitHub: <a href="https://github.com/bbeckley-hub" target="_blank">bbeckley-hub</a></p>
        <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
        <p class="timestamp">Summary generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
</div>
</body>
</html>'''
        with open(out_dir / "mlst_summary.html", 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"✅ Summary HTML: {out_dir / 'mlst_summary.html'}")


def main():
    parser = argparse.ArgumentParser(description="PseudoScope MLST for P. aeruginosa (single scheme)")
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file, directory, or glob pattern')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('--batch', action='store_true', help='Process multiple files (use with glob/directory)')
    args = parser.parse_args()

    base_dir = Path(__file__).parent.resolve()
    try:
        analyzer = PseudoMLSTAnalyzer(base_dir)
    except FileNotFoundError as e:
        print(f"❌ {e}")
        sys.exit(1)

    analyzer.process(args.input, Path(args.output), batch=args.batch)

if __name__ == "__main__":
    main()