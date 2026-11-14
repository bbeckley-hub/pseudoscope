#!/usr/bin/env python3
"""
PseudoScope - MLST Module for Pseudomonas aeruginosa
Author: Brown Beckley <brownbeckley94@gmail.com>
GitHub: bbeckley-hub
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-10-02

Note on invocation:
    The bundled `bin/mlst` is a Perl script. We invoke it via the `perl`
    interpreter rather than executing the file directly. This makes the
    pipeline robust against zip/tar round-trips that strip the Unix
    executable bit — perl only needs read access. 
"""

import os
import sys
import json
import glob
import argparse
import subprocess
import random
import shutil
import tempfile
import pandas as pd
from pathlib import Path
from typing import Dict, List, Any
from datetime import datetime


class PseudoMLSTAnalyzer:
    """Single-scheme MLST analyzer for P. aeruginosa (Pasteur scheme)."""

    def __init__(self, base_dir: Path):
        self.base_dir = base_dir
        self.mlst_bin = self._find_mlst_binary()
        # Defensive: restore exec bit on anything in bin/ (zip round-trips strip it)
        self._restore_exec_bits(self.base_dir / "bin")
        self.scheme = "paeruginosa"
        self.scheme_display = "PASTEUR"
        self.db_path = base_dir / "db" / "pubmlst" / self.scheme
        self._check_db()

        # ASCII art
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

    # ------------------------------------------------------------------
    # Binary / exec helpers
    # ------------------------------------------------------------------
    def _find_mlst_binary(self) -> Path:
        """Locate mlst binary in bin/ directory or system PATH."""
        local_bin = self.base_dir / "bin" / "mlst"
        if local_bin.exists() and local_bin.is_file():
            return local_bin
        system_bin = shutil.which("mlst")
        if system_bin:
            return Path(system_bin)
        raise FileNotFoundError("mlst binary not found in bin/ or in PATH")

    @staticmethod
    def _restore_exec_bits(root: Path) -> None:
        """Restore the +x bit on anything under `root` that should be runnable.

        zipfile.extractall() and some NFS configurations strip Unix mode bits.
        Perl scripts and shell scripts only need read access when invoked via
        their interpreter, but restoring +x keeps the tree clean and matches
        upstream conventions.
        """
        if not root.exists():
            return
        executable_exts = {".pl", ".sh", ".py"}
        for f in root.rglob("*"):
            if not f.is_file():
                continue
            try:
                needs_chmod = (
                    # Anything directly in bin/ or scripts/ is expected to run
                    f.parent.name in ("bin", "scripts")
                    or f.suffix.lower() in executable_exts
                )
                if needs_chmod:
                    f.chmod(f.stat().st_mode | 0o755)
            except Exception:
                # Non-fatal — we invoke via perl anyway
                pass

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

    # ------------------------------------------------------------------
    # FASTA discovery
    # ------------------------------------------------------------------
    def find_fasta_files(self, input_path: str) -> List[Path]:
        if os.path.isfile(input_path):
            return [Path(input_path)]

        if os.path.isdir(input_path):
            patterns = ["**/*.fna", "**/*.fasta", "**/*.fa", "**/*.fn",
                        "**/*.gb", "**/*.gbk", "**/*.gbff"]
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
        try:
            with open(path, 'r') as f:
                return f.readline().strip().startswith('>')
        except Exception:
            return False

    # ------------------------------------------------------------------
    # Core MLST runner
    # ------------------------------------------------------------------
    def run_mlst(self, input_file: Path, out_dir: Path) -> Dict[str, Any]:
        """Run mlst on a single FASTA file and return parsed results."""
        print(f"🔬 Processing: {input_file.name}")

        sample_dir = out_dir / input_file.stem
        sample_dir.mkdir(parents=True, exist_ok=True)

        # NOTE: invoke via `perl` explicitly. The bundled mlst is a perl script
        # and this makes the pipeline robust against stripped exec bits (zip
        # extraction, NFS, some container filesystems).
        cmd = [
            "perl", str(self.mlst_bin),
            str(input_file),
            "--scheme", self.scheme,
            "--csv",
            "--nopath",
        ]

        # The perl mlst script reads its database location from MLST_DATADIR.
        # Set it explicitly so an updated user-local db is always picked up.
        env = os.environ.copy()
        env["MLST_DATADIR"] = str(self.base_dir / "db")

        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, check=True, env=env,
            )

            mlst_data = self._parse_mlst_output(result.stdout, input_file.name)
            self._generate_reports(mlst_data, sample_dir)

            st_display = mlst_data.get('st', 'UNKNOWN')
            print(f"✅ {input_file.name} -> ST{st_display}")
            return mlst_data

        except FileNotFoundError as e:
            # `perl` itself not on PATH
            print(f"❌ MLST failed for {input_file.name}: perl not found ({e})")
            return self._empty_result(input_file.name)

        except subprocess.CalledProcessError as e:
            print(f"❌ MLST failed for {input_file.name}")
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
        if st_raw in ('-', '', '0'):
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
            "pubmlst_link": (
                f"https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_isolates"
                f"&page=query&scheme_id=1&ST={st}"
                if st.isdigit()
                else "https://pubmlst.org/organisms/pseudomonas-aeruginosa"
            ),
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
            "pubmlst_link": "https://pubmlst.org/organisms/pseudomonas-aeruginosa",
        }

    # ------------------------------------------------------------------
    # Report writers
    # ------------------------------------------------------------------
    def _generate_reports(self, data: Dict[str, Any], out_dir: Path):
        self._write_tsv_report(data, out_dir)
        self._write_html_report(data, out_dir)

    def _write_tsv_report(self, data: Dict, out_dir: Path):
        tsv_line = (
            f"{data['sample']}\t{data['original_filename']}\t"
            f"{data['scheme_display'].lower()}\t{data['scheme']}\t{data['st']}\t"
            f"{'Assigned' if data['mlst_assigned'] else 'Not Assigned'}\t"
            f"{data['confidence']}\t{data['allele_profile']}\n"
        )
        with open(out_dir / "mlst_report.tsv", 'w') as f:
            f.write("Sample\tOriginal_File\tScheme\tMLST_Database\tST\t"
                    "MLST_Status\tConfidence\tAllele_Profile\n")
            f.write(tsv_line)

    def _write_html_report(self, data: Dict, out_dir: Path):
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
            background: linear-gradient(135deg, #2A1408 0%, #5C2E0B 40%, #8B4513 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #fff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .ascii-container {{
            background: rgba(0,0,0,0.75);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            border: 2px solid rgba(210,180,140,0.4);
            overflow-x: auto;
            text-align: center;
            box-shadow: 0 8px 32px rgba(0,0,0,0.4);
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #E8C39E;
            text-shadow: 0 0 10px rgba(232,195,158,0.4);
            display: inline-block;
            text-align: left;
        }}
        .quote-container {{
            background: rgba(255,255,255,0.08);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            border-left: 4px solid #8B4513;
        }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 5px; }}
        .quote-author {{ font-size: 14px; color: #E8C39E; }}
        .card {{
            background: rgba(255,255,255,0.96);
            color: #2c2c2c;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
        }}
        h2 {{ color: #5C2E0B; border-bottom: 3px solid #8B4513; padding-bottom: 10px; margin-bottom: 20px; }}
        h3 {{ color: #5C2E0B; margin: 20px 0 10px; }}
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
        }}
        .metric {{
            background: linear-gradient(135deg, #A0522D 0%, #5C2E0B 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}
        .metric .label {{ font-size: 14px; opacity: 0.9; }}
        .metric .value {{ font-size: 28px; font-weight: bold; word-wrap: break-word; }}
        .profile-box {{
            background: #fdf6ec;
            padding: 15px;
            border-radius: 8px;
            border-left: 4px solid #8B4513;
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
            background: linear-gradient(135deg, #A0522D 0%, #5C2E0B 100%);
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
            background: rgba(0,0,0,0.35);
            border-radius: 10px;
            font-size: 14px;
        }}
        .footer a {{ color: #E8C39E; text-decoration: none; }}
        .footer a:hover {{ text-decoration: underline; }}
        .timestamp {{ color: #E8C39E; }}
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
                a.innerText = `\\u2014 ${{quotes[idx].author}}`;
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
            <div class="metric"><span class="label">Confidence</span><div class="value" style="color: {'#28a745' if data['confidence']=='HIGH' else '#5C2E0B'};">{data['confidence']}</div></div>
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
        <p style="text-align:center; margin-top:20px;"><a href="{data['pubmlst_link']}" target="_blank" style="background:#8B4513; color:white; padding:8px 16px; border-radius:4px; text-decoration:none;">🔗 View on PubMLST</a></p>
    </div>
    <div class="footer">
<p><strong>PSEUDOSCOPE</strong> - MLST Sequence Typing</p>
<p class="timestamp">Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
<div class="authorship">
<p><strong>Technical Support &amp; Inquiries:</strong></p>
<p>Maintainer: Brown Beckley</p>
<p>GitHub: <a href="https://github.com/bbeckley-hub" style="color:#E8C39E;">bbeckley-hub</a></p>
<p>Email: <a href="mailto:brownbeckley94@gmail.com" style="color:#E8C39E;">brownbeckley94@gmail.com</a></p>
<p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
<p style="margin-top:10px;font-size:11px;opacity:0.8;">
Powered by mlst (Seemann) and the PubMLST database
</p>
</div>
</div>
</div>
</body>
</html>'''
        with open(out_dir / "mlst_report.html", 'w', encoding='utf-8') as f:
            f.write(html)

    # ------------------------------------------------------------------
    # Batch processing
    # ------------------------------------------------------------------
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

        if all_results:
            self._create_summary(all_results, out_dir)

        print(f"\n✅ Analysis complete. Results in: {out_dir}")

    def _create_summary(self, results: List[Dict], out_dir: Path):
        df = pd.DataFrame([
            {
                "Sample": r['sample'],
                "ST": r['st'],
                "Confidence": r['confidence'],
                "Allele_Profile": r['allele_profile'],
            }
            for r in results
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
                    "samples_analyzed": len(results),
                },
                "results": results,
            }, f, indent=2)
        print(f"✅ Summary JSON: {json_path}")

        self._write_batch_html(results, out_dir)

    def _write_batch_html(self, results: List[Dict], out_dir: Path):
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
            background: linear-gradient(135deg, #2A1408 0%, #5C2E0B 40%, #8B4513 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #fff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1600px; margin: 0 auto; }}
        .ascii-container {{
            background: rgba(0,0,0,0.75);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            border: 2px solid rgba(210,180,140,0.4);
            overflow-x: auto;
            text-align: center;
            box-shadow: 0 8px 32px rgba(0,0,0,0.4);
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #E8C39E;
            text-shadow: 0 0 10px rgba(232,195,158,0.4);
            display: inline-block;
            text-align: left;
        }}
        .quote-container {{
            background: rgba(255,255,255,0.08);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            border-left: 4px solid #8B4513;
        }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 5px; }}
        .quote-author {{ font-size: 14px; color: #E8C39E; }}
        .stats {{
            display: flex;
            justify-content: space-around;
            flex-wrap: wrap;
            gap: 15px;
            margin-bottom: 20px;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #A0522D 0%, #5C2E0B 100%);
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
            background: #8B4513;
            color: white;
            border: none;
            border-radius: 5px;
            cursor: pointer;
            font-size: 16px;
        }}
        button:hover {{ background: #A0522D; }}
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
            color: #2c2c2c;
            font-size: 14px;
        }}
        th {{
            background: #5C2E0B;
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
        .unknown {{ color: #5C2E0B; font-weight: bold; }}
        .conf-high {{ color: #28a745; }}
        .conf-low {{ color: #5C2E0B; }}
        .profile {{ white-space: pre-wrap; word-wrap: break-word; max-width: 300px; }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0,0,0,0.35);
            border-radius: 10px;
        }}
        .footer a {{ color: #E8C39E; text-decoration: none; }}
        .footer a:hover {{ text-decoration: underline; }}
        .timestamp {{ color: #E8C39E; }}
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
                a.innerText = `\\u2014 ${{quotes[idx].author}}`;
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
<p><strong>PSEUDOSCOPE</strong> - MLST Sequence Typing</p>
<p class="timestamp">Summary generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
<div class="authorship">
<p><strong>Technical Support &amp; Inquiries:</strong></p>
<p>Maintainer: Brown Beckley</p>
<p>GitHub: <a href="https://github.com/bbeckley-hub" style="color:#E8C39E;">bbeckley-hub</a></p>
<p>Email: <a href="mailto:brownbeckley94@gmail.com" style="color:#E8C39E;">brownbeckley94@gmail.com</a></p>
<p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
<p style="margin-top:10px;font-size:11px;opacity:0.8;">
Powered by mlst (Seemann) and the PubMLST database
</p>
</div>
</div>
</div>
</body>
</html>'''
        with open(out_dir / "mlst_summary.html", 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"✅ Summary HTML: {out_dir / 'mlst_summary.html'}")


# ---------------------------------------------------------------------------
# Database update / pull logic
# ---------------------------------------------------------------------------

def update_mlst_database(db_dir: Path, script_dir: Path, scheme: str = "paeruginosa") -> bool:
    """
    Download all MLST schemes via mlstdb, keep only the P. aeruginosa (Pasteur)
    scheme inside db_dir, prune the rest, and rebuild the BLAST index.

    Args:
        db_dir: Target database root; will contain pubmlst/ and blast/.
        script_dir: Directory containing the mlst binary (i.e. <base>/bin).
        scheme: Scheme folder name to keep (default: paeruginosa).

    Returns:
        True on success, False otherwise.
    """
    if not shutil.which('mlstdb'):
        print("❌ 'mlstdb' not found. Install with: pip install mlstdb")
        print("   Then set up API key: mlstdb connect --db pubmlst --api-key")
        return False

    db_dir.mkdir(parents=True, exist_ok=True)

    temp_dir = Path(tempfile.mkdtemp(prefix="pseudoscope_mlst_update_"))
    print(f"📂 Temporary download directory: {temp_dir}")

    try:
        print("🔄 Running mlstdb update (full download, may take 10-20 minutes)...")
        print("   Progress will be shown below:\n")

        result = subprocess.run(
            ['mlstdb', 'update'],
            input='y\n',
            text=True,
            cwd=temp_dir,
            env=os.environ.copy(),
        )

        if result.returncode != 0:
            print(f"\n❌ mlstdb update failed with exit code {result.returncode}")
            return False

        print("\n✅ Download completed.")

        scheme_src = temp_dir / "pubmlst" / scheme
        if not scheme_src.exists():
            print(f"❌ Scheme '{scheme}' not found in downloaded data.")
            return False

        target_pubmlst = db_dir / "pubmlst"
        target_pubmlst.mkdir(parents=True, exist_ok=True)
        target_scheme = target_pubmlst / scheme
        if target_scheme.exists():
            shutil.rmtree(target_scheme)
        shutil.copytree(scheme_src, target_scheme)
        print(f"📁 Copied scheme '{scheme}' to {target_scheme}")

        blast_src = temp_dir / "blast"
        if blast_src.exists():
            target_blast = db_dir / "blast"
            if target_blast.exists():
                shutil.rmtree(target_blast)
            shutil.copytree(blast_src, target_blast)
            print(f"📁 Copied BLAST database to {target_blast}")

        removed = 0
        for item in target_pubmlst.iterdir():
            if item.is_dir() and item.name != scheme:
                shutil.rmtree(item)
                removed += 1
                print(f"   🗑️  Removed scheme: {item.name}")
        if removed:
            print(f"✅ Removed {removed} non-target schemes.")
        else:
            print("✅ No extra schemes to remove.")

        rebuild_script = script_dir.parent / "scripts" / "mlst-make_blast_db"
        if not rebuild_script.exists():
            alt = script_dir / "mlst-make_blast_db"
            if alt.exists():
                rebuild_script = alt

        if rebuild_script.exists():
            print(f"🔄 Rebuilding BLAST database for {scheme}...")
            rebuild_env = os.environ.copy()
            rebuild_env["MLST_DATADIR"] = str(db_dir)
            try:
                subprocess.run(
                    ['perl', str(rebuild_script)],
                    check=True,
                    cwd=str(script_dir.parent),
                    env=rebuild_env,
                )
                print("✅ BLAST database rebuilt successfully.")
            except subprocess.CalledProcessError as e:
                print(f"⚠️  BLAST rebuild failed ({e}); the mlstdb download already")
                print("   produced a blast/ directory, so analysis can still proceed.")
        else:
            print(f"⚠️  mlst-make_blast_db not found at {rebuild_script}")
            print("    The mlstdb download already produced a blast/ directory,")
            print("    so analysis can proceed; rebuild manually if you see errors.")

        print(f"✅ MLST database updated and pruned to '{scheme}' only in {db_dir}")
        return True

    except subprocess.CalledProcessError as e:
        print(f"❌ Update failed: {e}")
        return False
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        return False
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)
        print(f"🧹 Removed temporary directory: {temp_dir}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="PseudoScope MLST for P. aeruginosa (Pasteur scheme)")
    parser.add_argument('-i', '--input', help='Input FASTA file, directory, or glob pattern')
    parser.add_argument('-o', '--output', help='Output directory')
    parser.add_argument('--batch', action='store_true', help='Process multiple files (use with glob/directory)')
    parser.add_argument('--update-db', action='store_true', help='Download and refresh the P. aeruginosa MLST database, then exit')
    args = parser.parse_args()

    base_dir = Path(__file__).parent.resolve()
    db_dir = base_dir / "db"
    bin_dir = base_dir / "bin"

    if args.update_db:
        success = update_mlst_database(db_dir, bin_dir, scheme="paeruginosa")
        sys.exit(0 if success else 1)

    if not args.input or not args.output:
        parser.error("When not using --update-db, -i and -o are required.")

    try:
        analyzer = PseudoMLSTAnalyzer(base_dir)
    except FileNotFoundError as e:
        print(f"❌ {e}")
        sys.exit(1)

    analyzer.process(args.input, Path(args.output), batch=args.batch)


if __name__ == "__main__":
    main()