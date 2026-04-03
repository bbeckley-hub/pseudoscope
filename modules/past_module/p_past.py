#!/usr/bin/env python3
"""
PseudoScope - PAST Module for Pseudomonas aeruginosa serotyping
Author: Brown Beckley <brownbeckley94@gmail.com>
GitHub: bbeckley-hub
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-02-20
"""

import os
import sys
import json
import glob
import argparse
import subprocess
import random
import shutil
import csv
from pathlib import Path
from typing import Dict, List, Any, Optional
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

class PseudoPASTAnalyzer:
    """PAST (Pseudomonas aeruginosa serotyping) module using pasty (camlhmp-blast-regions)."""

    def __init__(self):
        self.pasty_bin = self._find_pasty()
        self.schema = "pasty"
        self.schema_version = "2.2.1"  # from test output
        self.camlhmp_version = "1.1.3"

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

    def _find_pasty(self) -> str:
        """Locate pasty executable in PATH."""
        pasty_path = shutil.which("pasty")
        if not pasty_path:
            raise FileNotFoundError("pasty not found in PATH. Please install pasty (conda install -c bioconda pasty).")
        return pasty_path

    def get_random_quote(self):
        return random.choice(self.science_quotes)

    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find FASTA files from a file, directory, or glob pattern."""
        if os.path.isfile(input_path):
            return [Path(input_path)]

        if os.path.isdir(input_path):
            patterns = ["**/*.fna", "**/*.fasta", "**/*.fa", "**/*.fn"]
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

    def run_pasty(self, input_file: Path, out_dir: Path, prefix: Optional[str] = None) -> Dict[str, Any]:
        """Run pasty on a single FASTA file and return parsed results."""
        if prefix is None:
            prefix = input_file.stem

        print(f"🔬 Processing: {input_file.name}")

        # Create sample output directory
        sample_dir = out_dir / prefix
        sample_dir.mkdir(parents=True, exist_ok=True)

        # Run pasty
        cmd = [
            self.pasty_bin,
            "-i", str(input_file),
            "-o", str(sample_dir),
            "--prefix", prefix,
            "--min-pident", "90",
            "--min-coverage", "90"
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            # Save stdout/stderr for debugging
            with open(sample_dir / "pasty_run.log", 'w') as f:
                f.write(f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}")

            # Parse the main result file (prefix.tsv)
            main_tsv = sample_dir / f"{prefix}.tsv"
            if not main_tsv.exists():
                raise FileNotFoundError(f"Expected main output file not found: {main_tsv}")

            result_data = self._parse_tsv(main_tsv, prefix, input_file.name)

            # Parse details file if exists
            details_tsv = sample_dir / f"{prefix}.details.tsv"
            if details_tsv.exists():
                result_data["details"] = self._parse_details_tsv(details_tsv)
            else:
                result_data["details"] = []

            # Generate reports
            self._generate_reports(result_data, sample_dir)

            print(f"✅ {input_file.name} -> {result_data['type']}")
            return result_data

        except subprocess.CalledProcessError as e:
            print(f"❌ pasty failed for {input_file.name}")
            with open(sample_dir / "pasty_run.log", 'w') as f:
                f.write(f"ERROR: {e.stderr}\n{e.stdout}")
            return self._empty_result(prefix, f"pasty error: {e.stderr}")

    def _parse_tsv(self, tsv_path: Path, sample_name: str, original_filename: str) -> Dict[str, Any]:
        """Parse the main TSV output (single row)."""
        with open(tsv_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            rows = list(reader)
        if not rows:
            return self._empty_result(sample_name, "No data in main TSV")
        row = rows[0]
        # Convert numeric fields
        try:
            coverage = float(row.get("coverage", 0))
        except:
            coverage = 0.0
        try:
            hits = int(row.get("hits", 0))
        except:
            hits = 0
        return {
            "sample": sample_name,
            "original_filename": original_filename,
            "type": row.get("type", "UNKNOWN"),
            "targets": row.get("targets", ""),
            "coverage": coverage,
            "hits": hits,
            "schema": row.get("schema", self.schema),
            "schema_version": row.get("schema_version", self.schema_version),
            "camlhmp_version": row.get("camlhmp_version", self.camlhmp_version),
            "params": row.get("params", ""),
            "comment": row.get("comment", "")
        }

    def _parse_details_tsv(self, tsv_path: Path) -> List[Dict]:
        """Parse the details TSV (multiple rows)."""
        details = []
        with open(tsv_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                # Convert boolean-like status
                status = row.get("status", "False")
                if status.lower() == "true":
                    status_bool = True
                else:
                    status_bool = False
                # Convert numeric fields
                try:
                    coverage = float(row.get("coverage", 0))
                except:
                    coverage = 0.0
                try:
                    hits = int(row.get("hits", 0))
                except:
                    hits = 0
                details.append({
                    "type": row.get("type", ""),
                    "status": status_bool,
                    "targets": row.get("targets", ""),
                    "missing": row.get("missing", ""),
                    "coverage": coverage,
                    "hits": hits,
                    "comment": row.get("comment", "")
                })
        return details

    def _empty_result(self, sample: str, error: str = "") -> Dict[str, Any]:
        return {
            "sample": sample,
            "original_filename": "",
            "type": "UNKNOWN",
            "targets": "",
            "coverage": 0.0,
            "hits": 0,
            "schema": self.schema,
            "schema_version": self.schema_version,
            "camlhmp_version": self.camlhmp_version,
            "params": "",
            "comment": error,
            "details": []
        }

    def _generate_reports(self, data: Dict[str, Any], out_dir: Path):
        """Create HTML, text, TSV, and JSON reports."""
        self._write_text_report(data, out_dir)
        self._write_tsv_report(data, out_dir)
        self._write_json_report(data, out_dir)
        self._write_html_report(data, out_dir)

    def _write_text_report(self, data: Dict, out_dir: Path):
        lines = [
            "PSEUDOSCOPE - PAST (Pseudomonas aeruginosa serotyping) Report",
            "==============================================================",
            f"Sample: {data['sample']}",
            f"Original File: {data.get('original_filename', 'N/A')}",
            f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"PAST Schema: {data['schema']} v{data['schema_version']}",
            "",
            "SEROTYPING RESULT:",
            "------------------",
            f"Predicted Type: {data['type']}",
            f"Targets matched: {data['targets']}",
            f"Coverage: {data['coverage']:.2f}%",
            f"Hits: {data['hits']}",
            f"Parameters: {data['params']}",
            f"Comment: {data.get('comment', '')}",
            "",
            "DETAILED RESULTS (all tested types):",
            "-----------------------------------"
        ]
        if data['details']:
            for d in data['details']:
                lines.append(f"  {d['type']}: {d['status']} (coverage={d['coverage']:.2f}%, hits={d['hits']})")
        else:
            lines.append("  No detailed data available.")

        with open(out_dir / "past_report.txt", 'w') as f:
            f.write('\n'.join(lines))

    def _write_tsv_report(self, data: Dict, out_dir: Path):
        # One-line summary TSV
        tsv_line = f"{data['sample']}\t{data.get('original_filename','')}\t{data['type']}\t{data['targets']}\t{data['coverage']:.2f}\t{data['hits']}\t{data['schema']}\t{data['schema_version']}\t{data['params']}\t{data.get('comment','')}\n"
        with open(out_dir / "past_report.tsv", 'w') as f:
            f.write("Sample\tOriginal_File\tType\tTargets\tCoverage\tHits\tSchema\tSchema_Version\tParams\tComment\n")
            f.write(tsv_line)

    def _write_json_report(self, data: Dict, out_dir: Path):
        with open(out_dir / "past_report.json", 'w') as f:
            json.dump({
                "metadata": {
                    "sample": data['sample'],
                    "original_filename": data.get('original_filename', ''),
                    "analysis_date": datetime.now().isoformat(),
                    "schema": data['schema'],
                    "schema_version": data['schema_version'],
                    "version": "1.0.0",
                    "tool": "PseudoScope PAST"
                },
                "result": {
                    "type": data['type'],
                    "targets": data['targets'],
                    "coverage": data['coverage'],
                    "hits": data['hits'],
                    "params": data['params'],
                    "comment": data.get('comment', '')
                },
                "details": data.get('details', [])
            }, f, indent=2)

    def _write_html_report(self, data: Dict, out_dir: Path):
        """Individual HTML report with ASCII art, quote, and footer."""
        quote = self.get_random_quote()

        # Build details table (collapsible)
        details_rows = ""
        if data['details']:
            for d in data['details']:
                status_icon = "✅" if d.get('status') else "❌"
                details_rows += f'''<tr>
                    <td>{d.get('type','')}</td>
                    <td>{status_icon}</td>
                    <td>{d.get('coverage','')}%</td>
                    <td>{d.get('hits','')}</td>
                </tr>\n'''
        else:
            details_rows = "<tr><td colspan='4' style='text-align:center;'>No detailed data available.</td></tr>"

        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PseudoScope PAST Report - {data['sample']}</title>
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
        .type-badge {{
            font-size: 36px;
            font-weight: bold;
            color: #dc3545;
            text-align: center;
            padding: 20px;
            background: #f8f9fa;
            border-radius: 10px;
            margin: 15px 0;
        }}
        .table-container {{
            max-height: 400px;
            overflow-y: auto;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            margin-top: 15px;
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
            padding: 10px;
            position: sticky;
            top: 0;
        }}
        td, th {{
            padding: 8px 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        tr:hover {{ background-color: #f8f9fa; }}
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
        .collapsible {{
            background-color: #f8f9fa;
            color: #a71d2a;
            cursor: pointer;
            padding: 15px;
            width: 100%;
            border: none;
            text-align: left;
            outline: none;
            font-size: 16px;
            font-weight: bold;
            border-radius: 5px;
            margin: 10px 0;
        }}
        .active, .collapsible:hover {{
            background-color: #e9ecef;
        }}
        .content {{
            padding: 0 18px;
            display: none;
            overflow: hidden;
            background-color: #f9f9f9;
            border-radius: 5px;
        }}
    </style>
    <script>
        // Rotate quotes
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

        // Toggle details
        function toggleDetails() {{
            var content = document.getElementById("detailsContent");
            var btn = document.getElementById("detailsBtn");
            if (content.style.display === "block") {{
                content.style.display = "none";
                btn.innerHTML = "🔽 Show Detailed Results";
            }} else {{
                content.style.display = "block";
                btn.innerHTML = "🔼 Hide Detailed Results";
            }}
        }}
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
            <div class="metric"><span class="label">PAST Schema</span><div class="value">{data['schema']} v{data['schema_version']}</div></div>
        </div>
    </div>
    <div class="card">
        <h2>🧬 Serotyping Result</h2>
        <div class="type-badge">{data['type']}</div>
        <div class="grid">
            <div class="metric"><span class="label">Targets matched</span><div class="value">{data['targets']}</div></div>
            <div class="metric"><span class="label">Coverage</span><div class="value">{data['coverage']:.2f}%</div></div>
            <div class="metric"><span class="label">Hits</span><div class="value">{data['hits']}</div></div>
        </div>
        <p><strong>Parameters:</strong> {data['params']}</p>
        <p><strong>Comment:</strong> {data.get('comment','')}</p>
    </div>
    <button id="detailsBtn" class="collapsible" onclick="toggleDetails()">🔽 Show Detailed Results</button>
    <div id="detailsContent" class="content">
        <div class="card" style="margin-top:10px;">
            <h3>Detailed Results (all tested types)</h3>
            <div class="table-container">
                <table>
                    <thead>
                        <tr><th>Type</th><th>Status</th><th>Coverage (%)</th><th>Hits</th></tr>
                    </thead>
                    <tbody>
                        {details_rows}
                    </tbody>
                </table>
            </div>
        </div>
    </div>
    <div class="footer">
        <p><strong>PSEUDOSCOPE</strong> — P. aeruginosa PAST Module</p>
        <p>Author: <a href="mailto:brownbeckley94@gmail.com">Brown Beckley</a> | GitHub: <a href="https://github.com/bbeckley-hub" target="_blank">bbeckley-hub</a></p>
        <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
        <p class="timestamp">Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
</div>
</body>
</html>'''
        with open(out_dir / "past_report.html", 'w', encoding='utf-8') as f:
            f.write(html)

    def process(self, input_path: str, out_dir: Path, cpus: int = 1):
        """Main entry: process one or multiple files (automatically batch if >1)."""
        fasta_files = self.find_fasta_files(input_path)
        if not fasta_files:
            print("❌ No FASTA files found.")
            return

        out_dir.mkdir(parents=True, exist_ok=True)

        all_results = []
        if len(fasta_files) == 1:
            # Single file, just process it directly (no thread overhead)
            res = self.run_pasty(fasta_files[0], out_dir, fasta_files[0].stem)
            all_results.append(res)
        else:
            print(f"📁 Found {len(fasta_files)} files. Processing...")
            if cpus > 1:
                print(f"🚀 Using {cpus} threads.")
                with ThreadPoolExecutor(max_workers=cpus) as executor:
                    future_to_file = {
                        executor.submit(self.run_pasty, f, out_dir, f.stem): f
                        for f in fasta_files
                    }
                    for future in as_completed(future_to_file):
                        try:
                            res = future.result()
                            all_results.append(res)
                        except Exception as e:
                            print(f"❌ Error processing {future_to_file[future]}: {e}")
            else:
                for f in fasta_files:
                    res = self.run_pasty(f, out_dir, f.stem)
                    all_results.append(res)

        if len(all_results) > 1:
            self._create_batch_summary(all_results, out_dir)

        print(f"\n✅ Analysis complete. Results in: {out_dir}")

    def _create_batch_summary(self, results: List[Dict], out_dir: Path):
        """Generate summary CSV/JSON/HTML for multiple samples."""
        # Prepare summary data as list of dicts
        summary_rows = []
        for r in results:
            summary_rows.append({
                "Sample": r['sample'],
                "Type": r['type'],
                "Targets": r['targets'],
                "Coverage": f"{r['coverage']:.2f}",
                "Hits": r['hits'],
                "Comment": r.get('comment', '')
            })

        # Write CSV
        csv_path = out_dir / "past_summary.csv"
        with open(csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=["Sample", "Type", "Targets", "Coverage", "Hits", "Comment"])
            writer.writeheader()
            writer.writerows(summary_rows)
        print(f"✅ Summary CSV: {csv_path}")

        # Write JSON
        json_path = out_dir / "past_summary.json"
        with open(json_path, 'w') as f:
            json.dump({
                "metadata": {
                    "analysis_date": datetime.now().isoformat(),
                    "schema": self.schema,
                    "schema_version": self.schema_version,
                    "samples_analyzed": len(results)
                },
                "results": results
            }, f, indent=2)
        print(f"✅ Summary JSON: {json_path}")

        # Write HTML
        self._write_batch_html(results, out_dir)

    def _write_batch_html(self, results: List[Dict], out_dir: Path):
        """Batch summary HTML with ASCII art, quotes, search, scrollable table."""
        quote = self.get_random_quote()
        rows = ''
        for i, r in enumerate(results, 1):
            rows += f'''<tr>
                <td>{i}</td>
                <td>{r['sample']}</td>
                <td><strong>{r['type']}</strong></td>
                <td>{r['targets']}</td>
                <td>{r['coverage']:.2f}%</td>
                <td>{r['hits']}</td>
                <td>{r.get('comment', '')}</td>
            </tr>\n'''

        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PseudoScope PAST Batch Summary</title>
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
            const rows = document.querySelectorAll('#pastTable tbody tr');
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
        <div class="stat-card"><div class="value">{len([r for r in results if r['type'] != 'UNKNOWN'])}</div><div class="label">Typed</div></div>
        <div class="stat-card"><div class="value">{len([r for r in results if r['type'] == 'UNKNOWN'])}</div><div class="label">Unknown</div></div>
    </div>
    <div class="controls">
        <input type="text" id="searchInput" placeholder="🔍 Search by sample, type, targets..." onkeyup="filterTable()">
        <button onclick="document.getElementById('searchInput').value=''; filterTable();">Clear</button>
    </div>
    <div class="table-container">
        <table id="pastTable">
            <thead>
                <tr><th>#</th><th>Sample</th><th>Type</th><th>Targets</th><th>Coverage</th><th>Hits</th><th>Comment</th></tr>
            </thead>
            <tbody>
                {rows}
            </tbody>
        </table>
    </div>
    <div class="footer">
        <p><strong>PSEUDOSCOPE</strong> — P. aeruginosa PAST Module</p>
        <p>Author: <a href="mailto:brownbeckley94@gmail.com">Brown Beckley</a> | GitHub: <a href="https://github.com/bbeckley-hub" target="_blank">bbeckley-hub</a></p>
        <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
        <p class="timestamp">Summary generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
</div>
</body>
</html>'''
        with open(out_dir / "past_summary.html", 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"✅ Summary HTML: {out_dir / 'past_summary.html'}")


def main():
    parser = argparse.ArgumentParser(
        description="PseudoScope PAST - Pseudomonas aeruginosa serotyping (using pasty)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single file
  python p_past.py -i sample.fasta -o results/

  # Multiple files with glob pattern 
  python p_past.py -i "*.fna" -o results/

  # Use multiple threads
  python p_past.py -i "genomes/*.fasta" -o results/ --cpus 4
        """
    )
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file, directory, or glob pattern')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('--cpus', type=int, default=1, help='Number of CPU cores for parallel processing (default: 1)')
    args = parser.parse_args()

    try:
        analyzer = PseudoPASTAnalyzer()
    except FileNotFoundError as e:
        print(f"❌ {e}")
        sys.exit(1)

    analyzer.process(args.input, Path(args.output), cpus=args.cpus)


if __name__ == "__main__":
    main()
