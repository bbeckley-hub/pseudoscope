from setuptools import setup

setup(
    name="pseudoscope",
    version="1.0.0",
    author="Brown Beckley & Riccardo Bolini",
    author_email="brownbeckley94@gmail.com",
    description="Advanced P. aeruginosa Genomic Analysis Platform",
    long_description="""Complete P. aeruginosa genomic analysis pipeline: MLST | K/O Locus | AMR | Plasmid | Virulence | Quality Control | Summary Reports

Critical Genes Tracked:
🔴 Carbapenemases (KPC, NDM, OXA-48)
🟠 Colistin (mcr)
🟡 Tigecycline (tetX)
🟢 ICEKp Markers (ybt, clb, iro, rmp)
🔵 Virulence Plasmid Markers (iro, iuc, rmp, rmpA2)
🟣 Biocides & Heavy Metals (qac, sil, mer, ars, pco)
⚪ Adhesins (fim, mrk, ecp)
🔶 Secretion Systems (tss)
💧 Siderophores
🧪 Toxins""",
    long_description_content_type="text/markdown",
    url="https://github.com/bbeckley-hub/pseudoscope",
    packages=['kleboscope'],
    include_package_data=True,
    python_requires=">=3.9",
    install_requires=[
        "pandas>=1.5.0",
        "biopython>=1.85",  
        "psutil>=5.9.0",
        "requests>=2.28.0",
        "tqdm>=4.64.0",
        "click>=8.0.0",
        "kaptive>=3.1.0", 
        "beautifulsoup4>=4.11.0",
        "lxml>=4.9.0",
        "matplotlib>=3.5.0",
        "seaborn>=0.12.0",
    ],
    entry_points={
        "console_scripts": [
            "pseudoscope=pseudoscope.pseudoscope:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
    ],
)
