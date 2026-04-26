from setuptools import setup, find_packages
from pathlib import Path

# Read long description from README if exists
readme_path = Path(__file__).parent / "README.md"
long_description = readme_path.read_text(encoding="utf-8") if readme_path.exists() else ""

setup(
    name="pseudoscope",
    version="1.0.0",
    author="Brown Beckley",
    author_email="brownbeckley94@gmail.com",
    description="Advanced P. aeruginosa Genomic Analysis Platform",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bbeckley-hub/pseudoscope",
    packages=find_packages(include=["pseudoscope", "pseudoscope.*"]),
    include_package_data=True,
    package_data={
        "pseudoscope": [
            "modules/*/*.py",
            "modules/*/*.sh",
            "modules/*/data/**/*",
            "modules/*/db/**/*",
            "modules/*/bin/**/*",
        ],
    },
    python_requires=">=3.9",
    install_requires=[
        "pandas>=1.5.0",
        "biopython>=1.85",
        "psutil>=5.9.0",
        "requests>=2.28.0",
        "tqdm>=4.64.0",
        "click>=8.0.0",
        "beautifulsoup4>=4.11.0",
        "lxml>=4.9.0",
        "matplotlib>=3.5.0",
        "seaborn>=0.12.0",
        "networkx>=2.8",
        "plotly>=5.10.0",
        "scikit-learn>=1.0.0",
        "scipy>=1.9.0",
    ],
    entry_points={
        "console_scripts": [
            "pseudoscope = pseudoscope.pseudoscope:main",
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