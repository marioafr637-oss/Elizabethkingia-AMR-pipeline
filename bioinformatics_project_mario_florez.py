#"""Project Introduction: Identification and analysis of antimicrobial resistance (AMR) genes Organism: Elizabethkingia anophelis Field: Red biotechnology / Bioinformatics Author: [Mario Andres Florez Ramirez] Date: [Date] Objective: To identify antimicrobial resistance genes in public genomes of Elizabethkingia anophelis and characterize their AMR profile."""

#Importing libraries
import os
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt

  
BASE_DIR = "."

print(os.path.exists("data/Elizabethkingia_Anophelis_Ag1.fna"))


fasta_path = "data/Elizabethkingia_Anophelis_Ag1.fna"
records = list(SeqIO.parse(fasta_path, "fasta"))
print(len(records))
print(records[0].id)
print(len(records[0]))

fasta_path = "data/Elizabethkingia_Anophelis_Ag1.fna"

BASE_DIR = "."
DATA_DIR = os.path.join(BASE_DIR, "data")
GENOMES_DIR = os.path.join(BASE_DIR, "genomes")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# Project structure
BASE_DIR = "."

DATA_DIR  = os.path.join(BASE_DIR, "data")
GENOMES_DIR = os.path.join(DATA_DIR, "genomes")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

os.makedirs(GENOMES_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs (RESULTS_DIR, exist_ok=True)

GENOMES_DIR, RESULTS_DIR

genome_filename = "Elizabethkingia_anophelis_Ag1.fna"
genome_path = fasta_path # Corrected to use the existing correct path

record = list(SeqIO.parse(genome_path, "fasta"))
 
print("Number of sequences:", len(record))
print("ID:", record[0].id)
print("Genome length:", len(record[0].seq))

genome_stat = {
    "Organism" : ["Elizabethkingia anophelis"],
    "Accesion" : [records [0].id] ,
    "Genome_length_bp" : [len(records[0].seq)]
}

df_genome_stat = pd.DataFrame(genome_stat)
df_genome_stat

sequence = records[0].seq.upper()

counts = {
    "A": sequence.count("A"),
    "T": sequence.count("T"),
    "G": sequence.count("G"),
    "C": sequence.count("C")
}

plt.figure()
plt.bar(counts.keys(), counts.values())
plt.title ("Nucleotide composition of the genome")
plt.xlabel("Nucleotide")
plt.ylabel("Frequency")
plt.show()

#"""DETECTION OF ANTIMICROBIAL RESISTENCE GENES""""

known_amr_genes = [
     "blaB",
     "blaGOB",
     "blaCME",
     "tetX",
     "cat",
     "aac",
     "aad",
     "erm",
]

from Bio.SeqUtils import gc_fraction

# Re-initialize records to ensure it contains SeqRecord objects
records = list(SeqIO.parse(fasta_path, "fasta"))

amr_candidates = []

for record in records:
  seq = record.seq.upper()
  seq_len = len(seq)

  if seq_len > 300:
    gc = gc_fraction(seq) * 100

    amr_candidates.append ({
        "id" : record.id,
        "length": seq_len,
        "gc_content": gc,
    })

# Create the DataFrame outside the loop after all candidates have been collected
amr_df = pd.DataFrame(amr_candidates)
amr_df.head()

filtered_amr = amr_df[(amr_df["length"] >= 250) &
 (amr_df["length"] <= 1200)]
filtered_amr.head()
amr_df.describe()

import subprocess

fasta_path = os.path.join (DATA_DIR, "Elizabethkingia_Anophelis_Ag1.fna")

proteins_path = os.path.join (RESULTS_DIR, "proteins.faa")
genes_path = os.path.join (RESULTS_DIR, "genes.fna")
gff_path = os.path.join (RESULTS_DIR,"genes.gff")
amr_output = os.path.join (RESULTS_DIR, "amrfinder_results.tsv")

print ("Running prodigal...")

subprocess.run([
  "prodigal",
  "-i", fasta_path,
  "-a", proteins_path,
  "-d", genes_path,
  "-f", "gff",
  "-o", gff_path
], check=True)

print ("prodigal finished...")

for path in [genes_path, proteins_path]:
    if not os.path.exists(path):
        raise FileNotFoundError(f"The expected file was not generated: {path}")

print ("Running AMRFinder...")

subprocess.run ([
  "amrfinder",
  "-p", proteins_path,
  "-o", amr_output
], check=True)

print ("AMRFinder finished...")

if not os.path.exists(proteins_path):
  raise FileNotFoundError ("proteins.faa was not generated")

record = list(SeqIO.parse(genes_path, "fasta"))

if not os.path.exists(fasta_path):
  raise FileNotFoundError(f"The file was not found: {fasta_path}")

data = []
for r in record:
  seq = str(r.seq)
  data.append({
      "id": r.id,
      "length": len(seq),
      "GC_content": gc_fraction(seq) * 100
  })

amr_df = pd.DataFrame(data)
amr_df.head()

amr_df.describe()

filtered_genes_for_plot = amr_df[(amr_df["length"] >= 250) & (amr_df["length"] <= 1200)]

plt.figure(figsize=(6, 4))
plt.hist(filtered_genes_for_plot["length"], bins=30)
plt.title("Distribution of candidate lenghts AMR")
plt.xlabel("Length (bp)")
plt.ylabel("Frequency")
plt.show()

GENOME = fasta_path

subprocess.run (["amrfinder","--version"], check=True)

subprocess.run([
"amrfinder",
"-p", proteins_path,
"--threads", "2",
"--output", amr_output
], check=True)

df = pd.read_csv (amr_output, sep= "\t")
print(df.head())