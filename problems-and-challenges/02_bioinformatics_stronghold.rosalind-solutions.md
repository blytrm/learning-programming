# Rosalind Problems + Answers

Solutions to bioinformatics problems from [Rosalind](http://rosalind.info/), implemented in both Python and R.

# ROSALIND Solutions – Python & R

## 1. Counting DNA Nucleotides

**Problem**  
Given: A DNA string `s` of length at most 1000 bp.  
Return: Four integers (separated by spaces) counting the respective number of occurrences of `A`, `C`, `G`, and `T` in `s`.

**Python**

```
with open("Counting DNA Nucleotides (1).txt", "r") as file:
    S = file.readline().strip()

bases = ("A", "C", "G", "T")
chr_count = {}

for base in S:
    chr_count[base] = chr_count.get(base, 0) + 1

counts = [chr_count.get(base, 0) for base in bases]
print(" ".join(str(count) for count in counts))
```

**R**

```
txt_path <- "Counting DNA Nucleotides (1).txt"
S <- readLines(txt_path)
bases <- strsplit(S, "")[][1]
count <- table(bases)
print(count)
```

---

## 2. Transcribing DNA into RNA

**Problem**  
Given: A DNA string `t` corresponding to a coding strand.  
Return: The transcribed RNA string `u` formed by replacing all `T` with `U`.

**Python**

```
with open("Transcribing DNA to RNA.txt", "r") as file:
    t = file.readline().strip()

transc = t.replace("T", "U")
print(transc)
```

**R**

```
txt <- "Transcribing DNA to RNA.txt"
t <- readLines(txt)
transc <- gsub(pattern = "T", replacement = "U", x = t)
print(transc)

# or
transc2 <- chartr(old = "T", new = "U", x = t)
print(transc2)
```

---

## 3. Complementing a Strand of DNA

**Problem**  
Given: A DNA string `s`.  
Return: The reverse complement of `s`.

**Python**

```
from Bio.Seq import Seq

with open("DNA Complement Strand.txt", "r") as file:
    s = file.readline().strip()

seq = Seq(s)
rev_comp = seq.reverse_complement()
print(str(rev_comp))
```

**R**

```
txt <- "DNA Complement Strand.txt"
s <- readLines(txt)[1]

rc <- function(x) {
  ss    <- unlist(strsplit(x, NULL))
  ssrev <- rev(ss)

  for (i in seq_along(ssrev)) {
    if (ssrev[i] == "T") {
      ssrev[i] <- "A"
    } else if (ssrev[i] == "A") {
      ssrev[i] <- "T"
    } else if (ssrev[i] == "G") {
      ssrev[i] <- "C"
    } else if (ssrev[i] == "C") {
      ssrev[i] <- "G"
    }
  }

  paste(ssrev, collapse = "")
}

rc_seq <- rc(s)
cat(rc_seq, "\n")
```

---

## 4. Rabbits and Recurrence Relations

**Problem**  
Given: Positive integers `n` and `k`. Start with 1 pair of rabbits; each reproduction‑age pair produces `k` new pairs each month.  
Return: The total number of rabbit pairs after `n` months.

**R**

```
rabbit_pairs_month <- c()
n <- 28  # months
k <- 2   # new pairs per reproducing pair

for (i in 1:n) {
  if (i == 1) {
    rabbit_pairs_month[i] <- 1
  } else if (i == 2) {
    rabbit_pairs_month[i] <- 1
  } else {
    rabbit_pairs_month[i] <- rabbit_pairs_month[i - 1] + k * rabbit_pairs_month[i - 2]
  }
}

rabbit_pairs_month[length(rabbit_pairs_month)]
```

**Python**

```
n, k = 28, 2

def mortal_r(n, k):
    if n == 1 or n == 2:
        return 1

    F_prev2 = 1  # F1
    F_prev1 = 1  # F2

    for _ in range(3, n + 1):
        F_curr = F_prev1 + k * F_prev2
        F_prev2, F_prev1 = F_prev1, F_curr

    return F_curr

print(mortal_r(n, k))


def fib(n, k):
    previous1, previous2 = 1, 1
    for i in range(2, n):
        current = previous1 + k * previous2
        previous2 = previous1
        previous1 = current
    return current
```

---

## 5. GC Content

**Problem**  
Given: At most 10 DNA strings in FASTA format (length ≤ 1 kbp).  
Return: The ID of the string with the highest GC content, followed by its GC percentage (to within 0.001).

**Python**

```
def readFile(filePath):
    with open(filePath, "r") as f:
        return [l.strip() for l in f.readlines()]

ffile = readFile("GC Content Analysis.txt")

fdict = {}
flabel = ""

for line in ffile:
    if ">" in line:
        flabel = line
        fdict[flabel] = ""
    else:
        fdict[flabel] += line

def gc_content(seq):
    return round(((seq.count("C") + seq.count("G")) / len(seq) * 100), 6)

resultsDict = {key: gc_content(value) for (key, value) in fdict.items()}

maxGC_key = max(resultsDict, key=resultsDict.get)
print(maxGC_key)
print(f"{maxGC_key[1:]}\n{resultsDict[maxGC_key]}")
```

**R**

```
file_path <- "GC Content Analysis.txt"
lines <- readLines(file_path)

seq_list  <- list()
current_id <- NULL

for (line in lines) {
  if (startsWith(line, ">")) {
    current_id <- substring(line, 2)
    seq_list[[current_id]] <- ""
  } else {
    seq_list[[current_id]] <- paste0(seq_list[[current_id]], line)
  }
}

gc_content <- function(seq) {
  chars <- strsplit(seq, "")[][1]
  g <- sum(chars == "G")
  c <- sum(chars == "C")
  gc <- (g + c) / length(chars) * 100
  round(gc, 6)
}

gc_values <- sapply(seq_list, gc_content)

max_id <- names(gc_values)[which.max(gc_values)]
max_gc <- gc_values[max_id]

cat(max_id, "\n", max_gc, "\n", sep = "")
```

---

## 6. Counting Point Mutations

**Problem**  
Given: Two DNA strings `s` and `t` of equal length.  
Return: The Hamming distance between them (number of differing positions).

**Python**

```
def readFile(filePath):
    with open(filePath, "r") as file:
        return [l.strip() for l in file.readlines()]

ffile = readFile("Counting Point Mutations (1).txt")
s = ffile
t = ffile[1]

difrncs = 0

if len(s) == len(t):
    for char1, char2 in zip(s, t):
        if char1 != char2:
            difrncs += 1

print(difrncs)
```

**R**

```
file_path <- "Counting Point Mutations (1).txt"
lines <- readLines(file_path)
s <- lines[1]
t <- lines[2]

hamming_distance <- function(a, b) {
  count <- 0
  for (i in seq_len(nchar(a))) {
    if (substr(a, i, i) != substr(b, i, i)) {
      count <- count + 1
    }
  }
  count
}

d <- hamming_distance(s, t)
cat(d, "\n")
```

---

## 7. Mendel’s First Law

**Problem**  
Given: Three positive integers `k`, `m`, `n` representing a population of `k + m + n` organisms:  
- `k` homozygous dominant (`AA`),  
- `m` heterozygous (`Aa`),  
- `n` homozygous recessive (`aa`).  

Return: The probability that two randomly selected organisms will produce an offspring with at least one dominant allele.

**Python**

```
with open("Mendel's First Law (1).txt", "r") as f:
    line = f.readline().strip()

values = line.split()
k, m, n = map(int, values)
N = k + m + n

P_AA_AA = (k / N) * ((k - 1) / (N - 1))
P_AA_Aa = (k / N) * (m / (N - 1)) + (m / N) * (k / (N - 1))
P_AA_aa = (k / N) * (n / (N - 1)) + (n / N) * (k / (N - 1))
P_Aa_Aa = (m / N) * ((m - 1) / (N - 1))
P_Aa_aa = (m / N) * (n / (N - 1)) + (n / N) * (m / (N - 1))
P_aa_aa = (n / N) * ((n - 1) / (N - 1))

Pr_dom = (
    P_AA_AA * 1.0 +
    P_AA_Aa * 1.0 +
    P_AA_aa * 1.0 +
    P_Aa_Aa * 0.75 +
    P_Aa_aa * 0.5 +
    P_aa_aa * 0.0
)

print(f"{Pr_dom:.5f}")
```

**R**

```
line <- readLines("Mendel's First Law (1).txt")
line <- trimws(line)

values <- strsplit(line, " ")[][1]
nums <- as.numeric(values)
k <- nums[1]
m <- nums[2]
n <- nums[3]
N <- k + m + n

P_AA_AA <- (k / N) * ((k - 1) / (N - 1))
P_AA_Aa <- (k / N) * (m / (N - 1)) + (m / N) * (k / (N - 1))
P_AA_aa <- (k / N) * (n / (N - 1)) + (n / N) * (k / (N - 1))
P_Aa_Aa <- (m / N) * ((m - 1) / (N - 1))
P_Aa_aa <- (m / N) * (n / (N - 1)) + (n / N) * (m / (N - 1))
P_aa_aa <- (n / N) * ((n - 1) / (N - 1))

Pr_dom <-
  P_AA_AA * 1.0 +
  P_AA_Aa * 1.0 +
  P_AA_aa * 1.0 +
  P_Aa_Aa * 0.75 +
  P_Aa_aa * 0.5 +
  P_aa_aa * 0.0

cat(sprintf("%.5f\n", Pr_dom))
```

---

## 8. Translating RNA into Protein

**Problem**  
Given: An RNA string `s` corresponding to an mRNA strand.  
Return: The protein string encoded by `s`.

**Python**

```
from Bio.Seq import Seq

def readFile(filePath):
    with open(filePath, "r") as f:
        content = f.readline().strip()
    return content

s = readFile("RNA to Protein (2).txt")

rna = Seq(s)
protein = rna.translate(to_stop=True)
print(protein)
```

**R**

```
library(Biostrings)
library(bioseq)

file_path <- "RNA to Protein.txt"
s <- trimws(readLines(file_path))[1]

# Biostrings
rna_bs <- RNAString(s)
protein_bs <- translate(rna_bs)
print(protein_bs)

# bioseq
rna_bq <- rna(s)
protein_bq <- seq_translate(rna_bq)
print(protein_bq)
```

---

## 9. Finding a Motif in DNA

**Problem**  
Given: Two strings `s` and `t`, with `t` a substring of `s`.  
Return: All locations (1‑based) where `t` appears in `s`.

**Python – simple double loop**

```
def read_file(path):
    with open(path) as f:
        lines = [l.strip() for l in f]
    return lines, lines[1]

s, t = read_file("Finding a Motif in DNA.txt")

n = len(s)
m = len(t)
positions = []

for i in range(n - m + 1):
    has_match = True
    for j in range(m):
        if s[i + j] != t[j]:
            has_match = False
            break
    if has_match:
        positions.append(i + 1)

print(" ".join(map(str, positions)))
```

**R (Biostrings)**

```
library(Biostrings)

lines <- readLines("Finding a Motif in DNA.txt")
s <- lines[1]
t <- lines[2]

dna   <- DNAString(s)
motif <- DNAString(t)

hits <- matchPattern(motif, dna)
pos  <- start(hits)

cat(paste(pos, collapse = " "), "\n")
```

---

## 10. Consensus and Profile

**Problem**  
Given: A collection of DNA strings in FASTA format, all of the same length.  
Return: A consensus string and a profile matrix (counts of A, C, G, T at each position).

**Python**

```
from Bio import motifs
from Bio.Seq import Seq
from Bio import SeqIO
from io import StringIO

def readFile(filePath):
    with open(filePath, "r") as f:
        return [l.strip() for l in f.readlines()]

lines = readFile("Consensus Profile Rosalind.txt")
fasta_content = "\n".join(lines)

seq_strings = []
for record in SeqIO.parse(StringIO(fasta_content), "fasta"):
    seq_strings.append(str(record.seq))

def profile_matrix(sequences):
    seq_length = len(sequences)
    profile = {
        "A":  * seq_length,
        "C":  * seq_length,
        "G":  * seq_length,
        "T":  * seq_length,
    }
    for seq in sequences:
        for i, n in enumerate(seq):
            if n in profile:
                profile[n][i] += 1
    return profile

def consensus(profile):
    cols = len(next(iter(profile.values())))
    con = []
    for i in range(cols):
        col_count = {nt: profile[nt][i] for nt in "ACGT"}
        con.append(max(col_count, key=col_count.get))
    return "".join(con)

profile = profile_matrix(seq_strings)
con = consensus(profile)

print(con)
for nt in "ACGT":
    print(f"{nt}: {' '.join(map(str, profile[nt]))}")

seq_obj = [Seq(s) for s in seq_strings]
m = motifs.create(seq_obj)
print(m.consensus)
print(m.counts)
```

**R**

```
library(Biostrings)
library(seqinr)

fasta_file <- "Consensus Proile Rosalind (2).txt"
seq_list   <- read.fasta(file = fasta_file, seqtype = "DNA", as.string = FALSE)

seq_strings <- vapply(seq_list, paste0, collapse = "", FUN.VALUE = character(1))

prof_matrix <- function(sequences) {
  seq_length <- nchar(sequences)[1]
  profile <- list(
    A = integer(seq_length),
    C = integer(seq_length),
    G = integer(seq_length),
    T = integer(seq_length)
  )
  for (seq in sequences) {
    nts <- strsplit(seq, "")[][1]
    for (i in seq_along(nts)) {
      n <- nts[i]
      if (n %in% names(profile)) {
        profile[[n]][i] <- profile[[n]][i] + 1
      }
    }
  }
  profile
}

profile <- prof_matrix(seq_strings)

consensus_matrix_acgt <- function(profile) {
  cols <- length(profile[])[1]
  consensus <- character(cols)
  for (i in seq_len(cols)) {
    col_counts <- sapply(c("A", "C", "G", "T"), function(nt) profile[[nt]][i])
    consensus[i] <- names(which.max(col_counts))
  }
  paste0(consensus, collapse = "")
}

cons <- consensus_matrix_acgt(profile)
cat(cons, "\n")

for (nt in c("A", "C", "G", "T")) {
  cat(nt, ":", paste(profile[[nt]], collapse = " "), "\n")
}
```
---
