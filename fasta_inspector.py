# FASTA Inspector: the Python logic on its own. Run with: python fasta_inspector.py
EXAMPLE = """>geneA first sample
ATGGCCTGGAAAGGCTTTTGGATCGATGCCCGGGAAATTTGGG
CCCATGGCCTGA
>geneB copy of geneA
ATGGCCTGGAAAGGCTTTTGGATCGATGCCCGGGAAATTTGGGCCCATGGCCTGA
>geneC has unknown bases
ATGNNNNCCTGGAAAGGCTTTNNGATCG
>geneD too short
ATGG
>geneE has a typing error
ATGGCCTGGAXAGGCTTTTGGATCGATG
>geneF
>geneG first copy
ATATATATGCGCGCGCATATATATGCGC
>geneG second copy with the same ID
GGGCCCGGGCCCGGGCCCATATATGCGC
"""


def parse_fasta(text):
    sequences = {}
    repeated_ids = []
    header = None
    for line in text.split("\n"):
        line = line.strip()
        if line == "":
            continue
        if line.startswith(">"):
            words = line[1:].split()
            if len(words) == 0:
                header = "unnamed"
            else:
                header = words[0]
            if header in sequences:
                repeated_ids.append(header)
            sequences[header] = ""
        elif header is not None:
            sequences[header] += line.upper()
    return sequences, repeated_ids


def gc_content(seq):
    if len(seq) == 0:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq) * 100


def count_invalid(seq):
    invalid = 0
    for base in seq:
        if base not in "ATGCN":
            invalid += 1
    return invalid


def find_duplicates(sequences):
    by_sequence = {}
    for seq_id, seq in sequences.items():
        by_sequence.setdefault(seq, []).append(seq_id)
    return by_sequence


def analyse(sequences, min_length):
    by_sequence = find_duplicates(sequences)
    rows = []
    for seq_id, seq in sequences.items():
        notes = []
        if len(seq) == 0:
            notes.append("empty")
        elif len(seq) < min_length:
            notes.append("short")
        if seq.count("N") > 0:
            notes.append("contains N")
        if count_invalid(seq) > 0:
            notes.append("invalid characters")
        if len(seq) > 0 and len(by_sequence[seq]) > 1:
            notes.append("duplicate sequence")
        if len(notes) == 0:
            notes.append("OK")
        rows.append((seq_id, len(seq), gc_content(seq), ", ".join(notes)))
    return rows


def summarise(sequences):
    lengths = []
    for seq in sequences.values():
        lengths.append(len(seq))

    all_bases = "".join(sequences.values())
    mean_length = sum(lengths) / len(lengths)
    return len(lengths), min(lengths), max(lengths), mean_length, gc_content(all_bases)


# Try it on the example
sequences, repeated_ids = parse_fasta(EXAMPLE)
count, shortest, longest, mean_length, gc = summarise(sequences)
rows = analyse(sequences, 20)
print(count, "sequences. Shortest", shortest, "longest", longest)
print("Mean length {:.1f}, overall GC {:.1f}%".format(mean_length, gc))
print("Repeated IDs:", repeated_ids)
print()
print("ID\tLength\tGC %\tNotes")
for seq_id, length, gc, notes in rows:
    print(seq_id + "\t" + str(length) + "\t" + "{:.1f}".format(gc) + "\t" + notes)
