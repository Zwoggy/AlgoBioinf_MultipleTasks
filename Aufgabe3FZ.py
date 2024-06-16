
"""
Algorithmische Bioinformatik Aufgabe 3:
Erstellen Sie ein Programm, mit dem in t DNA-Sequenzen mit jeweils N Nukleotiden das beste Motiv  der Länge k=10 gefunden wird.

Author: Florian Zwicker

"""

def calculate_score(seqs):
    """
    Berechnet den Score für eine Liste von Sequenzen.

    Args:
        seqs (list): Liste von Sequenzen als Strings.

    Returns:
        int: Der berechnete Score.
    """
    score_value = 0

    for pos in range(len(seqs[0])):
        nucleotide_counts = {"a": 0, "c": 0, "t": 0, "g": 0}

        for seq in seqs:
            nucleotide_counts[seq[pos]] += 1

        score_value += max(nucleotide_counts.values())

    return score_value

def generate_subsequences(seqs, k):
    """
    Generiert Teilsequenzen einer bestimmten Länge k aus einer Liste von Sequenzen.

    Args:
        seqs (list): Liste von Sequenzen als Strings.
        k (int): Länge der Teilsequenzen.

    Returns:
        list: Liste der generierten Teilsequenzen.
    """
    subsequences = set()

    for seq in seqs:
        for i in range(len(seq) - k + 1):
            subsequences.add(seq[i:i+k])

    return list(subsequences)

def find_regulatory_motif(dna, subsequences):
    """
    Findet das regulatorische Motiv in einer Liste von DNA-Sequenzen.

    Args:
        dna (list): Liste von DNA-Sequenzen als Strings.
        subsequences (list): Liste der generierten Teilsequenzen.

    Prints:
        Der höchste erreichte Score und das regulatorische Motiv.
    """
    max_scores = []
    try:
        for subseq in subsequences:
            total_score = 0

            for seq in dna:
                best_score = 0

                for i in range(len(seq) - len(subseq) + 1):
                    current_seqs = [subseq, seq[i:i+len(subseq)]]
                    current_score = calculate_score(current_seqs)

                    if current_score > best_score:
                        best_score = current_score

                total_score += best_score - 10

            max_scores.append(total_score)

        max_score = max(max_scores)
        max_index = max_scores.index(max_score)
        motif = subsequences[max_index]
    except:
        print("Verify user input to be strings containing solely a, c, t, g!")
    print(f"Maxscore: {max_score}")
    print(f"Motiv: {motif}")

# Eingabedaten
dna_sequences = [
    'tagtggtcttttgagtgtagatctgaagggaaagtatttccaccagttcggggtcacccagcagggcagggtgacttaat',
    'cgcgactcggcgctcacagttatcgcacgtttagaccaaaacggagttggatccgaaactggagtttaatcggagtcctt',
    'gttacttgtgagcctggttagacccgaaatataattgttggctgcatagcggagctgacatacgagtaggggaaatgcgt',
    'aacatcaggctttgattaaacaatttaagcacgtaaatccgaattgacctgatgacaatacggaacatgccggctccggg',
    'accaccggataggctgcttattaggtccaaaaggtagtatcgtaataatggctcagccatgtcaatgtgcggcattccac',
    'tagattcgaatcgatcgtgtttctccctctgtgggttaacgaggggtccgaccttgctcgcatgtgccgaacttgtaccc',
    'gaaatggttcggtgcgatatcaggccgttctcttaacttggcggtgcagatccgaacgtctctggaggggtcgtgcgcta',
    'atgtatactagacattctaacgctcgcttattggcggagaccatttgctccactacaagaggctactgtgtagatccgta',
    'ttcttacacccttctttagatccaaacctgttggcgccatcttcttttcgagtccttgtacctccatttgctctgatgac',
    'ctacctatgtaaaacaacatctactaacgtagtccggtctttcctgatctgccctaacctacaggtcgatccgaaattcg'
]

k_length = 10
subsequences = generate_subsequences(dna_sequences, k_length)
find_regulatory_motif(dna_sequences, subsequences)


