import numpy as np

def needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-1):
    # Initialize the scoring matrix
    n = len(seq1) + 1
    m = len(seq2) + 1
    scoring_matrix = np.zeros((n, m), dtype=int)
    traceback_matrix = np.zeros((n, m), dtype=int)

    # Fill the first row and column with gap penalties
    for i in range(1, n):
        scoring_matrix[i][0] = i * gap_penalty
    for j in range(1, m):
        scoring_matrix[0][j] = j * gap_penalty

    # Fill the scoring matrix and the traceback matrix
    for i in range(1, n):
        for j in range(1, m):
            match = scoring_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            delete = scoring_matrix[i - 1][j] + gap_penalty
            insert = scoring_matrix[i][j - 1] + gap_penalty
            scoring_matrix[i][j] = max(match, delete, insert)

            if scoring_matrix[i][j] == match:
                traceback_matrix[i][j] = 1
            elif scoring_matrix[i][j] == delete:
                traceback_matrix[i][j] = 2
            else:
                traceback_matrix[i][j] = 3

    # Traceback to get the aligned sequences
    aligned_seq1 = ""
    aligned_seq2 = ""
    i = n - 1
    j = m - 1

    while i > 0 or j > 0:
        if traceback_matrix[i][j] == 1:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 2:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1

    return aligned_seq1, aligned_seq2, scoring_matrix[-1][-1]

# Beispielsequenzen
seq1 = "GATTACA"
seq2 = "GCATGCU"

aligned_seq1, aligned_seq2, score = needleman_wunsch(seq1, seq2)
print("Seq1: ", aligned_seq1)
print("Seq2: ", aligned_seq2)
print("Score: ", score)
