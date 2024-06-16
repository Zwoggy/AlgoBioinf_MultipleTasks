"""
Aufgabe 9 Algorithmische Bioinformaik:
Implementieren Sie einen Algorithmus zur Sekundärstrukturvorhersage von RNA-Molekülen, welche auf Basenpaarungsmaximierung beruht.
Z.B. Nussinov Algo
"""


def Nussinov(seq):
    """
    This is the main function that implements the Nussinov algorithm.

    :param seq: The sequence of RNA molecule.
    :return: The secondary structure prediction of RNA molecule as a string.
    """
    S = seq
    N = len(seq)
    R = [[0 for col in range(N)] for row in range(N)]
    bp = ['.' for b in range(N)]

    for k in range(1, N):
        for i in range(N - k):
            j = i + k

            if j - i < 3:
                R[i][j] = 0
            else:
                R[i][j] = max([R[i + 1][j - 1] + 1 if is_pair(S[i], S[j]) else 0, R[i][j - 1], R[i + 1][j],
                               max([R[i][t] + R[t + 1][j] for t in range(i + 1, j - 1)])])

    traceback(R, S, bp, 0, N - 1)
    return "".join(bp)


def is_pair(bp1, bp2):
    """
    This function checks if two bases can be paired together.

    :param bp1: The first base.
    :param bp2: The second base.
    :return: True if the bases can be paired together, False otherwise.
    """
    if (bp1 == 'A' and bp2 == 'U') or (bp1 == 'U' and bp2 == 'A') or (bp1 == 'C' and bp2 == 'G') or (
            bp1 == 'G' and bp2 == 'C'):
        return True
    else:
        return False


def traceback(R, S, bp, i, j):
    """
    This function checks if two bases can be paired together.

    :param bp1: The first base.
    :param bp2: The second base.
    :return: True if the bases can be paired together, False otherwise.
    """
    if i < j:
        if R[i][j] == R[i + 1][j]:
            traceback(R, S, bp, i + 1, j)
        elif R[i][j] == R[i][j - 1]:
            traceback(R, S, bp, i, j - 1)
        elif R[i][j] == R[i + 1][j - 1] + 1 and is_pair(S[i], S[j]):
            bp[i] = "("
            bp[j] = ")"
            traceback(R, S, bp, i + 1, j - 1)
        else:
            for k in range(i + 1, j):
                if R[i][j] == R[i][k] + R[k + 1][j]:
                    traceback(R, S, bp, i, k)
                    traceback(R, S, bp, k + 1, j)
                    break




if __name__ == "__main__":
    """
    The main entry point of the code. 
    Uses the Nussinov algorithm on a sample RNA sequence.
    Verify the printed result on http://nibiru.tbi.univie.ac.at/forna/forna.html
    """
    sequence = "GGGAAAUCC"

    print(Nussinov(sequence))