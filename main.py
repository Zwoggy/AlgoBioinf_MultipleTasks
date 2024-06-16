"""
Algorithmische Bioinformatik Aufgabe 3:
Erstellen Sie ein Programm, mit dem in t DNA-Sequenzen mit jeweils N Nukleotiden das beste Motiv  der Länge k=10 gefunden wird.

Author: Florian Zwicker
"""

dna2 = [
    'tagtggtcttttgagtgtagatctgaagggaaagtatttccaccagttcggggtcacccagcagggcagggtgacttaat',
    'cgcgactcggcgctcacagttatcgcacgtttagaccaaaacggagttggatccgaaactggagtttaatcggagtcctt',
    'gttacttgtgagcctggttagacccgaaatataattgttggctgcatagcggagctgacatacgagtaggggaaatgcgt',
    'aacatcaggctttgattaaacaatttaagcacgtaaatccgaattgacctgatgacaatacggaacatgccggctccggg',
    'accaccggataggctgcttattaggtccaaaaggtagtatcgtaataatggctcagccatgtcaatgtgcggcattccac',
    'tagattcgaatcgatcgtgtttctccctctgtgggttaacgaggggtccgaccttgctcgcatgtgccgaacttgtaccc',
    'gaaatggttcggtgcgatatcaggccgttctcttaacttggcggtgcagatccgaacgtctctggaggggtcgtgcgcta',
    'atgtatactagacattctaacgctcgcttattggcggagaccatttgctccactacaagaggctactgtgtagatccgta',
    'ttcttacacccttctttagatccaaacctgttggcgccatcttcttttcgagtccttgtacctccatttgctctgatgac',
    'ctacctatgtaaaacaacatctactaacgtagtccggtctttcctgatctgccctaacctacaggtcgatccgaaattcg']




def create_subsequences(sequences):
    all_subsequences = []
    try:
        for seq in sequences:
            print(seq)
            reg_sequences = []
            for i in range((len(seq) - 9)):
                reg_sequences.append(seq[i:(i + 10)])
            all_subsequences.append(reg_sequences)

    except:
        print("Verify user Input to be a list of strings that contain a, t, g, c only!")
    return all_subsequences


def search_for_motif(all_subsequences):
    max_score = 0
    targets = 0
    for i in range(len(all_subsequences[0])):
        """for all subsequences of a sequence"""
        score = 0
        current_target = all_subsequences[0][i]
        for j in range(len(all_subsequences)):
            """for all 10 sequences"""

            for k in range(len(all_subsequences[0][0])):
                """for every base in a subsequence"""
                if current_target[j] != all_subsequences[j][i][k]:
                    score +=1

        if score > max_score:
            max_score = score
            motif = current_target
    print(motif)
    print(score)








if __name__ == "__main__":
    x = create_subsequences(dna2)
    print(x)
    search_for_motif(x)