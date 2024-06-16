"""
Aufgabe: Erstellen Sie eine Programm zur Rekonstruktion einer Sequenz anhand gegebener k-mere (Hamiltonpfad-Ansatz), wobei die k-mere jeweils um k-1

Beispielsatz = {'ATG', 'AGG', 'TGC', 'TCC', 'GTC', 'GGT', 'GCA', 'CAG'}

Author: Florian Zwicker
"""
k_meres = ['ATG', 'AGG', 'TGC', 'TCC', 'GTC', 'GGT', 'GCA', 'CAG']


def find_start(k_meres):
    """
    This function finds the k_mere such that the first two characters of the k_mere
    do not appear as the second and the third characters in any other k_mere

    :param k_meres: List of kmers
    :return: The starting k_mere
    """
    for i, k_mere in enumerate(k_meres):
        k_mere_list = k_meres[:]
        k_mere_list.pop(i)
        flag = 0
        for entry in k_mere_list:
            if entry[1:3] == k_mere[:2]:
                flag = 1
                break
        if flag == 0:
            return k_mere


def find_end(k_meres):
    """
    This function finds the k_mere such that the last two characters of the k_mere
    do not appear as the first two characters in any other k_mere

    :param k_meres: List of kmers
    :return: The ending k_mere
    """
    for i, k_mere in enumerate(k_meres):
        k_mere_list = k_meres[:]
        k_mere_list.pop(i)
        flag = 0
        for entry in k_mere_list:
            if entry[:2] == k_mere[1:]:
                flag = 1
                break
        if flag == 0:
            return k_mere


def reconstruct_sequence(k_meres):
    """
    This function reconstructs the sequence from the given kmers.

    :param k_meres: List of kmers
    :return: The reconstructed sequence
    """
    start = find_start(k_meres)
    end = find_end(k_meres)
    sequence = start
    k_meres.remove(start)
    k_meres.remove(end)
    while len(k_meres) > 0:
        for k_mere in k_meres:
            if sequence[-2:] == k_mere[:2]:
                sequence += k_mere[2:]
                k_meres.remove(k_mere)
                break
    sequence += end[2:]
    return sequence



if __name__ == '__main__':
    """
    The main method that invokes reconstruct_sequence function and prints 
    the reconstructed sequence from the given kmers.
    """
    reconstructed_sequence = reconstruct_sequence(k_meres)
    print(reconstructed_sequence)