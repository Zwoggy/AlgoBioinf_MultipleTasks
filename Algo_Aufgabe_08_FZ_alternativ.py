"""
Algorithmische Bioinformatik Aufgabe 8

Aufgabe: Entwerfen Sie einen effizienten Algorithmus, der die längste darin vorkommende Tandemwiederholung in einem Text findet.

Author: Florian Zwicker BBM23
"""

import time
import csv

def open_csv(filename):
    """
   Read data from a .csv file and append it to a string.

   :param filename: path and filename of a .csv file.
   :return: a string representing file content.
   """
    with open(filename, 'r') as file:
        text = ""
        reader = csv.reader(file)
        for i, row in enumerate(reader):
            if len(row) > 0 and i > 0:
                text += row[0].replace("\n", "")
    return text



class DNA(object):
    """
    A class to represent a DNA sequence.
    """
    def __init__(self, dna_sequence):
        """
        Construct a new 'DNA' object.

        :param dna_sequence: A string that represents the DNA sequence.
        :return: returns nothing.
        """
        self.dna_sequence = dna_sequence
        self.kmeres = []

    def get_kmeres(self, length):
        """
        Get 3-character long substrings (k-mers) from the DNA sequence
        and store them in the class' attribute: self.kmeres.
        """

        for i in range(len(self.dna_sequence) - (length - 1)):
            self.kmeres.append(self.dna_sequence[i:i + length])

    def find_tandem_repeats(self):
        """
        Find the longest tandem repeat in the DNA sequence.

        :return: a string that represents the longest tandem repeat in the DNA sequence.
        """
        matching_slices = ""
        for kmere in list(set(self.kmeres)):
            print(kmere)
            sequence = self.dna_sequence
            dna_slices = []
            repeat_indices = []
            possible_repeat_index = self.dna_sequence.find(kmere)
            while possible_repeat_index != -1:
                repeat_indices.append(possible_repeat_index)
                possible_repeat_index = self.dna_sequence.find(kmere, possible_repeat_index + 1)

            for index in repeat_indices[::-1]:
                dna_slices.append(sequence[index:])
                sequence = sequence[:index]

            #TCGTCGAA
            for i in range(len(dna_slices) - 1):

                #if i != j:  # don't compare a slice with itself
                    slice_1 = dna_slices[-i]
                    slice_2 = dna_slices[-(i+1)]

                    # check if slice_2 starts with slice_1
                    if slice_2.startswith(slice_1):
                        #print("TRUE")
                        matching_slice = slice_1
                        matching_slice_2 = slice_2[:len(slice_1)]
                        # update matching_slices if the new slice is longer
                        if len(matching_slice) > len(matching_slices):
                            #print("ALSO TRUE")

                            matching_slices = matching_slice + matching_slice_2
                            print(matching_slices)

        return matching_slices, self.dna_sequence.find((matching_slices))


if __name__ == '__main__':
    start_time = time.time()
    dna_sequence = open_csv('Homo_sapiens_ENST00000621226_2_sequence.fa')
    dna = DNA(dna_sequence)
    dna.get_kmeres(20)
    tandem_repeat, index = dna.find_tandem_repeats()
    end_time = time.time()
    print("Index: " + str(index))
    print(str(len(tandem_repeat)) + "\n" +  tandem_repeat + "\n" + tandem_repeat)
    elapsed_time = end_time - start_time
    print(f"The program took {elapsed_time} seconds to run.")