import argparse
import sys
import re
import os


class Protein_Seq_Features:
    amino_acids = [
        "A",
        "V",
        "I",
        "L",
        "M",
        "P",
        "G",
        "C",
        "F",
        "Y",
        "W",
        "S",
        "T",
        "N",
        "Q",
        "H",
        "K",
        "R",
        "D",
        "E",
    ]

    def __init__(self, seq):
        self.seq = seq.upper()

    def seq_length(self):
        """Returns the length of the protein sequence"""

        return len(self.seq)

    def aa_composition_num(self):
        """Returns a dictionary containing the number of each of 20 amino acids present in the protein sequence"""

        return {aa: self.seq.count(aa) for aa in Protein_Seq_Features.amino_acids}

    def aa_composition_percentage(self):
        """Returns a dictionary containing relative percentage of each of 20 amino acids present in the protein sequence"""

        return {
            aa: round(self.seq.count(aa) / len(self.seq) * 100, 2)
            for aa in Protein_Seq_Features.amino_acids
        }

    def aa_group_percentage(self):
        """Returns a dictionary containing percentage of nonpolar,aromatic,polar,all_charged, positive_charged
        and negative_charged amino acids present in the protein sequence """

        self.nonpolar = round(
            sum([self.seq.count(aa) for aa in ["A", "V", "I", "L", "M", "P", "G", "C"]])
            / len(self.seq)
            * 100,
            2,
        )
        self.aromatic = round(
            sum([self.seq.count(aa) for aa in ["F", "Y", "W"]]) / len(self.seq) * 100, 2
        )
        self.polar = round(
            sum([self.seq.count(aa) for aa in ["S", "T", "N", "Q"]])
            / len(self.seq)
            * 100,
            2,
        )
        self.all_charged = round(
            sum([self.seq.count(aa) for aa in ["H", "K", "R", "D", "E"]])
            / len(self.seq)
            * 100,
            2,
        )
        self.positive_charged = round(
            sum([self.seq.count(aa) for aa in ["H", "K"]]) / len(self.seq) * 100, 2
        )
        self.negative_charged = round(
            sum([self.seq.count(aa) for aa in ["R", "D", "E"]]) / len(self.seq) * 100, 2
        )
        return {
            "nonpolar": self.nonpolar,
            "aromatic": self.aromatic,
            "polar": self.polar,
            "all_charged": self.all_charged,
            "positive_charged": self.positive_charged,
            "negative_charged": self.negative_charged,
        }

    def Mw(self, use_monoisotopic_mass=False):
        """Returns the  average or monoisotopic molecular weight of the protein, depending on the argument of use_monotopic_mass;
        by default, returns the average molecular weight"""

        monoisotopic_mass = {
            "A": 71.03711,
            "R": 156.10111,
            "N": 114.04293,
            "D": 115.02694,
            "C": 103.00919,
            "E": 129.04259,
            "Q": 128.05858,
            "G": 57.02146,
            "H": 137.05891,
            "I": 113.08406,
            "L": 113.08406,
            "K": 128.09496,
            "M": 131.04049,
            "F": 147.06841,
            "P": 97.05276,
            "S": 87.03203,
            "T": 101.04768,
            "W": 186.07931,
            "Y": 163.06333,
            "V": 99.06841,
        }
        average_mass = {
            "A": 71.0788,
            "R": 156.1875,
            "N": 114.1038,
            "D": 115.0886,
            "C": 103.1388,
            "E": 129.1155,
            "Q": 128.1307,
            "G": 57.0519,
            "H": 137.1411,
            "I": 113.1594,
            "L": 113.1594,
            "K": 128.1741,
            "M": 131.1926,
            "F": 147.1766,
            "P": 97.1167,
            "S": 87.0782,
            "T": 101.1051,
            "W": 186.2132,
            "Y": 163.1760,
            "V": 99.1326,
        }
        Mw = 0
        if use_monoisotopic_mass:
            return round(
                sum(
                    [
                        self.seq.count(aa) * monoisotopic_mass[aa]
                        for aa in monoisotopic_mass
                    ]
                )
                + 18.010565,
                3,
            )
        else:
            return round(
                sum([self.seq.count(aa) * average_mass[aa] for aa in average_mass])
                + 18.01528,
                3,
            )

    def motif_finder(self, motif):
        """Takes a arguemt of the motif as a regular experssion, and returns the motifs matched and their position,
        if there is any match."""

        if motif:
            matched = [
                {"motif_match": match.group(), "position": (match.start()+1,match.end())}
                for match in re.finditer(motif.strip(), self.seq, re.I)
            ]
            if matched:
                return matched
            else:
                return f"No such motif exits in sequence"
        else:
            return f"The input is empty. Please provide a valid motif."

    def gravy_score(self):
        """Returns the GRAVY (grand average of hydropathy) values of a protein sequence,
        which is a measure of the hydrophobicity or hydrophilicity of the protein sequence."""

        hydropathy_values = {

            "A": 1.800,
            "R": -4.500,
            "N": -3.500,
            "D": -3.500,
            "C": 2.500,
            "E": -3.500,
            "Q": -3.500,
            "G": -0.400,
            "H": -3.200,
            "I": 4.500,
            "L": 3.800,
            "K": -3.900,
            "M": 1.900,
            "F": 2.800,
            "P": -1.600,
            "S": -0.800,
            "T": -0.700,
            "W": -0.900,
            "Y": -1.300,
            "V": 4.200,
        }

        return round(
            sum(
                [self.seq.count(aa) * hydropathy_values[aa] for aa in hydropathy_values]
            )
            / len(self.seq),
            3,
        )

    @property
    def seq(self):
        return self.__seq

    @seq.setter
    def seq(self, seq):
        if not seq.isalpha():
            raise ValueError("No sequence given")
        self.__seq = seq


def main():

    args = get_args()
    if args.m:
        motif = get_motif()

    if args.o:
        if os.path.exists(args.o):
            sys.exit("Given output file already exists!")
        else:
            output_file = open(args.o, "w")
    else:
        output_file = None

    for seq_id, seq in get_input_seq(args.i):
        protein = Protein_Seq_Features(seq)
        results = [
            f"ID:{seq_id}",
            seq,
            f"Length:{protein.seq_length()}",
            f"Number of amino acids:{protein.aa_composition_num()}",
            f"Percentage of amino acids:{protein.aa_composition_percentage()}",
            f"Percentage of amino acids in respective groups:{protein.aa_group_percentage()}",
            f"Molecular Weight:{protein.Mw(args.w)} Da",
            f"GRAVY Score:{protein.gravy_score()}",
        ]
        if args.m:
            results.append(f"Motifs matched:{protein.motif_finder(motif)}")

        for result in results:
            if output_file:
                output_file.write(result + "\n")
            else:
                print(result)

    if output_file:
        output_file.close()


def get_args():
    """Obtains command-line arguments from the user and parses them"""
    parser = argparse.ArgumentParser(
        description="Arguments to determine properties of a protein sequence"
    )
    parser.add_argument(
        "-i",
        required=False,
        type=str,
        help="Enter input filepath of a protein fasta file",
    )
    parser.add_argument("-o", required=False, type=str, help="Enter output filepath")
    parser.add_argument("-m", action="store_true", help="Search for motifs")
    parser.add_argument(
        "-w", action="store_true", help="get monoisoptopic molecular weight"
    )
    args = parser.parse_args()
    return args


def get_input_seq(args):
    """Obtains input from the user (a fasta file or a sequnce id and sequence in the command-line) and
    returns a generator of the sequence id(s) and sequence(s)."""

    if args:
        try:
            if ".fasta" in args:
                with open(args, "r") as file:
                    seq=""
                    seq_id=""
                    for line in file:
                        line=line.strip()
                        if line.startswith(">"):
                            if seq:
                                yield seq_id,seq
                                seq=""
                            seq_id=line.lstrip(">")
                        else:
                            seq+=line
                    yield seq_id,seq
            else:
                raise ValueError
        except FileNotFoundError:
            sys.exit("File not found")
        except ValueError:
            sys.exit("Not a fasta file")
    else:
        seq_id = input("Enter protein sequence identifier: ")
        seq = input("Enter protein sequence: ").strip()
        yield seq_id, seq


def get_motif():
    """Obtains a motif from the user as a regular expression and
    returns the motif"""

    print("{n} matches number of repeating units")
    print("{m,n} matches m to n repeating units")
    print(". matches any amino acid", "a|b indicates a or b")
    print("[] matches sets of charectars")
    print("[^] indicates charectars not in set")
    print("^ at beginning matches from N terminal and $ at end matches at C terminal")
    print("Example of C2H2 Zinc Finger motif: C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H")
    motif = input("Enter pattern of motif in regular expression: ").strip()
    return motif


if __name__ == "__main__":
    main()
