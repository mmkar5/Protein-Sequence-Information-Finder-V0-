# Protein Sequence Property Information Finder (V0)
#### Video Demo: <https://youtu.be/OXrPTR6PThM->
#### Description:
The program project.py takes an id & protein sequence or an input fasta file containing several sequences & provides information regarding the following features of the given protein sequence:
- Length of the sequence
- Number of each amino acids
- Percentage of each amino acids
- Percentage of amino acids in each groups: nonpolar, aromatic, polar, all_charged, positive_charged, negative_charged
- Protein molecular weight
- Protein sequence motif presence and location
- GRAVY (grand average of hydropathy) score

This program considers the standard 20 common amino acids and does not take into acoount any modified amino acids.

This project also contains test_project.py designed to test all the functions implemented in project.py

`usage: python project.py -h -i I -o O -m -w`

The program takes the following command line arguments to determine properties of a protein sequence:
```
-h, --help  show this help message and exit

-i I        Enter input filepath of a protein fasta file (optional)

-o O        Enter output filepath (optional)

-m          Search for motifs (optional)

-w          get monoisoptopic molecular weight (optional)
```
If no input filepath for a fasta file is provides, will propmpt user to enter "Enter protein sequence identifier: " and subsequently in next line, "Enter protein sequence: ". If the input filename is not a ".fasta" file or does not exist, the program will exit, indicating that "File not found" or "Not a fasta file" respectively.

If output filepath is provided, then the outputs will be saved in the provided file, otherwise it will will be printed on the command line interface. If the output file already exists, then the program will exit, indicating "Given output file already exists!", preventing any overwriting on an existing file.

if -w is provided as a command line argument, then the monoisotopic molecular weight will be calculated for the protein sequence. Otherwise, by default, the averal molecular weight will be calculated.

if -m is provided as a command line argument, then user will be prompted, to provide a regular expression pattern to be matched as a protein sequence motif for each of the given protein sequence. The format of regular expression is printed to the user with an example to make it easier for the user to enter the correct input.

Apart from the main function, there are the following functions defined at the same indentation level as main:
-   get_args():
        Obtains command-line arguments from the user and parses them.
-   get_input_seq(args):
        Obtains input from the user (a fasta file or a sequnce id and sequence in the command-line) and
    returns a generator of the sequence id(s) and sequence(s).
-   get_motif():
        Obtains a motif from the user as a regular expression and
    returns the motif

All the methods used for obtaining the protein sequence features are implemented within a class called Protein_Seq_Features. The object of the class is initialised as, `protein = Protein_Seq_Features(seq)`,where seq is an instance attribute, indicaing the protein sequence. The value of seq is checked to ensure it is made up of alphbets, otherwise there will be a ValueError indicating "No valid sequence given". The methods of this class are:
-   seq_length(self):
        Returns the length of the protein sequence
-   aa_composition_num(self):
        Returns a dictionary containing the number of each of 20 amino acids present in the protein sequence
-   aa_composition_percentage(self):
        Returns a dictionary containing relative percentage of each of 20 amino acids present in the protein sequence
-   aa_group_percentage(self):
        Returns a dictionary containing percentage of nonpolar,aromatic,polar,all_charged, positive_charged
        and negative_charged amino acids present in the protein sequence
-   Mw(self, use_monoisotopic_mass=False):
        Returns the  average or monoisotopic molecular weight of the protein, depending on the argument of use_monotopic_mass;
        by default, returns the average molecular weight
-   motif_finder(self, motif):
        Takes a arguemt of the motif as a regular experssion, and returns the motifs matched and their position,
        if there is any match.
-   gravy_score(self):
        Returns the GRAVY (grand average of hydropathy) values of a protein sequence,
        which is a measure of the hydrophobicity or hydrophilicity of the protein sequence.

Example-

input.fasta:
```
>sp|P02748|CO9_HUMAN Complement component C9 OS=Homo sapiens OX=9606 GN=C9 PE=1 SV=2
MSACRSFAVAICILEISILTAQYTTSYDPELTESSGSASHIDCRMSPWSEWSQCDPCLRQ
MFRSRSIEVFGQFNGKRCTDAVGDRRQCVPTEPCEDAEDDCGNDFQCSTGRCIKMRLRCN
GDNDCGDFSDEDDCESEPRPPCRDRVVEESELARTAGYGINILGMDPLSTPFDNEFYNGL
CNRDRDGNTLTYYRRPWNVASLIYETKGEKNFRTEHYEEQIEAFKSIIQEKTSNFNAAIS
LKFTPTETNKAEQCCEETASSISLHGKGSFRFSYSKNETYQLFLSYSSKKEKMFLHVKGE
IHLGRFVMRNRDVVLTTTFVDDIKALPTTYEKGEYFAFLETYGTHYSSSGSLGGLYELIY
VLDKASMKRKGVELKDIKRCLGYHLDVSLAFSEISVGAEFNKDDCVKRGEGRAVNITSEN
LIDDVVSLIRGGTRKYAFELKEKLLRGTVIDVTDFVNWASSINDAPVLISQKLSPIYNLV
PVKMKNAHLKKQNLERAIEDYINEFSVRKCHTCQNGGTVILMDGKCLCACPFKFEGIACE
ISKQKISEGLPALEFPNEK
```

Terminal:
```
python.py -i input.fasta -o output.txt -m
{n} matches number of repeating units
{m,n} matches m to n repeating units
. matches any amino acid a|b indicates a or b
[] matches sets of charectars
[^] indicates charectars not in set
^ at beginning matches from N terminal and $ at end matches at C terminal
Example of C2H2 Zinc Finger motif: C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H
Enter pattern of motif in regular expression: N[^P][ST][^P]
```
output.txt:
```
ID:sp|P02748|CO9_HUMAN Complement component C9 OS=Homo sapiens OX=9606 GN=C9 PE=1 SV=2
MSACRSFAVAICILEISILTAQYTTSYDPELTESSGSASHIDCRMSPWSEWSQCDPCLRQMFRSRSIEVFGQFNGKRCTDAVGDRRQCVPTEPCEDAEDDCGNDFQCSTGRCIKMRLRCNGDNDCGDFSDEDDCESEPRPPCRDRVVEESELARTAGYGINILGMDPLSTPFDNEFYNGLCNRDRDGNTLTYYRRPWNVASLIYETKGEKNFRTEHYEEQIEAFKSIIQEKTSNFNAAISLKFTPTETNKAEQCCEETASSISLHGKGSFRFSYSKNETYQLFLSYSSKKEKMFLHVKGEIHLGRFVMRNRDVVLTTTFVDDIKALPTTYEKGEYFAFLETYGTHYSSSGSLGGLYELIYVLDKASMKRKGVELKDIKRCLGYHLDVSLAFSEISVGAEFNKDDCVKRGEGRAVNITSENLIDDVVSLIRGGTRKYAFELKEKLLRGTVIDVTDFVNWASSINDAPVLISQKLSPIYNLVPVKMKNAHLKKQNLERAIEDYINEFSVRKCHTCQNGGTVILMDGKCLCACPFKFEGIACEISKQKISEGLPALEFPNEK
Length:559
Number of amino acids:{'A': 29, 'V': 28, 'I': 33, 'L': 42, 'M': 10, 'P': 19, 'G': 36, 'C': 26, 'F': 28, 'Y': 21, 'W': 4, 'S': 48, 'T': 32, 'N': 27, 'Q': 14, 'H': 9, 'K': 37, 'R': 33, 'D': 35, 'E': 48}
Percentage of amino acids:{'A': 5.19, 'V': 5.01, 'I': 5.9, 'L': 7.51, 'M': 1.79, 'P': 3.4, 'G': 6.44, 'C': 4.65, 'F': 5.01, 'Y': 3.76, 'W': 0.72, 'S': 8.59, 'T': 5.72, 'N': 4.83, 'Q': 2.5, 'H': 1.61, 'K': 6.62, 'R': 5.9, 'D': 6.26, 'E': 8.59}
Percentage of amino acids in respective groups:{'nonpolar': 39.89, 'aromatic': 9.48, 'polar': 21.65, 'all_charged': 28.98, 'positive_charged': 8.23, 'negative_charged': 20.75}
Molecular Weight:63173.376 Da
GRAVY Score:-0.45
Motifs matched:[{'motif_match': 'NETY', 'position': (277, 280)}, {'motif_match': 'NITS', 'position': (415, 418)}]
```
