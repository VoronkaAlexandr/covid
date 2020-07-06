import random

data = open('/home/bufnita/PycharmProjects/covid/data/covid_proteins.txt', 'r')
covid_proteins = data.read()
# открываю файл со скопированной последовательностью аминокислот слепленных в один файл
codon_dict = {'mRNA': ['UUU','UUC','UUA','UUG','UCU','UCC','UCA','UCG','UAU','UAC','UAA','UAG','UGU','UGC','UGA','UGG','CUU','CUC','CUA','CUG','CCU','CCC','CCA','CCG','CAU','CAC','CAA','CAG','CGU','CGC','CGA','CGG','AUU','AUC','AUA','AUG','ACU','ACC','ACA','ACG','AAU','AAC','AAA','AAG','AGU','AGC','AGA','AGG','GUU','GUC','GUA','GUG','GCU','GCC','GCA','GCG','GAU','GAC','GAA','GAG','GGU','GGC','GGA','GGG'],
       'AA': ['F','F','L','L','S','S','S','S','Y','Y','_','_','C','C','_','W','L','L','L','L','P','P','P','P','H','H','Q','Q','R','R','R','R','I','I','I','M','T','T','T','T','N','N','K','K','S','S','R','R','V','V','V','V','A','A','A','A','D','D','E','E','G','G','G','G']}

r2a = dict(zip(codon_dict['mRNA'],codon_dict['AA']))

covid_dna = []
codon = []

for i in covid_proteins:
       for j in i:
              for k, v in r2a.items():
                     if j == v:
                            codon.append(k)
              if codon != []:
                     covid_dna.append(random.choice(codon))
              codon = []

print(covid_dna)