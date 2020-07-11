import random
import pandas as pd
import json
from collections import Counter
import seaborn as sns

data = open('data/covid_proteins.txt', 'r')
covid_proteins = data.read()
# открываю файл со скопированной последовательностью аминокислот слепленных в один файл
codon_dict = {'DNA': ['TTT','TTC','TTA','TTG','TCT','TCC','TCA','TCG','TAT','TAC','TAA','TAG','TGT','TGC','TGA','TGG','CTT','CTC','CTA','CTG','CCT','CCC','CCA','CCG','CAT','CAC','CAA','CAG','CGT','CGC','CGA','CGG','ATT','ATC','ATA','ATG','ACT','ACC','ACA','ACG','AAT','AAC','AAA','AAG','AGT','AGC','AGA','AGG','GTT','GTC','GTA','GTG','GCT','GCC','GCA','GCG','GAT','GAC','GAA','GAG','GGT','GGC','GGA','GGG'],
       'AA': ['F','F','L','L','S','S','S','S','Y','Y','_','_','C','C','_','W','L','L','L','L','P','P','P','P','H','H','Q','Q','R','R','R','R','I','I','I','M','T','T','T','T','N','N','K','K','S','S','R','R','V','V','V','V','A','A','A','A','D','D','E','E','G','G','G','G']}
r2a = dict(zip(codon_dict['DNA'],codon_dict['AA']))
# Создал библиотеку соотношений кодонов и аминокислот
covid_dna = []
codon = []
# В цикле происходит обратная трансляция, аминоксилоты из переменной covid_proteins заменяются кодонами
# Если подходящих кодонов несколько, то рандомным из них
for i in covid_proteins:
       for j in i:
              for k, v in r2a.items():
                     if j == v:
                            codon.append(k)
              if codon != []:
                     covid_dna.append(random.choice(codon))
              codon = []

#print(covid_rna) # Тут хранится последовательность ДНК

# Далее считаю вероятность мутаций

covid_data = pd.read_csv('data/covid_next_strain.tsv', sep = '\t')
only_mutations = covid_data['mutations'].dropna()
only_mutations = only_mutations.str.replace("'", '"')
nuclear_mutations = []
for i in only_mutations:
       mutation = json.loads(i)
       for k, v in mutation.items():
              if k == 'nuc':
                     nuclear_mutations.append((v[0][0],v[0][-1]))

first = []
second = []
probability = []
summ_A = 0
summ_T = 0
summ_G = 0
summ_C = 0
for k,v in Counter(nuclear_mutations).items(): # Нормализую данные по колличеству мутаций по нуклеотидам
       if k[0] == 'A':
              summ_A += v
       if k[0] == 'T':
              summ_T += v
       if k[0] == 'G':
              summ_G += v
       if k[0] == 'C':
              summ_C += v
for k,v in Counter(nuclear_mutations).items():
       first.append(k[0])
       second.append(k[1])
       if k[0] == 'A':
              probability.append(v/summ_A)
       if k[0] == 'T':
              probability.append(v/summ_T)
       if k[0] == 'G':
              probability.append(v/summ_G)
       if k[0] == 'C':
              probability.append(v/summ_C)

mutation_probability = pd.DataFrame({'From' : first, 'To' : second, 'probability' : probability})
mutation_probability.to_csv('data/mutation_probability.tsv', sep = '\t')
# Сохранил вероятности мутаций в файл

# Далее на основании полученных данных мутируем нашу ДНК и подсчитываем количество мутаций

new_DNA = covid_dna
mutation_counter = 0
mutation_list = []
for i in range(100000):
       random_codon = random.choice(new_DNA) # Беру случайный кодон из ДНК
       AA = r2a[random_codon] # Присваиваю ему аминокислоту
       index_random_nucl = random.choice(range(len(random_codon)))
       random_nucl = random_codon[index_random_nucl] # Беру по индексу рандомный нуклеотид
       for k in range(12):
              if random_nucl == mutation_probability.From[k]: # Проверяю совпадает ли нуклеотид с таковыми из таблички с мутациями
                     my_dice = random.random() # Кидаю кубик на удачу чтоб выпала мутация
                     if my_dice <= mutation_probability.probability[k]:
                            random_codon = random_codon[:index_random_nucl]+mutation_probability.To[k]+random_codon[index_random_nucl+1:]
              # Если мутация выпала - заменил в кодоне нуклеотид на мутированный
       AA_new = r2a[random_codon] # Присвоил новому кодону аминокислоту
       if AA == AA_new: # Сравниваю аминокислоты
              mutation_counter += 1  # Если они совпадают - добавляю к счётчику мутаций еденичку
       else:
              mutation_list.append(mutation_counter)
              mutation_counter = 0

#print(mutation_list) Тут хранятся данные по количеству мутаций без изменения аминокислоты
data = pd.DataFrame(mutation_list, columns = ['number of mutation without changing amino acid'])
# Рисую график распределения
sns.set()
sns_plot = sns.distplot(data['number of mutation without changing amino acid'])
fig = sns_plot.get_figure()
fig.savefig('data/covid_mutation_distribution.png') # Сохранил график в файл