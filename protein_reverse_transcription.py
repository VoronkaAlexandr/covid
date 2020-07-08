import random
import pandas as pd
import json
from collections import Counter

data = open('/home/bufnita/PycharmProjects/covid/data/covid_proteins.txt', 'r')
covid_proteins = data.read()
# открываю файл со скопированной последовательностью аминокислот слепленных в один файл
codon_dict = {'mRNA': ['UUU','UUC','UUA','UUG','UCU','UCC','UCA','UCG','UAU','UAC','UAA','UAG','UGU','UGC','UGA','UGG','CUU','CUC','CUA','CUG','CCU','CCC','CCA','CCG','CAU','CAC','CAA','CAG','CGU','CGC','CGA','CGG','AUU','AUC','AUA','AUG','ACU','ACC','ACA','ACG','AAU','AAC','AAA','AAG','AGU','AGC','AGA','AGG','GUU','GUC','GUA','GUG','GCU','GCC','GCA','GCG','GAU','GAC','GAA','GAG','GGU','GGC','GGA','GGG'],
       'AA': ['F','F','L','L','S','S','S','S','Y','Y','_','_','C','C','_','W','L','L','L','L','P','P','P','P','H','H','Q','Q','R','R','R','R','I','I','I','M','T','T','T','T','N','N','K','K','S','S','R','R','V','V','V','V','A','A','A','A','D','D','E','E','G','G','G','G']}
r2a = dict(zip(codon_dict['mRNA'],codon_dict['AA']))
# Создал библиотеку соотношений кодонов и аминокислот
covid_rna = []
codon = []
# В цикле происходит обратная трансляция, аминоксилоты из переменной covid_proteins заменяются кодонами
# Если подходящих кодонов несколько, то рандомным из них
for i in covid_proteins:
       for j in i:
              for k, v in r2a.items():
                     if j == v:
                            codon.append(k)
              if codon != []:
                     covid_rna.append(random.choice(codon))
              codon = []

#print(covid_rna) # Тут хранится последовательность РНК

# Далее считаю вероятность мутаций

covid_data = pd.read_csv('/home/bufnita/PycharmProjects/covid/data/covid_next_strain.tsv', sep = '\t')
only_mutations = covid_data['mutations'].dropna()
data_lenght = len(only_mutations) # Колличество наблюдений для рассчёта вероятности
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
for k,v in Counter(nuclear_mutations).items():
       first.append(k[0])
       second.append(k[1])
       probability.append(v/data_lenght)

mutation_probability = pd.DataFrame({'frome' : first, 'to' : second, 'probability' : probability})
mutation_probability = mutation_probability.replace('T', 'U')
mutation_probability.to_csv('/home/bufnita/PycharmProjects/covid/data/mutation_probability.tsv', sep = '\t')
# Сохранил вероятности мутаций в файл

# Далее на основании полученных данных мутируем нашу РНК и подсчитываем количество мутаций

new_RNA = covid_rna
mutation_counter = 0
for i in range(100000):
       random_codon = random.choice(new_RNA) # Беру случайный кодон из РНК
       AA = r2a[random_codon] # Присваиваю ему аминокислоту
       index_random_nucl = random.choice(range(len(random_codon)))
       random_nucl = random_codon[index_random_nucl] # Беру по индексу рандомный нуклеотид
       for k in range(12):
              if random_nucl == mutation_probability.frome[k]: # Проверяю совпадает ли нуклеотид с таковыми из таблички с мутациями
                     my_dice = random.random() # Кидаю кубик на удачу чтоб выпала мутация
                     if my_dice <= mutation_probability.probability[k]:
                            random_codon = random_codon[:index_random_nucl]+mutation_probability.to[k]+random_codon[index_random_nucl+1:]
              # Если мутация выпала - заменил в кодоне нуклеотид на мутированный
       AA_new = r2a[random_codon] # Присвоил новому кодону аминокислоту
       if AA != AA_new: # Сравниваю аминокислоты
              break
       else:
              mutation_counter += 1 # Если они совпадают - добавляю к счётчику мутаций еденичку
print(mutation_counter)