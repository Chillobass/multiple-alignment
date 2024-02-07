import pandas as pd
import numpy as np
import textwrap as tw
from tkinter import *
from tkinter import filedialog

"""                                                     ГЛОБАЛЬНЫЕ ПЕРЕМЕННЫЕ                                                              """
num_A = 0
num_T = 0
num_G = 0
num_C = 0
num_gap = 0
pos = 0
                # num_x - счетчики нуклеотидов

f = []          # шаблон матрицы выравнивания
maxt = 0        # максимум выравнивания
sc1 = ''        # переменная, содержащая текущий скэффолд, изменяется функцией smena_scaff()
ii = 0          # координата maxt
jj = 0          # координата maxt
res = []        # список со всеми координатами пути от maxt до нуля
result = []     # список с координатами начала и конца каждой выравненной части
final = ''
l = []
consensus = ''



"""                                                         ФУНКЦИИ                                                         """
def ok(pos, d):          # функция "ok" считает встречавемость нуклеотидов в переданной позиции(pos)
  global num_A
  global num_T
  global num_G
  global num_C
  global num_gap
  for i in d:
    if len(i) <= pos:
      continue
    if i[pos] == '-':
      num_gap += 1
    elif i[pos] == 'A':
      num_A += 1
    elif i[pos] == 'T':
      num_T += 1
    elif i[pos] == 'G':
      num_G += 1
    elif i[pos] == 'C':
      num_C += 1




def smena_scaff(number_of_sc, aaa): # функция, кт отдаёт скаффолд (по индексу number_of_sc)
  global sc1
  sc1 = aaa[number_of_sc]



def mmax(consensus, sc1, maxt):  # находит максимальный score в матрице выравнивания
    global f, l
    global ii, jj

    l = [0, 0]
    m = len(consensus) + 1
    n = len(sc1) + 1
    f = np.zeros((n, m), dtype=int)
    for i in range(1, n):
        for j in range(1, m):
            if consensus[j-1] == sc1[i-1]:
                score = 1
            else:
                score = -2
            f[i][j] = max(f[i - 1][j] - 2, f[i][j - 1] - 2, f[i - 1][j - 1] + score, 0)

    for k in range(1, n):
        for t in range(1, m):
            if maxt <= f[k][t]:
                maxt = f[k][t]
                l[0] = k
                l[1] = t

    print(maxt)
    # ii = np.where(f == maxt)[0][0]
    # jj = np.where(f == maxt)[1][0]
    #np.savetxt('text.txt', f, fmt='%.0f')
    ii = l[0]
    jj = l[1]
    print(ii, jj)


    # return maxt, ii, jj  # получаем maxt - максимум матрицы, от кт мы восстановим нужную нам часть последовательности



def nucleotide_score(seq1, seq2):
  key=0
  if seq1==seq2:
    key=1
  if seq1!=seq2 and seq1!="N":
    key=-2
  if seq1!=seq2 and seq1=="N":
    key=0
  return key


def path(consensus):  # путь от максимума до начала
    global jj, ii
    global res, result
    global f
    global sc1
    res = []

    i = ii
    j = jj
    # print(i, j)
    # for k in range(len(f)):
    #     for t in range(len(f[0])):
    #         print(f[k][t], end=" ")
    #     print("\n")


    while True:
        if f[i][j] == 0:
            #print("end of cycle bruh...")
            break
        else:
            if f[i][j - 1] + 2 >= f[i - 1][j] + 2 and f[i][j - 1] + 2 >= f[i - 1][j - 1] - nucleotide_score(
                    consensus[j - 1], sc1[i - 1]):
                j -= 1
            elif f[i - 1][j] + 2 >= f[i][j - 1] + 2 and f[i - 1][j] + 2 >= f[i - 1][j - 1] - nucleotide_score(
                    consensus[j - 1], sc1[i - 1]):
                i -= 1
            elif [i - 1][j - 1] - nucleotide_score(consensus[j - 1], sc1[i - 1]) >= [i - 1][j] + 2 and f[i - 1][j - 1] - nucleotide_score(consensus[j - 1], sc1[i - 1]) >= f[i][j - 1] + 2:
                i -= 1
                j -= 1
        res.append([i, j])  # в res содержится список с индексами всего пути по матрице
    result.append([res[-1][1]-2, res[0][1]+1, sc1])  # добавляет в result индексы последовательности, на которую мы выровняли
    #print((consensus[res[-1][1]-2:res[0][1]+1]))
    #(result)


def doin_file(res, maxt):                      # функция, создающая файл с заголовком и полной итоговой последовательностью
  text = ''                               # эта переменная - всё, что пойдет в файл, то есть заголовок (zagolov) + final
  zagolov = '>'                    # делаем заголовок
  #zagolov += id.strip('\n')
  zagolov += '/'
  zagolov += str(res[-1][1])
  zagolov += '/'
  zagolov += str(res[0][1])
  zagolov += '/'
  zagolov +=str(res[-1][0]-1)
  zagolov += '/'
  zagolov += str(res[0][0]-1)
  zagolov += '/'
  zagolov += str(maxt)
  text+=zagolov
  text+='\n'

  wrapped_final = tw.wrap(final, 200)     # перезаписываем и изменяем final для добавления в файл по 200 символов в строке
  for i in range(len(wrapped_final)):
    text+=wrapped_final[i]
    text += '\n'
  return text





"""                                                                 ЗАПУСК ВСЕГО ЭТОГО БЕЗОБРАЗИЯ                                                                                       """
def clicked_consensus():
    s1 = filedialog.askopenfilename()
    global consensus
    gg = ''
    a = {}
    id = ''
    with open(s1) as file:
        for row in file:
            if row[0] == '>':
                id = row.split()[0][1:]
                gg = ''
            else:
                gg += row.strip()
            a[id] = gg

    d = list(a.values())  # соответственно словарь последовательностей
    print(d)
    global num_A
    global num_T
    global num_G
    global num_C
    global num_gap
    global maxt
    global pos
    for pos in range(len(d[0])):
        num_A = 0  # задаем счетчики для каждого нуклеотида и счетчик для гэпов
        num_T = 0  #
        num_G = 0
        num_C = 0
        num_gap = 0
        ok(pos, d)
        # print(num_A, num_T, num_G, num_C, num_gap)

        if num_gap == max(num_A, num_T, num_G, num_C, num_gap):   # если гэпов больше всего - пропускаем
            continue
        elif num_A == max(num_A, num_T, num_G, num_C):            # для каждой позиции поочередно узнаем какой нуклеотид(или гэп) встретился чаще остальных
            consensus += 'A'                                      # и его заносим в консунсус
        elif num_T == max(num_A, num_T, num_G, num_C):
            consensus += 'T'
        elif num_G == max(num_A, num_T, num_G, num_C):
            consensus += 'G'
        elif num_C == max(num_A, num_T, num_G, num_C):
            consensus += 'C'
    print(consensus)




def clicked_scaff():
    global consensus
    # скэффолды
    s = filedialog.askopenfilename()  # s - путь к файлу со скэффолдами
    aaa = []
    with open(s) as file:  # открываем скаффолды, добавляем в aaa всё кроме заголовков
        for row in file:
            if row[0] == '>':
                continue
            else:
                aaa.append(row.strip())
    aaa.pop(0)
    aaa = [i for i in aaa if i != '']  # форматируем aaa
    #print(aaa)  # aaa - наш список со всеми скэффолдами



    global final
    global result
    global res
    for i in range(0, len(aaa)):                                    # основной цикл запуска всех функций поочередно
        smena_scaff(i, aaa)                                         # берем скаффолд(sc1)
        mmax(consensus, sc1, maxt)                                  # находим максимум его выравнивания на consensus
        path(consensus)                                             # восстанавливаем путь, то есть получаем координаты отрезка, на кт мы выровняли sc1



    result = sorted(result)
    print(result)
    index_konca = len(final)
    for i in range(0, len(result)):
        final += consensus[index_konca:result[i][0]]
        final += result[i][2]
        index_konca = len(final)
    final += consensus[index_konca:]


    text = doin_file(res, maxt)
    with open(r'final.fasta', 'w') as fp:  # теперь там же, где скрипт, появится наш fasta файл
        fp.write(text)

def main():
    window = Tk()
    window.title("local_align")
    window.geometry("700x400")
    window.config(bg='White smoke')
    lbl = Label(window, text="1) Choose file with genomes",font=("Arial Bold", 18))
    lbl.grid(column=0, row=0, sticky='w')
    lbl1 = Label(window, text="2) Choose file with scaffolds", font=("Arial Bold", 18))
    lbl1.grid(column=0, row=3, sticky='w')

    btn1 = Button(window, text="genomes", command=clicked_consensus, bg="black", fg="white", width=10, height=3)
    btn1.grid(column=0, row=2, sticky='w')

    btn2 = Button(window, text="scaffolds", command=clicked_scaff, bg="black", fg="white", width=10, height=3)
    btn2.grid(column=0, row=4, sticky='w')
    window.mainloop()

main()