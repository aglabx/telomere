#import subprocess
import re
import os
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from collections import defaultdict
import statistics
import pandas as pd
from tabulate import tabulate
from PIL import Image, ImageDraw, ImageFont

###############################################################################

# Шаг 1.0, парсинг координат из .dat
# Подходит если trf запускался вместе с файлом впервые(повторы из .dat идут 
# в той же последовательности, что и контиги из .fasta и не требуют дополнительной
# сортировки).

# ВАЖНО: нули это не баг, а костыль, чтобы сохранять одинаковое количество контигов
# на результат не влияет, так как в статистике и визуализации они не рисуются

###############################################################################

def get_coord_from_dat(trf_path):
    print("################################")
    print("Шаг 1.0, парсинг координат из .dat")
    print("################################")
    print("Starting...")
    
    with open(trf_path) as input_file, open('output.dat', 'w') as output_file:
        lines = input_file.readlines()[6:]
        for line in lines:
            if not re.match(r'^\s*$', line) and not 'Parameters:' in line:
                output_file.write(line)
    
    repeat_input = open('output.dat', "r")
    repeat_output = open('output.txt', 'a+')
    print("Proceed, please wait...")

    repeat = repeat_input.readline()
    repeat = repeat.split()
    while True:
        if not repeat:
            break
        if "Sequence:" == repeat[0]:
            repeat_output.write((str(repeat[0]) + " " + str(repeat[1]) + "\n"))
            while True:
                repeat = repeat_input.readline()
                repeat = repeat.split()
                if not repeat:
                    break
                if ("Sequence:" != repeat[0]) and ("TTAGGG"*4 in repeat[-1] or "ttaggg"*4 in repeat[-1]):
                    repeat_output.write((str(repeat[0]) + " " + str(repeat[1]) + " forward" + "\n"))
                if ("Sequence:" != repeat[0]) and ("CCCTAA"*4 in repeat[-1] or "ccctaa"*4 in repeat[-1]):
                    repeat_output.write((str(repeat[0]) + " " + str(repeat[1]) + " reverse" + "\n"))
                if "Sequence:" == repeat[0]:
                    repeat_output.write("0 0\n")
                    break
            
    print("Done!")        
    repeat_input.close()
    repeat_output.close()
    os.remove('output.dat')

###############################################################################

# Шаг 1.1, парсинг координат из .trf
# Подоходит только для сборок, в которых контиги из .fna в порядке возрастания
# и необходимо их отсортировывать из случайного порядка .trf

###############################################################################

#забрать координаты из trf в порядке их нахождения в документе

def get_coord_from_trf(trf_path):
    print("################################")
    print("Шаг 1.0, парсинг координат из .dat")
    print("################################")
    print("Starting...")
    repeat_input = open(trf_path, "r")
    repeat_output = open('coord.txt', 'a+')
    print("Proceed, please wait...")


    while True:
        repeat = repeat_input.readline().split()
        if not repeat:
            break
        if "ttaggg"*4 in repeat[14]:
            repeat_output.write("Sequence: " + repeat[18]+"\n")
            repeat_output.write(repeat[6] + " " + repeat[7] + " forward" + "\n")
        if "ccctaa"*4 in repeat[14]:
            repeat_output.write("Sequence: " + repeat[18]+"\n")
            repeat_output.write(repeat[6] + " " + repeat[7] + " reverse" + "\n")
            
    print("Done!")        
    repeat_input.close()
    repeat_output.close()

    #Вывод текстового файла в словарь и очистка от пустых нулевых координат
    #ВАЖНО: Если на какой-то хромосоме не будет найдено повторов, визуализация поплывет вообще прям совсем,
    #поэтому нужно проверять вывод на соответствие

    data = open("coord.txt", "r")

    dic = {}
    coordinate = data.readline().split()
    while True:
        if not coordinate:
                break
        if coordinate[0] == "Sequence:":
            key = coordinate[1]
            dic.setdefault(key,[])
            while True:
                coordinate = data.readline().split()
                print(coordinate)
                dic[key].append(coordinate)
                coordinate = data.readline().split()
                if not coordinate:
                    break
                if coordinate[0] == "Sequence:":
                    break

            
    clear_dic = (dict(filter(lambda x:x[1], dic.items())))
    print(clear_dic)
    with open("sorted.txt", "w") as f:
        for key, value in clear_dic.items():
            f.write(key + ": " + str(value) + "\n")

    data.close()

    # из файла output.txt создастся "sorted.txt"
    # Выглядит вот так: CP100568.1: [['1', '10482', 'reverse'], ['15824212', '16059225', 'reverse']]

    # его нужно отсортировать по названию в порядке возрастания:
    # awk '{print $0}' sorted.txt | sort -k1 > sorted_file.txt

    #with open('sorted_file.txt', 'w') as output_file:
    #    subprocess.run(['awk', '{print $0}', 'sorted.txt'], stdout=subprocess.PIPE, text=True) \
    #        .stdout.pipe(['sort', '-k1'], stdout=output_file)
            
    # потом из sorted_file.txt нужно убрать скобки:
    # sed "s/\[//g; s/\]//g; s/'//g; s/,//g" sorted_file.txt > output_file.txt

    with open('sorted.txt') as input_file, open('clear.txt', 'w') as output_file:
        for line in input_file:
            line = re.sub(r'[\[\]\',]', '', line)
            output_file.write(line)
            
    with open('clear.txt') as input_file, open('output_file.txt', 'w') as output_file:
        lines = input_file.readlines()
        sorted_lines = sorted(lines, key=lambda x: x.split()[0])
        output_file.writelines(sorted_lines)

    # в итоге получается вот так:
    # CP100559.1: 20 4686 reverse 59453345 59473045 forward

    # Дальше идет парсинг output_file.txt в fin_res.txt, который имеет блочный вид типа:
    # Sequence: CP100555.1:
    # 15 7467 reverse
    # 9382219 9382281 forward
    # 11001358 11031560 forward
    # 19330916 19449732 forward

    with open('output_file.txt', 'r') as input_file, open('output.txt', 'w') as output_file:
        for line in input_file:
            parts = line.strip().split(' ')
            sequence_id = parts[0]
            num_triplets = (len(parts) - 1) // 3
            output_lines = []
            for i in range(num_triplets):
                start = parts[1 + 3*i]
                end = parts[2 + 3*i]
                direction = parts[3 + 3*i]
                output_lines.append(f"{start} {end} {direction}")
                if i == 0:
                    output_lines.insert(0, f"Sequence: {sequence_id}")

            output_file.write('\n'.join(output_lines) + '\n')

    # В итоге имеем на выходе аутпут блочного вида(output.txt), который используется дальше

    # ВРЕМЕННО: очистка всех насоздававшихся здесь txtшников:
    os.remove('clear.txt')
    os.remove('coord.txt')
    os.remove('output_file.txt')
    os.remove('sorted.txt')

###############################################################################

# Шаг 2, КОСТЫЛЬ
# Слепить контиги/хромосомы/етц в одну большую строку на каждый контиг
# желательно бы выпилить его, так как жрет очень много места на жд, но пока как есть 

###############################################################################

def contigs_join(contig_path):
    print("################################")
    print("Шаг 2, костылирование контигов")
    print("################################")
    
    contig = open(contig_path, "r")
    output = open('contig_output.txt', 'a+')
    print("Starting contigs filtering...")
    i = 0
    while True:
        i+=1
        print("In process: ", i)
        sequence = contig.readline()
        if not sequence:
            break
        if  sequence[0] != ">":
            while True:
                sequence = sequence.replace("\r","")
                sequence = sequence.replace("\n","")
                output.write("" + sequence)
                sequence = contig.readline()
                if not sequence:
                    break
                if sequence[0] == ">":
                    output.write("\n")
                    break

    print("Done!")
    contig.close()
    output.close()
    
    
    
###############################################################################

# Шаг 2.1, КОСТЫЛЬ ВЫСОКОПРОИЗВОДИТЕЛЬНЫЙ
# Вместо того чтобы слеплять просто берет длину строки
# 

###############################################################################

def contigs_len(contig_path):
    print("################################")
    print("Шаг 2.1, ВЫСОКОПРОИЗВОДИТЕЛЬНОЕ костылирование контигов")
    print("Берет только длину контигов")
    print("################################")
    
    contig = open(contig_path, "r")
    output = open('contig_output.txt', 'a+')
    print("Starting contigs filtering...")
    
    length = 0
    while True:
        sequence = contig.readline()
        if not sequence:
            if length != 0:
                output.write(str(length) + "\n")
            break
        if sequence[0] == ">":
            if length != 0:
                output.write(str(length) + "\n")
            length = 0
        if  sequence[0] != ">": 
            length += len(sequence)-1



    print("Done!")
    contig.close()
    output.close()   
    
###############################################################################

# Шаг 3, налепить на контиги координаты
# Инпут создается автоматически, поэтому аргументы не принимает
# на выход выдает result.txt, в котором каждый повтор помечен маленькими буквами
# В ПРИНЦИПЕ можно удалить и не использовать если переделать output.tx чтобы в нем все было сразу

###############################################################################

def make_result():
    print("################################")
    print("Шаг 3, вывод результата")
    print("################################")
    
    coord = open("output.txt", "r")
    contig = open("contig_output.txt", "r")
    result = open("result.txt", "a+")

    print("Starting...")
    coordinate = coord.readline().split()
    while True:
        print("In process..")
        if not coordinate:
            break
        if coordinate[0] == "Sequence:":
            result.write(coordinate[1]+"\n")
            coordinate = coord.readline().split()
            seq = contig.readline()
            if not seq:
                break
            if not coordinate:
                break
            if coordinate[0] != "0" and coordinate[1] !="0":
                while True:
                    if coordinate[0] == "Sequence:":
                        break           
                    c1 = int(coordinate[0])
                    c2 = int(coordinate[1])
                    seq = seq.replace(seq[c1:c2], seq[c1:c2].lower())
                    coordinate = coord.readline().split()
                    if not coordinate:
                        break
                result.write(""+seq)
            else:
               coordinate = coord.readline().split() 
    print("Done!")

    coord.close()
    contig.close()
    result.close()
    os.remove('contig_output.txt')
 
#в результате трех шагов остается только output.txt с координатами и result.txt с результатами    
 
    
###############################################################################


#Шаг 4 визуализация

#сделать выбор для разного масштаба в виде цифр
###############################################################################

def visualization(table_title):
    print("################################")
    print("Шаг 4, отрисовка")
    print("################################")
    def get_origin_data():
        x_temp = []
        y_temp = []
        file = open("result.txt", "r")#всегда result.txt
        while True:
            inp = file.readline()
            if not inp:
                break
            y_temp.append(inp)
            inp = file.readline()
            if not inp:
                break
            x_temp.append(int(inp)) #если перестало работать заменить inp на len и использовать шаг 2 вместо 2.1
        file.close()
        print(x_temp)
        print(y_temp)
        return(x_temp, y_temp)


    def get_marker_data():
        dic = {}
        with open("output.txt", "r") as coord:
            key = None
            for line in coord:
                coordinate = line.split()
                if coordinate[0] == "Sequence:":
                    key = coordinate[1]
                    dic[key] = []
                elif key is not None and len(coordinate) >= 2 and all(c.isdigit() for c in coordinate[:2]):
                    if coordinate[0] != "0" and coordinate[1] != "0":
                        dic[key].append(coordinate)
        clear_dic = (dict(filter(lambda x: x[1], dic.items())))
        return clear_dic

    def change_scale(x_mar, x0_mar, contig_size):
        for i in range(len(x_mar)):
            x_mar[i]-=x0_mar[i]
            x_mar[i] = x_mar[i] + 3 * (contig_size / 100) #изменение масштаба отрисовки маркера
            x0_mar[i] = x0_mar[i] - 3 * (contig_size / 100)
            if x0_mar[i] < 0:
                x0_mar[i] = 0
            if x_mar[i] > contig_size:
                x_mar[i] = contig_size
        return(x_mar, x0_mar)

    #visualize origin data
    origin_data = get_origin_data()
    x_temp = origin_data[0]
    y_temp = origin_data[1]
    fig = go.Figure()
    fig.add_trace(
    go.Bar(
        name = "Contigs",
        x = x_temp,
        y = y_temp,
        marker=go.bar.Marker(
            color="rgb(253, 240, 54)",
            line=dict(color="rgb(0, 0, 0)",
                        width=2)
        ),
        orientation="h"
        )
    )

    #visualize marker data
    marker_data = get_marker_data()
    j = 0
    for key in marker_data:
        print(key)
        print(int(x_temp[j]))
        y_mar = [y_temp[j]]
        coord = marker_data[key]
        i = 0
        for item in marker_data[key]:
            print(coord[i])
            c = coord[i]
            scale = change_scale([int(c[1])], [int(c[0])], int(x_temp[j]))
            x0_mar = scale[1]
            x_mar = scale[0]
            fig.add_trace(
                go.Bar(
                    name = key + " marker",
                    base = x0_mar,
                    x = x_mar,
                    y = y_mar,
                    marker=go.bar.Marker(
                        color="rgb(0, 176, 243)",
                        line=dict(color="rgb(0, 0, 0)",
                              width=0)
                    ),
                    orientation="h",
                    #showlegend = False
                    )
                )
            i += 1 
        j += 1
        
    # update layout properties
    fig.update_layout(
        autosize=True,
        height=1600,
        width=1400,
        bargap=0.15,
        bargroupgap=0.1,
        barmode="stack",
        hovermode="x",
        margin=dict(r=20, l=300, b=75, t=125),
        title=(table_title),########
    )

    fig.show()
    fig.write_html(table_title + '.html')#######

###############################################################################


#Шаг 5 статистика


###############################################################################

def statistika():
    def get_origin_data():
        x_temp = []
        y_temp = []
        file = open("result.txt", "r")
        while True:
            inp = file.readline()
            if not inp:
                break
            y_temp.append(inp)
            inp = file.readline()
            if not inp:
                break
            x_temp.append(int(inp)) #если перестало работать заменить inp на len и использовать шаг 2 вместо 2.1
        file.close()
        return(x_temp, y_temp)

    def get_coordinate_data():
        coord = open("output.txt", "r")
        dic = {}
        coordinate = coord.readline().split()
        while True:
            if not coordinate:
                break
            if coordinate[0] == "Sequence:":
                key = coordinate[1]
                dic.setdefault(key,[])
                coordinate = coord.readline().split()
                if not coordinate:
                    break
                if coordinate[0] != "0" and coordinate[1] !="0":
                    while True:
                        if coordinate[0] == "Sequence:":
                            break
                        if coordinate[0] != "0" and coordinate[1] !="0":
                            dic[key].append(coordinate)
                        coordinate = coord.readline().split()
                        if not coordinate:
                            break
                else:
                   coordinate = coord.readline().split()
        coord.close()
        clear_dic = (dict(filter(lambda x:x[1], dic.items())))
        return(clear_dic)

    origin_data = get_origin_data()
    contig_len = origin_data[0] # [248387329, 242696753, 201105949]
    contig_name = origin_data[1] # ['NC_060945.1\n', 'NC_060926.1\n', 'NC_060944.1\n']
    coordinate_data = get_coordinate_data() #value = coordinate_data[key]
    machine = open("Machine.txt", "a+")
    
    print("Статистика повторов: ")

    info = {}
    contig_count = 0
    locus_count = 0
    locus_sum = 0
    max_locus = 0
    delta_arr = []
    edge = 0
    non_edge = 0
    forward = 0
    reverse = 0
    i = 0
    for key in coordinate_data:
        print("Sequence ", key)
        machine.write("Sequence " + key + "\n")
        info.setdefault(key,[])
        data = coordinate_data[key]
        length = contig_len[i]
        contig_count += 1
        for value in data:
            if str(value[2]) == "wtf":
                contig_count -= 1
                continue
            locus_count += 1
            delta = int(value[1]) - int(value[0])
            delta_arr.append(delta)
            locus_sum += delta
            if delta > max_locus:
                max_locus = delta
            nach = int(value[0])
            nach_per = (nach/length)*100
            kon = (length - int(value[1]))
            kon_per = (kon/length)*100
            if kon_per <= 1 or nach_per <= 1:
                edge += 1
            else:
                non_edge +=1
            if str(value[2]) == "reverse":
                reverse += 1
            elif str(value[2]) == "forward":
                forward +=1
            #print(f"Длина повтора: {delta}, направление {value[2]}, расстояние от начала {nach}/{nach_per:.2f}%, расстояние от конца {kon}/{kon_per:.2f}% ")
            #print(f"Repeat length: {delta}, direction {value[2]}, distance from start {nach}/{nach_per:.2f}%, distance from end {kon}/{kon_per:.2f}%")
            line =f"{delta} {value[2]} {nach} {nach_per:.2f} {kon} {kon_per:.2f}%\n"
            machine.write(line)
            print(f"{delta} {value[2]} {nach} {nach_per:.2f} {kon} {kon_per:.2f}%")
            app = [delta, value[2], nach, nach_per, kon, kon_per]
            info[key].append(app)
        
        i += 1
    mean_locus = statistics.mean(delta_arr)
    stdev = statistics.stdev(delta_arr)
    median_locus = statistics.median(delta_arr)
    mad = statistics.median([abs(x - median_locus) for x in delta_arr])
    print("Количество контигов с теломерными повторами: ", contig_count)
    print("Количество найденных локусов: ", locus_count)
    print("Суммарная длина всех локусов: ", locus_sum)
    print("Максимальная длина всех локусов: ", max_locus)
    print("Средняя длина локуса: ", mean_locus, " среднеквадратическое отклонение: ", stdev)
    print("Медиана длины: ", median_locus, " медианное абсолютное отклонение: ", mad)
    print("Располагаются на концах хромосом:", edge)
    print("Располагаются в середине хромосом:", non_edge)
    print("Прямые чтения: ", forward, " обратные чтения: ", reverse)
    machine.close()

###############################################################################


#Шаг 6 таблицы

#сделать определение размеров таблицы в зависимости от количества строк
###############################################################################

def tables():
    data = []
    with open('Machine.txt', 'r') as f:
        for line in f:
            row = line.strip().split()
            if row[0] != "Sequence":
                data.append(row)
            if row[0] == "Sequence":
                data.append([row[0], row[1], '', '', '', ''])
            
    row_template = '{:<8} {:<9} {:<10} {:<9} {:<10} {:<10}'
    # Вывод заголовка таблицы
    header = row_template.format('length', 'direction', 'start_dist', 'start_per', 'end_dist', 'end_per')

    # Вывод строк таблицы
    for row in data:
        print(row_template.format(*row))
        
    table_html = tabulate(data, headers='firstrow', tablefmt='html', stralign="center")
    table_html = table_html.replace('<table>', '<table border="1" cellpadding="5" style="border: 2px solid black;">')
    print(table_html)
    with open('table.html', 'w') as f:
        f.write(table_html)


    # устанавливаем размер шрифта
    font_size = 16

    # создаем шрифт
    font = ImageFont.truetype('arial.ttf', font_size)

    # определяем размеры таблицы
    table_width = 1000
    table_height = 3100

    # создаем белый фон
    background_color = (255, 255, 255)

    # создаем изображение с белым фоном
    img = Image.new('RGB', (table_width, table_height), background_color)

    # создаем объект ImageDraw для рисования на изображении
    draw = ImageDraw.Draw(img)

    # устанавливаем координаты верхнего левого угла таблицы
    x = 5
    y = 5

    # рисуем заголовок таблицы
    headers = data[0]
    header_x = x
    header_y = y
    for i, header in enumerate(headers):
        draw.text((header_x, header_y), header, font=font, fill=(0, 0, 0))
        header_x += table_width / len(headers)

    # рисуем горизонтальные линии
    line_y = y + 2 * font_size
    for i in range(1, len(data)):
        draw.line((x, line_y, x + table_width, line_y), fill=(0, 0, 0), width=1)
        line_y += font_size

    # рисуем вертикальные линии
    line_x = x
    for i in range(len(headers) + 1):
        draw.line((line_x, y, line_x, line_y), fill=(0, 0, 0), width=1)
        line_x += table_width / len(headers)

    # выводим данные в ячейки
    text_y = y + 2 * font_size
    for row in data[1:]:
        text_x = x
        for i, cell in enumerate(row):
            draw.text((text_x, text_y), cell, font=font, fill=(0, 0, 0))
            text_x += table_width / len(headers)
        text_y += font_size

    # Сохраняем изображение
    img.save("table.png")

###############################################################################
#отладка
# примечание: на курином геноме для корректной отрисовки нужно добавить недостающие хромосомы и заново получить результат

#get_coord_from_trf("GCA_024206055.trf")
#contigs_len("GCA_024206055.fna")
#make_result()
#visualization("Chicken_genome")
statistika()
tables()
