import re
import os

def join_tel_prev_version_dont_work_correctly():
    with open('test_output_trf.txt', 'r') as f:
        content = f.read()
        sequences = content.split('>')[1:]
        with open('final_trf.txt', 'a+') as output_file:
            for seq in sequences:
                seq_name = seq.split('\n')[0]
                seq_lines = seq.split('\n')[1:]
                output_file.write(f">{seq_name}\n")
                print(f">{seq_name}\n Proceed...\n")
                prev_start = None
                prev_end = None
                saved_start = None
                saved_end = None
                saved_line = ''
                #for i in range(len(seq_lines)):
                i = 0
                while i < len(seq_lines):
                    line = seq_lines[i]
                    parts = seq_lines[i].split()
                    if len(parts) > 0:
                        start = int(parts[0])
                        end = int(parts[1])
                        seq_type = parts[2]
                        if prev_end is not None and start - prev_end < 1000:
                            saved_start = start
                            while start - prev_end < 1000:
                                prev_start = start
                                prev_end = end
                                saved_end = end
                                i += 1
                                parts = seq_lines[i].split()
                                #print(parts)
                                if parts == []:
                                    break
                                start = int(parts[0])
                                end = int(parts[1])
                            length = saved_end - saved_start
                            saved_line = f"{saved_start} {saved_end} {seq_type} {length} new\n"
                            output_file.write(saved_line)
                        else:
                            output_file.write(line + "\n")
                    prev_start = start
                    prev_end = end
                    i += 1

def find_tel():
    with open('GCA_024206055.fna', 'r') as f:
        content = f.read()

    #pattern = re.compile('(?:TTAGGG|CCCTAA|ttaggg|ccctaa){5,}')
    pattern = re.compile('(?:C{2,4}T{1,2}A{1,3}|T{1,3}A{1,2}G{2,4}|TTAGGG|CCCTAA|ttaggg|ccctaa){4,}')

    sequences = content.split('>')[1:]
    with open('test_output_trf.txt', 'a+') as output_file:
        for seq in sequences:
            seq_lines = seq.split('\n')[1:]  # Пропускаем первую строку с именем
            seq_name = seq.split('\n')[0]
            seq_string = ''.join(seq_lines)
            output_file.write(f">{seq_name}\n")
            print(f">{seq_name}\n Proceed...\n")
            matches = pattern.finditer(seq_string)
            for match in matches:
                direction = match.group()[:6]
                #if match.group()[:6] == "ttaggg":
                #    direction = "forward"
                #if match.group()[:6] == "ccctaa":
                #    direction = "reverse"
                start = match.start()
                end = match.end() 
                length = end - start
                #print(match.group()[:6], direction)
                #print(f"Паттерн найден в последовательности {seq_name} на позиции {start} - {end}")
                output_file.write(f"{start} {end} {direction} {length}\n")
                #print(f" {seq_name} на позиции {start} - {end} {direction}")
                    
def join_tel():
    with open('test_output_trf.txt', 'r') as f:
        content = f.read()
        sequences = content.split('>')[1:]
        with open('final_trf.txt', 'w') as output_file:
            for seq in sequences:
                seq_name = seq.split('\n')[0]
                seq_lines = seq.split('\n')[1:]
                output_file.write(f">{seq_name}\n")
                print(f">{seq_name}\n Proceed...\n")
                prev_end = 0
                saved_start = None
                saved_end = None
                saved_line = None
                for line in seq_lines:
                    parts = line.split()
                    if len(parts) > 0:
                        start = int(parts[0])
                        end = int(parts[1])
                        seq_type = parts[2]
                        if saved_start is None:
                            saved_start = start
                            saved_end = end
                            saved_line = line
                        elif start - prev_end < 1000:
                            saved_end = end
                            saved_line = f"{saved_start} {saved_end} {seq_type} {saved_end - saved_start} new"
                        else:
                            output_file.write(saved_line + "\n")
                            saved_start = start
                            saved_end = end
                            saved_line = line
                    prev_end = end

                if saved_line is not None:
                    output_file.write(saved_line + "\n")


                            

#find_tel()
join_tel()







