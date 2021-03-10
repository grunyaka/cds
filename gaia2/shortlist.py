infile1 = 'list.txt'
infile2 = 'list2.txt'
outfile = 'shortlist.txt'
with open(infile1, 'rt') as f:#это читаем входной файл по каждой звезде
    file_content = f.read()
lines = file_content.split('\n')
with open(infile2, 'rt') as f:#это читаем входной файл по каждой звезде
    file_content = f.read()
stars = file_content.split('\n')

indexdone = []
for line in lines:
    if len(line) < 1:
        continue
    RA = float(line.split(',')[5])
    indexes = [i for i in list(range(len(stars))) if not i in indexdone]
    for i in indexes:
        star = stars[i]
        if len(star) < 1:
            indexdone.append(i)
            continue
        if star[0] == 'G' or star[0] == 'R':
            indexdone.append(i)
            continue
        R = float(star.split('|')[0].split('\t')[0].strip())
        if abs(R - RA)<0.00000000000001:
            print(line)
            with open(outfile, 'a') as f:
                f.write(line + '\n')
            indexdone.append(i)
