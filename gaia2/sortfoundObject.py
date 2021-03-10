import os


def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)




filename = 'findObject.txt'
outfile = 'sortfoundObject.txt'
foldout = os.path.join('./','outSORT')
createFolder(foldout)
cats = ['CMC15', 'Pan-Starrs DR1', '2MASS All-Sky', 'SDSS DR12', 'URAT1', 'UCAC4']
logs = {}
for i in range(1, 64):
    logs['{:06d}'.format(int(bin(i).split('b')[1]))] = []
print(logs)
#print(logs.values())
#exit()
with open(filename, 'rt') as f:
    file_content = f.read()
lines = file_content.split('\n')
newlines = []
for i, line in enumerate(lines):
    if len(line) > 0 and (line[0] == 'G' or line[0] == 'R'):
        continue
    catsdata = line.split('|')
    code = 0
    newline = ''
    for index, a in enumerate(catsdata):
        #print(a[0])
        if index == 0:
            #print(a)
            newline = newline + a
        elif a[0] != 'n':
            #print(index, a)
            #print(a[0])
            newline = newline + '|' + a
            code = code + 10**(6-index)
        #print(index, a)

    for k, v in logs.items():
        if k == '{:06d}'.format(code):
            v.append(newline)
            logs[k] = v
    newline = newline + '|{:06d}'.format(code)
        #print(k,v)
    #logs['|s{:06d}'.format(code)] +=1
    #print(newline)
    #print(catsdata)
    #if i == 10:
    #    break
#print(logs)
for k, v in logs.items():
    if len(v) > 0:
        #print(k, len(v))
        headline = ''
        for index, a in enumerate(k):
            #print(a, cats[index])
            if a == '1':
             headline += cats[int(index)]+', '
        print(k)
        print(headline)
        with open(os.path.join(foldout, k+'.txt'), 'w') as f:
            f.write(headline+'\n')
            for line in v:
                f.write(line+'\n')

