import sys
[_, fai, window] = sys.argv
window = int(window)
with open(fai, 'r') as f:
    for line in f:
        dat = line.split("\t")
        CHR = dat[0]
        END = int(dat[1])
        i = 1
        while True:
            print(CHR + '\t' + str(i) + '\t' + str(i + window - 1))
            i += window
            if i + window - 1 >= END:
                break
