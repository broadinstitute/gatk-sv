import sys


with open(sys.argv[1]) as f1, open(sys.argv[2]) as f2, open(sys.argv[3]) as f3, open(sys.argv[4]) as f4, open(sys.argv[5], "w") as f5:
    for x, y, z, t in zip(f1, f2, f3,f4):
          f5.write(x.strip() + "\t" + y.strip() + '\t' + z.strip() + '\t' + t.strip() + '\n')
