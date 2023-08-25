filename = "/home/nfs_data/zhanggh/SparseWS/data/sizes0.txt"
file = open(filename, "r")
l = [int(i) for i in file.read().split(" ")[:-1]] 
keys = set(l)
for i in keys:
    cnt = l.count(i)
    print(i, cnt, round(cnt/len(l),4))