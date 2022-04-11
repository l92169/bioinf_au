file1 = open("C:/Users/liza_/Desktop/AU/2course/bioinf/output.fasta", "r")
file2 = open("C:/Users/liza_/Desktop/AU/2course/bioinf/output_startstop.fasta", "w")
lines = file1.readlines()
minn = []
maxx = []
minnn = []
maxxx = []
names = []
for line in lines:
    if len(line) > 3:
        result = line.split()
        i = -1
        if result[0] == "Sbjct":
            minnn.append(int(result[1]))
            maxxx.append(int(result[3]))
        if result[0][0] == ">":
            names.append((line[1:-1]).replace(" ", "_"))
            if minnn != []:
                minn.append(minnn)
                maxx.append(maxxx)
                i+=1
                minnn = []
                maxxx = []
minn.append(minnn)
maxx.append(maxxx)
i+=1
minnn = []
maxxx = []
print(len(names),len(minn))
for i in range(len(names)):
    if min(minn[i]) < max(maxx[i]):
        file2.write("{0}{1}{2}".format(str(names[i])+"\t", (str(min(minn[i])))+"\t", (str(max(maxx[i])))+"\t.\t.\t+")  + "\n")
    else:
        file2.write("{0}{1}{2}".format(str(names[i]) +"\t", (str(max(maxx[i])))+ "\t", (str(min(minn[i]))) + "\t.\t.\t-") + "\n")
file1.close
file2.close
