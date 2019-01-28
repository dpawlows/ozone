from numpy import array
def readOutputFile(file):
    f = open(file,'r')
    data = []
    started = False
    for line in f:
        if started:
            data.append([float(item) for item in line.split("\t")])
        if line[0] == '#':
            started = True
            header = line.strip('#').strip("\n").split("\t")
            nvars = len(header)

    data = array(data)
    f.close()

    return data,header
