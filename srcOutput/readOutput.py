from numpy import array
def readOutputFile(file):
    f = open(file,'r')
    data = []
    started = False
    for line in f:
        if started:
            data.append([float(item) for item in line.split("\t")])
        if line[0:4] == '#Alt':
            started = True
            header = line.strip('#').strip("\n").split("\t")
            nvars = len(header)
        if line[0:7] == "#Planet":
            t = line.split()
            planetDistance = float(t[1])
    data = array(data)
    f.close()

    header.append("r={}".format(planetDistance))
    return data,header
