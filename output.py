import settings as s

def output(type=None):
    '''Output to file.  Default is major species at all altitudes'''
    days, remainder = divmod(s.totaltime,s.nSecondsInDay)
    hours, remainder = divmod(remainder, s.nSecoundsInHour)
    minutes, seconds = divmod(remainder, s.nSecondsInMinute)
    cTime = \
     "output{:05d}_{:02d}{:02d}{:02d}".format(int(days),\
     int(hours),int(minutes),int(seconds))

    file = "data/{}.dat".format(cTime)
    try:
        outfile = open(file,'w')
    except:
        print("----Error: Issue opening file {}".format(file))
        print("----Does directory and file exist?")
        exit(1)

    outfile.write("#Alt\t[O2]\t[O]\t[O3]\n")

    for iAlt in range(s.nLayers):

        outfile.write('{:05.2f}\t{:09.3e}\t{:09.3e}\t{:09.3e}\n'\
            .format(s.Altitude[iAlt]/1e5,s.O2[iAlt],s.O[iAlt],\
            s.O3[iAlt]))


    outfile.close()
    return 0
