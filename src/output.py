"""Methods related to producing model output"""
import settings as s
import pdb
import inputs
def output(type=None):
    """Output to file.  Currently only option is to output
    all major species at all
    altitudes.
    """
    print("Writing output at {:6.1f}h".\
        format((s.runTime-inputs.startTime).total_seconds()/3600.))
    days, remainder = \
        divmod((s.runTime-inputs.startTime).total_seconds(),\
        s.nSecondsInDay)
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
    outfile.write("#PlanetDistance: {:6.3f}\n".\
        format(s.orbitalDistance))
    outfile.write("#Alt\t[O2]\t[O]\t[O3]\n")

    for iAlt in range(s.nLayers):
        if iAlt == 97:
            print(s.O3[iAlt])
        outfile.write('{:05.2f}\t{:09.7e}\t{:09.7e}\t{:09.7e}\n'\
            .format(s.Altitude[iAlt]/1e5,\
            s.O2[iAlt],s.O[iAlt],s.O3[iAlt]))


    outfile.close()
    return 0
