"""Methods related to producing model output"""
import settings as s
import inputs
import user

def output(path_input,fin):
    """Output to file based on output types specified in input file.
    """

    for type in inputs.cOutputType:
        print("Writing output type: {} at {:6.1f}h".\
            format(type,(s.runTime-inputs.startTime).total_seconds()/3600.))
        days, remainder = \
            divmod((s.runTime-inputs.startTime).total_seconds(),\
            s.nSecondsInDay)
        hours, remainder = divmod(remainder, s.nSecoundsInHour)
        minutes, seconds = divmod(remainder, s.nSecondsInMinute)
        cTime = \
         "{:05d}_{:02d}{:02d}{:02d}".format(int(days),\
         int(hours),int(minutes),int(seconds))
        if len(fin)>0:
            file = path_input+"/final.dat".format(type,cTime)
        else:
            file = path_input+"/{}_{}.dat".format(type,cTime)
        try:
            outfile = open(file,'w')
        except:
            print("----Error: Issue opening file {}".format(file))
            print("----Does directory and file exist?")
            exit(1)

        outfile.write("#PlanetDistance: {:6.3f}\n".\
            format(s.orbitalDistance))

        if type.lower() == "all":
            outfile.write("#Alt\t[O]\t[O2]\t[O3]\tOvmr\tO2vmr\tO3vmr\n")

            for iAlt in range(s.nLayers):

                outfile.write('{:05.2f}\t'.format(s.Altitude[iAlt]/1e5))
                outfile.write('\t'.join("{:09.7e}".format(d) for d \
                    in s.density[:,iAlt]))
                outfile.write('\t')
                outfile.write('\t'.join("{:09.7e}".format(d) for d \
                                            in s.density[:,iAlt]/s.N[iAlt]))
                outfile.write('\n')


        if type.lower() == "photo":

            outfile.write("#Alt\tJ(O3)\tJ(O2)\n")

            for iAlt in range(s.nLayers):

                outfile.write('{:05.2f}\t'.format(s.Altitude[iAlt]/1e5))
                outfile.write('\t'.join("{:09.7e}".format(d) for d \
                    in s.PhotoDissRate_Alt[:,iAlt]))
                outfile.write('\n')

        if type.lower() == "user":
            s.userdata[0,:] = s.Altitude/1.0e5
            iError = user.outputuser(outfile,s.userdata)

        outfile.close()
    return 0
