from numpy import zeros

def inituser(nlayers):
    '''This function is for quick output of non-standard variables

    Usage:
        Change nvars to the appropriate number and import the
        userdata variable in the file that contains the data
        for output.
        data[0,:] is already filled with altitude data, so don't use index 0.
        Other indices should be filled by user in the appropriate file.
        Note that nvars should include altitude.

    '''

    nvars = 5

    return zeros((nvars,nlayers))


def outputuser(file,data):
    '''Output the USER file type'''

    ### This shouldn't change unless you know what you are doing!
    file.write("#Alt\tr(O3+hv)\tr(O2+hv)\tk(O2_O)\tkO3_O\n")
    for iAlt in range(len(data[0,:])):
        file.write('{:05.2f}\t'.format(data[0,iAlt]))
        file.write('\t'.join("{:09.7e}".format(d) for d \
            in data[1:,iAlt]))
        file.write('\n')
    return 0
