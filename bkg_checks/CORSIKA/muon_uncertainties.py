
import numpy

from scipy import interpolate

def Muon_Primary_Spline(norm=False, kind='linear'):
    ''' 
    Constructs a muon primary uncertainty spline.

    There is a method to return this in a normalisation-preserving way. This 
    defaults to False because it's not exactly tested, but please try it! It is 
    expected that the .txt file is in a subdirectory called Uncertainties 
    wherever this python script is located with the name 
    Muon1SigmaUncertaintiesCosZenith.txt. This is how it was provided on SVN.

    There is a method to return this in a normalisation-preserving way. This 
    defaults to False because it's not exactly tested, but please try it!
    '''

    MuonUncs = open('Uncertainties/Muon1SigmaUncertaintiesCosZenith.txt')

    MuonZenithPoints = []
    MuonUncPoints = []

    for line in MuonUncs:
        MuonZenithPoints.append(float(line.rstrip().split(' ')[0]))
        MuonUncPoints.append(float(line.rstrip().split(' ')[1]))

    MuonZenithPoints = numpy.array(MuonZenithPoints)
    MuonUncPoints = numpy.array(MuonUncPoints)

    # Need to deal with zeroes. This may be a problem if this work ever
    # changes but is fine for the zenith spline here.
    while 0.0 in MuonUncPoints:
        ZeroIndices = numpy.where(MuonUncPoints == 0)[0]
        for ZeroIndex in ZeroIndices:
            MuonUncPoints[ZeroIndex] = MuonUncPoints[ZeroIndex+1]

    # Try to do this normalisation-preserving
    # Find the integral and then subtract an appropriate amount from each point
    if norm:
        BinWidth = MuonZenithPoints[1] - MuonZenithPoints[0]
        Integral = numpy.sum(MuonUncPoints*BinWidth)
        DiffFactor = Integral/(BinWidth*len(MuonUncPoints))
        MuonUncPoints = MuonUncPoints - DiffFactor

    # Add a dummy point for coszenith = 0
    MuonZenithPoints = numpy.insert(MuonZenithPoints,0,0.0)
    MuonUncPoints = numpy.insert(MuonUncPoints,0,MuonUncPoints[0])
    # Add a dummy poiny for coszenith = 1
    MuonZenithPoints = numpy.insert(MuonZenithPoints,-1,1.0)
    MuonUncPoints = numpy.insert(MuonUncPoints,-1,MuonUncPoints[-1])

    MuonUncF = interpolate.interp1d(MuonZenithPoints, MuonUncPoints, kind=kind)

    return MuonUncF

def Muon_Primary_Unc(cosZenith, MuonUncF=None, norm=False, kind='linear'):
    ''' 
    Takes the truth cosZenith of the muon and returns a percentage uncertainty 
    corresponding to that due to the uncertainty on the primary cosmic rays. 
    The construction of the spline is a different function so you only have to
    call it once. If you don't provide one then it will be made. It is expected
    that the .txt file is in a subdirectory called Uncertainties wherever this 
    python script is located with the name Muon1SigmaUncertaintiesCosZenith.txt.
    This is how it was provided on SVN.

    For details on the calculation please see 

    https://wiki.icecube.wisc.edu/index.php/DeepCore_Muon_Background_Systematics
    
    i.e. the section on Uncertainties from Primary Cosmic Rays

    There is a method to return this in a normalisation-preserving way. This 
    defaults to False because it's not exactly tested, but please try it!
    '''

    if MuonUncF is None:
        MuonUncF = Muon_Primary_Spline(norm=norm,kind=kind)

    return MuonUncF(cosZenith)

if __name__ == '__main__':
    '''
    This is a debug for the normalisation. 
    Make a plot see if it looks sensible.
    '''
    from matplotlib import pyplot
    pyplot.rcParams['text.usetex'] = True

    MuonUncs = open('Uncertainties/Muon1SigmaUncertaintiesCosZenith.txt')

    MuonZenithPoints = []
    MuonUncPoints = []

    for line in MuonUncs:
        MuonZenithPoints.append(float(line.rstrip().split(' ')[0]))
        MuonUncPoints.append(float(line.rstrip().split(' ')[1]))

    MuonZenithPoints = numpy.array(MuonZenithPoints)
    MuonUncPoints = numpy.array(MuonUncPoints)

    # Need to deal with zeroes. This may be a problem if this work ever
    # changes but is fine for the zenith spline here.
    while 0.0 in MuonUncPoints:
        ZeroIndices = numpy.where(MuonUncPoints == 0.0)[0]
        for ZeroIndex in ZeroIndices:
            MuonUncPoints[ZeroIndex] = MuonUncPoints[ZeroIndex+1]

    # Try to do this normalisation-preserving
    # Find the integral and then subtract an appropriate amount from each point
    BinWidth = MuonZenithPoints[1] - MuonZenithPoints[0]
    Integral = numpy.sum(MuonUncPoints*BinWidth)
    DiffFactor = Integral/(BinWidth*len(MuonUncPoints))
    NormMuonUncPoints = MuonUncPoints - DiffFactor

    # Add a dummy point for coszenith = 0
    MuonZenithPoints = numpy.insert(MuonZenithPoints,0,0.0)
    MuonUncPoints = numpy.insert(MuonUncPoints,0,MuonUncPoints[0])
    NormMuonUncPoints = numpy.insert(NormMuonUncPoints,0,NormMuonUncPoints[0])
    # Add a dummy poiny for coszenith = 1
    MuonZenithPoints = numpy.insert(MuonZenithPoints,-1,1.0)
    MuonUncPoints = numpy.insert(MuonUncPoints,-1,MuonUncPoints[-1])
    NormMuonUncPoints = numpy.insert(NormMuonUncPoints,-1,NormMuonUncPoints[-1])

    MuonUncF = interpolate.interp1d(MuonZenithPoints, MuonUncPoints)
    NormMuonUncF = interpolate.interp1d(MuonZenithPoints, NormMuonUncPoints)

    MuonZenithPointsNew = numpy.linspace(0,1,101)

    pyplot.plot(MuonZenithPoints,MuonUncPoints,'o',
                MuonZenithPointsNew,MuonUncF(MuonZenithPointsNew),'-',
                MuonZenithPointsNew,NormMuonUncF(MuonZenithPointsNew),'--')
    pyplot.legend(['data', 'standard', 'normalised'], loc='best')
    pyplot.xlabel(r'Muon $\cos\theta_Z$')
    pyplot.ylabel('Percentage Uncertainty')
    pyplot.savefig("splinesvisualised.png")
