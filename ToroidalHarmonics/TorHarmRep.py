# Written by Vassili Savinov on 06/11/2018
# class to handle repreentations via toroidal harmonics

# it is initialized from N*M*E array
# the first two indices are for theta and phi
# and must give a good grid
# the last index is for eta, and can have just few steps

# I will also allow to initialize it by providing call-back functions
# to functions in Carthesian and toroidal coordinates

class TorHarmRep:

    # torCoordRep=N*M*E array
    def __init__(self, torCoordRep):
        # check dimensions
        # do DFT2