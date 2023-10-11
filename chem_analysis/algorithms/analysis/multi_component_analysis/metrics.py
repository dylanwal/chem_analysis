
def mse(C, ST, D_actual, D_calculated):
    """ Mean square error """
    return ((D_actual - D_calculated)**2).sum()/D_actual.size
