from flint import acb

def string_to_acb(x):
    # Takes as argument a string x, where x is str(y) where y is an acb object
    # Returns the original acb number y
    # Does not conserve error bounds

    sp = x.split(" + ")
    y = acb(real=sp[0][1:-1], imag=sp[1][1:-1])

    return y