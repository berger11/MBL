from flint import acb
import numpy as np

def string_to_acb(x):
    # Takes as argument a string x, where x is str(y) where y is an acb object
    # Returns the original acb number y
    # Does not conserve error bounds

    sp = x.split(" + ")
    y = acb(real=sp[0][1:-1], imag=sp[1][1:-1])

    return y

def _spectrum_to_string(L, w, comp_vecs):
    # Converts the spectrum, and alternatively eigenvectors, to two arrays of strings such that it can be pickled
    # Returns a numpy array of the spectrum if comp_vecs=False
    # Returns a tuple of the spectrum and the eigenvector matrix if comp_vecs=True

    if not comp_vecs:
        spectrum = np.zeros((4 * L), dtype=object)  # There are 4*L eigenvalues of the shape matrix
        for i in range(4 * L):
            spectrum[i] = str(w[i])  # Saves eigenvalues as strings since acb objects cant be pickled
        eigens = spectrum

    else:
        spectrum = np.zeros((4 * L), dtype=object)
        eigenvectors = np.zeros((4 * L, 4 * L), dtype=object)
        for i in range(4 * L):
            spectrum[i] = str(w[0][i])
            for j in range(4 * L):
                eigenvectors[i, j] = str(w[1][i, j])
        eigens = (spectrum, eigenvectors)

    return eigens