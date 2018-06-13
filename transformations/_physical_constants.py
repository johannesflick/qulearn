"""Physical and mathematical constants.

   Contains
   --------

          SIGMA_X (2X2 array): Pauli matrix
          SIGMA_Y (2X2 array): Pauli matrix
          SIGMA_Z (2X2 array): Pauli matrix
          SIGMA_PLUS (2X2 array): Pauli matrix
          SIGMA_MINUS (2X2 array): Pauli matrix

"""
#author tjdim, Cambridge MA 2018
from scipy.sparse import diags
import numpy as np

# flip array vertically (axis=0)
SIGMA_X = np.flip(diags([1, 1]).toarray(), 0)
SIGMA_Y = np.flip(diags([1j, -1j]).toarray(), 0)
SIGMA_Z = diags([1, -1], 0).toarray()
SIGMA_PLUS = 0.5 * (SIGMA_X - 1j*SIGMA_Y)
SIGMA_MINUS = 0.5 * (SIGMA_X + 1j*SIGMA_Y)
