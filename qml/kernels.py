# MIT License
#
# Copyright (c) 2016 Anders Steen Christensen, Felix Faber
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import numpy as np
from numpy import empty, asfortranarray, ascontiguousarray, zeros

from .fkernels import fgaussian_kernel
from .fkernels import flaplacian_kernel
from .fkernels import fget_vector_kernels_gaussian
from .fkernels import fget_vector_kernels_laplacian


def laplacian_kernel(A, B, sigma):
    """ Calculates the Laplacian kernel matrix K, where K_ij:

            K_ij = exp(-1 * sigma**(-1) * || A_i - B_j ||_1)

        Where A_i and B_j are descriptor vectors.

        K is calculated using an OpenMP parallel Fortran routine.

        :param arg1: np.array of np.array of descriptors.
        :type arg1: 2D np.array (N, size), where N is the size of each representation.
        :param arg2: np.array of np.array of descriptors.
        :type arg1: 2D np.array (M, size), where M is the size of each representation.
        :param arg3: The value of sigma in the kernel matrix.

        :return: The Laplacian kernel matrix.
        :rtype: 2D np.array (N, M)
    """

    na = A.shape[0]
    nb = B.shape[0]

    K = empty((na, nb), order='F')

    # Note: Transposed for Fortran
    flaplacian_kernel(A.T, na, B.T, nb, K, sigma)

    return K


def gaussian_kernel(A, B, sigma):
    """ Calculates the Gaussian kernel matrix K, where K_ij:

            K_ij = exp(-0.5 * sigma**(-2) * || A_i - B_j ||_2)

        Where A_i and B_j are descriptor vectors.

        :param arg1: np.array of np.array of descriptors.
        :type arg1: 2D np.array (N, size), where N is the size of each representation.
        :param arg2: np.array of np.array of descriptors.
        :type arg1: 2D np.array (M, size), where M is the size of each representation.
        :param arg3: The value of sigma in the kernel matrix.

        :return: The Gaussian kernel matrix.
        :rtype: 2D np.array (N, M)
    """

    na = A.shape[0]
    nb = B.shape[0]

    K = empty((na, nb), order='F')

    # Note: Transposed for Fortran
    fgaussian_kernel(A.T, na, B.T, nb, K, sigma)

    return K
