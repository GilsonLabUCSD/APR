import numpy as np
import subprocess as sp
import sys, os
import apr_restraints as apr_restraints

def factors(n):
    # Return list of integer factors
    factor = []
    sqrt = int(round(np.sqrt(n) + 0.5))
    i = 1
    while i <= sqrt:
        if n % i == 0:
            factor.append(i)
            j = n / i
            if j != i:
                factor.append(j)
        i += 1
    return sorted(factor, key=int)


def nearest_integer_big_factor(n):
    # Nearest number with the most number of factors.
    # Check numbers down to 100 less than the input number and is divisible by two

    maxfac = 0
    if n % 2 == 0:
        low = n - 100
        high = n
    else:
        low = n - 101
        high = n - 1
    if low < 0:
        low = 0
    for i in range(low, high + 2, 2):
        numfac = len(factors(i))
        if numfac >= maxfac:
            maxfac = numfac
            mostfac = i
    return mostfac


def seom(N, arr):
    # Return the maximum value in the blocking curve
    Facs = factors(N)
    Bn = np.zeros([len(Facs)], np.int32)  # Number of blocks for a given block size
    Bmean = np.zeros([len(Facs), Facs[-1]], np.float64)  # Means for each block of each block size
    SEOM = np.zeros([len(Facs) - 2], np.float64)
    for i in range(len(Facs) - 2):  # Run over all block sizes except the final two: two blocks, one block.
        for j in range(Facs[-i - 1]):  # Run over all blocks in the data for a specific size.
            Bmean[i, j] = np.mean(arr[j * Facs[i]:(j + 1) * Facs[i]])
        Bn[i] = j + 1
        SEOM[i] = np.std(Bmean[i, 0:Bn[i]], ddof=0) / np.sqrt(Bn[i] - 1)
    return np.max(SEOM)  # Assume the blocking plateau corresponds to the max histogram count!


def interpolate(x, y, x_new, axis=-1, out=None):
    # Copyright (c) 2007-2015, Christoph Gohlke
    # Copyright (c) 2007-2015, The Regents of the University of California
    # Produced at the Laboratory for Fluorescence Dynamics
    # All rights reserved.
    #
    # Redistribution and use in source and binary forms, with or without
    # modification, are permitted provided that the following conditions are met:
    #
    # * Redistributions of source code must retain the above copyright
    #   notice, this list of conditions and the following disclaimer.
    # * Redistributions in binary form must reproduce the above copyright
    #   notice, this list of conditions and the following disclaimer in the
    #   documentation and/or other materials provided with the distribution.
    # * Neither the name of the copyright holders nor the names of any
    #   contributors may be used to endorse or promote products derived
    #   from this software without specific prior written permission.
    #
    # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    # AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    # ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
    # LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    # CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    # SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    # INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    # CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    # ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    # POSSIBILITY OF SUCH DAMAGE.
    x = np.array(x, dtype=np.float64, copy=True)
    y = np.array(y, dtype=np.float64, copy=True)
    xi = np.array(x_new, dtype=np.float64, copy=True)
    if axis != -1 or out is not None or y.ndim != 1:
        raise NotImplementedError("implemented in C extension module")
    if x.ndim != 1 or xi.ndim != 1:
        raise ValueError("x-arrays must be one dimensional")
    n = len(x)
    if n < 3:
        raise ValueError("array too small")
    if n != y.shape[axis]:
        raise ValueError("size of x-array must match data shape")
    dx = np.diff(x)
    if any(dx <= 0.0):
        raise ValueError("x-axis not valid")
    if any(xi < x[0]) or any(xi > x[-1]):
        raise ValueError("interpolation x-axis out of bounds")
    m = np.diff(y) / dx
    mm = 2.0 * m[0] - m[1]
    mmm = 2.0 * mm - m[0]
    mp = 2.0 * m[n - 2] - m[n - 3]
    mpp = 2.0 * mp - m[n - 2]
    m1 = np.concatenate(([mmm], [mm], m, [mp], [mpp]))
    dm = np.abs(np.diff(m1))
    f1 = dm[2:n + 2]
    f2 = dm[0:n]
    f12 = f1 + f2
    ids = np.nonzero(f12 > 1e-9 * np.max(f12))[0]
    b = m1[1:n + 1]
    b[ids] = (f1[ids] * m1[ids + 1] + f2[ids] * m1[ids + 2]) / f12[ids]
    c = (3.0 * m - 2.0 * b[0:n - 1] - b[1:n]) / dx
    d = (b[0:n - 1] + b[1:n] - 2.0 * m) / dx ** 2
    bins = np.digitize(xi, x)
    bins = np.minimum(bins, n - 1) - 1
    bb = bins[0:len(xi)]
    wj = xi - x[bb]
    return ((wj * d[bb] + c[bb]) * wj + b[bb]) * wj + y[bb]



def get_forces(method, restraint_list):
    ### Calculate the mean force/sem for a simulation window
    ### Determine Number Of Simulation Frames.  Assume a single line header.
    with open('./restraints.dat', 'r') as f:
        for i, line in enumerate(f):
            pass
        number_of_frames = i


    # We should have seven things whether jacks are used or not.
    # The last thing in the list is the number of jacks (this may be zero)
    # The 3 is for 1 distance and 2 angle restraint for the guest
    # The host translation and rotational restraint are not used here
    number_restraints = 3 + restraint_list[6]
    Reqs = np.zeros([number_restraints], np.float64)  # Restraint equilibrium (target) value
    Rfcs = np.zeros([number_restraints], np.float64)  # Restraint force constant
    Vals = np.zeros([number_of_frames, number_restraints], np.float64)  # Restraint coordinate values from simulation
    Forces = np.zeros([number_of_frames], np.float64)  # Forces

    # This is the distance restraint equilibrium target value
    Reqs[0] = restraint_list[0]
    # This is the distance restraint weight (force constant)
    Rfcs[0] = restraint_list[1]
    # This is the angle restraint equilibrium target value
    Reqs[1:3] = restraint_list[2]
    # This is the angle restraint weight (force constant)
    # This needs to be converted (rad**2 to degrees**2), but the target does not because of AMBER.
    Rfcs[1:3] = restraint_list[3] * (np.pi / 180.0) * (np.pi / 180.0)



    number_of_jacks = restraint_list[6]
    if number_of_jacks > 0:
        for index in range(3, number_of_jacks + 3):
            # This is the jacks distance equilibrium target value
            Reqs[index] = restraint_list[4]
            # This is the jacks distance equilibrium weight (force constant)
            Rfcs[index] = restraint_list[5]


    ### Read the restraint coordinate values (restraints.dat)
    with open('./restraints.dat', 'r') as f:
        frame = 0
        for line in f.readlines()[1:]:  # Skip header
            cols = line.split()
            for restraint in range(number_restraints):
                Vals[frame, restraint] = float(cols[restraint + 1])
            frame += 1

    ### Calculate everything in terms of TI "forces".
    if method == 'attachment':
        Forces = np.sum(Rfcs[0:number_restraints] * ((Vals[0:number_of_frames, 0:number_restraints] -
                                                      Reqs[0:number_restraints]) ** 2), axis=1)
    elif method == 'translation':
        Forces = 2.0 * Rfcs[0] * (Vals[0:number_of_frames, 0] - Reqs[0])
        # A hack, assumes the pulling restraint comes first in the array!
    elif method == 'release':
        # Since this only has to do with removing jacks restraints, we don't start at zero.
        Forces = np.sum(Rfcs[3:number_restraints] * ((Vals[0:number_of_frames, 3:number_restraints] -
                                                      Reqs[3:number_restraints]) ** 2), axis=1)


    ### Nearest integer to the total number of data points that has the maximum number of factors
    Nm = nearest_integer_big_factor(number_of_frames)

    ### Calculate Mean and SEM of Forces with the blocking approach
    mean = np.mean(Forces)
    blkmean = np.mean(Forces[0:Nm])
    blksem = seom(Nm, Forces[0:Nm])
    blksti = np.var(Forces[0:Nm]) / (blksem ** 2)  # Statistical Inefficiency

    return [number_of_frames, mean, Nm, blkmean, blksem, blksti]



def error_estimate(scale_w, method, restraint_list):
    frcvals = get_forces(method, restraint_list)
    return scale_w * frcvals[4]


