FCHL18 - Forces and Electric Fields
-----------------------------------

Representation
~~~~~~~~~~~~~~

Generate Representation
^^^^^^^^^^^^^^^^^^^^^^^

The representation for a molecule is a 3-d array of the dimensions (max_size, 5, max_neighbors).
It can only be used with the FCHL-specific kernel functions, as the distance is an analytically calculated integral.
There are several different representations available, which all share a few common keyword arguments.

The representation for a molecule can be generated with the following code:


.. code:: python

    import numpy as np
    import qml
    from qml.fchl import generate_representation


    # Dummy coordinates for a water molecule
    coordinates = np.array([[1.464, 0.707, 1.056],
                            [0.878, 1.218, 0.498],
                            [2.319, 1.126, 0.952]])

    # Oxygen, Hydrogen, Hydrogen
    nuclear_charges = np.array([8, 1, 1])

    rep = generate_representation(
        coordinates,        # np.array with coordinates
        nuclear_charges,    # np.array with nuclear charges as ints.
        max_size=23,        # Max number of atoms in any molecule
        cell=None,          # None for molecules, otherwise a 3x3 matrix with the unit cell vectors
        cut_distance=5.0,   # Cut-off distance
        neighbors=23,       # Max number of neighbours within the cutoff, only necessary when cell is not None
        )

    print(rep)
    print(rep.shape)


This code should print something like:

.. code:: python

    [[[ 0.00000000e+000  9.57016719e-001  9.57811046e-001 ...
        1.00000000e+100  1.00000000e+100  1.00000000e+100]
    ...
      [ 0.00000000e+000  0.00000000e+000  0.00000000e+000 ...
        0.00000000e+000  0.00000000e+000  0.00000000e+000]]]
    (23, 5, 23)


The representation for e.g. the first atom in the molecule can then be obtained as ``rep[0]``.


Generate Force Representation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To use the FCHL kernels for numerical forces, it is necessary to first generate representations for displaced geometries.
In addition the the keyword arguments for the normal FCHL representation, the representation takes the finite difference displacement (dx) as keyword. 0.005 should work well for most uses.

.. code:: python

    from qml.fchl import generate_displaced_representations

    rep = generate_displaced_representations(
        coordinates,        # np.array with coordinates
        nuclear_charges,    # np.array with nuclear charges as ints.
        dx=0.005,           # Displacement in Angstrom
        )


Since the numerical finite-difference method gives very small errors, it is possible to use 5-point finite differences to get derivatives very close to machine accuracy. To generate the representation for this kernel use the dedicated function:

.. code:: python

    from qml.fchl import generate_displaced_representations_5point

    rep = generate_displaced_representations_5point(
        coordinates,        # np.array with coordinates
        nuclear_charges,    # np.array with nuclear charges as ints.
        dx=0.005,           # Displacement in Angstrom
        )



Generate Electric Field Representation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Kernels that include the electric field needs a set of fictitious partial charges when the representation is generated.
Default behavior is to use the Gasteiger charge model as implemented in Open Babel.
This requires the ``pybel`` Python module to be available at runtime.
Alternatively, a list of fictitious partial charges can be supplied by the user:

.. code:: python

    from qml.fchl import generate_representation_electric_field

    rep1 = generate_representation_electric_field(
        coordinates, nuclear_charges,
        fictitious_charges="gasteiger", # Using the Gasteiger charge model
        )

    rep2 = generate_representation_electric_field(
        coordinates, nuclear_charges,
        fictitious_charges=[-0.5, 0.25, 0.25], # list of partial charges
        )


QML can use any charge model that is available from Open Babel via the keyword.
A list of valid models can be obtained via ``obabel -L charges``, e.g.:

.. code:: bash

    $ obabel -L charges
    eem    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Bultinck B3LYP/6-31G*/MPA
    eem2015ba    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf B3LYP/6-311G/AIM
    eem2015bm    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf B3LYP/6-311G/MPA
    eem2015bn    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf B3LYP/6-311G/NPA
    eem2015ha    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf HF/6-311G/AIM
    eem2015hm    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf HF/6-311G/MPA
    eem2015hn    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf HF/6-311G/NPA
    eqeq    Assign EQEq (charge equilibration) partial charges.
    fromfile    Assign charges from file containing {'atom-name', charge} pairs
    gasteiger    Assign Gasteiger-Marsili sigma partial charges
    mmff94       Assign MMFF94 partial charges
    none    Clear all partial charges
    qeq    Assign QEq (charge equilibration) partial charges (Rappe and Goddard, 1991)
    qtpie    Assign QTPIE (charge transfer, polarization and equilibration) partial charges (Chen and Martinez, 2007)


The valid names are the strings in the first column.


Standard Kernel
~~~~~~~~~~~~~~~

The standard kernels used for e.g. kernel ridge regression (KRR) come in three flavors for FCHL. These are described below.

Local Kernels
^^^^^^^^^^^^^

In the "local" kernel betweent two chemical compounds(e.g. `I` and :math:`J`), the kernel elements are sums over the pair-wise kernels between the atoms in the two molecules, that is:

    :math:`k(I,J) = \sum_{i \in I} \sum_{j \in J} k(i,j)`

This kernel is implemented in the `fchl` module. There are two version, one for kernels between two sets of compounds, and one for the symmetric case, for example the symmetric training kernel in KRR.

.. code:: python

    from qml.fchl import get_local_kernels
    from qml.fchl import get_local_symmetric_kernels

First function argument is the set of query compounds and the second argument is the set of compounds uses as basis (most often the training set).
In the symmetric case only one argument should be given.

Additionally, the kernels take a number of keyword arguments - these are described in a later section. The default values are carefully selected and should work well for most cases, but could be optimized for each specific dataset.

An example KRR program looks like:

.. code:: python

    # Generate representations
    reps = np.array([generate_representation(mol.coordinates, mol.nuclear_charges) for mol in mols])

    # Energies for each molecule
    energies = np.array([mol.energy for mol in mols])

    # Divide in training and test representations
    X  = reps[:100]
    Xs = reps[100:]

    # Divide in training and test energies
    U  = energies[:100]
    Us = energies[100:]

    # Generate training and test kernel
    K_training = get_local_symmetric_kernels(X)[0]
    K_test = get_local_kernels(Xs, X)[0]

    # Solve the regression using lambda=1e-7
    alphas = cho_solve(K_training, U, l2reg=1e-7)

    # Make predictions using the test kernel
    U_test = np.dot(K_test, alphas)


Note that since it is possible to get a number of kernels for different hyperparameters for free, the resulting kernel is a numpy array with three axis, where the first is the kernel index corresponding to the hyperparameters.
In the above example only the default hyperparameters are given, and thus only the 0'th kernel is used.

Global Kernels
^^^^^^^^^^^^^^

The in contrast to the local kernel, the global kernel uses a different summation so it is always a number between 1 and 0.
Besides the naming, these work similarly to the local kernels:

.. code:: python

    from qml.fchl import get_global_kernels
    from qml.fchl import get_global_symmetric_kernels

    # Generate training and test kernel
    K_training = get_global_symmetric_kernels(X)[0]
    K_test = get_global_kernels(Xs, X)[0]


Atomic Kernels
^^^^^^^^^^^^^^

In QML, atomic kernels are pairwise kernels between atoms or atomic environments in a chemical compound.


For example, to compare the atomic environments of the first atom in a set of molecules we could do the following:

.. code:: python

    from qml.fchl import get_atomic_kernels
    from qml.fchl import get_atomic_symmetric_kernels

    # Generate some molecular representations
    mol_reps = np.array([generate_representation(mol.coordinates, mol.nuclear_charges) for mol in mols])

    # Extract the representations of the first atom in each representation, and divide in training/test
    X  = np.array([rep[0] for rep in mol_reps[:100])
    Xs = np.array([rep[0] for rep in mol_reps[100:])

    # Generate training and test kernel
    K_training = get_atomic_symmetric_kernels(X)[0]
    K_test = get_atomic_kernels(Xs, X)[0]


The dimensions of the representations parsed to the atomic kernel functions should be (n_atoms,5,max_neighbors).


Common Kernel Arguments
~~~~~~~~~~~~~~~~~~~~~~~

A number of keyword arguments are shared amongst all the FCHL kernel functions. These control how two representations are compared.


.. code:: python

    K = get_local_kernels(Xs, X,

        # Weight of the two-body term, relative to the one-body term
        two_body_scaling=np.sqrt(8),

        # Weight of the three-body term, relative to the one-body term
        three_body_scaling=1.6,

        # With of the Gaussians used to compare two-body terms
        two_body_width=0.2,

        # With of the Gaussians used to compare three-body terms
        three_body_width=np.pi,

        # 1/R^n decay of two-body terms
        two_body_power=4.0,

        # 1/R^(3*n) decay of three-body terms
        three_body_power=2.0,

        # Cut-off distance in Angstrom
        cut_distance=5.0,

        # Fraction of the cut-off distance at which a switching function is turned on
        cut_start=1.0,

        # Truncation order of the Fourier expansion
        fourier_order=1,

        # Type of alchemical similarity
        alchemy='periodic-table',

        # Gaussian width for 'periodic-table" alchemy
        alchemy_period_width=1.6,
        alchemy_group_width=1.6,

        # See below
        kernel="gaussian",
        kernel_args=None,
    )


The default values are optimized to work well for molecules as well as periodic compounds, and rarely need to be changed.

However, changing to ``alchemy='off'`` can be beneficial in many situation, especially if the dataset does not have diverse stochiometries, for example different conformations of the same molecule.
``'off'`` may also reduce the computational time substantially, depending on the chemical composition of the dataset.

IMPORTANT: Remember to always use the same cut-off to generate the representation and call the kernel function!

Kernel Functions
~~~~~~~~~~~~~~~~

There are two keywords that control the kernel function (e.g. Gaussian, etc.). Read more about kernels here:
http://crsouza.com/2010/03/17/kernel-functions-for-machine-learning-applications/#kernel_functions

The first keyword ``kernel="gussian"`` controls which kernel function is used. The ``kernel_args`` keyword is a dictionary of hyperparameters for the kernel. The implemented kernel functions are:


Gaussian Kernel
^^^^^^^^^^^^^^^

This kernel contains the hyperparameters :math:`\sigma`.

    :math:`k(x,y) = \exp \left( -\frac{\|x-y\|_2^2}{2\sigma^2}\right)`

Example to explicitly call this kernel:

.. code:: python

    K = get_local_kernels(Xs, X,
        kernel="gaussian",
        kernel_args={
            "sigma": [2.5],
            }
    )


Note that this is the default kernel, and 2.5 is the default sigma.
For local kernels, 2.5 is close to optimal for most use cases.

To get the kernels for :math:`\sigma \in \{1, 2.5, 10\}`, call as:

.. code:: python

    K = get_local_kernels(Xs, X,
        kernel="gaussian",
        kernel_args={
            "sigma": [1.0, 2.5, 10.0],
            }
    )


Linear Kernel
^^^^^^^^^^^^^

This kernel contains the hyperparameters :math:`c`.

    :math:`k(x,y) = x^{T}y + c`

Example to explicitly call this kernel:

.. code:: python

    K = get_local_kernels(Xs, X,
        kernel="linear",
        kernel_args={
            "c": [0.0],
            }
    )

Ploynomial Kernel
^^^^^^^^^^^^^^^^^

This kernel contains the hyperparameters :math:`\alpha`, :math:`c`, and :math:`d`.

    :math:`k(x,y) = \left(\alpha x^{T}y + c\right)^d`

Example to explicitly call this kernel:

.. code:: python

    K = get_local_kernels(Xs, X,
        kernel="polynomial",
        kernel_args={
            "alpha": [1.0],
            "c": [0.0],
            "d": [1.0]
            }
    )


Sigmoid Kernel
^^^^^^^^^^^^^^

This kernel contains the hyperparameters :math:`\alpha` and :math:`c`.

    :math:`k(x,y) = \tanh\left(\alpha x^{T}y + c\right)`

Example to explicitly call this kernel:

.. code:: python

    K = get_local_kernels(Xs, X,
        kernel="sigmoid",
        kernel_args={
            "alpha": [1.0],
            "c": [0.0],
            }
    )


Multiquadratic Kernel
^^^^^^^^^^^^^^^^^^^^^

This kernel contains the hyperparameter :math:`c`.

    :math:`k(x,y) = \sqrt{\|x - y\|^2_2 + c^2}`

Example to explicitly call this kernel:

.. code:: python

    K = get_local_kernels(Xs, X,
        kernel="multiquadratic",
        kernel_args={
            "c": [0.0],
            }
    )


Inverse Multiquadratic Kernel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This kernel contains the hyperparameter :math:`c`.

    :math:`k(x,y) = \frac{1}{\sqrt{\|x - y\|^2_2 + c^2}}`

Example to explicitly call this kernel:

.. code:: python

    K = get_local_kernels(Xs, X,
        kernel="inverse-multiquadratic",
        kernel_args={
            "c": [0.0],
            }
    )


Bessel Kernel
^^^^^^^^^^^^^

This kernel contains the hyperparameters :math:`\sigma`, :math:`v`, and :math:`n`.

    :math:`k(x,y) = \frac{J_{v+1}\left(\sigma \|x - y\| \right)}{\|x - y\|^{-n\left(v+1\right)}}`

where :math:`J_{v+1}` is the Bessel function of the first kind.

Example to explicitly call this kernel:

.. code:: python

    K = get_local_kernels(Xs, X,
        kernel="bessel",
        kernel_args={
            "sigma": [1.0],
            "v": [1.0],
            "n": [1.0]
            }
    )

L2 Kernel
^^^^^^^^^

This kernel contains the hyperparameter :math:`c`.

    :math:`k(x,y) = \sqrt{\|x - y\| + c^2}`

This is simply the L2-distance between two representations if :math:`c=0`.

Example to explicitly call this kernel:

.. code:: python

    K = get_local_kernels(Xs, X,
        kernel="l2",
        kernel_args={
            "c": [0.0],
            }
    )


Matern Kernel
^^^^^^^^^^^^^

This kernel contains the hyperparameter :math:`\sigma` and :math:`n`.

    :math:`k(x,y)=C_n(d) = \sigma^2\frac{2^{1-n}}{\Gamma(n)}\Bigg(\sqrt{2n}\frac{\|x - y\|}{\rho}\Bigg)^n K_n\Bigg(\sqrt{2n}\frac{\|x - y\|}{\rho}\Bigg)`,

    with  :math:`rho` depending on n and :math:`\sigma`.


Example to explicitly call this kernel:

.. code:: python

    K = get_local_kernels(Xs, X,
        kernel="l2",
        kernel_args={
            "sigma": [10.0],
            "n": [1.0],
            }
    )


Cauchy Kernel
^^^^^^^^^^^^^

This kernel contains the hyperparameter :math:`\sigma`.

    :math:`k(x,y) = \frac{1}{1+\frac{\|x - y\|^2_2}{\sigma^2}}`

Example to explicitly call this kernel:

.. code:: python

    K = get_local_kernels(Xs, X,
        kernel="l2",
        kernel_args={
            "sigma": [2.0],
            }
    )


Force Kernels
~~~~~~~~~~~~~

aasf

Electric Field-Dependent Kernels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

asfasf

Dipole-Moment Kernels
~~~~~~~~~~~~~~~~~~~~~

sdgsdg
