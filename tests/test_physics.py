import importlib.resources as pkg_resources
import numpy as np
import os
import pytest
from pistachio import transfer_matrix as tm
from pistachio import default

yaml_config = None
with pkg_resources.path(default, 'fabry-perot.yaml') as yml:
        yaml_config = os.path.abspath(yml)

test_fp = tm.Structure()
test_fp.load_struct_from_config(yaml_config)

lmbda = 1.e-5
thickness = 3 * lmbda

num_layers = 100
layers = []

test_layer = tm.Layer(material='glass', thickness=thickness, wavelengths=[lmbda],
                      n_real=[1.0], n_imag=[0.j], num_points=1)

for i in range(num_layers):
    layers.append(test_layer)

test_structure = tm.Structure(layers=layers, wavelengths=[lmbda], theta=[0.])


def test_t_r_conservation():
    # No absorption
    M = test_structure.calculate_transfer_matrix(1. + 0.j, lmbda, 0., 's-wave')
    T, R = test_structure.calculate_t_r(M)

    M_fp = test_fp.calculate_transfer_matrix(1. + 0.j, lmbda, 0., 's-wave')
    T_fp, R_fp = test_fp.calculate_t_r(M_fp)

    assert (1 - T - R == 0.)
    assert (1 - T_fp - R_fp == 0.)
