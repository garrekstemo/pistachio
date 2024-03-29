# Test Structure class

import importlib.resources as pkg_resources
import numpy as np
import os
import pytest
from numpy import testing
from pistachio import transfer_matrix as tm
from pistachio import default


def make_structure(num_layers, lmbda, thickness):
    structure = tm.Structure()
    for i in range(num_layers):
        material = "Material_" + str(i)
        test_layer = tm.Layer(material=material, thickness=thickness)
        test_layer.wavelengths = [lmbda]
        structure.add_layer(test_layer)
    return structure


def test_add_layer():
    lmbda = 1.e-5
    thickness = 3 * lmbda
    new_layer = tm.Layer(material="New Material", thickness=1.e-6)
    new_layer.wavelengths = lmbda
    s = make_structure(3, lmbda, thickness)
    s.add_layer(new_layer)
    assert s.layers[-1].material == "New Material"

def test_delete_layer():
    s = make_structure(3, 1.e-6, 1.e-3)
    s.delete_layer(1)
    assert s.layers[1].material == "Material_2"

# Test yaml structure has three layers
#   Air substrate
#   Thin Au layer
#   Air superstrate

yaml_file = 'Au_thin_layer.yaml'
with pkg_resources.path(default, yaml_file) as yml:
    yaml_config = os.path.abspath(yml)
s_from_yaml = tm.Structure()
s_from_yaml.load_struct_from_config(yaml_config)


def test_load_struct_from_yaml():
    npoints = 100
    assert len(s_from_yaml.layers) == 3
    assert len(s_from_yaml.wavelengths) == npoints
    assert len(s_from_yaml.layers[0].wavelengths) == npoints
    assert len(s_from_yaml.layers[1].wavelengths) == npoints
    assert len(s_from_yaml.layers[2].wavelengths) == npoints


def test_radians_after_initialize_from_yaml():
    # Test whether angles provided in yaml config file (degrees)
    # are converted to radians after initialization.
    # Incidence angle > pi/2 is not expected.
    test_theta = np.full(len(s_from_yaml.theta), np.pi/2)
    testing.assert_array_less(s_from_yaml.theta, test_theta)


def test_calculate_t_r():
    M = np.array([[1.0 + 0.j, 0.0 + 0.j],
              [0.0 + 0.j,  1.0 + 0.j]])

    T, R = tm.Structure().calculate_t_r(M)
    assert T == 1.0
    assert R == 0.0



