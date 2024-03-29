# Test Layer class

import numpy as np
import pytest
from numpy import testing
from pistachio import transfer_matrix as tm


lmbda = 1.0  # μm
thickness = 3 * lmbda
num_layers = 100

test_layer = tm.Layer(material="Air", thickness=thickness, num_points=1)


def test_make_datapoints():
    layer = tm.Layer(material="Air")
    layer.refractive_index = [1.0]
    layer.extinction_coeff = [0.0]
    num_new_wavelengths = 3
    lmbda = np.linspace(1.0, 10.0, num_new_wavelengths)
    layer.make_datapoints(lmbda)
    assert len(layer.wavelengths) == num_new_wavelengths


def test_make_datapoints_interpolation():
    # Test how much interpolation alters precision compared to original data.
    layer = tm.Layer()
    layer.get_data_from_csv("Au.csv")
    old_wavelengths = layer.wavelengths
    old_refractive = layer.refractive_index
    old_extinction = layer.extinction_coeff

    new_wavelengths = np.linspace(0.2, 15.0, 10000)
    print(new_wavelengths)
    layer.make_datapoints(new_wavelengths)
    new_refractive = layer.refractive_index
    new_extinction = layer.extinction_coeff
    
    # Get the array indices for a particular wavelength
    r_old1 = np.where(np.isclose(old_wavelengths, 0.2))[0][0]
    r_old2 = np.where(np.isclose(old_wavelengths, 1.0))[0][0]
    r_old3 = np.where(np.isclose(old_wavelengths, 2.0))[0][0]
    r_old4 = np.where(np.isclose(old_wavelengths, 12.711))[0][0]

    r_new1 = np.where(np.isclose(new_wavelengths, 0.2))[0][0]
    r_new2 = np.where(np.isclose(new_wavelengths, 1.0, atol=1e-2))[0][0]
    r_new3 = np.where(np.isclose(new_wavelengths, 2.0, atol=1e-2))[0][0]
    r_new4 = np.where(np.isclose(new_wavelengths, 12.7, atol=1e-2))[0][0]

    
    testing.assert_almost_equal(old_wavelengths[0], 0.19077)
    testing.assert_almost_equal(old_wavelengths[264], 1.0164)

    testing.assert_almost_equal(new_refractive[r_new1], old_refractive[r_old1], decimal=3)
    testing.assert_almost_equal(new_refractive[r_new2], old_refractive[r_old2], decimal=2)
    testing.assert_almost_equal(new_refractive[r_new3], old_refractive[r_old3], decimal=2)
    testing.assert_almost_equal(new_refractive[r_new4], old_refractive[r_old4], decimal=1)

    testing.assert_almost_equal(new_extinction[r_new1], old_extinction[r_old1], decimal=2)
    testing.assert_almost_equal(new_extinction[r_new2], old_extinction[r_old2], decimal=1)
    testing.assert_almost_equal(new_extinction[r_new3], old_extinction[r_old3], decimal=1)
    testing.assert_almost_equal(new_extinction[r_new4], old_extinction[r_old4], decimal=1)



def test_dynamical_matrix():
    A = tm.Layer().dynamical_matrix(1. + 1.j, 0., "s-wave")
    B = tm.Layer().dynamical_matrix(1. + 1.j, 0., "p-wave")
    C = np.array([[1., 1.,], [1. + 1.j, -1. - 1.j]])

    testing.assert_allclose(A, B)
    testing.assert_allclose(A, C)
    testing.assert_allclose(B, C)

    #TODO: Test empty array

def test_propagation_matrix():
    kx = 2.e5 + 0.j
    d = np.pi / (2 * kx)
    P_actual = tm.Layer().propagation_matrix(kx, d)
    P_desired = np.array([[-1.j, 0.], [0., 1.j]])

    testing.assert_allclose(P_actual, P_desired)


def test_transfer_matrix():
    layer = test_layer
    M_actual = layer.calculate_transfer_matrix(1. + 0.j, lmbda, 0., "s-wave")[-1]

    testing.assert_allclose(M_actual, np.identity(2), atol=1e-4)
