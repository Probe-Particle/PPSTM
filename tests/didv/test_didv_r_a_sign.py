import math
from functools import partial
from typing import Callable
from abc import ABC, abstractmethod
from typing import Callable, Tuple

import numpy as np

from pyPPSTM import ProbeSTM
from pyPPSTM.probe_stm_numpy import ProbeStmNumpy
from pyPPSTM.probe_stm_opencl import ProbeSTMOpenCLSequential, ProbeSTMOpenCLParallel
from pyPPSTM.probe_stm_pytorch import ProbeStmPytorch


class _TestDidvRaSign(ABC):
    _TIP_ORBITAL: np.ndarray
    _BACKENDS_NAME_FN_RTOL: Tuple[Tuple[str, Callable, float]]
    _RTOL_REV_R_A_SIGN: float

    _AB = 1.889725989
    _EV = 0.036749034
    _N_P = math.sqrt(3)

    _WORK_FUNCTION = 1/(2*_EV)  # eV, decay = 1

    # sample DOS = 1
    _V = 0.0  # voltage
    _EIGENENERGIES = np.asarray([0.0])
    _ETA = 2 / math.pi  # eV, Lorentzian width of states in energy scale

    # ||r_a|| = ||tip_pos - sample_atom_pos|| = 1
    _SAMPLE_ATOMS_POSITIONS = np.asarray([[0.0, 0.0, 0.0]])
    _TIP_POSITIONS = np.asarray([[[[1.0, 0.5, 2.0]]]])
    _TIP_POSITIONS /= np.linalg.norm(_TIP_POSITIONS, ord=2, axis=-1)*_AB
    _R_A = (_TIP_POSITIONS.squeeze() - _SAMPLE_ATOMS_POSITIONS.squeeze())*_AB
    _R_A_NORM = np.linalg.norm(_R_A, ord=2, axis=-1)
    assert round(_R_A_NORM, 8) == 1.0

    _SAMPLE_ORBITAL_COUNT = 4  # s and p sample orbitals

    #  selection of LCAO coefficients that help detecting issues in sign of r_a
    #  by including the s sample orbital plus at least one p sample orbital
    _SAMPLE_ORBITALS_LCAO_COEFFICIENTS = (
        ("s,py",       np.array([[1., 1., 0., 0.]])),
        ("s,pz",       np.array([[1., 0., 1., 0.]])),
        ("s,px",       np.array([[1., 0., 0., 1.]])),
        ("s,py,pz",    np.array([[1., 1., 1., 0.]])),
        ("s,py,px",    np.array([[1., 1., 0., 1.]])),
        ("s,pz,px",    np.array([[1., 0., 1., 1.]])),
        ("s,py,pz,px", np.array([[1., 1., 1., 1.]])),
    )

    def test_didv_r_a_sign(self, lcao_coefficients: np.ndarray, didv_backend: Callable, rtol: float):
        """Verify the correctness of the sign of r_a.

        Args:
            lcao_coefficients (np.ndarray): Orbital coefficients for sample, shape (1, 4)
            didv_backend (Callable): dI/dV implementation
            rtol (float): Relative tolerance
        """
        actual_didv = self._actual_didv(lcao_coefficients=lcao_coefficients, backend=didv_backend)
        expected_didv = self._expected_didv(lcao_coefficients)
        np.testing.assert_allclose(actual_didv, expected_didv, rtol=rtol, atol=0.)

    def test_reversing_r_a_sign_changes_didv(self, lcao_coefficients: np.ndarray):
        """Verify that dI/dV differs when r_a sign is reversed.

        Args:
            lcao_coefficients (np.ndarray): Orbital coefficients for sample, shape (1, 4)
        """
        expected_didv = self._expected_didv(lcao_coefficients)
        reverse_r_a_sign_didv = self._expected_didv(lcao_coefficients, reverse_r_a_sign=True)
        assert not np.allclose(reverse_r_a_sign_didv, expected_didv, rtol=self._RTOL_REV_R_A_SIGN, atol=0.)

    def _expected_didv(self, lcao_coefficients: np.ndarray, reverse_r_a_sign: bool = False):
        """Compute dI/dV (differential conductance).

        Args:
            lcao_coefficients (np.ndarray): Orbital coefficients for sample, shape (1, 4)
            reverse_r_a_sign (bool): True when the correct sign of r_a should be reversed,
            False when the correct sign of r_a should be kept. Default is False.
        """
        if not reverse_r_a_sign:
            r_a = self._R_A
        else:
            r_a = -self._R_A
        radial = np.exp(-self._R_A_NORM)
        g = 16 * math.pi ** 3 * (self._g_tip_orb(lcao_coefficients, r_a) * radial) ** 2
        return g.reshape(1, 1, 1)

    @abstractmethod
    def _g_tip_orb(self, lcao_coefficients: np.ndarray, r_a: np.ndarray):
        """Orbital conductances for different sp orbitals of sample on an orbital tip.

        Args:
            lcao_coefficients (np.ndarray): Orbital coefficients for sample, shape (1, 4)
            r_a (np.ndarray): shape (1, 1, 1, 3)

        Returns:
            np.ndarray, shape (1, 1, 1)
        """
        pass

    def _actual_didv(self, lcao_coefficients: np.ndarray, backend: Callable) -> np.ndarray:
        """Compute differential conductance using the specified backend.

        Args:
            lcao_coefficients (np.ndarray): Orbital coefficients for sample, shape (1, 4)
            backend (Callable): dI/dV implementation

        Returns:
            Array of conductance values, shape: shape (n_z, n_y, n_x)
        """
        return backend(V=self._V,
                       WF=self._WORK_FUNCTION,
                       eta=self._ETA,
                       eig=self._EIGENENERGIES,
                       R=self._TIP_POSITIONS,
                       Rat=self._SAMPLE_ATOMS_POSITIONS,
                       coes=lcao_coefficients,
                       tip_coes=self._TIP_ORBITAL,
                       orb_t=self._SAMPLE_ORBITAL_COUNT)

    def pytest_generate_tests(self, metafunc):
        if "didv_backend" in metafunc.fixturenames and "rtol" in metafunc.fixturenames:
            metafunc.parametrize(("didv_backend", "rtol"),
                                 ((fn, rtol) for _, fn, rtol in self._BACKENDS_NAME_FN_RTOL),
                                 ids=(name for name, _, _ in self._BACKENDS_NAME_FN_RTOL))

        if "lcao_coefficients" in metafunc.fixturenames:
            metafunc.parametrize("lcao_coefficients",
                                 (lcao_coeffs for _, lcao_coeffs in self._SAMPLE_ORBITALS_LCAO_COEFFICIENTS),
                                 ids=(f"saorb={sample_orbitals}"
                                      for sample_orbitals, _ in self._SAMPLE_ORBITALS_LCAO_COEFFICIENTS))

class TestDidvRaSignSTipOrbital(_TestDidvRaSign):
    _TIP_ORBITAL = np.asarray([1.0] + [0.0] * 8)  # s tip orbital

    _BACKENDS_NAME_FN_RTOL = (
        ("C++",                        ProbeSTM.dIdV_sp_sp,                                     3e-8),
        ("OpenCL parallel",            ProbeSTMOpenCLParallel.didv,                             2e-7),
        ("OpenCL sequential",          ProbeSTMOpenCLSequential.didv,                           2e-7),
        ("NumPy (default chunking)",   ProbeStmNumpy.didv,                                      2e-7),
        ("PyTorch (default chunking)", ProbeStmPytorch.didv,                                    2e-7),
        ("NumPy (no chunking)",        partial(ProbeStmNumpy.didv,    n_tip_position_chunks=1), 2e-7),
        ("PyTorch (no chunking)",      partial(ProbeStmPytorch.didv,  n_tip_position_chunks=1), 2e-7),
        ("NumPy (2 chunks)",           partial(ProbeStmNumpy.didv,    n_tip_position_chunks=2), 2e-7),
        ("PyTorch (2 chunks)",         partial(ProbeStmPytorch.didv,  n_tip_position_chunks=2), 2e-7),
    )
    _RTOL_REV_R_A_SIGN = 7e-1

    def _g_tip_orb(self, lcao_coefficients: np.ndarray, r_a: np.ndarray):
        """Orbital conductances for different sp orbitals of sample on an s orbital tip.

        Args:
            lcao_coefficients (np.ndarray): Orbital coefficients for sample, shape (1, 4)
            r_a (np.ndarray): shape (1, 1, 1, 3)

        Returns:
            np.ndarray, shape (1, 1, 1)
        """
        t = lcao_coefficients[..., 0] / self._N_P              # s  orb. of sample
        t = t + lcao_coefficients[..., 1] * r_a[..., 1]  # py orb. of sample
        t = t + lcao_coefficients[..., 2] * r_a[..., 2]  # pz orb. of sample
        t = t + lcao_coefficients[..., 3] * r_a[..., 0]  # px orb. of sample
        return t * self._N_P


class TestDidvRaSignDz2TipOrbital(_TestDidvRaSign):
    _I_3 = 1 / 3

    _TIP_ORBITAL = np.asarray([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0])  # dz2 tip orbital

    _BACKENDS_NAME_FN_RTOL = (
        ("C++",                        ProbeSTM.dIdV_sp_sp,                                     2e-6),
        ("OpenCL parallel",            ProbeSTMOpenCLParallel.didv,                             4e-6),
        ("OpenCL sequential",          ProbeSTMOpenCLSequential.didv,                           4e-6),
        ("NumPy (default chunking)",   ProbeStmNumpy.didv,                                      4e-6),
        ("PyTorch (default chunking)", ProbeStmPytorch.didv,                                    4e-6),
        ("NumPy (no chunking)",        partial(ProbeStmNumpy.didv,    n_tip_position_chunks=1), 4e-6),
        ("PyTorch (no chunking)",      partial(ProbeStmPytorch.didv,  n_tip_position_chunks=1), 4e-6),
        ("NumPy (2 chunks)",           partial(ProbeStmNumpy.didv,    n_tip_position_chunks=2), 4e-6),
        ("PyTorch (2 chunks)",         partial(ProbeStmPytorch.didv,  n_tip_position_chunks=2), 4e-6),
    )
    _RTOL_REV_R_A_SIGN = 2e-1

    def _g_tip_orb(self, lcao_coefficients: np.ndarray, r_a: np.ndarray):
        """Orbital conductances for different sp orbitals of sample on a dz2 orbital tip.

        Args:
            lcao_coefficients (np.ndarray): Orbital coefficients for sample, shape (1, 4)
            r_a (np.ndarray): shape (1, 1, 1, 3)

        Returns:
            np.ndarray, shape (1, 1, 1)
        """
        t = lcao_coefficients[..., 0] * (2 * r_a[..., 2] ** 2 - 1. - self._I_3) / self._N_P              # s  orb. of sample
        t = t + lcao_coefficients[..., 1] * r_a[..., 1] * (7 * r_a[..., 2] ** 2 - 2. - self._I_3)  # py orb. of sample
        t = t + lcao_coefficients[..., 2] * r_a[..., 2] * (7 * r_a[..., 2] ** 2 - 6. - self._I_3)  # pz orb. of sample
        t = t + lcao_coefficients[..., 3] * r_a[..., 0] * (7 * r_a[..., 2] ** 2 - 2. - self._I_3)  # px orb. of sample
        return t * self._N_P
