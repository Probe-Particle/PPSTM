import math
from functools import partial
from abc import ABC, abstractmethod
from typing import Callable, Tuple

import numpy as np

from pyPPSTM import ProbeSTM
from pyPPSTM.probe_stm_numpy import ProbeStmNumpy
from pyPPSTM.probe_stm_opencl import ProbeSTMOpenCLSequential, ProbeSTMOpenCLParallel
from pyPPSTM.probe_stm_pytorch import ProbeStmPytorch


class _TestDidvSpSampleOrbital(ABC):
    _TIP_ORBITAL: np.ndarray
    _BACKENDS_NAME_FN_RTOL: Tuple[Tuple[str, Callable, float]]
    _RTOL_WRONG_SAMPLE_ORBITALS: float

    _AB = 1.889725989
    _EV = 0.036749034
    _N_P = math.sqrt(3)
    _N_D = math.sqrt(15)
    _N_D2 = math.sqrt(5) * 0.5

    _WORK_FUNCTION = 1/(2*_EV)  # eV, decay = 1

    # sample DOS = 1
    _V = 0.0  # voltage
    _EIGENENERGIES = np.asarray([0.0])
    _ETA = 2 / math.pi  # eV, Lorentzian width of states in energy scale

    # ||r_a|| = ||tip_pos - sample_atom_pos|| = 1
    _SAMPLE_ATOMS_POSITIONS = np.asarray([[0.0, 0.0, 0.0]]*3)
    _TIP_POSITIONS = np.asarray([[[[1.0, 0.5, 2.0]]]])
    _TIP_POSITIONS /= np.linalg.norm(_TIP_POSITIONS, ord=2, axis=-1)*_AB
    _N_TIP_POSITIONS = 1
    _R_A = (_TIP_POSITIONS.squeeze() - _SAMPLE_ATOMS_POSITIONS)*_AB
    _R_A_NORM = np.linalg.norm(_R_A, ord=2, axis=-1)
    np.testing.assert_allclose(_R_A_NORM, np.ones_like(_R_A_NORM), rtol=1e-8)
    _RADIAL = np.exp(-_R_A_NORM)

    _SAMPLE_ORBITAL_COUNT = 4  # s and p sample orbitals

    #  selection of LCAO coefficients that can help detecting
    #  the insertion of the computation of g specific for d sample orbitals
    #  in paths specific for s/p only sample orbitals
    _SAMPLE_ORBITALS_LCAO_COEFFICIENTS = (
        ("dxy",       "dxy",   np.array([[0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.]])),  # target: g for dxy sample orbital
        ("dyz",       "dyz",   np.array([[0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.]])),  # target: g for dyz sample orbital
        ("dz2",       "dz2",   np.array([[0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.]])),  # target: g for dz2 sample orbital
        ("dxz",       "dxz",   np.array([[0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.]])),  # target: g for dxz sample orbital
        ("dxy|dx2y2", "dx2y2", np.array([[0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.]])),  # target: g for dxy or dx2y2 sample orbitals
    )

    def test_correct_didv_sample_orbital(self, lcao_coefficients: np.ndarray, didv_backend: Callable, rtol: float):
        """Verify the correctness of dI/dV for a single sample orbital.

        Args:
            lcao_coefficients (np.ndarray): Orbital coefficients for sample, shape (1, 4)
            didv_backend (Callable): dI/dV implementation
            rtol (float): Relative tolerance
        """
        actual_didv = self._actual_didv(lcao_coefficients=lcao_coefficients, backend=didv_backend)
        expected_didv = self._expected_didv(lcao_coefficients)
        np.testing.assert_allclose(actual_didv, expected_didv, rtol=rtol, atol=0.)

    def test_wrong_didv_sample_orbital(self,
                                       lcao_coefficients: np.ndarray,
                                       didv_backend: Callable,
                                       sample_orbital_insertion: str):
        """Verify that dI/dV is far from incorrect ones.

        Args:
            lcao_coefficients (np.ndarray): Orbital coefficients for sample, shape (1, 4)
            didv_backend (Callable): dI/dV implementation
            sample_orbital_insertion (str): sample orbital for which to include g
        """
        actual_didv = self._actual_didv(lcao_coefficients=lcao_coefficients, backend=didv_backend)
        expected_didv = self._expected_didv(lcao_coefficients, sample_orbital_insertion)
        assert not np.allclose(actual_didv, expected_didv, rtol=self._RTOL_WRONG_SAMPLE_ORBITALS, atol=0.)

    def _expected_didv(self, lcao_coefficients: np.ndarray, sample_orbital_insertion: str|None = None):
        """Compute dI/dV (differential conductance).

        Args:
            lcao_coefficients (np.ndarray): Orbital coefficients for sample, shape (1, 4)
        """
        g = 0
        for a in range(len(self._SAMPLE_ATOMS_POSITIONS)):
            g += self._RADIAL[a] * self._g_tip_orb(lcao_coefficients[0, a*self._SAMPLE_ORBITAL_COUNT:],
                                                   self._R_A[a],
                                                   sample_orbital_insertion)
        g = 16 * math.pi ** 3 * g ** 2
        return np.asarray(g).reshape(*(self._N_TIP_POSITIONS,)*3)

    @abstractmethod
    def _g_tip_orb(self, lcao_coefficients: np.ndarray, r_a, sample_orbital_insertion: str|None = None):
        """Orbital conductances for different sp orbitals of sample on an orbital tip.

        Args:
            lcao_coefficients (np.ndarray): Orbital coefficients for sample, shape (1, 4)

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
        if "didv_backend" in metafunc.fixturenames:
            if "sample_orbital_insertion" in metafunc.fixturenames:
                metafunc.parametrize("didv_backend",
                                     (fn for _, fn, _ in self._BACKENDS_NAME_FN_RTOL),
                                     ids=(name for name, _, _ in self._BACKENDS_NAME_FN_RTOL))
            elif "rtol" in metafunc.fixturenames:
                metafunc.parametrize(("didv_backend", "rtol"),
                                     ((fn, rtol) for _, fn, rtol in self._BACKENDS_NAME_FN_RTOL),
                                     ids=(name for name, _, _ in self._BACKENDS_NAME_FN_RTOL))

        if "lcao_coefficients" in metafunc.fixturenames:
            if "sample_orbital_insertion" in metafunc.fixturenames:
                metafunc.parametrize(("lcao_coefficients", "sample_orbital_insertion"),
                                     ((lcao_coeffs, sample_orbital_insertion)
                                      for _, sample_orbital_insertion, lcao_coeffs
                                      in self._SAMPLE_ORBITALS_LCAO_COEFFICIENTS),
                                     ids=(f"saorb={sample_orbitals}-i:{sample_orbital_insertion}"
                                          for sample_orbitals, sample_orbital_insertion, _
                                          in self._SAMPLE_ORBITALS_LCAO_COEFFICIENTS))
            else:
                metafunc.parametrize("lcao_coefficients",
                                     (lcao_coeffs for _, _, lcao_coeffs in self._SAMPLE_ORBITALS_LCAO_COEFFICIENTS),
                                     ids=(f"saorb={sample_orbitals}"
                                          for sample_orbitals, _, _ in self._SAMPLE_ORBITALS_LCAO_COEFFICIENTS))

class TestDidvSpSampleOrbitalSTipOrbital(_TestDidvSpSampleOrbital):
    _TIP_ORBITAL = np.asarray([1.0] + [0.0] * 8)

    _BACKENDS_NAME_FN_RTOL = (
        ("C++",                        ProbeSTM.dIdV_sp_sp,                                             4e-8),
        ("OpenCL parallel",            ProbeSTMOpenCLParallel.didv,                                     2e-7),
        ("OpenCL sequential",          ProbeSTMOpenCLSequential.didv,                                   2e-7),
        ("NumPy (default chunking)",   ProbeStmNumpy.didv,                                              3e-7),
        ("PyTorch (default chunking)", ProbeStmPytorch.didv,                                            2e-7),
        ("NumPy (no chunking)",        partial(ProbeStmNumpy.didv,   n_tip_position_chunks=1), 3e-7),
        ("PyTorch (no chunking)",      partial(ProbeStmPytorch.didv, n_tip_position_chunks=1), 2e-7),
        ("NumPy (2 chunks)",           partial(ProbeStmNumpy.didv,   n_tip_position_chunks=2), 3e-7),
        ("PyTorch (2 chunks)",         partial(ProbeStmPytorch.didv, n_tip_position_chunks=2), 2e-7),
    )

    _RTOL_WRONG_SAMPLE_ORBITALS = 3e-1

    def _g_tip_orb(self, lcao_coefficients: np.ndarray, r_a: np.ndarray, sample_orbital_insertion: str|None = None):
        """Orbital conductances for different sp orbitals of sample on an s orbital tip.

        Args:
            lcao_coefficients (np.ndarray): Orbital coefficients for sample, shape (1, 4)

        Returns:
            np.ndarray, shape (1, 1, 1)
        """
        t = lcao_coefficients[0]                           # s  orb. of sample
        t = t + lcao_coefficients[1] * self._N_P * r_a[1]  # py orb. of sample
        t = t + lcao_coefficients[2] * self._N_P * r_a[2]  # pz orb. of sample
        t = t + lcao_coefficients[3] * self._N_P * r_a[0]  # px orb. of sample
        if len(lcao_coefficients) > 4 and sample_orbital_insertion == "dxy":
            t = t + lcao_coefficients[4] * self._N_D * r_a[0] * r_a[1]
        if len(lcao_coefficients) > 5 and sample_orbital_insertion == "dyz":
            t = t + lcao_coefficients[5] * self._N_D * r_a[1] * r_a[2]
        if len(lcao_coefficients) > 6 and sample_orbital_insertion == "dz2":
            t = t + lcao_coefficients[6] * self._N_D2 * (3 * r_a[2] ** 2 - 1)
        if len(lcao_coefficients) > 7 and sample_orbital_insertion == "dxz":
            t = t + lcao_coefficients[7] * self._N_D * r_a[0] * r_a[2]
        if len(lcao_coefficients) > 8 and sample_orbital_insertion == "dx2y2":
            t = t + lcao_coefficients[8] * self._N_D * .5 * (r_a[0] ** 2 - r_a[1] ** 2)
        return t

class TestDidvSpSampleOrbitalPzTipOrbital(_TestDidvSpSampleOrbital):
    _TIP_ORBITAL = np.asarray([0.0, 0.0, 1.0] + [0.0] * 6)

    _BACKENDS_NAME_FN_RTOL = (
        ("C++",                        ProbeSTM.dIdV_sp_sp,                                             4e-8),
        ("OpenCL parallel",            ProbeSTMOpenCLParallel.didv,                                     3e-7),
        ("OpenCL sequential",          ProbeSTMOpenCLSequential.didv,                                   3e-7),
        ("NumPy (default chunking)",   ProbeStmNumpy.didv,                                              1e-7),
        ("PyTorch (default chunking)", ProbeStmPytorch.didv,                                            3e-7),
        ("NumPy (no chunking)",        partial(ProbeStmNumpy.didv,   n_tip_position_chunks=1), 1e-7),
        ("PyTorch (no chunking)",      partial(ProbeStmPytorch.didv, n_tip_position_chunks=1), 3e-7),
        ("NumPy (2 chunks)",           partial(ProbeStmNumpy.didv,   n_tip_position_chunks=2), 1e-7),
        ("PyTorch (2 chunks)",         partial(ProbeStmPytorch.didv, n_tip_position_chunks=2), 3e-7),
    )

    _RTOL_WRONG_SAMPLE_ORBITALS = 8e-2

    def _g_tip_orb(self, lcao_coefficients: np.ndarray, r_a: np.ndarray, sample_orbital_insertion: str|None = None):
        """Orbital conductances for different sp orbitals of sample on an s orbital tip.

        Args:
            lcao_coefficients (np.ndarray): Orbital coefficients for sample, shape (1, 4)

        Returns:
            np.ndarray, shape (1, 1, 1)
        """
        t = lcao_coefficients[0] * r_a[2]                                  # s  orb. of sample
        t = t + lcao_coefficients[1] * self._N_P * r_a[2] * r_a[1] * 2     # py orb. of sample
        t = t + lcao_coefficients[2] * self._N_P * (-1 + 2 * r_a[2] ** 2)  # pz orb. of sample
        t = t + lcao_coefficients[3] * self._N_P * r_a[2] * r_a[0] * 2     # px orb. of sample
        if len(lcao_coefficients) > 4 and sample_orbital_insertion == "dxy":
             t = t + lcao_coefficients[4] * self._N_D * r_a[0] * r_a[1] * r_a[2] * 3
        if len(lcao_coefficients) > 5 and sample_orbital_insertion == "dyz":
             t = t + lcao_coefficients[5] * self._N_D * r_a[1] * (3 * r_a[2] ** 2 - 1)
        if len(lcao_coefficients) > 6 and sample_orbital_insertion == "dz2":
            t = t + lcao_coefficients[6] * self._N_D2 * (6 * (r_a[2] ** 3 - r_a[2]) + 3 * r_a[2] ** 2 - 1)
        if len(lcao_coefficients) > 7 and sample_orbital_insertion == "dxz":
            t = t + lcao_coefficients[7] * self._N_D * r_a[0] * (3 * r_a[2] ** 2 - 1)
        if len(lcao_coefficients) > 8 and sample_orbital_insertion == "dx2y2":
            t = t + lcao_coefficients[8] * self._N_D * r_a[2] * (r_a[0] ** 2 - r_a[1] ** 2) * 1.5
        return t

class TestDidvSpSampleOrbitalPyTipOrbital(_TestDidvSpSampleOrbital):
    _TIP_ORBITAL = np.asarray([0.0, 1.0] + [0.0] * 7)

    _BACKENDS_NAME_FN_RTOL = (
        ("C++",                        ProbeSTM.dIdV_sp_sp,                                             4e-8),
        ("OpenCL parallel",            ProbeSTMOpenCLParallel.didv,                                     2e-7),
        ("OpenCL sequential",          ProbeSTMOpenCLSequential.didv,                                   2e-7),
        ("NumPy (default chunking)",   ProbeStmNumpy.didv,                                              1e-7),
        ("PyTorch (default chunking)", ProbeStmPytorch.didv,                                            2e-7),
        ("NumPy (no chunking)",        partial(ProbeStmNumpy.didv,   n_tip_position_chunks=1), 1e-7),
        ("PyTorch (no chunking)",      partial(ProbeStmPytorch.didv, n_tip_position_chunks=1), 2e-7),
        ("NumPy (2 chunks)",           partial(ProbeStmNumpy.didv,   n_tip_position_chunks=2), 1e-7),
        ("PyTorch (2 chunks)",         partial(ProbeStmPytorch.didv, n_tip_position_chunks=2), 2e-7),
    )

    _RTOL_WRONG_SAMPLE_ORBITALS = 8e-1

    def _g_tip_orb(self, lcao_coefficients: np.ndarray, r_a: np.ndarray, sample_orbital_insertion: str|None = None):
        """Orbital conductances for different sp orbitals of sample on an s orbital tip.

        Args:
            lcao_coefficients (np.ndarray): Orbital coefficients for sample, shape (1, 4)

        Returns:
            np.ndarray, shape (1, 1, 1)
        """
        t = lcao_coefficients[0] * r_a[1]                                  # s  orb. of sample
        t = t + lcao_coefficients[1] * self._N_P * (-1 + 2 * r_a[1] ** 2)     # py orb. of sample
        t = t + lcao_coefficients[2] * self._N_P * r_a[1] * r_a[2] * 2  # pz orb. of sample
        t = t + lcao_coefficients[3] * self._N_P * r_a[1] * r_a[0] * 2     # px orb. of sample
        if len(lcao_coefficients) > 4 and sample_orbital_insertion == "dxy":
             t = t + lcao_coefficients[4] * self._N_D * r_a[0] * (3 * r_a[1] ** 2 - 1)
        if len(lcao_coefficients) > 5 and sample_orbital_insertion == "dyz":
             t = t + lcao_coefficients[5] * self._N_D * r_a[2] * (3 * r_a[1] ** 2 - 1)
        if len(lcao_coefficients) > 6 and sample_orbital_insertion == "dz2":
            t = t + lcao_coefficients[6] * self._N_D2 * r_a[1] * (6 * r_a[2] ** 2 + (3 * r_a[2] ** 2 - 1))
        if len(lcao_coefficients) > 7 and sample_orbital_insertion == "dxz":
            t = t + lcao_coefficients[7] * self._N_D * r_a[0] * r_a[1] * r_a[2] * 3
        if len(lcao_coefficients) > 8 and sample_orbital_insertion == "dx2y2":
            t = t + lcao_coefficients[8] * self._N_D * r_a[1] * ((r_a[0] ** 2 - r_a[1] ** 2) * 1.5 + 1)
        return t
