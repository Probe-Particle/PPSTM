import itertools
from functools import partial
from typing import Callable, List, Iterator, Tuple

import numpy as np

from ppstm import ProbeSTM
from ppstm.probe_stm_numpy import ProbeStmNumpy
from ppstm.probe_stm_opencl import ProbeSTMOpenCLSequential, ProbeSTMOpenCLParallel
from ppstm.probe_stm_pytorch import ProbeStmPytorch


class TestDidvBackendConsistencyPairwiseTipSampleOrbitals:
    _WORK_FUNCTION = 5.0  # eV
    _V = -2.0  # voltage
    _ETA = 1e-7  # eV, Lorentzian width of states in energy scale

    _EIGENENERGIES = np.asarray([-2.5])
    _SAMPLE_ATOMS_POSITIONS = np.asarray([[6.5, 0.001, 4.7]])
    _TIP_POSITIONS = np.asarray([[[[5.5, 0., 8.7]]]])

    _BACKENDS_NAME_FN_RTOL = (
        ("OpenCL parallel",            ProbeSTMOpenCLParallel.didv,                             2e-6),
        ("OpenCL sequential",          ProbeSTMOpenCLSequential.didv,                           2e-6),
        ("NumPy (default chunking)",   ProbeStmNumpy.didv,                                      9e-6),
        ("PyTorch (default chunking)", ProbeStmPytorch.didv,                                    9e-6),
        ("NumPy (no chunking)",        partial(ProbeStmNumpy.didv,    n_tip_position_chunks=1), 9e-6),
        ("PyTorch (no chunking)",      partial(ProbeStmPytorch.didv,  n_tip_position_chunks=1), 9e-6),
        ("NumPy (2 chunks)",           partial(ProbeStmNumpy.didv,    n_tip_position_chunks=2), 9e-6),
        ("PyTorch (2 chunks)",         partial(ProbeStmPytorch.didv,  n_tip_position_chunks=2), 9e-6),
    )

    def test_didv_backend_consistency_pairwise_tip_sample_orbitals(self,
                                                                   n_sample_orbitals: int,
                                                                   tip_orbitals: np.ndarray,
                                                                   lcao_coefficients: np.ndarray,
                                                                   didv_backend: Callable,
                                                                   rtol: float):
        """Check dI/dV consistency between the official C++ backend and alternative backends
        for single-orbital tip/sample pairs,
        using s/p tip orbitals with s/p/d sample orbitals and s/p/d tip orbitals with s/p sample orbitals."""
        actual_didv = self._didv(n_sample_orbitals=n_sample_orbitals,
                                 lcao_coefficients=lcao_coefficients,
                                 tip_orbitals=tip_orbitals,
                                 backend=didv_backend)
        expected_didv = self._didv(n_sample_orbitals=n_sample_orbitals,
                                   lcao_coefficients=lcao_coefficients,
                                   tip_orbitals=tip_orbitals,
                                   backend=ProbeSTM.dIdV_sp_sp)
        np.testing.assert_allclose(actual_didv, expected_didv, rtol=rtol)

    def _didv(self,
              n_sample_orbitals: int,
              tip_orbitals: np.ndarray,
              lcao_coefficients: np.ndarray,
              backend: Callable) -> np.ndarray:
        """Compute differential conductance using the specified backend.

        Args:
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
                       tip_coes=tip_orbitals,
                       orb_t=n_sample_orbitals)

    def pytest_generate_tests(self, metafunc):
        if "didv_backend" in metafunc.fixturenames and \
           "rtol" in metafunc.fixturenames:
            metafunc.parametrize(("didv_backend", "rtol"),
                                 ((fn, rtol) for _, fn, rtol in self._BACKENDS_NAME_FN_RTOL),
                                 ids=(name for name, _, _ in self._BACKENDS_NAME_FN_RTOL))

        if "n_sample_orbitals" in metafunc.fixturenames and \
           "tip_orbitals" in metafunc.fixturenames and \
           "lcao_coefficients" in metafunc.fixturenames:
            metafunc.parametrize(
                ("n_sample_orbitals", "tip_orbitals", "lcao_coefficients"),
                (
                    *self._tip_sample_orbital_pairs(self._SPD_TIP_ORBITALS, self._N_SP_ORBITALS),
                    *self._tip_sample_orbital_pairs(self._SP_TIP_ORBITALS, self._N_SPD_ORBITALS),
                ),
                ids=(
                    *self._tip_sample_orbital_pair_ids(self._SPD_TIP_ORBITAL_NAMES, self._SP_SAMPLE_ORBITAL_NAMES),
                    *self._tip_sample_orbital_pair_ids(self._SP_TIP_ORBITAL_NAMES,  self._SPD_SAMPLE_ORBITAL_NAMES),
                )
            )

    @classmethod
    def _tip_sample_orbital_pairs(
            cls,
            tip_orbitals: List[np.ndarray],
            n_sample_orbitals: int
    ) -> Iterator[Tuple[int, np.ndarray, np.ndarray]]:
        return itertools.product(
            [n_sample_orbitals],
            tip_orbitals,
            cls._lcao_coefficients(n_sample_orbitals),
        )

    @staticmethod
    def _tip_sample_orbital_pair_ids(
            tip_orbital_names: List[str],
            sample_orbital_names: List[str]
    ) -> List[str]:
        return [
            f"#sOrbs={n_sample_orbitals}_tOrb={tip_orb}_sOrb={sample_orb}"
            for n_sample_orbitals, tip_orb, sample_orb in
                itertools.product([len(sample_orbital_names)], tip_orbital_names, sample_orbital_names)
        ]

    @classmethod
    def _lcao_coefficients(cls, n_sample_orbitals: int) -> List[np.ndarray]:
        return [
            np.eye(N=1, M=n_sample_orbitals, k=i)
            for i in range(n_sample_orbitals)
        ]

    _N_SP_ORBITALS = 4
    _N_SPD_ORBITALS = 9
    _SPD_SAMPLE_ORBITAL_NAMES = [
        "s",
        "py",
        "pz",
        "px",
        "dxy",
        "dyz",
        "dz2",
        "dxz",
        "dx2y2",
    ]
    _SP_SAMPLE_ORBITAL_NAMES = _SPD_SAMPLE_ORBITAL_NAMES[:_N_SP_ORBITALS]
    _SPD_TIP_ORBITAL_NAMES = [
        "s",
        "py",
        "pz",
        "px",
        "dyz",
        "dz2",
        "dxz",
    ]
    _SP_TIP_ORBITAL_NAMES = _SPD_TIP_ORBITAL_NAMES[:_N_SP_ORBITALS]
    _SPD_TIP_ORBITAL_INDICES = (0, 1, 2, 3, 5, 6, 7)
    _SP_TIP_ORBITAL_INDICES = _SPD_TIP_ORBITAL_INDICES[:_N_SP_ORBITALS]
    _SPD_TIP_ORBITALS = [
        np.eye(N=1, M=9, k=i)[0]
        for i in _SPD_TIP_ORBITAL_INDICES
    ]
    _SP_TIP_ORBITALS = [
        np.eye(N=1, M=9, k=i)[0]
        for i in _SP_TIP_ORBITAL_INDICES
    ]
