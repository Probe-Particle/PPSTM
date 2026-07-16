import warnings
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Tuple

import numpy as np

import pytest

from pyPPSTM.ReadSTM import read_FIREBALL_all


class _TestReadDFTDeprecated(ABC):
    _ORIG_FERMI: float
    _SAMPLE_ORBS = 'sp'
    _SAMPLE_ORBITAL_COUNT = 4  # sample orbs: 'sp'
    _CUT_ATOMS = -1
    _LOWER_ATOMS = []
    _LOWER_COEFFS = []
    _FERMI = 0.0
    _CUT_MIN = -0.045
    _CUT_MAX = -0.005
    _PBC = [0, 0]
    # LVS = False
    _CELL = np.zeros((3, 3))  # LVS is not {list, tuple, np.ndarray, str} and PBC = [0, 0]
    _IMAGINARY = False
    _EIGENENERGIES_ALL = np.arange(-0.04, 0.04 + 0.01, 0.01)
    _SAMPLE_ATOM_ELEMENTS = ["H"]
    _SAMPLE_ATOM_POSITIONS_EXPECTED = np.asarray([[0.0, 0.0, 0.0]])
    _LCAO_KS_COEFFICIENTS_EXPECTED = np.diag([1.0,]*_SAMPLE_ORBITAL_COUNT)
    _LCAO_KS_COEFFICIENTS_ALL = np.pad(_LCAO_KS_COEFFICIENTS_EXPECTED,
                                       ((0, len(_EIGENENERGIES_ALL) - len(_LCAO_KS_COEFFICIENTS_EXPECTED)),
                                        (0, 0),),
                                       constant_values=0)

    def test_n_eigenenergies_eq_n_lcao_rows(self,
                                            eigenenergies_actual: np.ndarray,
                                            lcao_ks_coefficients_actual: np.ndarray):
        assert len(eigenenergies_actual) == len(lcao_ks_coefficients_actual)

    def test_n_lcao_cols_eq_n_sample_positions_times_n_sample_orbitals(self,
                                                                       eigenenergies_actual: np.ndarray,
                                                                       lcao_ks_coefficients_actual: np.ndarray,
                                                                       sample_atom_positions_actual: np.ndarray):
        assert lcao_ks_coefficients_actual.shape[1] == self._SAMPLE_ORBITAL_COUNT * len(sample_atom_positions_actual)

    def test_sample_atom_positions(self, sample_atom_positions_actual: np.ndarray):
        np.testing.assert_equal(sample_atom_positions_actual, self._SAMPLE_ATOM_POSITIONS_EXPECTED)

    def test_lcao_ks_coefficients(self, lcao_ks_coefficients_actual: np.ndarray):
        np.testing.assert_equal(lcao_ks_coefficients_actual, self._LCAO_KS_COEFFICIENTS_EXPECTED)

    @pytest.fixture
    def eigenenergies_actual(self, dft_actual: Tuple[np.ndarray, np.ndarray, np.ndarray]):
        yield dft_actual[0]

    @pytest.fixture
    def lcao_ks_coefficients_actual(self, dft_actual: Tuple[np.ndarray, np.ndarray, np.ndarray]):
        yield dft_actual[1]

    @pytest.fixture
    def sample_atom_positions_actual(self, dft_actual: Tuple[np.ndarray, np.ndarray, np.ndarray]):
        yield dft_actual[2]

    @pytest.fixture(scope="class")
    @abstractmethod
    def dft_actual(self, geometry_path: Path, input_dir: Path):
        pass

    @pytest.fixture(scope="class")
    def input_dir(self, tmp_path_factory):
        yield tmp_path_factory.mktemp("input")


class TestReadDFTFireballDeprecated(_TestReadDFTDeprecated):
    _ORIG_FERMI = 0.0  # @TODO: switch to something different from zero to make the tests stronger
    _GEOMETRY_FILE_NAME = 'input_plot.xyz'

    @pytest.fixture(scope="class")
    def dft_actual(self, geometry_path: Path, input_dir: Path):
        eigenenergies, lcao_ks_coefficients, sample_atom_positions = read_FIREBALL_all(
            name=str(input_dir / "phik_0001_"),
            geom=str(geometry_path),
            lvs=self._CELL,
            fermi=self._FERMI,
            pbc=self._PBC,
            orbs=self._SAMPLE_ORBS,
            cut_min=self._CUT_MIN,
            cut_max=self._CUT_MAX,
            cut_at=self._CUT_ATOMS,
            lower_atoms=self._LOWER_ATOMS,
            lower_coefs=self._LOWER_COEFFS,
            imaginary=self._IMAGINARY,
        )
        yield eigenenergies, lcao_ks_coefficients, sample_atom_positions

    @pytest.fixture(scope="class")
    def geometry_path(self, input_dir: Path):
        geometry_path = input_dir / self._GEOMETRY_FILE_NAME
        with open(geometry_path, "w") as f:
            f.writelines([
                f"{len(self._SAMPLE_ATOM_POSITIONS_EXPECTED)}\n",
                *[f"{el} {str(pos).strip("[]")}\n" for el, pos
                  in zip(self._SAMPLE_ATOM_ELEMENTS, self._SAMPLE_ATOM_POSITIONS_EXPECTED)],
            ])
        yield geometry_path

    @pytest.fixture(scope="class", autouse=True)
    def phik_0001_s(self, input_dir: Path):
        data = np.stack([
            self._EIGENENERGIES_ALL,
            self._LCAO_KS_COEFFICIENTS_ALL[:,0],
        ], axis=1)
        np.savetxt(input_dir / "phik_0001_s.dat", data,
                   header=f"{len(self._SAMPLE_ATOM_POSITIONS_EXPECTED)} {len(self._EIGENENERGIES_ALL)} {self._ORIG_FERMI}",
                   comments='')

    @pytest.fixture(scope="class", autouse=True)
    def phik_0001_py(self, input_dir: Path):
        data = np.stack([
            self._EIGENENERGIES_ALL,
            self._LCAO_KS_COEFFICIENTS_ALL[:,1],
        ], axis=1)
        np.savetxt(input_dir / "phik_0001_py.dat", data,
                   header=f"{len(self._SAMPLE_ATOM_POSITIONS_EXPECTED)} {len(self._EIGENENERGIES_ALL)} {self._ORIG_FERMI}",
                   comments='')

    @pytest.fixture(scope="class", autouse=True)
    def phik_0001_pz(self, input_dir: Path):
        data = np.stack([
            self._EIGENENERGIES_ALL,
            self._LCAO_KS_COEFFICIENTS_ALL[:,2],
        ], axis=1)
        np.savetxt(input_dir / "phik_0001_pz.dat", data,
                   header=f"{len(self._SAMPLE_ATOM_POSITIONS_EXPECTED)} {len(self._EIGENENERGIES_ALL)} {self._ORIG_FERMI}",
                   comments='')

    @pytest.fixture(scope="class", autouse=True)
    def phik_0001_px(self, input_dir: Path):
        data = np.stack([
            self._EIGENENERGIES_ALL,
            self._LCAO_KS_COEFFICIENTS_ALL[:,3],
        ], axis=1)
        np.savetxt(input_dir / "phik_0001_px.dat", data,
                   header=f"{len(self._SAMPLE_ATOM_POSITIONS_EXPECTED)} {len(self._EIGENENERGIES_ALL)} {self._ORIG_FERMI}",
                   comments='')

# class TestReadDFTAims(_TestReadDFT):
#     _ORIG_FERMI = 0.0
#
#     _NAME =
#     # if ((spin == None) or (spin == False)):
#     #     name = 'KS_eigenvectors.band_1.kpt_1.out'
#     # elif spin in ['up', 'alpha']:
#     #     name = 'KS_eigenvectors_up.band_1.kpt_1.out'
#     # elif spin in ['down', 'beta', 'dn']:
#     #     name = 'KS_eigenvectors_dn.band_1.kpt_1.out'
#     # elif spin == 'both':
#     #     name_up = 'KS_eigenvectors_up.band_1.kpt_1.out'
#     #     name_dn = 'KS_eigenvectors_dn.band_1.kpt_1.out'
#     # else:
#     #     raise ValueError(f"Unknown spin: {spin}")
#
#     _GEOMETRY_FILE_NAME = 'input.xyz'
#
#     @pytest.fixture(scope="class")
#     def dft_actual(self, geometry_path: Path, input_dir: Path):
#         eigenenergies, lcao_ks_coefficients, sample_atom_positions = read_AIMS_all(
#             name=str(input_dir / self._NAME),
#             geom=str(geometry_path),
#             fermi=self._FERMI,
#             pbc=self._PBC,
#             orbs=self._SAMPLE_ORBS,
#             cut_min=self._CUT_MIN,
#             cut_max=self._CUT_MAX,
#             cut_at=self._CUT_ATOMS,
#             lower_atoms=self._LOWER_ATOMS,
#             lower_coefs=self._LOWER_COEFFS,
#             imaginary=self._IMAGINARY,
#         )
#         yield eigenenergies, lcao_ks_coefficients, sample_atom_positions
