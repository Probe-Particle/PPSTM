from pathlib import Path
from typing import Tuple

import numpy as np
import pytest

from pyPPSTM.ReadSTM import read_FIREBALL_all


class TestReadDFTEigenenergiesFireball:
    _CUT_ATOMS = -1
    _LOWER_ATOMS = []
    _LOWER_COEFFS = []
    _CUT_MIN = -0.045
    _CUT_MAX = -0.005
    _PBC = [0, 0]
    # LVS = False
    _CELL = np.zeros((3, 3))  # LVS is not {list, tuple, np.ndarray, str} and PBC = [0, 0]
    _IMAGINARY = False
    _EIGENENERGIES_ALL = np.arange(-0.04, 0.04 + 0.01, 0.01)
    _SAMPLE_ATOM_ELEMENTS = ["H"]
    _SAMPLE_ATOM_POSITIONS_EXPECTED = np.asarray([[0.0, 0.0, 0.0]])
    _SPD_SAMPLE_ORBITAL_COUNT = 9
    _LCAO_KS_COEFFICIENTS_SPD_FULL = np.pad(np.diag([1.0,]*_SPD_SAMPLE_ORBITAL_COUNT),
                                            ((0, len(_EIGENENERGIES_ALL) - _SPD_SAMPLE_ORBITAL_COUNT),
                                            (0, 0),),
                                            constant_values=0)
    _FILE_NAME = "phik_0001_"
    _GEOMETRY_FILE_NAME = 'input_plot.xyz'

    @pytest.mark.parametrize(
        (    "dft_files",          "dft_actual"), [
            ({"orig_fermi":  0.0}, {"sample_orb":"sp",  "cut_min": -0.045, "cut_max": -0.005, "fermi":  0.0}),
            ({"orig_fermi": -0.1}, {"sample_orb":"sp",  "cut_min": -0.045, "cut_max": -0.005, "fermi":  0.0}),
            ({"orig_fermi":  0.0}, {"sample_orb":"spd", "cut_min": -2.5,   "cut_max":  2.5,   "fermi":  0.0}),
            ({"orig_fermi": -0.1}, {"sample_orb":"spd", "cut_min": -2.5,   "cut_max":  2.5,   "fermi":  0.0}),
            ({"orig_fermi":  0.0}, {"sample_orb":"sp",  "cut_min": -0.045, "cut_max": -0.005, "fermi": -0.035}),
            ({"orig_fermi": -0.1}, {"sample_orb":"sp",  "cut_min": -0.045, "cut_max": -0.005, "fermi": -0.035}),
            ({"orig_fermi":  0.0}, {"sample_orb":"spd", "cut_min": -2.5,   "cut_max":  2.5,   "fermi": -0.035}),
            ({"orig_fermi": -0.1}, {"sample_orb":"spd", "cut_min": -2.5,   "cut_max":  2.5,   "fermi": -0.035}),
        ], indirect=True, ids=[
            "sp,cm-0.045,cM-0.005,f0,of0",
            "sp,cm-0.045,cM-0.005,f0,of-.1",
            "spd,cm-2.5,cM2.5,f0,of0",
            "spd,cm-2.5,cM2.5,f0,of-.1",
            "sp,cm-0.045,cM-0.005,f-.035,of0",
            "sp,cm-0.045,cM-0.005,f-.035,of-.1",
            "spd,cm-2.5,cM2.5,f-.035,of0",
            "spd,cm-2.5,cM2.5,f-.035,of-.1",
        ])
    def test_n_eigenenergies_eq_n_lcao_rows(self,
                                            eigenenergies_actual: np.ndarray,
                                            lcao_ks_coefficients_actual: np.ndarray,
                                            dft_files):
        assert len(eigenenergies_actual) == len(lcao_ks_coefficients_actual)

    @pytest.mark.parametrize(
        (    "dft_files",          "dft_actual",                                                                "cut_min"), [
            ({"orig_fermi":  0.0}, {"sample_orb":"sp",  "cut_min": -0.045, "cut_max": -0.005, "fermi":  0.0},   -0.045),
            ({"orig_fermi": -0.1}, {"sample_orb":"sp",  "cut_min": -0.045, "cut_max": -0.005, "fermi":  0.0},   -0.045),
            ({"orig_fermi":  0.0}, {"sample_orb":"spd", "cut_min": -2.5,   "cut_max":  2.5,   "fermi":  0.0},   -2.5),
            ({"orig_fermi": -0.1}, {"sample_orb":"spd", "cut_min": -2.5,   "cut_max":  2.5,   "fermi":  0.0},   -2.5),
            ({"orig_fermi":  0.0}, {"sample_orb":"sp",  "cut_min": -0.045, "cut_max": -0.005, "fermi": -0.035}, -0.045),
            ({"orig_fermi": -0.1}, {"sample_orb":"sp",  "cut_min": -0.045, "cut_max": -0.005, "fermi": -0.035}, -0.045),
            ({"orig_fermi":  0.0}, {"sample_orb":"spd", "cut_min": -2.5,   "cut_max":  2.5,   "fermi": -0.035}, -2.5),
            ({"orig_fermi": -0.1}, {"sample_orb":"spd", "cut_min": -2.5,   "cut_max":  2.5,   "fermi": -0.035}, -2.5),
        ], indirect=["dft_actual"], ids=[
            "sp,cm-0.045,cM-0.005,f0,of0",
            "sp,cm-0.045,cM-0.005,f0,of-.1",
            "spd,cm-2.5,cM2.5,f0,of0",
            "spd,cm-2.5,cM2.5,f0,of-.1",
            "sp,cm-0.045,cM-0.005,f-.035,of0",
            "sp,cm-0.045,cM-0.005,f-.035,of-.1",
            "spd,cm-2.5,cM2.5,f-.035,of0",
            "spd,cm-2.5,cM2.5,f-.035,of-.1",
        ])
    def test_eigenenergies_gt_cut_min(self,
                                      eigenenergies_actual: np.ndarray,
                                      cut_min: float,
                                      dft_files):
        assert (eigenenergies_actual > cut_min).all()

    @pytest.mark.parametrize(
        (    "dft_files",          "dft_actual",                                                                "cut_max"), [
            ({"orig_fermi":  0.0}, {"sample_orb":"sp",  "cut_min": -0.045, "cut_max": -0.005, "fermi":  0.0},   -0.005),
            ({"orig_fermi": -0.1}, {"sample_orb":"sp",  "cut_min": -0.045, "cut_max": -0.005, "fermi":  0.0},   -0.005),
            ({"orig_fermi":  0.0}, {"sample_orb":"spd", "cut_min": -2.5,   "cut_max":  2.5,   "fermi":  0.0},    2.5),
            ({"orig_fermi": -0.1}, {"sample_orb":"spd", "cut_min": -2.5,   "cut_max":  2.5,   "fermi":  0.0},    2.5),
            ({"orig_fermi":  0.0}, {"sample_orb":"sp",  "cut_min": -0.045, "cut_max": -0.005, "fermi": -0.035}, -0.005),
            ({"orig_fermi": -0.1}, {"sample_orb":"sp",  "cut_min": -0.045, "cut_max": -0.005, "fermi": -0.035}, -0.005),
            ({"orig_fermi":  0.0}, {"sample_orb":"spd", "cut_min": -2.5,   "cut_max":  2.5,   "fermi": -0.035},  2.5),
            ({"orig_fermi": -0.1}, {"sample_orb":"spd", "cut_min": -2.5,   "cut_max":  2.5,   "fermi": -0.035},  2.5),
    ], indirect=["dft_actual"], ids=[
        "sp,cm-0.045,cM-0.005,f0,of0",
        "sp,cm-0.045,cM-0.005,f0,of-.1",
        "spd,cm-2.5,cM2.5,f0,of0",
        "spd,cm-2.5,cM2.5,f0,of-.1",
        "sp,cm-0.045,cM-0.005,f-.035,of0",
        "sp,cm-0.045,cM-0.005,f-.035,of-.1",
        "spd,cm-2.5,cM2.5,f-.035,of0",
        "spd,cm-2.5,cM2.5,f-.035,of-.1",
    ])
    def test_eigenenergies_lt_cut_max(self,
                                      eigenenergies_actual: np.ndarray,
                                      cut_max: float,
                                      dft_files):
        assert (eigenenergies_actual < cut_max).all()

    @pytest.mark.parametrize(
        (    "dft_files", "dft_actual",                                                                         "cut_min", "cut_max", "fermi", "orig_fermi"), [
            ({"orig_fermi":  0.0}, {"sample_orb":"sp",  "cut_min": -0.045, "cut_max": -0.005, "fermi":  0.0},   -0.045,    -0.005,     0.0,     0.0),
            ({"orig_fermi": -0.1}, {"sample_orb":"sp",  "cut_min": -0.045, "cut_max": -0.005, "fermi":  0.0},   -0.045,    -0.005,     0.0,    -0.1),
            ({"orig_fermi":  0.0}, {"sample_orb":"spd", "cut_min": -2.5,   "cut_max":  2.5,   "fermi":  0.0},   -2.5,       2.5,       0.0,     0.0),
            ({"orig_fermi": -0.1}, {"sample_orb":"spd", "cut_min": -2.5,   "cut_max":  2.5,   "fermi":  0.0},   -2.5,       2.5,       0.0,    -0.1),
            ({"orig_fermi":  0.0}, {"sample_orb":"sp",  "cut_min": -0.045, "cut_max": -0.005, "fermi": -0.035}, -0.045,    -0.005,    -0.035,   0.0),
            ({"orig_fermi": -0.1}, {"sample_orb":"sp",  "cut_min": -0.045, "cut_max": -0.005, "fermi": -0.035}, -0.045,    -0.005,    -0.035,  -0.1),
            ({"orig_fermi":  0.0}, {"sample_orb":"spd", "cut_min": -2.5,   "cut_max":  2.5,   "fermi": -0.035}, -2.5,       2.5,      -0.035,   0.0),
            ({"orig_fermi": -0.1}, {"sample_orb":"spd", "cut_min": -2.5,   "cut_max":  2.5,   "fermi": -0.035}, -2.5,       2.5,      -0.035,  -0.1),
    ], indirect=["dft_files", "dft_actual"], ids=[
        "sp,cm-0.045,cM-0.005,f0,of0",
        "sp,cm-0.045,cM-0.005,f0,of-.1",
        "spd,cm-2.5,cM2.5,f0,of0",
        "spd,cm-2.5,cM2.5,f0,of-.1",
        "sp,cm-0.045,cM-0.005,f-.035,of0",
        "sp,cm-0.045,cM-0.005,f-.035,of-.1",
        "spd,cm-2.5,cM2.5,f-.035,of0",
        "spd,cm-2.5,cM2.5,f-.035,of-.1",
    ])
    def test_eigenenergies(self,
                           eigenenergies_actual: np.ndarray,
                           cut_min: float, cut_max: float, fermi: float, orig_fermi: float,
                           dft_files):
        eigenenergies = self._EIGENENERGIES_ALL - (fermi + orig_fermi)
        eigenenergies_expected = eigenenergies[
            (cut_min < eigenenergies) & (eigenenergies < cut_max)
        ]
        np.testing.assert_equal(eigenenergies_actual, eigenenergies_expected)

    @pytest.fixture
    def eigenenergies_actual(self, dft_actual: Tuple[np.ndarray, np.ndarray, np.ndarray]):
        yield dft_actual[0]

    @pytest.fixture
    def lcao_ks_coefficients_actual(self, dft_actual: Tuple[np.ndarray, np.ndarray, np.ndarray]):
        yield dft_actual[1]

    @pytest.fixture(scope="class")
    def dft_actual(self, geometry_path: Path, input_dir: Path, request):
        eigenenergies, lcao_ks_coefficients, sample_atom_positions = read_FIREBALL_all(
            name=str(input_dir / self._FILE_NAME),
            geom=str(geometry_path),
            lvs=self._CELL,
            fermi=request.param["fermi"],
            pbc=self._PBC,
            orbs=request.param["sample_orb"],
            cut_min=request.param["cut_min"],
            cut_max=request.param["cut_max"],
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

    @pytest.fixture(scope="class")
    def dft_files(self, input_dir: Path, request):
        for i, sample_orbital in enumerate(["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2y2"]):
            data = np.stack([
                self._EIGENENERGIES_ALL,
                self._LCAO_KS_COEFFICIENTS_SPD_FULL[:, i],
            ], axis=1)
            np.savetxt(input_dir / f"{self._FILE_NAME}{sample_orbital}.dat", data,
                       header=f"{len(self._SAMPLE_ATOM_POSITIONS_EXPECTED)}"
                              f" {len(self._EIGENENERGIES_ALL)}"
                              f" {request.param["orig_fermi"]}",
                       comments='')

    @pytest.fixture(scope="class")
    def input_dir(self, tmp_path_factory):
        yield tmp_path_factory.mktemp("input")
