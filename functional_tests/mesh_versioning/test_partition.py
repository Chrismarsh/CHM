"""Run with python3 test_partition.py /path/to/partition [--mpiexec /path/to/mpirun]."""
import argparse
from pathlib import Path
import json
import os
import shutil
import subprocess
import tempfile

import h5py
import numpy as np

FIXTURE = Path(__file__).resolve().parent


def launch(binary, directory, args, workers=1, launcher=None, expect_success=True):
    cmd = [str(binary), "-m", "mesh.h5", "-p", "params.h5", "--mpi-ranks", "31", *args]
    if workers != 1:
        cmd = [launcher, "-n", str(workers), *cmd]
    run = subprocess.run(cmd, cwd=directory, capture_output=True, text=True, timeout=120)
    assert (run.returncode == 0) == expect_success, (cmd, run.returncode, run.stdout, run.stderr)
    return run.stdout


def verify(directory, kept, radius, standalone=False):
    """Independent whole-mesh reference using CHM's original per-seed graph search."""
    suffix = f".standalone.{kept[0]}" if standalone else ""
    doc = json.loads((directory / f"mesh.np31{suffix}.partition").read_text())
    with h5py.File(directory / "mesh.h5") as src:
        elem = src["mesh/elem"][:]
        xyz = src["mesh/vertex"][:]
        nbr = src["mesh/neighbor"][:]
        sizes = src["mesh/local_sizes"][:]
        permutation = src["mesh/cell_global_id"][:]
        geographic = src.attrs["/mesh/is_geographic"][0]
    with h5py.File(directory / "params.h5") as src:
        fields = {k: v[:] for k, v in src["parameters"].items()}
    if (directory / "extra.h5").exists():
        with h5py.File(directory / "extra.h5") as src:
            fields.update({k: v[:] for k, v in src["parameters"].items()})
    center = xyz[elem].mean(axis=1)[:, :2]
    owner = np.empty(len(elem), dtype=int)
    owner[permutation] = np.repeat(np.arange(len(sizes)), sizes)
    assert int(doc["ranks"]) == len(kept)
    assert int(doc["num_global_faces"]) == sum(sizes[kept])
    assert float(doc["max_ghost_distance"]) == radius
    np.testing.assert_array_equal(list(map(int, doc["local_sizes"])), sizes[kept])
    assert len(doc["meshes"]) == len(doc["parameters"]) == len(kept)
    ids = np.arange(len(elem))
    if len(kept) < len(sizes):
        selected = permutation[np.isin(owner[permutation], kept)]
        ids[selected] = np.arange(len(selected))

    for index, rank in enumerate(kept):
        owned = set(np.flatnonzero(owner == rank))
        labels = {f: 0 for f in owned}
        seeds = [f for f in owned if any(n not in owned for n in nbr[f])]
        if not standalone:
            for f in seeds:
                for n in nbr[f]:
                    if n >= 0 and n not in owned:
                        labels[n] = 1
            for seed in seeds:
                a = center[seed]
                if geographic:
                    x, y = np.radians(center).T
                    dx, dy = x-x[seed], y-y[seed]
                    h = np.sin(dy/2)**2 + np.cos(y[seed])*np.cos(y)*np.sin(dx/2)**2
                    within = 2*6378137*np.arctan2(np.sqrt(h), np.sqrt(1-h)) <= radius
                else:
                    within = np.linalg.norm(center-a, axis=1) <= radius
                pending, visited = [seed], {seed}
                while pending:
                    f = pending.pop()
                    labels.setdefault(f, 2)
                    for n in nbr[f]:
                        if n >= 0 and n not in visited and within[n]:
                            visited.add(n)
                            pending.append(n)
        rows = sorted((r for r in labels if owner[r] in kept), key=ids.__getitem__)
        expected_nbr = np.where(np.isin(nbr[rows], rows), ids[nbr[rows]], -1)
        local_verts, local_elem = np.unique(elem[rows].reshape(-1), return_inverse=True)
        with h5py.File(directory / doc["meshes"][index]) as f:
            g = f["mesh"]
            assert f.attrs["/mesh/version"][0] == b"3.0.0"
            assert bool(f.attrs["/mesh/is_partition"][0]) == (not standalone)
            assert bool(f.attrs["/mesh/is_geographic"][0]) == bool(geographic)
            np.testing.assert_array_equal(g["cell_global_id"][:], ids[rows])
            np.testing.assert_array_equal(g["ghost_type"][:], [labels[r] for r in rows])
            np.testing.assert_array_equal(g["owner"][:], [kept.index(owner[r]) for r in rows])
            np.testing.assert_array_equal(g["local_sizes"][:], sizes[kept])
            np.testing.assert_array_equal(g["neighbor"][:], expected_nbr)
            np.testing.assert_array_equal(g["elem"][:], local_elem.reshape(-1, 3))
            np.testing.assert_array_equal(g["vertex"][:], xyz[local_verts])
        with h5py.File(directory / doc["parameters"][index]) as f:
            assert set(f["parameters"]) == set(fields)
            for name, column in fields.items():
                np.testing.assert_array_equal(f["parameters"][name][:], column[rows])


def exercise(binary, launcher):
    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)
        shutil.copyfile(FIXTURE / "slope.metis_mesh.h5", directory / "mesh.h5")
        shutil.copyfile(FIXTURE / "slope.metis_param.h5", directory / "params.h5")
        workers = 3 if launcher else 1
        all_ranks = list(range(31))
        for radius in [0, 0.25, 100, 400]:
            launch(binary, directory, ["-g", str(radius)], workers, launcher)
            verify(directory, all_ranks, radius)

        with h5py.File(directory / "extra.h5", "w") as extra:
            extra.create_dataset("parameters/face_tag", data=np.arange(2618) + 0.25)
        args = ["-p", "extra.h5", "-g", "100"]
        launch(binary, directory, [*args, "-v", "9", "-v", "2"], workers, launcher)
        verify(directory, [2, 9], 100)
        launch(binary, directory, [*args, "-s", "12", "--write-vtu"], workers, launcher)
        verify(directory, [12], 100, standalone=True)
        assert (directory / "rank.12.vtu").is_file()

        with h5py.File(directory / "mesh.h5", "r+") as src:
            src["mesh/cell_global_id"][:] = np.arange(2618)[::-1]
            p = src["mesh/vertex"][:]
            p[:, :2] = (p[:, :2] - p[0, :2]) / 100000 + [-115, 60]
            src["mesh/vertex"][:] = p
            src.attrs.modify("/mesh/is_geographic", [1])
        launch(binary, directory, args, workers, launcher)
        verify(directory, all_ranks, 100)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("partition", type=Path)
    parser.add_argument("--mpiexec")
    options = parser.parse_args()
    exercise(options.partition.resolve(), options.mpiexec)
    print("Partition equivalence checks passed.")
