"""Run CHM through a snowfall and rain event and compare its NetCDF restarts.

Usage: python check_restart.py /path/to/CHM
Add --partition /path/to/partition --mpiexec /path/to/mpiexec for two-rank parity.
"""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import subprocess
import tempfile

import netCDF4
import numpy as np


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")


def run(command, directory, threads=1):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    completed = subprocess.run(command, cwd=directory, env=env, capture_output=True, text=True)
    if completed.returncode:
        raise RuntimeError(f"{command} failed\n{completed.stdout[-4000:]}\n{completed.stderr[-4000:]}")


def inputs(directory):
    mesh = {
        "version": "2.0.0", "nvertex": 4, "nelem": 2, "is_geographic": 0,
        "proj4": "+proj=utm +zone=11 +datum=WGS84 +units=m +no_defs",
        "vertex": [[500000, 5538600, 1000], [500100, 5538600, 1000],
                   [500000, 5538700, 1000], [500200, 5538700, 1020]],
        "elem": [[0, 1, 2], [1, 3, 2]], "neigh": [[1, -1, -1], [-1, 0, -1]],
        "partition_method": "metis", "local_size": [1, 1], "cell_global_id": [0, 1],
    }
    write_json(directory / "test.mesh", {"mesh": mesh})
    write_json(directory / "test.param", {"area": [5000, 10000]})
    start = dt.datetime(2020, 1, 2)
    with (directory / "forcing.txt").open("w", encoding="utf-8") as stream:
        stream.write("datetime t rh u p Qli iswr_direct iswr_diffuse\n")
        for i in range(30):
            date = (start + dt.timedelta(hours=i)).strftime("%Y%m%dT%H%M%S")
            cold = i < 13
            snow = 2 if i < 8 else .00005 if i == 9 else 0
            stream.write(f"{date} {-6 if cold else 4} 80 2 {snow if cold else int(i == 15)} "
                         f"{280 if cold else 340} {0 if cold else 150} 25\n")


def checkpoints(directory):
    meta = sorted(directory.glob("checkpoint/checkpoint_*.json"))[-1]
    paths = [meta.parent / row for row in json.loads(meta.read_text())["files"]]
    fields = {}
    ids = []
    for path in paths:
        with netCDF4.Dataset(path) as f:
            ids.extend(f["global_id"][:])
            for name in ("crocus_state", "crocus_pending_snow", "crocus_fluxes", "crocus_cumulative_subl"):
                fields.setdefault(name, []).append(np.asarray(f[name][:]))
            assert f["crocus_state"].ground_temperature_k == 268.0
    order = np.argsort(ids)
    return meta, {key: np.concatenate(rows)[order] for key, rows in fields.items()}


def scenario(directory, chm, name, mesh, end, restart=None, threads=1, launcher=()):
    out = directory / name
    config = {
        "option": {"enddate": end,
                   "per_triangle_timeseries": False, "ui": False},
        "modules": ["point_mode", "solar", "threshold_p_phase", "Crocus"],
        "config": {"Crocus": {"ground_temperature_k": 268.0},
                   "point_mode": {"provide": {"U_R": False, "vw_dir": False, "iswr": False,
                                              "iswr_diffuse": True, "iswr_direct": True}}},
        "meshes": {"mesh": str(mesh.relative_to(directory))},
        "forcing": {"UTC_offset": 0, "interpolant": "nearest", "num_stations_to_use": 1, "station": {"file": "forcing.txt", "longitude": -117.,
                    "latitude": 50., "elevation": 1000}},
        "output": {"output_dir": name},
        "checkpoint": {"save_checkpoint": True, "on_last": True},
    }
    if mesh.suffix == ".mesh":
        config["meshes"]["parameters"] = {"file": "test.param"}
    if restart:
        config["checkpoint"]["load_checkpoint_path"] = str(restart.relative_to(directory))
    path = directory / f"{name}.json"
    write_json(path, config)
    run([*launcher, chm, "-f", path], directory, threads)
    return checkpoints(out)


def compare(expected, actual):
    for name in expected:
        np.testing.assert_array_equal(expected[name], actual[name], err_msg=name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chm", type=Path)
    parser.add_argument("--partition", type=Path)
    parser.add_argument("--mpiexec", default="mpiexec")
    args = parser.parse_args()
    chm = args.chm.resolve()
    with tempfile.TemporaryDirectory() as folder:
        directory = Path(folder)
        inputs(directory)
        mesh = directory / "test.mesh"
        end = "20200103T040000"
        _, serial = scenario(directory, chm, "continuous", mesh, end)
        split, first = scenario(directory, chm, "first", mesh, "20200102T090000")
        assert np.all(first["crocus_state"][:, 0].sum(axis=1) > 0)
        _, resumed = scenario(directory, chm, "resumed", mesh, end, restart=split)
        compare(serial, resumed)
        _, threaded = scenario(directory, chm, "threaded", mesh, end, threads=2)
        compare(serial, threaded)
        if args.partition:
            run([args.partition.resolve(), "--mesh-file", mesh, "--param-file", "test.param", "--mpi-ranks", "2"], directory)
            partition = next(directory.rglob("*.np2.partition"))
            launcher = (args.mpiexec, "-n", "2")
            _, parallel = scenario(directory, chm, "parallel", partition, end, launcher=launcher)
            compare(serial, parallel)
            split, _ = scenario(directory, chm, "parallel_first", partition, "20200102T090000", launcher=launcher)
            _, resumed = scenario(directory, chm, "parallel_resumed", partition, end, restart=split, launcher=launcher)
            compare(parallel, resumed)
    print("Crocus CHM restart and decomposition comparisons passed")


if __name__ == "__main__":
    main()
