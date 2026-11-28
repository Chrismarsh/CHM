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
    write_json(directory / "test.param", {"area": [5000, 10000], "CanopyHeight": [6., 0.], "LAI": [2.5, 0.]})
    write_forcing(directory / "forcing.txt", initial_snow=2)
    write_forcing(directory / "pde_forcing.txt", initial_snow=.00005)


def write_forcing(path, initial_snow):
    start = dt.datetime(2020, 1, 2)
    with path.open("w", encoding="utf-8") as stream:
        stream.write("datetime t rh u p Qli iswr_direct iswr_diffuse U_R vw_dir Qsi\n")
        for i in range(30):
            date = (start + dt.timedelta(hours=i)).strftime("%Y%m%dT%H%M%S")
            cold = i < 13
            snow = 2 if i < 8 else .00005 if i == 9 else 0
            if i < 2:
                snow = initial_snow
            stream.write(f"{date} {-6 if cold else 4} 80 2 {snow if cold else int(i == 15)} "
                         f"{280 if cold else 340} {0 if cold else 150} 25 {12 if 5 <= i < 11 else 1} 270 {25 if cold else 175}\n")


def checkpoint_fields(meta):
    paths = [meta.parent / row for row in json.loads(meta.read_text())["files"]]
    fields = {}
    ids = []
    for path in paths:
        with netCDF4.Dataset(path) as f:
            ids.extend(f["global_id"][:])
            names = ["crocus_state", "crocus_pending_snow", "crocus_fluxes", "crocus_cumulative_subl"]
            transport_fields = ("PBSM3D:sum_drift", "snow_slide:delta_avalanche_mass_sum",
                                "snowslide_pde:delta_avalanche_mass", "snowslide_pde:delta_avalanche_snowdepth",
                                "snowslide_pde:delta_avalanche_mass_sum", "snowslide_pde:delta_avalanche_snowdepth_sum",
                                "Simple_Canopy:cum_net_snow", "Simple_Canopy:cum_net_rain")
            names.extend(name for name in transport_fields if name in f.variables)
            for name in names:
                fields.setdefault(name, []).append(np.asarray(f[name][:]))
            assert f["crocus_state"].ground_temperature_k == 268.0
    order = np.argsort(ids)
    return {key: np.concatenate(rows)[order] for key, rows in fields.items()}


def scenario(directory, chm, name, mesh, end, restart=None, threads=1, launcher=(), avalanche=None, every_step=False):
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
    if avalanche:
        config["modules"].remove("threshold_p_phase")
        config["modules"].append("Harder_precip_phase")
        config["modules"].extend(["iswr_from_obs", "iswr", "scale_wind_vert", "Simple_Canopy", "PBSM3D", avalanche])
        config["config"]["point_mode"]["provide"].update(
            {"U_2m_above_srf": False, "U_R": True, "vw_dir": True, "iswr": False,
             "iswr_diffuse": False, "iswr_direct": False})
        config["remove_depency"] = {name: "Crocus" for name in ("scale_wind_vert", "Simple_Canopy", "PBSM3D", avalanche)}
        config["config"]["PBSM3D"] = {"nLayer": 3, "use_tanh_fetch": False, "enable_veg": False,
                                       "min_sd_trans": .001, "use_R94_lambda": False}
        # Lower the holding depth on this short test to exercise routing during a small storm.
        config["config"][avalanche] = {"reserve_drift_erosion": True, "avalache_mult": 1.5}
        if avalanche == "snowslide_pde":
            config["config"][avalanche]["release_profile"] = "chm_vionnet_2021"
            # Expose the avalanche donor to the gust so both modules erode the same face.
            config["config"]["scale_wind_vert"] = {"ignore_canopy": True}
            config["forcing"]["station"]["file"] = "pde_forcing.txt"
    if every_step:
        config["checkpoint"]["on_frequency"] = 1
    if mesh.suffix == ".mesh":
        config["meshes"]["parameters"] = {"file": "test.param"}
    if restart:
        config["checkpoint"]["load_checkpoint_path"] = str(restart.relative_to(directory))
    path = directory / f"{name}.json"
    write_json(path, config)
    run([*launcher, chm, "-f", path], directory, threads)
    meta = sorted(out.glob("checkpoint/checkpoint_*.json"))[-1]
    return meta, checkpoint_fields(meta)


def compare(expected, actual, repartitioned=False):
    for name in expected:
        if repartitioned and (name == "PBSM3D:sum_drift" or name.startswith("snowslide_pde:")):
            # Distributed transport solves can differ at double roundoff.
            np.testing.assert_allclose(expected[name], actual[name], rtol=1e-11, atol=1e-10, err_msg=name)
        else:
            np.testing.assert_array_equal(expected[name], actual[name], err_msg=name)


def check_pde_transfers(directory, areas):
    drift_sum = np.zeros_like(areas)
    precip_sum = np.zeros_like(areas)
    swe = np.zeros_like(areas)
    trace_only = np.zeros(len(areas), dtype=bool)
    saw_trace = saw_shared_donor = False
    for meta in sorted(directory.glob("checkpoint/checkpoint_*.json")):
        fields = checkpoint_fields(meta)
        volume = fields["snowslide_pde:delta_avalanche_mass"]  # m³ water per face; exterior edges are closed
        np.testing.assert_allclose(volume.sum(), 0., atol=1e-7, err_msg=meta.name)
        avalanche = 1000. * volume / areas
        drift = fields["PBSM3D:sum_drift"] - drift_sum
        demand = np.maximum(-drift, 0.) + np.maximum(-avalanche, 0.)
        assert np.all(demand <= swe + 1e-6), f"Transport overdrew a donor at {meta.name}"
        saw_shared_donor |= np.any((drift < -1e-7) & (avalanche < -1e-7))
        if np.any(trace_only):
            saw_trace = True
            assert np.all(avalanche[trace_only] >= 0.), f"Routed pending Crocus snow at {meta.name}"
        pending = fields["crocus_pending_snow"]
        layers = fields["crocus_state"][:, 0].sum(axis=1, dtype=np.float64)
        precip = fields["Simple_Canopy:cum_net_snow"] + fields["Simple_Canopy:cum_net_rain"]
        # The hourly Crocus fluxes store runoff and sublimation in kg/m²/s.
        loss = 3600. * fields["crocus_fluxes"][:, 4:6].sum(axis=1, dtype=np.float64)
        expected_swe = swe + (precip - precip_sum) + drift + avalanche - loss
        swe = layers + pending
        np.testing.assert_allclose(swe, expected_swe, atol=2e-4, rtol=2e-5, err_msg=meta.name)
        drift_sum = fields["PBSM3D:sum_drift"]
        precip_sum = precip
        trace_only = (pending > 0.) & (layers == 0.)
    assert saw_trace, "Expected a PDE step starting with unlayered snowfall"
    assert saw_shared_donor, "Expected simultaneous wind and avalanche erosion on one face"


def coupled_event(directory, chm, mesh, end, avalanche, launcher=()):
    name = avalanche + ("_mpi" if launcher else "_serial")
    pde = avalanche == "snowslide_pde"
    _, continuous = scenario(directory, chm, name, mesh, end, avalanche=avalanche, launcher=launcher, every_step=pde)
    assert np.any(continuous[f"{avalanche}:delta_avalanche_mass_sum"] != 0), "expected avalanche transport"
    assert np.any(continuous["PBSM3D:sum_drift"] != 0), "expected blowing snow"
    assert np.ptp(continuous["Simple_Canopy:cum_net_snow"]) > 0, "expected canopy interception"
    if pde:
        areas = np.array(json.loads((directory / "test.param").read_text())["area"], dtype=float)
        check_pde_transfers(directory / name, areas)
    split, _ = scenario(directory, chm, name + "_first", mesh, "20200102T090000",
                        avalanche=avalanche, launcher=launcher)
    _, resumed = scenario(directory, chm, name + "_resumed", mesh, end,
                          avalanche=avalanche, launcher=launcher, restart=split)
    compare(continuous, resumed)
    return continuous


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
        for avalanche in ("snow_slide", "snowslide_pde"):
            coupled = coupled_event(directory, chm, mesh, end, avalanche)
            if args.partition:
                parallel = coupled_event(directory, chm, partition, end, avalanche, launcher)
                compare(coupled, parallel, repartitioned=True)
    print("Crocus CHM restart and decomposition comparisons passed")


if __name__ == "__main__":
    main()
