# Copyright (C) 2023 Luigi Pertoldi <gipert@pm.me>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ruff: noqa: F821, T201

from __future__ import annotations

import argparse
import json

import ROOT
import uproot
from legendmeta import LegendMetadata


def process_mage_id(mage_id):
    m_id = str(mage_id)
    is_ged = bool(int(m_id[0]))
    if not is_ged:
        return False

    string = int(m_id[3:5])
    pos = int(m_id[5:7])

    for key, value in chmap.items():
        if isinstance(value, dict) and "location" in value:
            location = value["location"]
            usable = value["analysis"]["usability"] == "on"
            if (
                location.get("string") == string
                and location.get("position") == pos
                and usable
            ):
                return {"name": key, "ch": value["daq"]["rawid"], "mage_id": mage_id}

    return False


parser = argparse.ArgumentParser(
    prog="build_pdf", description="build LEGEND pdf files from evt tier files"
)
parser.add_argument("--config", "-c", help="configuration file")
parser.add_argument("--output", "-o", help="output file name")
parser.add_argument("--metadata", "-m", help="path to legend-metadata")
parser.add_argument("input_files", nargs="+", help="evt tier files")

args = parser.parse_args()

if not isinstance(args.input_files, list):
    args.input_files = [args.input_files]

# with Path(args.config).open() as f:
with open(args.config, "r") as f:
    rconfig = json.load(f)

meta = LegendMetadata(args.metadata)
chmap = meta.channelmap(rconfig["timestamp"])

for file_name in args.input_files:

    # load in all the data into a pandas dataframe
    with uproot.open(f"{file_name}:simTree") as pytree:
        df_data = pytree.arrays(['energy', 'npe_tot', 'mage_id'], library="pd")
    df = df_data[df_data["energy"] > rconfig['energy_threshold']]
    subentry_counts = df.index.get_level_values('entry').value_counts()
    n_primaries = len(df_data)

    uniq_mage_ids = df_data.mage_id.unique()
    mage_names = {
        mage_id: process_mage_id(mage_id)
        for mage_id in uniq_mage_ids
        if process_mage_id(mage_id)
    }

    out_file_name = file_name.split("/")[-1].split(".")[0].replace("tier_evt", "tier_pdf") + ".root"
    out_file = uproot.recreate(args.output + out_file_name)
    out_file["number_of_primaries"] = str(n_primaries)

    for _cut_name, _cut_string in rconfig["cuts"].items():
        dir = out_file.mkdir(_cut_name)
        # Define the grouped hists to be filled and store them in memory
        cut_hists = {
            f"{_type}": ROOT.TH1F( 
                f"{_type}", f"All {_type} energy deposit",
                rconfig['hist']['nbins'], rconfig['hist']['emin'], rconfig['hist']['emax']
            ) for _type in ["bege", "coax", "icpc", "ppc"]
        } 
        cut_hists['all'] = ROOT.TH1F( 
            f"all", f"All energy deposit",
            rconfig['hist']['nbins'], rconfig['hist']['emin'], rconfig['hist']['emax']
        ) 

        # We want to cut on multiplicity for all detectors >25keV
        # Include them in the dataset then apply cuts - then filter them out
        # Don;t store AC detectors
        exec(_cut_string)

        for _mage_id, _mage_names in mage_names.items():

            if chmap[_mage_names["name"]]["analysis"]["usability"] != 'on': continue
            hist = ROOT.TH1F(
                f"ch{_mage_names['ch']}", f"{_mage_names['name']} energy deposition",
                rconfig['hist']['nbins'], rconfig['hist']['emin'], rconfig['hist']['emax']
            )
            df_channel = df_cut[df_cut.mage_id == _mage_id]

            for energy in df_channel.energy.values:
                hist.Fill(energy * 1000) # energy in keV

            dir[hist.GetName()] = hist
            cut_hists[chmap[_mage_names["name"]]["type"]].Add(hist)
            cut_hists['all'].Add(hist)
            del hist
        
        for _type, _type_hist in cut_hists.items():
            dir[_type_hist.GetName()] = _type_hist
            del _type_hist
        del cut_hists

    out_file.close()
