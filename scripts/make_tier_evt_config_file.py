# Script to generate the json file that contains the FCCD, energy smearing parameters and the usability mask.
# ---
# Contribution: Aparajita Mazumdar (LANL), Elisabetta Bossio (TUM) [Aug 2, 2023]
# Updated for compatibility with Snakemake.
# ---
#
# Notes about this version:
# This version uses the new par_hit_results.json files, processed after the issues with the calibration were discovered. These files were grabbed from the LNGS server from the path `/data2/public/prodenv/prod-blind/ref/v01.06/generated/par/hit/cal`. These have stored in the local folder new-par-hit-results. This path will have to be changed once the files have been transferred to the project directory on NERSC. A gaussian smearing has been considered for the peaks, i.e., no tails have been implemented. The reason for this is that par_hit_results.json provides FWHM fits while mage-post-proc needs sigma and tau fits separately. It would be good if the par_hit_results.json started including the errors on sigma and tau fits as well, to incorporate tails.
# 
# 1. Uses the keylists to search for the relevant par_hit_results.json file.
# 2. par_hit_results.json file is parsed. The channel list, fitted peak list and eres parameters are grabbed.
# 3. Some parameters are being pulled out of the legend-metadata, i.e., detector name, usability mask, string and position. The string number and the detector position are later used to generate the mage-id, which is used as the index in the mage-pars.json. Conditions for the golden and silver dataset have been implemented. Golden dataset condition: Detectors marked as "no_psd" have been set to "ac" for simulations. Silver dataset condition: Detectors marked as "no_psd" have been set to "on" for simulations. These points are where usability masks for legend-dataflow-config and legend-simflow-config diverge.
# 4. Added Elisabetta's code to include the full charge collection depth parameter. This relies on reading the excel file 'fccd-reviewed.xlsx'.
# 5. "off" detectors have null energy smearing.
#
# Notes about legend-metadata:
# The config folder in legend-metadata is pulled manually from dataflow-config (v1.06) [updated on Aug 9, 2023].
# ---
#
# Running the script: Keep 'fccd-reviewed.xlsx' in the same directory. Open an instance in legend-base and run the script
# $ cenv legend-base
# $ python3 make_tier_evt_config_file.py -d golden -w l200-p04-r002-phy
# ---
#
################################################################################################

# Importing the packages
import json 
import numpy as np 
import pandas as pd 
import math 
from legendmeta import LegendMetadata 
from datetime import datetime 
import os 
import glob 
import argparse
import sys
from pathlib import Path

# Brief description accessed through --help.
parser = argparse.ArgumentParser(description="""Generates config files for mage-post-proc. Code depends on legend-metadata, and looks for input excel sheet "fccd-reviewed.xlsx" in the same folder.""")
parser.add_argument("-d", "--dataset", type=str, default='golden', help="golden or silver dataset. Default:golden")
parser.add_argument("-w", "--wildcard", type=str, help="Passing on the snakemake wildcard.")
args = parser.parse_args()

# If an invalid --dataset flag is set, abort code.
if (args.dataset == "golden"):
    print("Detectors marked as no_psd have been set to ac for simulations.")
elif (args.dataset == "silver"):
    print("Detectors marked as no_psd have been set to on for simulations.")
else:
    print("Error: Invalid --dataset flag. Abort code.")
    sys.exit(1)

# Define a null variable.
null=None

# Path to legend-metadata.
# If pointing to a user specific path, remember to define it inside the paranthesis
lmeta = LegendMetadata('/tmp/legend-metadata-mazumdar/') # LNGS

# Reading in the 'fccd-reviewed.xlsx' excel sheet with the nplus-fccd parameters.
# Operation performed outside the loop because for now, since we have only constant, i.e., time independent values.
# This might change if time-dependent dead layers are considered.
xls = pd.ExcelFile('fccd-reviewed.xlsx', engine='openpyxl')

# Reading in the various sheets from the excel, into dataframes
df1 = pd.read_excel(xls, 'l140-det')
df2 = pd.read_excel(xls, 'gerda-det')
df3 = pd.read_excel(xls, 'mirion-icpc')
df4 = pd.read_excel(xls, 'ortec-icpc')
df5 = pd.read_excel(xls, 'majorana-ppc')

# List of L140 detectors
l140_det = df1['name'].to_numpy()
# Drop unnecessary columns
mirion_icpc = df3.drop(columns=['fccd-mm (Ba)','fccd-mm (Am)'])
# Drop detectors that were not used in the L140 configuration.
mirion_icpc = mirion_icpc.loc[mirion_icpc['det name'].isin(l140_det)].reset_index(drop=True)
gerda_det = df2.loc[df2['det name'].isin(l140_det)].reset_index(drop=True)
ortec_icpc = df4.drop(columns=['fccd-mm (manufacturer)'])
ortec_icpc = ortec_icpc.loc[ortec_icpc['det name'].isin(l140_det)].reset_index(drop=True)
majorana_ppc = df5.loc[df5['det name'].isin(l140_det)].reset_index(drop=True)

# Make one large frame
frames = [mirion_icpc,gerda_det,ortec_icpc,majorana_ppc]
all_det = pd.concat(frames)
all_det = all_det.set_index('det name')

# Read in from the snakemake wildcard and drop phy for the wildcard search
try:
    wcard = args.wildcard[:-4]
except NoInput:
    wcard = snakemake.wildcards.runid[:-4]
except WildCardError:
    print("Error: Invalid wildcard. Abort code.")
    sys.exit(2)

# Find the relavant par_hit_results.json file
# This is a temporary path on my NERSC account, which should be changed once the files are transferred to NERSC
parfile=glob.glob(f'/pscratch/sd/j/jita/jita/rushcalV3-mpp/pargen-parser/new-par-hit-results/cal/*/*/{wcard}*par_hit_results.json')[0]

#LNGS
#parfile=glob.glob(f'/data2/public/prodenv/prod-blind/ref/v01.06/generated/par/hit/cal/*/*/{wcard}*par_hit_results.json')[0]
print(parfile)

# Extract the datetimestamp and get the channel map as on that timestamp
start_ind = parfile.index("-cal-") + len("-cal-")
stop_ind = parfile.index("-par_")
datetimestamp=parfile[start_ind:stop_ind]
chmap = lmeta.hardware.configuration.channelmaps.on(datetimestamp)

# Parse the par_hit_results.json
fo = open(parfile) # io object to open the file
data=json.load(fo) # load json file
channels=list(data.keys()) # Grab processed channels
fname = parfile[parfile.index("l200-"):parfile.index("-cal-")]
outstr={}
# Detectors in L140
# Could have used l140_det from the FCCD section directly, but preferred to use metadata
detlist=[]
for det, val in chmap.items():
    if val.system == "geds":
        detlist.append(det)
for gedname in detlist:
    ged = lmeta.channelmap(datetimestamp)[gedname]
    mask = ged.analysis.usability
    if (args.dataset == "golden"):
        # set "no_psd" detectors as "ac" for simulations
        if (mask == "no_psd"):
            mask = "ac"
    else:
        # set "no_psd" detectors as "on" for simulations
        if (mask == "no_psd"):
            mask = "on"
    ch='ch'+str(ged["daq"]["rawid"])
    fccd=all_det.loc[gedname]['fccd-mm']
    # generate the mage-id 10C0P0D, which will be the index
    cryostat=1
    string=chmap[gedname].location.string
    posn=chmap[gedname].location.position
    mageid=1000000+(cryostat*10000)+(string*100)+posn
    if ch in channels:
        eres=data[ch]['ecal']['cuspEmax_ctc_cal']['eres_pars']# grab eres_pars
        outstr[mageid] = {"name":gedname, "nplus-fccd-mm":fccd, "energy":{"sig0":math.sqrt(eres[0])/2.355, "sig1":math.sqrt(eres[1])/2.355, "sig2":0}, "usability":mask}
    else:
        outstr[mageid] = {"name":gedname, "nplus-fccd-mm":fccd, "energy":{"sig0":null, "sig1":null, "sig2":null}, "usability":mask}

# Writing out the JSON file
json_obj = json.dumps(outstr, indent=4)
jsonfile = f'{wcard}-phy-build_evt.json'
#jsonfile = (Path(snakemake.input[0])/ "simprod"/ "config"/ "tier"/ "evt"/ snakemake.config["experiment"]/ f"{wcard}-phy-build_evt.json")
with open(jsonfile, "w") as outfile:
    outfile.write(json_obj)
