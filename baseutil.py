import platform
import sys
import sys
import json
from json import dump
import os
import re
from os.path import exists
from pathlib import Path
import logging
import shutil

print("python version " + platform.python_version())

from configparser import ConfigParser

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

config = ConfigParser()
if not exists(str(Path.home()) + '/.kbase/config'):    
    if exists("/scratch/shared/code/sharedconfig.cfg"):
        shutil.copyfile("/scratch/shared/code/sharedconfig.cfg",str(Path.home()) + '/.kbase/config')
    else:
        print("You much create a config file in ~/.kbase/config before running this notebook. See instructions: https://docs.google.com/document/d/1fQ6iS_uaaZKbjWtw1MgzqilklttIibNO9XIIJWgxWKo/edit")
        sys.exit(1)
config.read(str(Path.home()) + '/.kbase/config')
paths = config.get("DevEnv","syspaths").split(";")
codebase = config.get("DevEnv","codebase",fallback="")
for i,filepath in enumerate(paths):
    if filepath[0:1] != "/":
        paths[i] = codebase+"/"+filepath
sys.path = paths + sys.path

from chenry_utility_module.kbdevutils import KBDevUtils

class BaseUtil(KBDevUtils):
    def __init__(self,name):
        KBDevUtils.__init__(self,name,output_root=os.path.dirname(os.path.realpath(__file__)))
        self.kbdevutil = self
        self.annoapi = self.anno_client(native_python_api=True)
        self.obs_ec = None
        self.msseedrecon()
        self.media={"auxo":self.msrecon.get_media("94026/Auxotrophy_media")}
        self.msrecon.util = self

    def simulate_biolog_phenotypes(
            self,
            model,
            objective="bio1",
            phenotypeset="157564/Full_Biolog_Dataset",
            base_media="auxo",
            base_uptake=100,
            base_excretion=1000,
            global_atom_limits={"C":1},
            add_missing_exchanges=True,
            save_fluxes=False,
            save_reaction_list=False,
            ignore_experimental_data=False,
            flux_coefficients=None):
        if base_media and type(base_media) == str:
            base_media = self.media[base_media]
        if type(phenotypeset) == str:
            phenotypeset = self.msrecon.get_phenotypeset(phenotypeset, base_media=base_media, base_uptake=base_uptake, base_excretion=base_excretion, global_atom_limits=global_atom_limits)
        model.model.objective = objective
        output = phenotypeset.simulate_phenotypes(model,multiplier=3,add_missing_exchanges=add_missing_exchanges,save_fluxes=save_fluxes,save_reaction_list=save_reaction_list,ignore_experimental_data=ignore_experimental_data,flux_coefficients=flux_coefficients)
        output["phenotypeset"] = phenotypeset
        return output