import sys
import os
from os import path
from zipfile import ZipFile

# Add the parent directory to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from baseutil import *

import math
import json
import hashlib
import pandas as pd
from pandas import DataFrame, read_csv, concat, set_option
from cobrakbase.core.kbasefba import FBAModel
import cobra
import optlang
import cobrakbase
from IPython.core.display import HTML
from modelseedpy import AnnotationOntology, MSPackageManager, MSMedia, MSModelUtil, MSBuilder, MSATPCorrection, MSGapfill, MSGrowthPhenotype, MSGrowthPhenotypes, ModelSEEDBiochem
from modelseedpy.core.msprobability import MSProbability
from modelseedpy.core.annotationontology import convert_to_search_role, split_role
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.core.msgenome import normalize_role
from modelseedpy.core.msensemble import MSEnsemble
from modelseedpy.community.mscommunity import MSCommunity
from modelseedpy.helpers import get_template
from cobra.flux_analysis import flux_variability_analysis


class PLPUtil(BaseUtil):
    def __init__(self):
        BaseUtil.__init__(self,"PLP")
        self.model = cobra.io.load_json_model("newdata/iML1515_PLP_Jun27_JTB.json")
        media_dict = {
            "glc__D":100,#22.5
            "nh4":100,
            "k":100,
            "pydxn":0.00002,
            'so4':100,
            'pi':100,
            'mn2': 100,
            'fe2': 100,
            'na1': 100,
            'mg2': 100,
            'o2': 20,
            'h2o': 100,
            'h': 100,
            'fe3': 100.0,
            'zn2': 100.0,
            'ca2': 100.0,
            'ni2': 100.0,
            'cu2': 100.0,
            'sel': 100.0,
            'cobalt2': 100.0,
            'mobd': 100.0,
            'cl': 100.0,
            'cbl1':100,
            'tungs': 100.0,
            'slnt': 100.0
        }
        #self.media = MSMedia.from_dict(media_dict)
        #self.pkgmgr = MSPackageManager.get_pkg_mgr(self.model)
        #self.pkgmgr.getpkg("KBaseMediaPkg").build_package(self.media)

    def compute_expression_fluxes(self,model,threshold,exp_data,growth,save_model=False):
        newmodel = cobra.io.json.from_json(cobra.io.json.to_json(model))
        #Reaction bounds
        rxn_bounds = {}
        for rxn in newmodel.reactions:
            rxn_bounds[rxn.id] = [rxn.lower_bound,rxn.upper_bound]
        newmodel.objective = "BIOMASS_Ec_iML1515_WT_75p37M"
        newmodel.objective.direction = "max"
        max_bio = newmodel.slim_optimize()
        killed = 0
        rev = 0
        nogene = 0
        above_threshold = 0
        output = {"threshold":threshold,"kos":[],"killed":0,"rev":0,"nogene":0,"above_threshold":0}
        if max_bio > growth:
            for rxn in newmodel.reactions:
                highest = None
                for gene in rxn.genes:
                    if gene.id in exp_data:
                        if not highest or exp_data[gene.id] > highest:
                            highest = exp_data[gene.id]
                if not highest:
                    output["nogene"] += 1
                elif highest < threshold:
                    rxn.lower_bound = 0
                    rxn.upper_bound = 0
                    new_bio = newmodel.slim_optimize()
                    if new_bio < growth:
                        output["rev"] += 1
                        rxn.lower_bound = rxn_bounds[rxn.id][0]
                        rxn.upper_bound = rxn_bounds[rxn.id][1]
                    else:
                        output["killed"] += 1
                        output["kos"].append(rxn.id)
                else:
                    output["above_threshold"] +=1
        newmodel.reactions.BIOMASS_Ec_iML1515_WT_75p37M.lower_bound = growth
        newmodel.reactions.BIOMASS_Ec_iML1515_WT_75p37M.upper_bound = growth
        newmodel.objective = "EX_glc__D_e"
        newmodel.objective.direction = "min"
        if save_model:
            output["model"] = cobra.io.json.to_json(newmodel)
        return output
    
    def process_reaction_class(self,sign,clshash,rxnid,nrxnid):
        if rxnid in clshash:
            for currcls in clshash[rxnid]:
                basecls = currcls[-1:]
                newval = 0
                if len(currcls) > 1:
                    newval = int(currcls[0:-1])
                if newval <= 0 and sign < 0:
                    newval += -1
                elif newval > 0 and sign > 0:
                    newval += 1
                newcls = str(newval)+basecls
                if nrxnid not in clshash:
                    clshash[nrxnid] = {}
                if newcls not in clshash[nrxnid]:
                    clshash[nrxnid][newcls] = 0
                clshash[nrxnid][newcls] += 1
    
util = PLPUtil() 