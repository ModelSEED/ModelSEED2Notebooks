import sys
import os
from os import path
from zipfile import ZipFile

# Add the parent directory to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from baseutil import *

import hashlib
import pandas as pd
from pandas import DataFrame, read_csv, concat, set_option
from cobrakbase.core.kbasefba import FBAModel
from cobra.io import write_sbml_model, read_sbml_model
from modelseedpy import AnnotationOntology, MSPackageManager, MSMedia, MSModelUtil, MSBuilder, MSATPCorrection, MSGapfill, MSGrowthPhenotype, MSGrowthPhenotypes, ModelSEEDBiochem
from modelseedpy.core.msprobability import MSProbability
from modelseedpy.core.annotationontology import convert_to_search_role, split_role
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.core.msgenome import normalize_role
from modelseedpy.core.msensemble import MSEnsemble
from modelseedpy.community.mscommunity import MSCommunity
from modelseedpy.helpers import get_template

class ModelAnalysisUtil(BaseUtil):
    def __init__(self):
        BaseUtil.__init__(self,"Analysis")
        self.com_media = self.msrecon.get_media("KBaseMedia/Complete")
        self.gmm_media = self.msrecon.get_media("KBaseMedia/Carbon-D-Glucose")
        self.auxo_media = self.msrecon.get_media("94026/Auxotrophy_media")
    
util = ModelAnalysisUtil() 