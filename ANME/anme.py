import sys
import os
from zipfile import ZipFile

# Add the parent directory to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)
folder_name = os.path.basename(script_dir)

import matplotlib.pyplot as plt
import math
import numpy as np
from baseutil import *

import hashlib
import pandas as pd
from pandas import DataFrame, read_csv, concat, set_option
from cobrakbase.core.kbasefba import FBAModel
import cobra
from cobra.flux_analysis import pfba
from modelseedpy import AnnotationOntology, MSPackageManager, MSMedia, MSModelUtil, MSBuilder, MSATPCorrection, MSGapfill, MSGrowthPhenotype, MSGrowthPhenotypes, ModelSEEDBiochem, MSExpression
from modelseedpy.core.msprobability import MSProbability
from modelseedpy.core.annotationontology import convert_to_search_role, split_role
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.core.msgenome import normalize_role
from modelseedpy.core.msensemble import MSEnsemble
from modelseedpy.community.mscommunity import MSCommunity
from modelseedpy.helpers import get_template

class NotebookUtil(BaseUtil):
    def __init__(self):
        BaseUtil.__init__(self,folder_name)
        
util = NotebookUtil() 