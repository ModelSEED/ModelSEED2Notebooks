import sys
import os

# Add the parent directory to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from baseutil import *

import pandas as pd
from pandas import DataFrame, read_csv, concat, set_option
from cobrakbase.core.kbasefba import FBAModel
from modelseedpy import AnnotationOntology, MSPackageManager, MSMedia, MSModelUtil, MSBuilder, MSATPCorrection, MSGapfill, MSGrowthPhenotype, MSGrowthPhenotypes, ModelSEEDBiochem
from modelseedpy.core.annotationontology import convert_to_search_role, split_role
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.core.msgenome import normalize_role
from modelseedpy.core.msensemble import MSEnsemble
from modelseedpy.helpers import get_template

class EnsembleModelingUtil(BaseUtil):
    def __init__(self):
        BaseUtil.__init__(self,"EnsembleModeling")

    def build_ensemble_model(self,params):
        self.initialize_call("build_ensemble_model",params,True)
        self.validate_args(params,["workspace","genome_ref"],{
            "ontology_events":[],
            "all_events":True,
            "run_gapfilling":True,
            "atp_safe":True,
            "forced_atp_list":[],
            "base_media":"157893/R2A_M.media",
            "gapfilling_media_list":[
                "KBaseMedia/Complete"
            ],
            "suffix":".EnsembleMdl",
            "core_template":"auto",
            "gs_template":"auto",
            "template_reactions_only":True,
            "output_core_models":True,
            "save_report_to_kbase":False,
            "return_data":True,
            "return_model_objects":True,
            "change_to_complete":False,
            "gapfilling_mode":"Independent",
            "max_gapfilling":0,
            "gapfilling_delta":0,
            "sample_count":100,
            "gs_template_ref":None,
            "core_template_ref":None
        })
        msoutput = self.msrecon.build_metabolic_models({
            "workspace":params["workspace"],
            "genome_refs":[params["genome_ref"]],
            "ontology_events":params["ontology_events"],
            "extend_model_with_ontology":True,
            "all_events":params["all_events"],
            "run_gapfilling":False,
            "atp_safe":params["atp_safe"],
            "forced_atp_list":params["forced_atp_list"],
            "base_media":params["base_media"],
            "gapfilling_media_list":[],
            "suffix":params["suffix"],
            "core_template":params["core_template"],
            "gs_template":params["gs_template"],
            "template_reactions_only":params["template_reactions_only"],
            "output_core_models":params["output_core_models"],
            "save_report_to_kbase":False,
            "return_data":params["return_data"],
            "save_models_to_kbase":False,
            "return_model_objects":True,
            "change_to_complete":params["change_to_complete"],
            "gapfilling_mode":params["gapfilling_mode"],
            "max_gapfilling":params["max_gapfilling"],
            "gapfilling_delta":params["gapfilling_delta"],
            "sample_count":params["sample_count"],
            "gs_template_ref":params["gs_template_ref"],
            "core_template_ref":params["core_template_ref"]
        }) 
        ensemble = MSEnsemble(msoutput["model_obj"])
        ensemblemdl = ensemble.sample_from_probabilities(from_reaction_probabilities=True,sample_count=1000)
        if params["run_gapfilling"]:
            self.gapfill_ensemble_model({
                "media_objs":params["gapfilling_media_objs"],#
                "gapfilling_media_list":params["gapfilling_media_list"],#
                "atp_safe":params["atp_safe"],#
                "workspace":params["workspace"],#
                "suffix":params["suffix"],#
                "default_objective":"bio1",#
                "output_data":msoutput["data"],#
                "forced_atp_list":params["forced_atp_list"],
                "templates":[self.msrecon.gs_template],
                "internal_call":True,
                "gapfilling_mode":params["gapfilling_mode"],
                "base_media":params["base_media"],
                "compound_list":params["compound_list"],
                "base_media_target_element":params["base_media_target_element"]
            })
        else:
            self.msrecon.save_model(ensemblemdl,params["workspace"])
        self.msrecon.build_report(msoutput["data"])
        if params["save_report_to_kbase"]:
            output = self.msrecon.save_report_to_kbase()
        if params["return_data"]:
            output["data"] = msoutput["data"]
        if params["return_model_objects"]:
            output["model_obj"] = ensemblemdl
        return output
        
    def gapfill_ensemble_model(self,params):
        self.initialize_call("gapfill_ensemble_model",params,True)
        self.validate_args(params,["workspace","model_ref"],{
            "change_to_complete":False,
            "atp_safe":True,
            "forced_atp_list":[],
            "suffix":".gf",
            "media_list":["KBaseMedia/Complete"],
            "limit_medias":[],
            "limit_thresholds":[],
            "is_max_limits":[],
            "templates":None,
            "source_models":None,
            "target_reaction":"bio1",
            "minimum_objective":0.01,
            "reaction_exlusion_list":[],
            "default_objective":"bio1",
            "save_report_to_kbase":False,
            "gapfilling_mode":"Cumulative",
            "base_media":None,
            "compound_list":[],
            "base_media_target_element":"C",
            "output_data":None,
            "model_obj":None,
            "media_objs":None,
            "core_template_ref":None
        })
        base_comments = []
        default_media = "KBaseMedia/AuxoMedia"
        if params["change_to_complete"]:
            base_comments.append("Changing default to complete.")
            default_media = "KBaseMedia/Complete"
        current_output = {"Model":None,"Genome":None,"Genes":None,"Class":None,
            "Model genes":None,"Reactions":None,"Core GF":None,"GS GF":None,"Growth":None,"Comments":[]}
        if params["output_data"]:
            current_output = params["output_data"]
        if not params["model_obj"]:
            params["model_obj"] = self.get_model(params["model_ref"])
        mdlutl = params["model_obj"]
        #Processing media
        if not params["media_objs"]:
            params["media_objs"] = self.process_media_list(params["media_list"],default_media,params["workspace"])
        #Processing compound list
        if params["compound_list"]:
            if not params["base_media"]:
                base_comments.append("No base media provided. Ignoring compound list.")
            else:
                for cpd in params["compound_list"]:
                    newmedia = MSMedia.from_dict({cpd:100})
                    newmedia.merge(params["base_media"])
                    params["media_objs"].append(newmedia)
        #Compiling additional tests
        additional_tests = []
        for i,limit_media in enumerate(params["limit_medias"]):
            additional_tests.append({
                "objective":params["limit_objectives"][i],
                "media":self.get_media(limit_media,None),
                "is_max_threshold":params["is_max_limits"][i],
                "threshold":params["limit_thresholds"][i]
            })
        #Getting core template
        if params["core_template_ref"]:
            self.core_template = self.get_template(params["core_template_ref"],None)
        else:
            self.core_template = self.get_template(self.templates["core"],None)
        #Iterating over each model and running gapfilling
        current_output["Comments"] = base_comments
        current_output["Model"] = mdlutl.wsid+params["suffix"]+'<br><a href="'+mdlutl.wsid+params["suffix"]+'-recon.html" target="_blank">(see reconstruction report)</a><br><a href="'+mdlutl.wsid+params["suffix"]+'-full.html" target="_blank">(see full view)</a>'
        if params["output_data"] and mdlutl in params["output_data"]:
            current_output = params["output_data"][mdlutl]
        #Setting the objective
        
        #Computing tests for ATP safe gapfilling
        if params["atp_safe"]:
            tests = mdlutl.get_atp_tests(core_template=self.core_template,atp_media_filename=self.module_dir+"/data/atp_medias.tsv",recompute=False)
            additional_tests.extend(tests)
        #Creating gapfilling object and configuring solver
        #mdlutl.model.solver = config["solver"]
        if not params["templates"]:
            params["templates"] = [self.get_template(mdlutl.model.template_ref)]
        msgapfill = MSGapfill(
            mdlutl,
            params["templates"],
            params["source_models"],
            additional_tests,
            blacklist=params["reaction_exlusion_list"],
            default_target=params["target_reaction"],
            minimum_obj=params["minimum_objective"],
            base_media=params["base_media"],
            base_media_target_element=params["base_media_target_element"]
        )
        #Running gapfilling in all conditions
        mdlutl.gfutl.cumulative_gapfilling = []
        growth_array = []
        solutions = msgapfill.run_multi_gapfill(
            params["media_objs"],
            target=params["target_reaction"],
            default_minimum_objective=params["minimum_objective"],
            binary_check=False,
            prefilter=True,
            check_for_growth=True,
            gapfilling_mode=params["gapfilling_mode"],
            run_sensitivity_analysis=True,
            integrate_solutions=True
        )
        for media in params["media_objs"]:
            if media.id in solutions and "growth" in solutions[media]:
                growth_array.append(media.id+":"+str(solutions[media]["growth"]))
        current_output["Growth"] = "<br>".join(growth_array)
        current_output["GS GF"] = len(mdlutl.gfutl.cumulative_gapfilling)
        current_output["Reactions"] = mdlutl.nonexchange_reaction_count()
        current_output["Model genes"] = len(mdlutl.model.genes)
        #Saving completely gapfilled model
        self.save_model(mdlutl,params["workspace"],None,params["suffix"])
        output = {}
        if not params["internal_call"]:
            self.msrecon.build_report(current_output)
            if params["save_report_to_kbase"]:
                output = self.msrecon.save_report_to_kbase()
            if params["return_data"]:
                output["data"] = current_output
            if params["return_model_objects"]:
                output["model_obj"] = mdlutl
        return output
    
    def run_ensemble_fba(self,params):
        self.msrecon.initialize_call("gapfill_ensemble_model",params,True)
        self.msrecon.validate_args(params,["workspace","model_ref","fba_output_id"],{
            "media":"KBaseMedia/Complete",
            "maximize":True,
            "target_reactions":["bio1",1],
            "thermodynamic_constraints":False,
            "fva":True,
            "pfba":True,
            "simulate_all_ko":False,
            "feature_ko_list":[],
            "reaction_ko_list":[],
            "media_supplement_list":[],
            "max_c_uptake":None,
            "max_n_uptake":None,
            "max_p_uptake":None,
            "max_s_uptake":None,
            "max_o_uptake":None
        })
        mdlutl = self.msrecon.get_model(params["model_ref"])
        ensemble = MSEnsemble(mdlutl)
        if not ensemble.is_ensemble():
            logger.critical("Input model is not an ensemble model!")
        media = self.msrecon.get_media(params["media"])
        fbaobj = ensemble.run_fba(media,params["target_reactions"],params["maximize"],params["gene_ko_list"],params["reaction_ko_list"],params["pfba"],params["fva"])
        self.msrecon.save_fba(fbaobj,params["workspace"],params["fba_output_id"])
        output = {}
        if not params["internal_call"]:
            self.build_dataframe_report(result_table,params["model_objs"])
            if params["save_report_to_kbase"]:
                output = self.save_report_to_kbase()
            if params["return_data"]:
                output["data"] = result_table.to_json()
            if params["return_model_objects"]:
                output["model_objs"] = params["model_objs"]
        return output

    def extract_models_from_ensemble(self,params):
        self.msrecon.initialize_call("gapfill_ensemble_model",params,True)
        self.msrecon.validate_args(params,["workspace","model_ref"],{
            "model_list":None
        })
        mdlutl = self.msrecon.get_model(params["model_ref"])
        ensemble = MSEnsemble(mdlutl)
        models = ensemble.unpack_models(model_list=params["model_list"])
        for (i,mdl) in enumerate(models):
            self.save_model(mdl,params["workspace"],None,"."+str(i))
        if not params["internal_call"]:
            self.build_dataframe_report(result_table,params["model_objs"])
            if params["save_report_to_kbase"]:
                output = self.save_report_to_kbase()
            if params["return_data"]:
                output["data"] = result_table.to_json()
            if params["return_model_objects"]:
                output["model_objs"] = params["model_objs"]
        return output
    
util = EnsembleModelingUtil()