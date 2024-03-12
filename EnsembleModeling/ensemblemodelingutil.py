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
        BaseUtil.__init__(self)
        self.get_kbdevutil("EnsembleModeling")
        self.get_msrecon("EnsembleModeling")

    def build_ensemble_model(self,params):
        self.msrecon.initialize_call("build_ensemble_model",params,True)
        self.msrecon.validate_args(params,["workspace","genome_ref"],{
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
        default_media = "KBaseMedia/AuxoMedia"
        if params["change_to_complete"]:
            default_media = "KBaseMedia/Complete"
        #Processing media
        params["gapfilling_media_objs"] = self.process_media_list(params["gapfilling_media_list"],default_media,params["workspace"])
        #Preloading core and preselected template
        self.gs_template = None
        if params["gs_template_ref"]:
            self.gs_template = self.get_template(params["gs_template_ref"],None)
        if params["core_template_ref"]:
            self.core_template = self.get_template(params["core_template_ref"],None)
        else:
            self.core_template = self.get_template(self.templates["core"],None)  
        #Initializing classifier
        genome_classifier = self.get_classifier()
        #Initializing output data tables
        current_output = {"Model":None,"Genome":None,"Genes":None,"Class":None,
                          "Model genes":None,"Reactions":None,
                          "Core GF":None,"GS GF":None,"Growth":None,"Comments":[]}
        #Retrieving genomes and building models one by one
        template_type = params["gs_template"]
        genome = self.get_msgenome(params["genome_ref"])
        #Initializing output row
        current_output["Comments"] = []
        gid = genome.id
        current_output["Model"] = gid+params["suffix"]+'<br><a href="'+gid+params["suffix"]+'-recon.html" target="_blank">(see reconstruction report)</a><br><a href="'+gid+params["suffix"]+'-full.html" target="_blank">(see full view)</a>'
        current_output["Genome"] = genome.info.metadata["Name"]
        current_output["Genes"] = genome.info.metadata["Number of Protein Encoding Genes"]
        #Pulling annotation priority
        current_output["Comments"].append("Other annotation priorities not supported by this app yet. Using RAST.")
        if template_type == "auto":
            current_output["Class"] = genome_classifier.classify(genome)
            if current_output["Class"] == "P":
                current_output["Class"] = "Gram Positive"
                template_type = "gp"
            elif current_output["Class"] == "N" or current_output["Class"] == "--":
                current_output["Class"] = "Gram Negative"
                template_type = "gn"
            elif current_output["Class"] == "A":
                current_output["Class"] = "Archaea"
                template_type = "ar"
            elif current_output["Class"] == "C":
                current_output["Class"] = "Cyanobacteria"
                template_type = "cyano"
                current_output["Comments"].append("Cyanobacteria not yet supported. Skipping genome.")
            else:
                current_output["Comments"].append("Unrecognized genome class "+current_output["Class"]+". Skipping genome.")
        if not self.gs_template:
            self.gs_template = self.get_template(self.templates[template_type],None)
        #Building model            
        base_model = FBAModel({'id':gid+params["suffix"], 'name':genome.scientific_name})
        builder = MSBuilder(genome, self.gs_template)
        annoapi = self.msrecon.anno_client(native_python_api=True)
        annoont = AnnotationOntology.from_kbase_data(annoapi.get_annotation_ontology_events({
            "input_ref" : params["genome_ref"]
        }),params["genome_ref"],self.module_dir+"/data/")
        gene_term_hash = annoont.get_gene_term_hash(None,params["ontology_events"],True,False)
        
        for gene in gene_term_hash:
            for term in gene_term_hash[gene]:
                if term.ontology.id == "SSO":
                    name = annoont.get_term_name(term)
                    f_norm = normalize_role(name)
                    if f_norm not in builder.search_name_to_genes:
                        builder.search_name_to_genes[f_norm] = set()
                        builder.search_name_to_original[f_norm] = set()
                    builder.search_name_to_original[f_norm].add(name)
                    builder.search_name_to_genes[f_norm].add(gene.id)
        mdl = builder.build(base_model, '0', False, False)
        mdl.genome = genome
        mdl.template = self.gs_template
        mdl.core_template_ref = str(self.core_template.info)
        mdl.genome_ref = str(genome.info)
        mdl.template_ref = str(self.gs_template.info)
        current_output["Core GF"] = "NA" 
        mdlutl = MSModelUtil.get(mdl)
        if params["atp_safe"]:
            atpcorrection = MSATPCorrection(mdlutl,self.core_template,params["atp_medias"],load_default_medias=params["load_default_medias"],max_gapfilling=params["max_gapfilling"],gapfilling_delta=params["gapfilling_delta"],forced_media=params["forced_atp_list"],default_media_path=self.module_dir+"/data/atp_medias.tsv")
            tests = atpcorrection.run_atp_correction()
            current_output["Core GF"] = len(atpcorrection.cumulative_core_gapfilling)
        #Setting the model ID so the model is saved with the correct name in KBase
        mdlutl.get_attributes()["class"] = current_output["Class"]
        mdlutl.wsid = gid+params["suffix"]
        #Running gapfilling
        current_output["GS GF"] = "NA"
        if params["run_gapfilling"]:
            self.gapfill_ensemble_model({
                "media_objs":params["gapfilling_media_objs"],#
                "model_obj":mdlutl,#
                "atp_safe":params["atp_safe"],#
                "workspace":params["workspace"],#
                "suffix":"",#
                "default_objective":"bio1",#
                "output_data":current_output,#
                "forced_atp_list":params["forced_atp_list"],
                "templates":[self.gs_template],
                "internal_call":True,
                "gapfilling_mode":params["gapfilling_mode"],
                "base_media":params["base_media"],
                "compound_list":params["compound_list"],
                "base_media_target_element":params["base_media_target_element"]
            })
        else:
            self.msrecon.save_model(mdlutl,params["workspace"],None)
            mdlutl.model.objective = "bio1"
            mdlutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(None)
            current_output["Growth"] = "Complete:"+str(mdlutl.model.slim_optimize())
        current_output["Reactions"] = mdlutl.nonexchange_reaction_count()
        current_output["Model genes"] = len(mdlutl.model.genes)
        output = {}
        self.msrecon.build_report(current_output)
        if params["save_report_to_kbase"]:
            output = self.msrecon.save_report_to_kbase()
        if params["return_data"]:
            output["data"] = current_output
        if params["return_model_objects"]:
            output["model_obj"] = mdlutl
        return output
        
    def gapfill_ensemble_model(self,params):
        self.msrecon.initialize_call("gapfill_ensemble_model",params,True)
        self.msrecon.validate_args(params,["workspace","model"],{
            "atp_safe":True,
            "forced_atp_list":[],
            "suffix":".gf",
            "media_list":["KBaseMedia/Complete"],
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
            "base_media_target_element":"C"
        })
        base_comments = []
        if params["change_to_complete"]:
            base_comments.append("Changing default to complete.")
            default_media = "KBaseMedia/Complete"
        current_output = {"Model":None,"Genome":None,"Genes":None,"Class":None,
            "Model genes":None,"Reactions":None,"Core GF":None,"GS GF":None,"Growth":None,"Comments":[]}
        mdlutl = self.get_model(params["model"])
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
        if not params["internal_call"]:
            result_table = result_table.append(current_output, ignore_index = True)
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