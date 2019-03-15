#                                                            _
# z2labelmap ds app
#
# (c) 2016 Fetal-Neonatal Neuroimaging & Developmental Science Center
#                   Boston Children's Hospital
#
#              http://childrenshospital.org/FNNDSC/
#                        dev@babyMRI.org
#

import os
import  os
from    os          import listdir, sep
from    os.path     import abspath, basename, isdir
import  shutil
import  pudb
import  sys
import  time
import  random

# import the Chris app superclass
from chrisapp.base import ChrisApp


class Z2labelmap(ChrisApp):
    """
    Convert a file of per-structure z-scores to a FreeSurfer labelmap..
    """
    AUTHORS         = 'FNNDSC (dev@babyMRI.org)'
    SELFPATH        = os.path.dirname(os.path.abspath(__file__))
    SELFEXEC        = os.path.basename(__file__)
    EXECSHELL       = 'python3'
    TITLE           = 'z-score to FreeSurfer label map'
    CATEGORY        = 'FreeSurfer'
    TYPE            = 'ds'
    DESCRIPTION     = 'Convert a file of per-structure z-scores to a FreeSurfer labelmap.'
    DOCUMENTATION   = 'http://wiki'
    VERSION         = '0.1'
    ICON            = '' # url of an icon image
    LICENSE         = 'Opensource (MIT)'
    MAX_NUMBER_OF_WORKERS = 1  # Override with integer value
    MIN_NUMBER_OF_WORKERS = 1  # Override with integer value
    MAX_CPU_LIMIT         = '' # Override with millicore value as string, e.g. '2000m'
    MIN_CPU_LIMIT         = '' # Override with millicore value as string, e.g. '2000m'
    MAX_MEMORY_LIMIT      = '' # Override with string, e.g. '1Gi', '2000Mi'
    MIN_MEMORY_LIMIT      = '' # Override with string, e.g. '1Gi', '2000Mi'
    MIN_GPU_LIMIT         = 0  # Override with the minimum number of GPUs, as an integer, for your plugin
    MAX_GPU_LIMIT         = 0  # Override with the maximum number of GPUs, as an integer, for your plugin

    # Fill out this with key-value output descriptive info (such as an output file path
    # relative to the output dir) that you want to save to the output meta file when
    # called with the --saveoutputmeta flag
    OUTPUT_META_DICT = {}
 
    def structList_define(self):
        """
        The list of structures in the a2009s cortical parcellation
        """
        self.l_a2009s = [
            'G_and_S_frontomargin',
            'G_and_S_occipital_inf',
            'G_and_S_paracentral',
            'G_and_S_subcentral',
            'G_and_S_transv_frontopol',
            'G_and_S_cingul-Ant',
            'G_and_S_cingul-Mid-Ant',
            'G_and_S_cingul-Mid-Post',
            'G_cingul-Post-dorsal',
            'G_cingul-Post-ventral',
            'G_cuneus',
            'G_front_inf-Opercular',
            'G_front_inf-Orbital',
            'G_front_inf-Triangul',
            'G_front_middle',
            'G_front_sup',
            'G_Ins_lg_and_S_cent_ins',
            'G_insular_short',
            'G_occipital_middle',
            'G_occipital_sup',
            'G_oc-temp_lat-fusifor',
            'G_oc-temp_med-Lingual',
            'G_oc-temp_med-Parahip',
            'G_orbital',
            'G_pariet_inf-Angular',
            'G_pariet_inf-Supramar',
            'G_parietal_sup',
            'G_postcentral',
            'G_precentral',
            'G_precuneus',
            'G_rectus',
            'G_subcallosal',
            'G_temp_sup-G_T_transv',
            'G_temp_sup-Lateral',
            'G_temp_sup-Plan_polar',
            'G_temp_sup-Plan_tempo',
            'G_temporal_inf',
            'G_temporal_middle',
            'Lat_Fis-ant-Horizont',
            'Lat_Fis-ant-Vertical',
            'Lat_Fis-post',
            'Pole_occipital',
            'Pole_temporal',
            'S_calcarine',
            'S_central',
            'S_cingul-Marginalis',
            'S_circular_insula_ant',
            'S_circular_insula_inf',
            'S_circular_insula_sup',
            'S_collat_transv_ant',
            'S_collat_transv_post',
            'S_front_inf',
            'S_front_middle',
            'S_front_sup',
            'S_interm_prim-Jensen',
            'S_intrapariet_and_P_trans',
            'S_oc_middle_and_Lunatus',
            'S_oc_sup_and_transversal',
            'S_occipital_ant',
            'S_oc-temp_lat',
            'S_oc-temp_med_and_Lingual',
            'S_orbital_lateral',
            'S_orbital_med-olfact',
            'S_orbital-H_Shaped',
            'S_parieto_occipital',
            'S_pericallosal',
            'S_postcentral',
            'S_precentral-inf-part',
            'S_precentral-sup-part',
            'S_suborbital',
            'S_subparietal',
            'S_temporal_inf',
            'S_temporal_sup',
            'S_temporal_transverse'
        ]

    def randomZscoreFile_generate(self):
        """
        Generate a "random" z-score file, based on the range given in the
        --random <range> command line flag.
        """
        sampleLen   = len(self.l_a2009s)

    def define_parameters(self):
        """
        Define the CLI arguments accepted by this plugin app.
        """

        self.add_argument("-v", "--verbosity",
                            help        = "verbosity level for app",
                            type        = str,
                            dest        = 'verbosity',
                            optional    = True,
                            default     = "0")
        self.add_argument("-p", "--posRange",
                            help        = "positive range for max deviation",
                            type        = float,
                            dest        = 'f_posRange',
                            optional    = True,
                            default     = 1.0)
        self.add_argument("-n", "--negRange",
                            help        = "negative range for max deviation",
                            type        = float,
                            dest        = 'f_negRange',
                            optional    = True,
                            default     = -1.0)
        self.add_argument('--random',
                            help        = 'if specified, generate a z-score file',
                            type        = bool,
                            dest        = 'b_random',
                            action      = 'store_true',
                            optional    = True,
                            default     = False)
        self.add_argument('--version',
                            help        = 'if specified, print version number',
                            type        = bool,
                            dest        = 'b_version',
                            action      = 'store_true',
                            optional    = True,
                            default     = False)

    def run(self, options):
        """
        Define the code to be run by this plugin app.
        """
        if options.b_version:
            print('Plugin Version: %s' % Z2labelmap.VERSION)
            sys.exit(0)



# ENTRYPOINT
if __name__ == "__main__":
    app = Z2labelmap()
    app.launch()
