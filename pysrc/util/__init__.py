#SETTINGS FOR CLUSTER JOB ARRAYS, TO WORK WITH CELL COGNITION MASTER - not especially to work with R

jobSize = 10
progFolder = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/pysrc'
scriptFolder = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/scripts'
path_command = """setenv PATH /cbio/donnees/nvaroquaux/.local/bin:${PATH}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/share/apps/libxml2/lib:/share/apps/libxslt/lib:/cbio/donnees/nvaroquaux/.local/lib
setenv LIBRARY_PATH /cbio/donnees/nvaroquaux/.local/lib
setenv PYTHONPATH /cbio/donnees/twalter/workspace/cellh5/pysrc:/cbio/donnees/aschoenauer/software/lib/python2.7/site-packages:/cbio/donnees/aschoenauer/workspace2/cecog_1.5/cecog/:/cbio/donnees/aschoenauer/workspace2/Xb_screen/pysrc
setenv DRMAA_LIBRARY_PATH /opt/gridengine/lib/lx26-amd64/libdrmaa.so
"""
pbsOutDir = '/cbio/donnees/aschoenauer/PBS/OUT'
pbsErrDir = '/cbio/donnees/aschoenauer/PBS/ERR'
pbsArrayEnvVar = 'SGE_TASK_ID'


#ctrl wells for labteks, id<50
typeD={}
#CONTROLES NEGATIFS
typeD["scrambled"] = ["015", "026", "063", "074", "304", "315", "352"]

#CONTROLES POSITIFS
#typeD["bcop"] = ["363"]
#typeD["eg5"] = ["290", "301", "338", "349"]
#typeD["incenp"] = ["004", "049", "052"]

#CONTROLES QUI ONT ETE SUPPRIMES OU SONT INUTILISABLES
#typeD["empty"] = []
#typeD["marker"] = ["001"]

typeD2={}
#CONTROLES NEGATIFS
typeD2["scrambled"] = ["015", "026", "063", "074", "304", "315", "352", "363"]

#CONTROLES POSITIFS
#typeD2["bcop"] = ["333", "336", "381", "384"]
#typeD2["eg5"] = ["290", "301", "338", "349"]
#typeD2["incenp"] = ["004", "049", "052"]

#CONTROLES QUI ONT ETE SUPPRIMES OU SONT INUTILISABLES
#typeD2["empty"] = ["311", "322", "359", "370"]
#typeD2["marker"] = ["001"]

ctrlWell1 = []; ctrlWell2=[]
for el in typeD:
    ctrlWell1.extend(typeD[el])
for el in typeD2:
    ctrlWell2.extend(typeD2[el])


"""
NB: we only want to use negative controls for migration studies
"""