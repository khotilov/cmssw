import FWCore.ParameterSet.Config as cms
import copy

#--------------------------------------------------------------------------------
# generic utility functions to modify analysis sequence,
# usable for all channels
#--------------------------------------------------------------------------------

def switchAnalyzers(analysisSequence, analyzers):
    for pset in analysisSequence:
        if hasattr(pset, "analyzers") : setattr(pset, "analyzers", analyzers)

def addAnalyzer(genAnalyzerModule, analyzer, afterCutName, replaceStatement = None):
    
    # add configuration parameterset of analyzer to genAnalyzerModule
    # in case it is not alread added
    analyzerName = getattr(analyzer, "pluginName").value()
    analyzer_isDefined = False
    
    genAnalyzerModule_analyzers = getattr(genAnalyzerModule, "analyzers")
    for genAnalyzerModule_analyzer in genAnalyzerModule_analyzers:
        genAnalyzerModule_analyzerName = getattr(genAnalyzerModule_analyzer, "pluginName").value()
        if genAnalyzerModule_analyzerName == analyzerName:
            analyzer_isDefined = True

    if not analyzer_isDefined:
        analyzers = getattr(genAnalyzerModule, "analyzers").value()
        analyzers.append(analyzer)
        setattr(genAnalyzerModule, "analyzers", analyzers)

    # insert analyzer into analysis sequence
    # at position after cut passed as function argument has been applied
    genAnalyzerModule_analysisSequence = getattr(genAnalyzerModule, "analysisSequence")
    analysisSequence_mod = []
    insert = False
    for pset in genAnalyzerModule_analysisSequence:
        if insert:
            if hasattr(pset, "analyzers"):
                pset_analyzers = getattr(pset, "analyzers").value()
                pset_analyzers.append(analyzerName)
                setattr(pset, "analyzers", cms.vstring(pset_analyzers))

                if replaceStatement is not None:
                    pset_replace = getattr(pset, "replace").value()
                    pset_replace.append(replaceStatement)
                    setattr(pset, "replace", cms.vstring(pset_replace))
            else:
                pset_analyzers = cms.PSet(
                    analyzers = cms.vstring(analyzerName)
                )
                
                if replaceStatement is not None:
                    setattr(pset_analyzers, "replace", cms.vstring(replaceStatement))

                analysisSequence_mod.append(pset_analyzers)

            insert = False
        
        if hasattr(pset, "filter"):
            pset_filter = getattr(pset, "filter").value()
            if pset_filter == afterCutName:
                insert = True

        analysisSequence_mod.append(pset)

    setattr(genAnalyzerModule, "analysisSequence", cms.VPSet(analysisSequence_mod))

def replaceAnalyzerInputTags(analysisSequence, replacements):
    # check that replacement is a list of type [ replace_1, replace_2 .. replace_N ]
    # with each replace_i being a list with exactly two entries
    # (first entry = string to be replace, second entry = substitute string)
    for replacement in replacements:
        assert isinstance(replacement,list) and len(replacement) == 2, \
               "Error in <replaceAnalyzerInputTags>: Invalid argument replacements = " + replacements

    # substitute in all "replacement" attributes of analysisSequence object given as function argument
    # all occurences of replace_i[first entry] by replace_i[second entry] (for all replace_i)
    for pset in analysisSequence:
        if hasattr(pset, "replace"):
            entries_orig = getattr(pset, "replace")
            for iEntry in range(len(entries_orig)):
                entry_orig = entries_orig[iEntry]
                entry_mod = entry_orig
                
                for replacement in replacements:
                    replacement_orig = replacement[0]
                    replacement_mod = replacement[1]

                    if entry_mod.find(replacement_orig) != -1 : entry_mod = entry_mod.replace(replacement_orig, replacement_mod)

                entries_orig[iEntry] = entry_mod

def replaceAnalyzerTitles(analysisSequence, replacements):
    # check that replacement is a list of type [ replace_1, replace_2 .. replace_N ]
    # with each replace_i being a list with exactly two entries
    # (first entry = string to be replace, second entry = substitute string)
    for replacement in replacements:
        assert isinstance(replacement,list) and len(replacement) == 2, \
               "Error in <replaceTitles>: Invalid argument replacements = " + replacements

    # substitute in all "title" attributes of analysisSequence object given as function argument
    # all occurences of replace_i[first entry] by replace_i[second entry] (for all replace_i)
    for pset in analysisSequence:
        if hasattr(pset, "title"):
            entries_orig = getattr(pset, "title")

            for replacement in replacements:
                replacement_orig = replacement[0]
                replacement_mod = replacement[1]    

                if str(entries_orig) == str(replacement_orig):
                   setattr(pset, 'title', replacement_mod)

def removeAnalyzer(analysisSequence, analyzerName):
    # remove all analyzers with name given as function argument from analysisSequence object
    for pset in analysisSequence:
        if hasattr(pset, "analyzers"):
            analyzers = getattr(pset, "analyzers")

            if analyzers.count(analyzerName) > 0:
                analyzers.remove(analyzerName)

