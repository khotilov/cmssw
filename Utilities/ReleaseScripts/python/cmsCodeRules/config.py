__author__="Aurelija"
__date__ ="$2010-08-06 14.27.51$"

import os
ordering = ['1', '2', '3', '4', '5', '6']

# default values for directories

checkPath = os.getcwd()
picklePath = os.getcwd()
txtPath = os.getcwd()
htmlPath = os.getcwd()

# exception for directories and files 

exceptPaths = []

# --------------------------------------------------------------------------------
# configuration info for each rule ...

rulesNames = []
Configuration = {}

# --------------------------------------------------------------------------------

# configuration for rule 1

ruleName = '1'
rulesNames.append(ruleName)
Configuration[ruleName] = {}

Configuration[ruleName]['description'] = 'Search for "using namespace" or "using std::" in header files'
Configuration[ruleName]['filesToMatch'] = ['*.h']
Configuration[ruleName]['exceptPaths'] = []
Configuration[ruleName]['skipComments']  = True
Configuration[ruleName]['filter'] = '(\susing|\Ausing)\s+(namespace|std::)' #should be regular expression
Configuration[ruleName]['exceptFilter'] = []

# --------------------------------------------------------------------------------

# configuration for rule 2

ruleName = '2'
rulesNames.append(ruleName)
Configuration[ruleName] = {}

Configuration[ruleName]['description'] = 'Search for CXXFLAGS flags that are set to -g or -O0 in BuildFile'
Configuration[ruleName]['filesToMatch'] = ['BuildFile', 'BuildFile.xml']
Configuration[ruleName]['exceptPaths'] = []
Configuration[ruleName]['skipComments']  = True
Configuration[ruleName]['filter'] = '\s(CXXFLAGS|CPPFLAGS)(\+|=|\w|\"|\'|-|\s)*(-g|-O0)(\s|\'|\")' #should be regular expression
Configuration[ruleName]['exceptFilter'] = []

# --------------------------------------------------------------------------------

# configuration for rule 3

ruleName = '3'
rulesNames.append(ruleName)
Configuration[ruleName] = {}

Configuration[ruleName]['description'] = 'Search for "catch(...)" statements in *.cc, *.cxx files'
Configuration[ruleName]['filesToMatch'] = ['*.cc', '*.cxx']
Configuration[ruleName]['exceptPaths'] = []
Configuration[ruleName]['skipComments']  = True
Configuration[ruleName]['filter'] = 'catch\(\s*\.\.\.\s*\)' #should be regular expression
Configuration[ruleName]['exceptFilter'] = []
# --------------------------------------------------------------------------------

# configuration for rule 4

ruleName = '4'
rulesNames.append(ruleName)
Configuration[ruleName] = {}

Configuration[ruleName]['description'] = 'Search for "copyright" declaration in *.c, *.cc, *.cxx, *.h files'
Configuration[ruleName]['filesToMatch'] = ['*.h', '*.c', '*.cc', '*.cxx']
Configuration[ruleName]['exceptPaths'] = []#could be file name, dir, fileName:line. But path should be only from that directory in which we are searching
Configuration[ruleName]['skipComments']  = False
Configuration[ruleName]['filter'] = '(\A|\W)(c|C)(o|O)(p|P)(y|Y)(r|R)(i|I)(g|G)(h|H)(t|T)\W(\+|=|\w|\"|\'|-|\s)*(\((c|C)\)|\d{4})' #should be regular expression
Configuration[ruleName]['exceptFilter'] = ["\sLineo,"]
# --------------------------------------------------------------------------------

# configuration for rule 5

ruleName = '5'
rulesNames.append(ruleName)
Configuration[ruleName] = {}

Configuration[ruleName]['description'] = 'Search for "pragma" statement in *.c, *.cc, *.cxx, *.h files'
Configuration[ruleName]['filesToMatch'] = ['*.h', '*.c', '*.cc', '*.cxx']
Configuration[ruleName]['exceptPaths'] = []#could be file name, dir, fileName:line. Path should be only from that directory in which we are searching
Configuration[ruleName]['skipComments']  = True
Configuration[ruleName]['filter'] = '#\s*pragma\s' #should be regular expression
Configuration[ruleName]['exceptFilter'] = []
# --------------------------------------------------------------------------------
# configuration for rule 5

ruleName = '6'
rulesNames.append(ruleName)
Configuration[ruleName] = {}

Configuration[ruleName]['description'] = 'Search for "flags" statements in BuildFile'
Configuration[ruleName]['filesToMatch'] = ['BuildFile', 'BuildFile.xml']
Configuration[ruleName]['exceptPaths'] = []#could be file name, dir, fileName:line. Path should be only from that directory in which we are searching
Configuration[ruleName]['skipComments']  = True
Configuration[ruleName]['filter'] = '<\s*(f|F)(l|L)(a|A)(g|G)(s|S)\s+' #should be regular expression
Configuration[ruleName]['exceptFilter'] = ['EDM_PLUGIN', 'GENREFLEX_ARGS', 'TEST_RUNNER_ARGS', 'INSTALL_SCRIPTS']
# --------------------------------------------------------------------------------

rulesDescription  = "Rule number    Description\n"
rulesDescription += "----------------------------------------------------------------------------------------\n"
for key, value in Configuration.items():
    rulesDescription += "     %s         %s\n" %(key, value['description'])

# --------------------------------------------------------------------------------
helpMsg  = "-----------------------------------------------------------HELP-----------------------------------------------------------\n"
helpMsg += "cmsCodeRulesChecker.py [-h] [-html] [-s [DIRECTORY]] [-S [DIRECTORY]] [-p] [-r ruleNumber[,ruleNumber[, ...]]] [-d DIRECTORY]\n\n"
helpMsg += "-r     Specifies rule or rules to be checked. After this parameter should\n       be at least one rule given.\n"
helpMsg += "-d     Specifies that rules should be checked in DIRECTORY. Default \n       directory - current directory\n"
helpMsg += "-S     Specifies to save results in python pickle files. Directory specifies\n       where to store these files. Default directory - current directory\n"
helpMsg += "-s     Specifies to save results in .txt files. Directory specifies where to\n       store these files. Default directory - current directory\n"
helpMsg += "-p     Specifies to print results into a screen\n"
helpMsg += "-h     Prints help message\n"
helpMsg += "-html  Reads pickle files and creates cmsCRPage.html\n\n"
helpMsg += "By default cmsCodeRulesChecker.py checks all rules in current directory and prints results into screen.\n\n"
helpMsg += rulesDescription
