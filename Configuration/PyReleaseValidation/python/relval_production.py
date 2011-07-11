
# import the definition of the steps and input files:
from  Configuration.PyReleaseValidation.relval_steps import *

# here only define the workflows as a combination of the steps defined above:
workflows = {}

# each workflow defines a name and a list of steps to be done. 
# if no explicit name/label given for the workflow (first arg),
# the name of step1 will be used

## data production test
workflows[1000] = [ '',['RunMinBias2011A','TIER0','HARVESTD','ALCASPLIT']]
#workflows[1001] = [ '',['RunMinBias2011A','TIER0EXP','ALCAPROMPT','ALCAHARVD']]
workflows[1002]=['RRD',['RunMinBias2011A','RECOD','RECODFROMRAWRECO']]

## MC production test
#workflows[1100] = [ '',[]]

workflows[1102]=['RR', ['TTbar','DIGI','RECO','RECOFROMRECO']]


