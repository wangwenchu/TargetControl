# TargetControl
Codes for control capacity in target control


**MainTargetSampling**

c++ code for sampling the control capacity of all nodes
- input file: the edgelists in file netfile, each line represents a directed link
- output file: the control conpacity and fractions of redundant,cirtical,intermittent nodes

if the target set is randomly picked, one should set the random_pick_target to "true" and give the fraction of target nodes, otherwise, the target set should be read from file(target_file)

**MainEnumerateAllMatchedSetStates**

c++ code for enumerating all the matched sets in out-set and in-set, and the all matched 2-triple (MMSout,MMSin)

