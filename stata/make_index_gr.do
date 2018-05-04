/***********************************************
*INVERSE COVARIANCE WEIGHTED INDEX V.0.2
*Cyrus Samii, NYU, December 2017
*(cds2083@nyu.edu)
*
*This ado file creates an index using inverse covariance weighting, 
*a la Anderson (2008, p. 1485).
*
*The command is 'make_index_gr' and it takes four arguments in a required 
*sequence: 
*	(1) Suffix for the name of the new index
*	(2) Treatment/sampling probability weight variable.
*	(4) Binary (0/1) variable to indicate set of units to use for standardizing.
*	(4) Name of a local macro with the variables to incorporate in the index.
*	
*An example is as follows.  Suppose we have variables "x", "y", and "z",
*along with sampling weights in a variable "wgt" and then an indicator
*variable stdgroup that takes the value 1 for units that we want to use
*as the basis of standardization. This is usually the "control group" in
*an experiment (as in Anderson, 2008). Then we could run,
*
*	. local local_foo x y z
*	. make_index_gr foo wgt stdgroup `local_foo' 
*	
*which adds a new variable called "index_foo" to your dataset. This variable 
*is the inverse-covariance weighted average of "x", "y", and "z", adjusted
*for the sampling weights in "wgt" and standardized against the means and
*standard deviations of units with "stdgroup==1".
*
*Note that I am not a professional programmer, so please use at your own
*risk! If you discover any errors, please write me to indicate. Known 
*issues:
* 	- If you construct your index "by hand" using Stata's "summarize" command,
*	  your results may differ slightly. This is because of differences in the
*	  way Stata's "summarize" normalizes weights when estimating a weighted
*	  standard deviation.  The implications should be negligible. 
*
*	- If you do not have any weights in your analysis, simply create a variable
*	  that is equal to 1 for all units (e.g., "generate wgt = 1") and use that
*	  in place of the weight variable. 
*
*	- If you wish to standaridized against the full-sample pooled means and
*	  standard and standard deviations, you can do similar (e.g., 
*	  "generate stdgroup=1", which standardizes against the pooled values).
*
*	- This replaces an older version of the program, which did not allow you
*	  to change the standardization group.
*
*	- I wrote this as a .do file that you can run or copy-paste into
*	  your own analysis .do file. Alternatively, you could save this as an
*	  .ado in your personal .ado directory, in which case you could just
*	  issue the "make_index_gr" command without having to run the .do.
*
*	- The final index may not have mean zero and standard deviation 1. This
*	  is because we are doing the first stage standardization with respect
*	  to group for which the standardization group indicator is 1. You may
*	  want to standardize the index again for the purposes of interpretability
*	  in your analysis. 
*
*References
*
*Anderson ML (2008) "Multiple Inference & Gender Differences..." JASA 103(484).
*
***********************************************/

capture program drop make_index_gr
capture mata: mata drop icwxmata()
program make_index_gr
version 11.1
    syntax anything [if]
    gettoken newname anything: anything
    gettoken wgt anything: anything
	gettoken sgroup anything: anything
    local Xvars `anything'
	marksample touse
  	mata: icwxmata(("`Xvars'"),"`wgt'","`sgroup'", "index")
	rename index index_`newname'
end
 
mata:
	function icwxmata(xvars, wgts, sgroups, indexname)
	{
		st_view(X=0,.,tokens(xvars))
		st_view(wgt=.,.,wgts)
		st_view(sgroup=.,.,sgroups)
		nr = rows(X)
		nc = cols(X)	
		sg_wgt = wgt:*sgroup
		sg_wgtst = sg_wgt/sum(sg_wgt)
		all_wgtst = wgt/sum(wgt)
		sg_wgtstdM = J(1,nc,1) # sg_wgtst
		all_wgtstdM = J(1,nc,1) # all_wgtst
		sg_wgtdmeans = colsum(X:*sg_wgtstdM)

		sgroup2 = sgroup
		sgroupdM = J(1,nc,1) # sgroup2
		sg_meandevs = ((X:*sgroupdM) - (J(nr,1,1) # sg_wgtdmeans):*sgroupdM)
		all_wgtdmeandevs = X - (J(nr,1,1) # sg_wgtdmeans)
		sg_wgtdstds = sqrt(colsum(sg_wgt:*(sg_meandevs:*sg_meandevs)):/(sum(sg_wgt)-1))
		Xs = all_wgtdmeandevs:/(J(nr,1,1) # sg_wgtdstds)

		S = variance(Xs, wgt)
		invS = invsym(S)

		ivec = J(nc,1,1)
		indexout_sc = (invsym(ivec'*invS*ivec)*ivec'*invS*Xs')'
		indexout = indexout_sc/sqrt(variance(indexout_sc, sg_wgt))
		st_addvar("float",indexname)
		st_store(.,indexname,indexout)
	}
end

