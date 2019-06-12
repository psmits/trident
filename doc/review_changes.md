Below are a list of changes I've made to our first submission and why.


- Saupe's minor textual suggestions
  - did not answer question about why out-of-sample predictions can be better than in-sample (answer is random chance and the uncertainty of our parameter estimates).
- Wolfgang's minor textual suggestions
  - did not deal with questions about justifying our choice of covariates

- measure of geographic range changed to minimum spanning tree instead of great circle distance
  - this measure is not scaled to the maximum distance in bin because that breaks the meaning of the lags/differences
  - if scaled, measures are a percent. the differences would then be change in percent, but because we've got a shifting denominator this percentage is not really coherent.

- increased the number of geographic range differences used as predictors
  - from 1 difference to 3 (range + range\_diff1 + range\_diff2 + range\_diff3)
  - I've not done this for temperature as that wasn't requested and i don't know if it is really necessary

- added explainer to text about how AIC, WAIC, etc are all approximations of cross-validation. Cross-validation naturally deals with model complexity because we're looking at out-of-sample estimates. Because this wasn't clear in the next, we've updated it to address the reviewers concerns.
  - brought up by reviewer 1

- improved figure captions to explain the alternating grey blocks which represent geological ages
  - requested by reviewer 2 and 3

- replaced all accidental references to dinoflaggelates with diatoms.
  - this was a mistake

- improved wording in figure captions point to table explaining our models
  - requested by reviewer 3.

- improved language explaining how Neptune has its own synonymy table, which makes identifications at least internally consistent.
  - requested by reviewer 3.

- added language describing our measure of global temperature more fully
  - reviewer 1 asked for discussion of the issues surrounding temperature estimates based on Mg/Ca
  - reviewer 3 asked for discussion of how our temperature estimate is a deep water estimate but our organisms are planktonic


**What remains**: 

- once the cross-validation procedure for models has completed, I will regenerate all figures and update our results sections to reflect these new figures.
- covariate effects figure



**What I need from you**: 

- rewrite the remaining sections where reviewers 2 and 3 requested clarification or additional discussion.
  - previous major events (PETM etc) where are model appears to maintain nominally fair predictive accuracy. the reviewer might not realize that we don't have out-of-sample estimates for some climatic events because they are in the first fold, and thus only used as training data.
  - 
- explain why we chose are covariates. reviewer 1 really harps on this and reviewer 3 asks why we choose to always include range and temperature
  - not including range or temp would mean we are assuming with 100% probability that they have no effect. i'm not comfortable making this assertion. 
  - our modeling strategy is to include predictors we suspect will matter and use regularizing priors to shrink their estimates towards 0.
- write/take lead on writing the actual response to reviewers. I'm not comfortable writing responses to reviewer 1.
