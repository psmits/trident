**Reviewer comment**

We are confronted with four models and even the simplest one is already quite complex. Is there a chance to provide for each model the relative importance of each variable for the model performance? Alternatively, as standing geographic range is among the best predictors of extinction risk in paleo studies, can the authors add a very simple model of just event ~ range? I’d be curious to see if the AUC is then “significantly” lower than in the other models. Of course, in a frequentist world, I would start with a complex model and then perform model selection with a step function. I’m not Bayesian enough to suggest a specific strategy for this paper, but somehow model performance needs to be balanced by model complexity. The authors seem to have compared just model performance.




**Author response**

We've included a figure comparing group-level effects for our covariates, along with the estimates of their variation over time for models that allow this variation. The estimates depicted in this figure are the *average* effects of those covariates over the Cenozoic. In models where covariate effects are allowed to vary over time, the actual effect for an interval may be different. The magnitude of the variance for that covariate is a good indicator of how consistent the effect is on average. Covariate effects with high variance can be positive in one interval and negative in the next.

This figure allows direct comparison in the relative magnitude of effect for each of the covariates. We consider this much more informative than fitting and comparing additional models because covariate effects are much more illustrative of *what* matters than measures of performance using simpler models. 

Additionally, simpler models (e.g. event ~ range) are statements that, with 100% probability, that covariates like temperature, age, and taxonomic group have 0 effect of extinction. I don't know about the reviewer, but I'm extremely uncomfortable making that assumption. Instead, we are using regularizing priors to shrink the effects of weak covariates towards 0. We consider this approach perfrable because we do not make strong assumptions about our covariates, we are able to include all covariates of substatntitive and scientific interest simultaneously, and our models better reflect our (lack) of knowledge. Additionally, covariate *effects* are the fundamental measure of variable importance and utility, not comparisons between different models.

The bias-variance tradeoff between complexity and performance for a set of models is estimated from how well that model is capable of predicting out-of-sample data. AIC is an approximation of out-of-sample model performance as measured by leave-one-out cross-validation, which makes AIC useful in situations where cross-validation is slow. We've used 5-fold cross-validation to estiamte out-of-sample model performance. This means are model comparisons in herently take complexity into account because we're measureing each models ability to predict out-of-sample data. A model that is too complex or not complex enough will have poor out-of-sample model performance. So, in summary, our measure of model performance are fundamentally taking into account model complexity because we are measuring our models ability to predict out-of-sample data. We've updated some of the language in the manuscript to better explain this distinction.

We've made a table comparing the WAIC and LOOIC estimates from our models as proof of their functionally identical meaning between these measures and using out-of-sample performance. We caution, however, as the WAIC/LOOIC estimates do not take into account the the temporal structure of our data. Those methods exist for LOOIC that accounts for temporal structure (CITATION), none of have been implemented in existing R packages.

Though see http://mc-stan.org/loo/articles/loo2-lfo.html ? It will take a while to implement fully and we're on a major time crunch. Fitting 4 models to single data takes ~ 2 day. Fitting each model 60+ times, even to smaller amounts of data, will take too long with our 6 week timeline.

Also, we caution to reviewer to avoid step-wise procedures in their own work because those procedures are fundamentally biased. Please see https://www.stata.com/support/faqs/statistics/stepwise-regression-problems/ and https://stats.stackexchange.com/questions/20836/algorithms-for-automatic-model-selection/20856?fbclid=IwAR26uuqd6bV_6d0ex_cqPuwu2dremR01c67j34jugxsaKCEnj65CDFSZbFU#20856 for an extended discussion as to the pitfalls behind these algorithms.


