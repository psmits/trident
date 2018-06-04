intro
=====

how does geographic range size influence species risk of extinction.

global geographic range is a property of the species (endogenous variable). in contrast, global temperature would be a exogenous variable beause it is not a property of the species but of the environment.

change to geographic range has two properties: current state (position) and change in state over time (slope/velocity). both of these may impact survival. we "know" that geographic range is an important and major predictor of species extinction risk, with small geographic range predicting greater extinction risk than large geographic range. however, we don't have a good idea of how **the change in geographic range over time** affects instantneous extinction risk. does extinction risk change as a function of changing geographic range?

all of this implies the possibility of non-constant hazard with species age. there is some evidence of increasing risk with age which is consistent with a few theories (mammals, forams). there is also some evidence of decreasing risk with age (foote). of course, we consider the default or null hypothesis to be Law of Constant Extinction which states that extinction risk is memoryless and cares not for taxon age.

can you have "wax and wane" effect of geographic range but still estimate a fundamentally monotonic baseline hazard function? this means that conditional on etavalue/etaslope, hazard is monotonic and could be modeled effectively by the expontial/weibul distribution(s).

if we assume range decreases towards the end of species duration, than this would imply increasing risk at old ages. however, this would also imply high risk in young species (there is evidence of this). this would actually mean non-monotonic hazard, something that has never been thoroughly investigated. 

both exponential and weibull hazard are monotonic functions and would be unable to estimate nonmonotonic hazard functions. 

log-normal and generalized gamma hazards an produce nonmonotonic hazard functions.

weird case is when we have nonmonotonic hazard due to changes to geographic range through time but that accounts for all the time-dependence and can yield ultimately "monotonic hazard" conditional on geographic range size. this has never been considered.

instantaneous value as well as derivative.

less concerned with overall shape of geographic range through time than how those changes map onto species instantaneous risk of extinction rate.

can also talk about effect of species age on extinction risk because of parameteric hazard.




theoretical issues in play
--------------------------

law of constant extinction: how does species age affect extinction risk?

rise-and-fall of geographic range: species geographic range is not constant through time and tends to rise then fall. note: we're less concerned with the shape of range over time so much as the effect those changes have on extinction risk.




analytical tools and concepts in play
------------------------

survival data. key details: survival function, hazard function, censoring.

longitudinal data. key details: value over time.

joint model of survival and longitudinal data. key details: non-random drop out (censoring) and endogenous time-varying covariates need joint model.

hazard function. instantaneous risk of failure. conditional statement.

survival function. probability of survival at time t given survival up till that point.

censoring. incomplete observation or event of interest has not occurred during observation period.

joint model allows man aspects of the longitudinal model to be used as a covariate in the survival process. e.g. expected value, derivative of function, area under function, etc.




discrete-time hazard model. this is pretty standard when we're dealing with naturally discrete time units. seth has the pieces that are assembled to form this model, but has never considered it this way. this model is built directly out of the discrete-time survival literature and represents something close to best practices for this type of data.

paleo data has been analyed as both discrete time and continuous time (including by peter). however, as we look deeper at the data is becomes clear to in most cases a discrete time approach makes more sense. the time units in paleontological data are brackets between t and t + dt. also, in continuous time there should never be "ties" or events occurring at identical moments; these are the hallmark of an actually discrete time occurrence. while there are certainly events in earth's history that could cause ties (e.g. mass extinctions) these are inherently an exception and not a rule. Additionally, these ties are a function of the time averaging of the fossil record, so what might be reasonable ties at the million year discrete time scale, if we were to actually observe these events in continuous time there would be no ties (species go extinct as a random process).
