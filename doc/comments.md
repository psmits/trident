Reviewer comments
==================


Reviewer 1
----------

Smits and Finnegan use the Neptune database to check which combination of potentially influencing variables has the best predictive power of extinction risk in the marine plankton. They use sophisticated logistic modelling of the time-binned raw data to suggest that considering changes of parameters only marginally improves prediction accuracy. Overall, this is a very well written, although difficult to comprehend, account, which is to the point and rich in insights.

While the main result is well founded by the sophisticated analyses (in fact too sophisticated to be understandable in all details by the reviewer), it is somewhat disappointing, as the relative importance of parameters in multivariate models remain vague. It is also disappointing to see that change of range had a negligible effect in the multivariate context, as this is one of the main arguments for using the fossil record in extinction risk assessments (besides taxon age and intrinsic turnover rates). I argue that, to really discard change of range, the authors would need to look at longer trajectories of range. Bin-to-bin changes at this temporal resolution (1 myr) are perhaps too noisy to inform models. I am not asking for additional analyses for this paper, but this point should be added to the discussion.

We are confronted with four models and even the simplest one is already quite complex. Is there a chance to provide for each model the relative importance of each variable for the model performance? Alternatively, as standing geographic range is among the best predictors of extinction risk in paleo studies, can the authors add a very simple model of just event ~ range? I’d be curious to see if the AUC is then “significantly” lower than in the other models. Of course, in a frequentist world, I would start with a complex model and then perform model selection with a step function. I’m not Bayesian enough to suggest a specific strategy for this paper, but somehow model performance needs to be balanced by model complexity. The authors seem to have compared just model performance.

The methods are very well explained (especially in the supplement), although uncertainties in parameter estimates merit some further prose. Using maximum great circle distance (GCD) as a measure of geographic range is ok but prone to sampling bias. Better sampled intervals will inevitably have greater mean ranges than more poorly sampled intervals. While the scaling applied by the authors takes different means and variance into account it may not be sufficient. If GCD is used, I’d recommend normalizing to the maximum observed geographic range in an interval prior to scaling. Grid cell occupancy or minimum spanning trees are alternative measures, which are worthy of exploration. 

Looking at change of range is good, but the authors could also look at changes beyond two intervals. Although models incorporating multiple bins may result in overfitting, the argument is that long-term decline will increase extinction risk. I’d argue that two bins (especially at the relatively fine temporal scale applied) can hardly measure long-term decline. Finally, the authors should provide a rationale why these and not other parameters were used in the study (see also below). Does latitude or skeletal mineralogy have little to add to model performance?

Presentation of results could be improved. Showing raw data of extinction rates in the four groups, and individual analyses how age and phyla actually influence extinction risk would be preferred over the wiggly lines in Figs. 2 and 3. The latter could be moved to the supplement. The text in the results section is very technical, describing too much detail. 

The discussion is good, but may have to be partly rewritten in light of my comments. 

Specific comments (line numbers refer to the smaller ones on the right):

- “the relative risk of extinction exhibited by different taxonomic groups” could be abbreviated as “intrinsic extinction risk”
- l. 68-69. Partial repetition of sentences before. “Incredible” should be omitted. If incredible why trust the age models?
- l. 79. “Ecologically unusual” should be omitted from reasoning. Not needed. 
- l. 88-91. Authors may want to add that Mg/Ca has its own issues, especially the steep rise of Mg/Ca in seawater through the Cenozoic. Cramer et al. have discussed and accounted for this trend, but with some assumptions, which result in considerable uncertainties.
- l. 132: paragraph starts with k-fold but is specified for 5-fold. I recommend rewriting this paragraph for the general application and then specify.
- l. 149-150: “so we interpret values between 0.7 and 0.8 could then be considered”  “so we consider values between   0.7 and 0.8  as” 
- Caption Fig. 1: Cut “a” from “with a higher AUC values”. It is not entirely clear what the error bars and bean plots represent. Are the AUC values derived from the 62 time bins or is there more to it? Please explain.
- l. 173: “Dinoflagellates” should be “diatoms” unless I have missed something. Same in labels of Figs. 2-4 and some incidences in the text. Please correct throughout.
- Xlab of Fig. 2 should be “Mya” or similar instead of “CI”
- l. 176-177. Repetition from methods.
- l. 183. “Comparison”  “Comparing”
- l. 192-193. No single-sentence paragraphs
- l. 197. “diatoms” instead of “Dinoflagellates”
- l. 214-216. To my understanding AUC does not simply measure a rank-order analysis but model performance. That the AUC’s of the eight models are statistically indistinguishable (are they?) does not necessarily imply a correct forecasting of rank-order extinction probability.  Please rephrase!
- l. 223-224. “Variation of extinction intensity”. Would be useful to either depict this variation or report some summary statistics (median, range, variance).
- l. 233. What does “scientifically significant” mean is this context? I am aware that simple model selection is not applicable here but it would be good to have a comparable metric.
- l. 236-237. Kiessling and Kocsis and used more than just corals. Their coarser temporal resolution (geological stages) is an obvious suspect for the greater role of change of occupancy than reported here. 
- l. 244-245. Poor reasoning. For example skeletal mineralogy (silica versus calcium carbonate), trophic level (phyto- vs. zooplankton), latitudinal preference (e.g., median paleolatitude) or additional physico-chemical parameters (e.g., sea-level, Mg/Ca, pCO2) are intuitive variables, whereas taxon-age is less intuitive in light of the recent literature. I’m not asking for additional analyses or more complex models but a sound rationale.
- l. 247. Add “for” between “not exist” and “all”. In any case, a poor defense given the obvious traits to look at (see above).
- l. 254. Modern “conservation determinations” are not continuous but ordinal. The difference between “critically endangered” and “endangered” may be very different from the difference between “endangered” and “vulnerable”. 
- l. 255-256. I argue that the authors’ results actually show that fossil data are not very relevant. If this impression is wrong, the authors should perhaps rephrase some parts of the text.

Also check use of “My” vs. “Mya” and cases in plankton groups. As a general rule these should be upper case in a formal context (e.g., Radiolaria) and lower case when used informally (e.g., radiolarians). 

The supplement is good but lacks sensitivity tests and a rationale of why the specific parameters and not others where chosen in the models. 

Review by Wolfgang Kiessling (wolfgang.kiessling@fau.de)




Reviewer 2
----------

This manuscript uses a high-resolution record of marine microfossils to test the predictive power of fossil record extinctions. A lot of paleontological research has focused on determining the importance of particular biological and environmental risk factors, so this study is novel in testing whether past history can predict future extinction. The methodology is explained well and the main result is well-supported.

I only have two broader thoughts and two more specific notes:

I find it interesting that the predictability of extinction doesn’t obviously differ between background intervals and intervals with more unusual environmental perturbations (such as the PETM or E-O transition). Perhaps that’s because those events weren’t mass extinctions in the vein of the P/T or K/Pg. There’s also a fair amount of volatility in the AUC time-series. Nevertheless, does this consistent predictive ability regardless of differences in conditions say anything interesting about extinction predictions?

The caveat that human impacts may dramatically alter extinction risk seems like an important one and, as it stands, seems like it could undercut the (rather terse) final conclusion. Would conservation decisions really be bolstered by including fossil data? I would agree with that (and of course I would, as a paleobiologist), but it might be useful to expand a bit on this topic. Do you mean that a model such as yours, incorporating geographic range and other such parameters, would be helpful? In what way could it be used? Or do you mean using models based on past intervals with environmental changes inferred to be similar – the PETM as an analogue for warming and acidification, for example? I wonder if it would be helpful to draw a stronger connection linking your finding to conservation biology, given the different stressors now.

I found the description on lines 171-174 to be a bit confusing. The phrase “this pattern” (line 171) seems to refer to the posterior distribution less than/equal to 0.5. The sentence states that the pattern is absent in absent in VP model for forams and radiolarians, but to me the radiolarian pattern doesn’t look too different from the V model (which also doesn’t dip below 0.5). The next sentence argue that there are fewer period of low performance for calcareous nannos and dinoflagellates in the VP model, but the differences seem minimal. Perhaps I’m not familiar with these graphs, but there only seems to be a difference in calcareous nannos around 50 Ma, and among dinoflagellates a small difference around 58 Ma and 30 Ma?

Figures 2-4. I guessed (and later confirmed) that the alternating gray and white bars represented stages of the geological time scale, but it would be helpful to indicate that in the captions.

Sincerely,

Matthew Clapham


Reviewer 3
----------

This is an interesting paper on extinction prediction using microfossil plankton data from deep-sea cores.
It’s very important that they found simple model with only a few parameters of geographic range, temperature etc can robustly predict extinction risk regardless of minor model difference.

I enjoyed my reading and hope the comments below help to improve the ms.

1. just in case, for foraminifera, benthic species are deleted from the database? Or they originally include planktonic foram only?

2. Neptune database has an internally consistent taxonomic identification strategy. Is this true? I guess it’s non-critical compilation of census data done by shipboard micropaleontologists over ~50 years (ie taxonomic sense of different generation micropaleontologists are very different). David Lazarus (in charge person of Neptune) is radiolarian person, and so radiolarian may be better in taxonomy? Planktonic foram diversity is low, so may be better straightforward. Anyway, I recommend deeper reference search and reading to see how is Neptune taxonomy and explain it in the method section more in details.

3. Fig 1, better to explain what are C, V, CP, VP, AUC in the caption

4. Fig 2 and 3, What do gray bars mean? Gray bay = envelope?

5. It will help our understanding if you explain the main point/message of each fig in the contest of discussion in each caption 

6 I am not sure if it’s possible, but I am curious where geographic range or temperature are the key? Or both are essential? It is possible to try geographic range only model and temperature only model to see such?

7. The authors use deep-water temperature (using geochemistry of benthic foram from deep-sea cores) data. And their biotic data is plankton. Their plankton can be mostly from ocean surface or photic zone (eg diatoms), but some have good bunch of deep-water species (radioralians).
The authors need some discussion and justification on this.

8. While geographical range consider spatial stuff, the authors use one temperature curve for all (I understand there will be no other way). Is there any difference in temperature sensitivity between narrow and wide distribution species; or low latitude and high latitude species? I understand this may be out of the scope. The authors may ignore.

Sincerely yours,






