treatmentResponseDiallel
===========

## Summary

This package is used for the analysis of treatment-response phenotypes in a forthcoming influenza diallel paper. It can be installed and used in conjunction with code and specific dependencies housed at [https://github.com/mauriziopaul/flu-diallel](https://github.com/mauriziopaul/flu-diallel).

The findings in this study have been submitted for publication as Maurizio et al., 2017, in *G3: Genes, Genomes, Genetics* (submitted 2017-01-05, in revision). A static version of the data, software, and scripts used to analyze this data upon submission is available at DOI: [10.5281/zenodo.293015](http://dx.doi.org/10.5281/zenodo.293015).

## Installation

You can install `treatmentResponseDiallel` using the following code.

First, in R:

```
install.packages(c('coda', 'corpcor','R.oo'))
```

Some additional dependencies are located here [https://github.com/mauriziopaul/flu-diallel/tree/master/packages](https://github.com/mauriziopaul/flu-diallel/tree/master/packages). Then, on the command line (for Mac, replacing * with version number):

```
R CMD install BayesDiallel_*.tar.gz;
R CMD install BayesSpike_*.tar.gz;
R CMD install cmdline_*.tar.gz;
R CMD install WVmisc_*.tar.gz;
R CMD install configfile_*.tar.gz;
```

Then, from R:

```
# install.packages("devtools")
devtools::install_github("mauriziopaul/treatmentResponseDiallel")
library("treatmentResponseDiallel")
```

## Related Software

1. **BayesDiallel** [http://valdarlab.unc.edu/software/bayesdiallel/BayesDiallel.html](http://valdarlab.unc.edu/software/bayesdiallel/BayesDiallel.html)

2. **Diploffect** [http://valdarlab.unc.edu/software/Diploffect/build/html/index.html](http://valdarlab.unc.edu/software/Diploffect/build/html/index.html)

## More Information

> For more information about my research interests, please visit [https://mauriziopaul.github.io/](https://mauriziopaul.github.io/).
> To learn more about research in the Heise and Valdar labs, please visit and [https://unclineberger.org/people/mark-heise](https://unclineberger.org/people/mark-heise) and [http://valdarlab.unc.edu](http://valdarlab.unc.edu).

## Key References

Maurizio PL, Ferris MT, Keele GR, Miller DR, Shaw GD, Whitmore AC, West A, Morrison CR, Noll KE, Plante KS, Cockrell AS, Threadgill DW, Pardo-Manuel de Villena F, Baric RS, Heise MT & Valdar W (2018) Bayesian diallel analysis reveals *Mx1*-dependent and *Mx1*-independent effects on response to influenza A virus in mice. *G3: Genes, Genomes, Genetics* 8(2):427-445. [10.1534/g3.117.300438](https://doi.org/10.1534/g3.117.300438). PMID:29187420.

Maurizio PL & Ferris MT (2017) “Chapter 28: The Collaborative Cross Resource for Systems Genetics Research of Infectious Diseases.” *Methods in Molecular Biology: Systems Genetics – Methods and Protocols*. Springer Science+Business Media, New York, NY. Klaus Schughart and Robert W. Williams (eds.) 1488:579-596. eBook ISBN: 978-1-4939-6427-7, hardcover ISBN: 978-1-4939-6425-3. doi: [10.1007/978-1-4939-6427-7_28](https://doi.org/10.1007/978-1-4939-6427-7_28). PMID:27933545.

Zhang Z, Wang W, Valdar W (2014) Bayesian modeling of haplotype effects in multiparent populations. *Genetics* 198(1):139-56. doi: [10.1534/genetics.114.166249](http://dx.doi.org/10.1534/genetics.114.166249)

Crowley JJ, Kim Y, Lenarcic AB, Quackenbush CR, Barrick C, Adkins DE, Shaw GS, Miller DR, Pardo Manuel de Villena F, Sullivan PF, Valdar W (2014) Genetics of adverse reactions to haloperidol in a mouse diallel: A drug-placebo experiment and Bayesian causal analysis. *Genetics* 196(1):321-47. doi: [10.1534/genetics.113.156901](http://dx.doi.org/10.1534/genetics.113.156901)

Ferris MT, Aylor DL, Bottomly D, Whitmore AC, Aicher LD, Bell TA, Bradel-Tretheway B, Bryan JT, Buus RJ, Gralinski LE, Haagmans BL, McMillan L, Miller DR, Rosenzweig E, Valdar W, Wang J, Churchill GA, Threadgill DW, McWeeney SK, Katze MG, Pardo-Manuel de Villena F, Baric RS, Heise MT (2013) Modeling host genetic regulation of influenza pathogenesis in the Collaborative Cross. *PLoS Pathogens* 9(2):e1003196. doi: [10.1371/journal.ppat.1003196](http://dx.doi.org/10.1371/journal.ppat.1003196)

Lenarcic AB, Svenson KL, Churchill GA, Valdar W (2012) A general Bayesian approach to analyzing diallel crosses of inbred strains. *Genetics* 190:413-435. doi: [10.1534/genetics.111.132563](http://dx.doi.org/10.1534/genetics.111.132563)
