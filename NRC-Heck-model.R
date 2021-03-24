#=====================================================================
#
#	SIMULATION PROGRAM FOR NON-RESPONSE MODEL TESTING
#
#	PS LC, Sept. 1st, 2020
#
#     5th version for Working paper
#
#=====================================================================



install.packages("sampleSelection")
install.packages("latex2exp")
install.packages("lubridate")

library(sampleSelection)
library(latex2exp)
library(lubridate)

rm(list=ls(all=TRUE))

#For randomness setting, corresponding to the results presented in the article
set.seed(10010) 


#------------------------------------------------------------------------------------------------
#Number of loops for bootstrap and simulations (i.e. draw in parent population)

#Set the parameters of the loops for the simulations (parameters used for the article)
#Computing time for the whole program on AUSV3: about 8 hours
NbBouclesbootstrap <- 20000
NbSimulations <- 2000

#To be faster (attention! the results will be different from the article)
#Divide by 10 the number of bootstrap loops and the number of simulations
#Computing time for the whole program on AUSV3: about 1 hour
#	NbBouclesbootstrap <- 2000
#	NbSimulations <- 200

#------------------------------------------------------------------------------------------------

tt0 <- now()



#==========================================================================
#1. Population generation --------------------------------------
#==========================================================================

taille_pop_men <- 10000  #Number of households in general population


generepopulation <- function(taille_pop_menF) {
	x1 <- runif(n=taille_pop_men,min=2,max=5)
	x2 <- runif(n=taille_pop_men,min=0,max=2)
	x3 <- runif(n=taille_pop_men,min=0,max=1)
	rev_men <- 2*x1+1*x2-0.5*x3 + rnorm(n=taille_pop_menF,mean=0,sd=2)
	return(data.frame(x1,x2,x3,rev_men))
	}

population <- generepopulation(taille_pop_men)
x1 <- population$x1
x2 <- population$x2
x3 <- population$x3
rev_men <- population$rev_men

sprintf("Average income of simulated households = %1.5f +/- %1.5f (1 ecty)",
        mean(rev_men),(var(rev_men)/length(rev_men))^0.5)


#Figure 3 or the article
hist(rev_men,main = "", xlab=TeX("$y_i$"),ylab="density")


#==========================================================================
#2. DGP for endogenous selection    --------------------------------
#==========================================================================


#Vector of response probabilities (depends on per capita household income)

#The latent participation variable of the first [nprobrenforce] households is boosted by 0.2:
nprobrenforce <- 3000  #Number of households with increased probability of selection


#Definition of the generating function of respondents 
#-> returns respondents TRUE/FALSE 
autoselection <- function (rev_menF,taille_pop_menF,nprobrenforceF,cte,dev_revenu,sdnorm,surcroitprob) {
	sel <- c()
	u <- participation_star <- cte+rnorm(n=taille_pop_menF,mean=0,sd=sdnorm)+(max(rev_menF)-rev_menF)/dev_revenu
	zz <- c(rep(surcroitprob,nprobrenforceF ),rep(0,taille_pop_menF-nprobrenforceF ))
	participation_star <- participation_star+zz
	select_men <- (participation_star>=0)
	sel$pi <- pnorm((cte+(max(rev_menF)-rev_menF)/dev_revenu+zz)/sdnorm)
	sel$pisszz <-  pnorm((cte+(max(rev_menF)-rev_menF)/dev_revenu)/sdnorm)
	sel$selectionne <- select_men
	sel$w <- as.numeric(zz>0)
	return(sel)  # list of selected households
}


#Application of self-selection on a draw of respondents:
sel <- autoselection(rev_men,taille_pop_men,nprobrenforce,-0.4,30,0.2,0.2)
selectionne <- sel$selectionne #binary variable of participation
probMen <- sel$pi #response probability
revenus_selectionnes <- rev_men[selectionne]
sprintf("Number of self-selected households in mean = %1.0f (/ %1.0f)",sum(probMen),length(probMen))
sprintf("Number of self-selected households = %1.0f (/ %1.0f)",length(revenus_selectionnes),length(probMen))
sprintf("Average income of self-selected households = %1.5f +/- %1.5f (1 ecty)",mean(revenus_selectionnes),
	(var(revenus_selectionnes)/length(revenus_selectionnes))^0.5)
sprintf("Hajek estimator with true probabilities of inclusion = %1.5f",
	sum(revenus_selectionnes/probMen[selectionne])/sum(1/probMen[selectionne]))

wt <- (1/probMen[selectionne])/sum(1/probMen[selectionne])
muH <- sum(wt*revenus_selectionnes)
varmuH <- sum(wt^2*(revenus_selectionnes-muH)^2)
sprintf("Hajek estimator with true probabilities of inclusion = %1.5f +/- %1.5f (1 ecty)",
	muH,varmuH^0.5)


#Figure 4 of the article
kdepi <- density(sel$pi[1:3000],  kernel = "gaussian")
kdesszz <- density(sel$pi[3001:10000],  kernel = "gaussian")

plot(kdepi, main = "",type="l",lwd=3, col=1,xlim=c(0,1.1),ylim=c(0,4),xlab=TeX("$\\pi_i$"),ylab="density")
lines(kdesszz,col=2,lwd=3)
legend(x=0.2,y=3.5,
	legend = c(TeX("$i<=3000$"),TeX("i>3000")),
	col = c(1,2), lty = 1,lwd=3)

#Figure 5 of the article
plot(x=rev_men,y=sel$pi,xlab=TeX("$y_i$"),ylab=TeX("$\\pi_i$"))




#==================================================================================
#3. Heckman estimate with sampleSelection and estimated corrected averages
#   1-step model
#==================================================================================

#Definition of the calculation function of the Heckman estimator
#We estimate in 1 step the model by maximum likelihood


#Construction of the instrument (ww) and data.frame of observables
ww <- sel$w
donnees <- data.frame(c(rev_men,selectionne,x1,x2,x3,ww))


#Heckman in 1 step
ttz0 <- now()
heck1Et <- selection(selectionne ~ x1+x2+x3+ww,
	rev_men ~  x1+x2+x3, 
	data=donnees,
	method="ml")
summary(heck1Et)
ttz1 <- now()
sprintf("Heckman 1-step model calculation time = %s", as.duration(round(ttz1-ttz0,digits=3)))


ess <- c()
#prediction of y, conditionnal to participation r for the whole population
ess$eqn1cond <- predict(heck1Et,part="outcome",type="conditional",
	newdata=data.frame(x1,x2,x3,ww))
#prediction of y, unconditionnal to participation r for the whole population
ess$eqn1Ncond <- predict(heck1Et,part="outcome",type="unconditional",
	newdata=data.frame(x1,x2,x3,ww))


#3.A- CORRECTION BY IMPUTATION

#the rev_men vector is corrected for non-response: we keep the income of those who responded 
#and we impute the one of those who did not answer with E[outcome|Z,R=1] 
#where outcome is the conditional predictor from the previous Heckman model
rev_men_cor <- rev_men
rev_men_cor[as.logical(1-selectionne)] <- 
	ess$eqn1cond[as.logical(1-selectionne),"E[yo|ys=0]"]


#Study of the difference in predicted outcome between self-non-selected and participants
vv <- ess$eqn1cond[as.logical(1-selectionne),"E[yo|ys=0]"]-
	ess$eqn1cond[as.logical(1-selectionne),"E[yo|ys=1]"]
hist(vv)
#Comment: vv>0, so the incomes of those who self-non-selected are, 
#everything else being equal, higher than the incomes of those who participated.

sprintf("E[YC=Corrected income] = %1.5f",mean(rev_men_cor))
sprintf("Recall: E[Y0=Original income] = %1.5f; E[Y0|R=1] = %1.5f )",mean(rev_men),mean(revenus_selectionnes))


#3.B- CORRECTION BY ENDOGENOUS WEIGHT CALCULATION

#Study of conditional weight
ess$eqn2cond <- predict(heck1Et,part="selection",type="link") 
ess$residout <- residuals(heck1Et,part="outcome")
proba <- pnorm((ess$eqn2cond+coef(heck1Et)["rho"]/coef(heck1Et)["sigma"]*ess$residout)/
	((1-coef(heck1Et)["rho"]^2)^(0.5)))
summary(proba)
sprintf("Nb of probabilities <0.1 : %1.f (=>weight threshold=10)",sum(proba<0.1,na.rm=TRUE))

#Winsorisation of endogenous weights and estimation
probac <- proba
probac[probac<0.1] <- 0.1

sprintf("E[YC|endogenous weights] = %1.5f",sum(rev_men/probac,na.rm=TRUE)/sum(1/probac,na.rm=TRUE))
sprintf("Recall: E[Y0=Original income] = %1.5f; E[Y0|R=1] = %1.5f )",mean(rev_men),mean(revenus_selectionnes))



#==================================================================================
#4. Heckman estimate with sampleSelection and estimated corrected averages
#   2-steps model
#==================================================================================

#Definition of the calculation function of the Heckman estimator
#We estimate in 2 steps the model
# (1) rev_men = c+a*(srev_men)+b*(Nbindiv_dans_men)+epsilon
# (2) D = c'+a'*(srev_men)+b'*(Nbindiv_dans_men)+gamma


#Construction of the instrument (ww) and data.frame of observables
ww <- sel$w
donnees <- data.frame(c(rev_men,selectionne,x1,x2,x3,ww))

#Heckman in 2 steps
ttz0 <- now()
heck2Et <- selection(selectionne ~ x1+x2+x3+ww,
	rev_men ~  x1+x2+x3, 
	data=donnees,
	method="2step")
summary(heck2Et)
ttz1 <- now()
sprintf("Heckman 2-steps model calculation time = %s", as.duration(round(ttz1-ttz0,digits=3)))


ess <- c()
ess$eqn1cond <- predict(heck2Et,part="outcome",type="conditional",
	newdata=data.frame(x1,x2,x3,ww))
ess$eqn1Ncond <- predict(heck2Et,part="outcome",type="unconditional",
	newdata=data.frame(x1,x2,x3,ww))


#4.A- CORRECTION BY IMPUTATION

#the rev_men vector is corrected for non-response: we keep the income of those who responded 
#and we impute the one of those who did not answer with E[outcome|Z,R=1] 
#where outcome is the conditional predictor from the previous Heckman 
rev_men_cor <- rev_men
rev_men_cor[as.logical(1-selectionne)] <- 
	ess$eqn1cond[as.logical(1-selectionne),"E[yo|ys=0]"]


#Study of the difference in predicted outcome between self-non-selected and participants
vv <- ess$eqn1cond[as.logical(1-selectionne),"E[yo|ys=0]"]-
	ess$eqn1cond[as.logical(1-selectionne),"E[yo|ys=1]"]
hist(vv)
#Comment: vv>0, so the incomes of those who self-non-selected are, 
#everything else being equal, higher than the incomes of those who participated.


sprintf("E[YC=Corrected income] = %1.5f",mean(rev_men_cor))
sprintf("Recall: E[Y0=Original income] = %1.5f; E[Y0|R=1] = %1.5f )",mean(rev_men),mean(revenus_selectionnes))

 

#4.B- CORRECTION BY ENOGENOUS WEIGHT COMPUTATION

#Here it is necessary to be careful. Indeed, contrary to the one-step case, 
#the residuals of the equation of outcome integrates all regressors, including 
#the Mills ratio. However, for the calculation of weights, residuals from the 
#outcome equation are needed without the Mills ratio. Therefore, the same procedure 
#as in the one-step case should not be used again, otherwise the weight 
#estimates would be biased. This is the correct procedure.


ess$eqn2cond <- predict(heck2Et,part="selection",type="link")
#ess$residout <- residuals(heck2Et,part="outcome") est changé en
ess$residout <- model.frame(heck2Et)$rev_men-predict(heck2Et,part="outcome",type="unconditional")

proba <- pnorm((ess$eqn2cond+coef(heck2Et)["rho"]/coef(heck2Et)["sigma"]*ess$residout)/
	((1-coef(heck2Et)["rho"]^2)^(0.5)))
summary(proba)
sprintf("Nb of probabilities <0.1 : %1.f (=>poids limite=10)",sum(proba<0.1,na.rm=TRUE))

#Winsorisation of endogenous weights l.e.to 1 and estimation
probac <- proba
probac[probac<0.1] <- 0.1

sprintf("E[YC|endogenous weights] = %1.5f",sum(rev_men/probac,na.rm=TRUE)/sum(1/probac,na.rm=TRUE))
sprintf("Recall: E[Y0=Original income] = %1.5f; E[Y0|R=1] = %1.5f )",mean(rev_men),mean(revenus_selectionnes))



#==========================================================================
#5. Multiple generation and selection of general and selected populations -
#==========================================================================

#5.A- SOME REFERENCE RESULTS

#sel <- autoselection(rev_men,taille_pop_men,nprobrenforce,-0.4,30,0.2,0.2)
#selectionne <- sel$selectionne 
#probMen <- sel$pi
revenus_selectionnes <- rev_men[selectionne]
sprintf("Average number of self-selected households = %1.0f (/ %1.0f)",sum(probMen),length(probMen))
sprintf("Number of households actually self-selected  = %1.0f (/ %1.0f)",length(revenus_selectionnes),length(probMen))
sprintf("Average income of self-selected households = %1.5f",mean(revenus_selectionnes))
sprintf("Hajek estimate with the real probabilities of inclusion = %1.5f",sum(revenus_selectionnes/probMen[selectionne])/sum(1/probMen[selectionne]))
sprintf("Average income for all households = %1.5f",mean(rev_men))


#5.B- VARIANCE ASSOCIATED TO THE ONLY SELF-SELECTION AND ESTIMATION BY HECKMAN CORRECTION IN 1 STEP

EstimeHeckman1EtapeSim1 <- function(rev_menF,taille_pop_menF,nprobrenforceF,x1F,x2F,x3F) {
	selF <- autoselection(rev_menF,taille_pop_menF,nprobrenforceF,-0.4,30,0.2,0.2)
	selectionneF <- selF$selection
	wwF <- selF$w
	heck1EtF <- selection(selectionneF ~ x1F+x2F+x3F+wwF,
			rev_menF ~  x1F+x2F+x3F, 
			data=data.frame(c(rev_menF,selectionneF,x1F,x2F,x3F,wwF)),
			method="ml")
	

	pp1 <- predict(heck1EtF,part="selection",type="link")
	rr1 <- residuals(heck1EtF,part="outcome")
	probaF <- pnorm((pp1+coef(heck1EtF)["rho"]/coef(heck1EtF)["sigma"]*rr1)/
	((1-coef(heck1EtF)["rho"]^2)^(0.5)))
	probacF <- probaF
	probacF[probacF<0.1] <- 0.1
	estim <- sum(rev_menF/probacF,na.rm=TRUE)/sum(1/probacF,na.rm=TRUE)
	return(estim)
	}

NbSim <- NbSimulations 

Heck_tirages <- as.numeric(lapply(1:NbSim,
	function (x) EstimeHeckman1EtapeSim1(rev_men,taille_pop_men,nprobrenforce,x1,x2,x3)))

sprintf("Heckman corrected mean  = %3.3f +/- %3.4f (1 ecty)",mean(Heck_tirages),sd(Heck_tirages))


#5.C VARIANCE ASSOCIATED WITH THE PROCESS OF POPULATION GENERATION + SELF-SELECTION ERROR, 
# AND HECKMAN CORRECTION IN 1 STEP BY RE-WEIGHTING

EstimeHeckman1EtapeSim2 <- function(taille_pop_menF,nprobrenforceF) {
	populationF <- generepopulation(taille_pop_menF)
	x1F <- populationF$x1
	x2F <- populationF$x2
	x3F <- populationF$x3
	rev_menF <- populationF$rev_men	
	selF <- autoselection(rev_menF,taille_pop_menF,nprobrenforceF,-0.4,30,0.2,0.2)
	selectionneF <- selF$selection
	wwF <- selF$w
	heck1EtF <- selection(selectionneF ~ x1F+x2F+x3F+wwF,
			rev_menF ~  x1F+x2F+x3F, 
			data=data.frame(c(rev_menF,selectionneF,x1F,x2F,x3F,wwF)),
			method="ml")
	
	pp1 <- predict(heck1EtF,part="selection",type="link")
	rr1 <- residuals(heck1EtF,part="outcome")
	probaF <- pnorm((pp1+coef(heck1EtF)["rho"]/coef(heck1EtF)["sigma"]*rr1)/
	((1-coef(heck1EtF)["rho"]^2)^(0.5)))
	probacF <- probaF
	probacF[probacF<0.1] <- 0.1
	estim <- sum(rev_menF/probacF,na.rm=TRUE)/sum(1/probacF,na.rm=TRUE)
	HT <- mean(rev_menF)
	return(c(estim, HT))
	}

NbSim <- NbSimulations 

Heck_tirages <- sapply(1:NbSim,
	function (x) EstimeHeckman1EtapeSim2(taille_pop_men,nprobrenforce))
sprintf("Heckman corrected mean = %3.3f +/- %3.4f (1 ecty)",mean(Heck_tirages[1,]),sd(Heck_tirages[1,]))
sprintf("Horvitz-Thompson mean  = %3.3f +/- %3.4f (1 ecty)",mean(Heck_tirages[2,]),sd(Heck_tirages[2,]))


#5.D VARIANCE ASSOCIATED WITH THE PROCESS OF POPULATION GENERATION + SELF-SELECTION ERROR, AND HECKMAN 
# CORRECTION IN 1 STEP BY IMPUTATION

EstimeHeckman1EtapeSim3 <- function(taille_pop_menF,nprobrenforceF) {
	populationF <- generepopulation(taille_pop_menF)
	x1F <- populationF$x1
	x2F <- populationF$x2
	x3F <- populationF$x3
	rev_menF <- populationF$rev_men	
	selF <- autoselection(rev_menF,taille_pop_menF,nprobrenforceF,-0.4,30,0.2,0.2)
	selectionneF <- selF$selection
	wwF <- selF$w
	heck1EtF <- selection(selectionneF ~ x1F+x2F+x3F+wwF,
			rev_menF ~  x1F+x2F+x3F, 
			data=data.frame(c(rev_menF,selectionneF,x1F,x2F,x3F,wwF)),
			method="ml")
	
	rev_men_corF <- rev_menF
	vvF <- predict(heck1EtF,part="outcome",type="conditional",
		newdata=data.frame(x1F,x2F,x3F,wwF))
	rev_men_corF[as.logical(1-selectionneF)] <- vvF[as.logical(1-selectionneF),"E[yo|ys=0]"]
	HT <- mean(rev_men_corF)
	return(HT)
	}

NbSim <- NbSimulations 

Heck_tirages <- sapply(1:NbSim,function (x) 
				EstimeHeckman1EtapeSim3(taille_pop_men,nprobrenforce))
sprintf("Heckman corrected mean by imputation = %3.3f +/- %3.4f (1 ecty)",mean(Heck_tirages),sd(Heck_tirages))


#5.E VARIANCE ASSOCIATED WITH THE POPULATION GENERATION PROCESS + SELF-SELECTION ALEA AND HECKMAN CORRECTION 
# IN 2 STEPS BY IMPUTATION

EstimeHeckman2EtapeSim1 <- function(taille_pop_menF,nprobrenforceF) {
	populationF <- generepopulation(taille_pop_menF)
	x1F <- populationF$x1
	x2F <- populationF$x2
	x3F <- populationF$x3
	rev_menF <- populationF$rev_men	
	selF <- autoselection(rev_menF,taille_pop_menF,nprobrenforceF,-0.4,30,0.2,0.2)
	selectionneF <- selF$selection
	wwF <- selF$w
	heck2EtF <- selection(selectionneF ~ x1F+x2F+x3F+wwF,
			rev_menF ~  x1F+x2F+x3F, 
			data=data.frame(c(rev_menF,selectionneF,x1F,x2F,x3F,wwF)),
			method="2step")
	
	rev_men_corF <- rev_menF
	vvF <- predict(heck2EtF,part="outcome",type="conditional",
		newdata=data.frame(x1F,x2F,x3F,wwF))
	rev_men_corF[as.logical(1-selectionneF)] <- vvF[as.logical(1-selectionneF),"E[yo|ys=0]"]
	HT <- mean(rev_men_corF)
	return(HT)
	}

NbSim <- NbSimulations 

Heck_tirages <- sapply(1:NbSim,function (x) 
				EstimeHeckman2EtapeSim1(taille_pop_men,nprobrenforce))
sprintf("Heckman corrected mean by imputation = %3.3f +/- %3.4f (1 ecty)",mean(Heck_tirages),sd(Heck_tirages))


#5.F VARIANCE ASSOCIATED WITH THE PROCESS OF POPULATION GENERATION + SELF-SELECTION ALEA, AND 
# HECKMAN CORRECTION IN 2 STEPS BY RE-WEIGHTING

EstimeHeckman2EtapeSim2 <- function(taille_pop_menF,nprobrenforceF) {
	populationF <- generepopulation(taille_pop_menF)
	x1F <- populationF$x1
	x2F <- populationF$x2
	x3F <- populationF$x3
	rev_menF <- populationF$rev_men	
	selF <- autoselection(rev_menF,taille_pop_menF,nprobrenforceF,-0.4,30,0.2,0.2)
	selectionneF <- selF$selection
	wwF <- selF$w
	heck2EtF <- selection(selectionneF ~ x1F+x2F+x3F+wwF,
			rev_menF ~  x1F+x2F+x3F, 
			data=data.frame(c(rev_menF,selectionneF,x1F,x2F,x3F,wwF)),
			method="2step")
	
	pp1 <- predict(heck2EtF,part="selection",type="link")
	rr1 <- model.frame(heck2EtF)$rev_menF-predict(heck2EtF,part="outcome",type="unconditional")
	probaF <- pnorm((pp1+coef(heck2EtF)["rho"]/coef(heck2EtF)["sigma"]*rr1)/
	((1-coef(heck2EtF)["rho"]^2)^(0.5)))
	probacF <- probaF
	probacF[probacF<0.1] <- 0.1
	estim <- sum(rev_menF/probacF,na.rm=TRUE)/sum(1/probacF,na.rm=TRUE)
	HT <- mean(rev_menF)
	return(c(estim, HT))
	}

NbSim <- NbSimulations 

Heck_tirages <- sapply(1:NbSim,
	function (x) EstimeHeckman2EtapeSim2(taille_pop_men,nprobrenforce))
sprintf("Heckman corrected mean  = %3.3f +/- %3.4f (1 ecty)",mean(Heck_tirages[1,]),sd(Heck_tirages[1,]))
sprintf("Horvitz-Thompson mean = %3.3f +/- %3.4f (1 ecty)",mean(Heck_tirages[2,]),sd(Heck_tirages[2,]))


                                                                                  
#==========================================================================
#6. Bootstrap on the selected population -------------------------------------
#==========================================================================


#population <- generepopulation(taille_pop_men)
#sel <- autoselection(population$rev_men,taille_pop_men,3000,-0.4,30,0.2,0.2)

#6.0- SOME REFERENCE RESULTS BY BOOTSTRA

#Boostrap loop to calibrate the variance of the mean estimator
SBootHV2 <- as.numeric(lapply(1:NbBouclesbootstrap ,function (x) mean(sample(rev_men,size=length(rev_men),replace=TRUE))))
sprintf("Boostrap estimate of HV  = %3.3f +/- %3.4f (1 ecty)",mean(SBootHV2),sd(SBootHV2))

#Boostrap loop to compute the variance of the mean DGP estimator
SBootRepondants <- as.numeric(lapply(1:NbBouclesbootstrap ,function (x) mean(sample(revenus_selectionnes,size=length(revenus_selectionnes),replace=TRUE))))
sprintf("Boostrap estimate for respondents  = %3.3f +/- %3.4f (1 ecty)",mean(SBootRepondants),sd(SBootRepondants))


#6.A BOOTSTRAP VARIANCE OF HECKMAN CORRECTION IN 1 STEP BY REWEIGHTING

EstimeHeckman1EtapeBoot3 <- function(taille_pop_menF,selectionFF,x1FF,x2FF,x3FF,wwFF,rev_menFF) {
	echB <- sample(1:taille_pop_menF,size=taille_pop_menF,replace=TRUE)

	x1F <- x1FF[echB]
	x2F <- x2FF[echB]
	x3F <- x3FF[echB]
	rev_menF <- rev_menFF[echB]
	selectionneF <- selectionFF[echB]
	wwF <- wwFF[echB]
	heck1EtF <- selection(selectionneF ~ x1F+x2F+x3F+wwF,
			rev_menF ~  x1F+x2F+x3F, 
			data=data.frame(c(rev_menF,selectionneF,x1F,x2F,x3F,wwF)),
			method="ml")
	
	pp1 <- predict(heck1EtF,part="selection",type="link")
	rr1 <- residuals(heck1EtF,part="outcome")
	probaF <- pnorm((pp1+coef(heck1EtF)["rho"]/coef(heck1EtF)["sigma"]*rr1)/
	((1-coef(heck1EtF)["rho"]^2)^(0.5)))
	probacF <- probaF
	probacF[probacF<0.1] <- 0.1
	estim <- sum(rev_menF/probacF,na.rm=TRUE)/sum(1/probacF,na.rm=TRUE)
	return(estim)
	}

Nbboot <- NbBouclesbootstrap 

ttz0 <- now()
Heck_tirages <- as.numeric(lapply(1:Nbboot,
	function (x) EstimeHeckman1EtapeBoot3(taille_pop_men,sel$selectionne,
	population$x1,population$x2,population$x3,sel$w,population$rev_men)))
ttz1 <- now()
sprintf("Computation time = %s", as.duration(round(ttz1-ttz0,digits=2)))
sprintf("Heckman 1S corrected  mean BOOT reweighting  = %3.3f +/- %3.4f (1 ecty)",mean(Heck_tirages),sd(Heck_tirages))


#6.B BOOTSTRAP VARIANCE OF HECKMAN CORRECTION IN 1 STEP BY IMPUTATION

EstimeHeckman1EtapeBoot4 <- function(taille_pop_menF,selectionFF,x1FF,x2FF,x3FF,wwFF,rev_menFF) {
	echB <- sample(1:taille_pop_menF,size=taille_pop_menF,replace=TRUE)

	x1F <- x1FF[echB]
	x2F <- x2FF[echB]
	x3F <- x3FF[echB]
	rev_menF <- rev_menFF[echB]
	selectionneF <- selectionFF[echB]
	wwF <- wwFF[echB]
	heck1EtF <- selection(selectionneF ~ x1F+x2F+x3F+wwF,
			rev_menF ~  x1F+x2F+x3F, 
			data=data.frame(c(rev_menF,selectionneF,x1F,x2F,x3F,wwF)),
			method="ml")
	
	rev_men_corF <- rev_menF
	vvF <- predict(heck1EtF,part="outcome",type="conditional",
		newdata=data.frame(x1F,x2F,x3F,wwF))
	rev_men_corF[as.logical(1-selectionneF)] <- vvF[as.logical(1-selectionneF),"E[yo|ys=0]"]
	HT <- mean(rev_men_corF)

	return(HT)
	}

Nbboot <- NbBouclesbootstrap 

ttz0 <- now()
Heck_tirages <- as.numeric(lapply(1:Nbboot,
	function (x) EstimeHeckman1EtapeBoot4(taille_pop_men,sel$selectionne,
	population$x1,population$x2,population$x3,sel$w,population$rev_men)))
ttz1 <- now()
sprintf("Computation time = %s", as.duration(round(ttz1-ttz0,digits=2)))
sprintf("Heckman 1S corrected  mean BOOT imputation  = %3.3f +/- %3.4f (1 ecty)",mean(Heck_tirages),sd(Heck_tirages))


#6.C BOOTSTRAP VARIANCE OF HECKMAN CORRECTION IN 2 STEPS BY REWEIGHTING

EstimeHeckman1EtapeBoot5 <- function(taille_pop_menF,selectionFF,x1FF,x2FF,x3FF,wwFF,rev_menFF) {
	echB <- sample(1:taille_pop_menF,size=taille_pop_menF,replace=TRUE)

	x1F <- x1FF[echB]
	x2F <- x2FF[echB]
	x3F <- x3FF[echB]
	rev_menF <- rev_menFF[echB]
	selectionneF <- selectionFF[echB]
	wwF <- wwFF[echB]
	heck2EtF <- selection(selectionneF ~ x1F+x2F+x3F+wwF,
			rev_menF ~  x1F+x2F+x3F, 
			data=data.frame(c(rev_menF,selectionneF,x1F,x2F,x3F,wwF)),
			method="2step")
	
	pp1 <- predict(heck2EtF,part="selection",type="link")
	#rr1 <- residuals(heck2EtF,part="outcome")
	rr1 <- model.frame(heck2EtF)$rev_menF-predict(heck2EtF,part="outcome",type="unconditional")

	probaF <- pnorm((pp1+coef(heck2EtF)["rho"]/coef(heck2EtF)["sigma"]*rr1)/
	((1-coef(heck2EtF)["rho"]^2)^(0.5)))
	probacF <- probaF
	probacF[probacF<0.1] <- 0.1
	estim <- sum(rev_menF/probacF,na.rm=TRUE)/sum(1/probacF,na.rm=TRUE)
	return(estim)
	}

Nbboot <- NbBouclesbootstrap 

ttz0 <- now()
Heck_tirages <- as.numeric(lapply(1:Nbboot,
	function (x) EstimeHeckman1EtapeBoot5(taille_pop_men,sel$selectionne,
	population$x1,population$x2,population$x3,sel$w,population$rev_men)))
ttz1 <- now()
sprintf("Computation time = %s", as.duration(round(ttz1-ttz0,digits=2)))
sprintf("Heckman 2S corrected  mean BOOT reweighting  = %3.3f +/- %3.4f (1 ecty)",mean(Heck_tirages),sd(Heck_tirages))


#6.B VARIANCE BOOTSTRAP DE CORRESTION D'HECKMAN EN 2 ETAPES PAR IMPUTATION

EstimeHeckman1EtapeBoot4 <- function(taille_pop_menF,selectionFF,x1FF,x2FF,x3FF,wwFF,rev_menFF) {
	echB <- sample(1:taille_pop_menF,size=taille_pop_menF,replace=TRUE)

	x1F <- x1FF[echB]
	x2F <- x2FF[echB]
	x3F <- x3FF[echB]
	rev_menF <- rev_menFF[echB]
	selectionneF <- selectionFF[echB]
	wwF <- wwFF[echB]
	heck2EtF <- selection(selectionneF ~ x1F+x2F+x3F+wwF,
			rev_menF ~  x1F+x2F+x3F, 
			data=data.frame(c(rev_menF,selectionneF,x1F,x2F,x3F,wwF)),
			method="2step")
	
	rev_men_corF <- rev_menF
	vvF <- predict(heck2EtF,part="outcome",type="conditional",
		newdata=data.frame(x1F,x2F,x3F,wwF))
	rev_men_corF[as.logical(1-selectionneF)] <- vvF[as.logical(1-selectionneF),"E[yo|ys=0]"]
	HT <- mean(rev_men_corF)

	return(HT)
	}

Nbboot <- NbBouclesbootstrap 

ttz0 <- now()
Heck_tirages <- as.numeric(lapply(1:Nbboot,
	function (x) EstimeHeckman1EtapeBoot4(taille_pop_men,sel$selectionne,
	population$x1,population$x2,population$x3,sel$w,population$rev_men)))
ttz1 <- now()
sprintf("Computation time = %s", as.duration(round(ttz1-ttz0,digits=2)))
sprintf("Heckman 2S corrected  mean BOOT imputation 1S  = %3.3f +/- %3.4f (1 ecty)",mean(Heck_tirages),sd(Heck_tirages))



tt1 <- now()
sprintf("Computation time = %s", as.duration(round(tt1-tt0,digits=2)))




