####################################
####MANDAI ANALYSIS FOR MS 1 #######
######Changes over 5 years##########
####################################

#KY: setwd("C:/Dropbox/Mandai for publication/1st MS/Data/")
#HR: setwd("C:/Users/A0135555/Dropbox/Botany Lab/Mandai MS/1st MS/Data")

####load libraries

library(Hmisc)
library(car) #for Type II ANOVA
library(lme4)
library(multcomp) #for generalized linear Tukey contrasts
library(aod) #for vcov-corrected wald test in underdispersed models
library(vegan)
library(FD)

##Alex's graphing function, modified by Kwek Yan

draw.plot.col=c("green", "darkgreen", "orange", "red")
draw.plot.pch=c(2, 17, 1, 16)
draw.plot.lty=c(2, 2, 1, 1)

draw.plot<-function(var, year, yt, labely,
	main=NULL, offset=0.05, ylim.drawplot="n"){
  std.e<-function(x) sqrt(var(x)/length(x))
  
  sum.mean<-aggregate(var~year+yt, FUN=mean)
  sum.se<-aggregate(var~year+yt, FUN=std.e)
  sum.sum<-data.frame(year=sum.mean$year, yt=sum.mean$yt, mean=sum.mean$var, se=sum.se$var)

  ylim.drawplot<-
	if(!is.numeric(ylim.drawplot)) with(sum.sum, c(min(mean)-max(se), max(mean)+max(se))) else ylim.drawplot

  with(sum.sum, plot(year, mean, type="n", main=main, xaxp  = c(2011, 2015, 4), ylim=ylim.drawplot, xlim=c(2011-3*offset, 2015+3*offset), ylab=labely, xlab="Year"))
  with(sum.sum[sum.sum$yt=="YU",], errbar(year-3*offset, mean, mean+se, mean-se, add=T, cap=.01, cex=1.5,  pch=draw.plot.pch[1], col=draw.plot.col[1]))
  with(sum.sum[sum.sum$yt=="YU",], lines(year-3*offset, mean, cex=1.5, lwd=1.4, pch=draw.plot.pch[1], col=draw.plot.col[1], lty=draw.plot.lty[1]))
  with(sum.sum[sum.sum$yt=="OU",], errbar(year-offset, mean, mean+se, mean-se, add=T, cap=.01, cex=1.5,  pch=draw.plot.pch[2], col=draw.plot.col[2]))
  with(sum.sum[sum.sum$yt=="OU",], lines(year-offset, mean, cex=1.5, lwd=1.4, pch=draw.plot.pch[2], col=draw.plot.col[2], lty=draw.plot.lty[2]))
  with(sum.sum[sum.sum$yt=="YA",], errbar(year+offset, mean, mean+se, mean-se, add=T, cap=.01, cex=1.5,  pch=draw.plot.pch[3], col=draw.plot.col[3]))
  with(sum.sum[sum.sum$yt=="YA",], lines(year+offset, mean, cex=1.5, lwd=1.4, pch=draw.plot.pch[3], col=draw.plot.col[3], lty=draw.plot.lty[3]))
  with(sum.sum[sum.sum$yt=="OA",], errbar(year+3*offset, mean, mean+se, mean-se, add=T, cap=.01, cex=1.5,  pch=draw.plot.pch[4], col=draw.plot.col[4]))
  with(sum.sum[sum.sum$yt=="OA",], lines(year+3*offset, mean, cex=1.5, lwd=1.4, pch=draw.plot.pch[4], col=draw.plot.col[4], lty=draw.plot.lty[4]))

  return(sum.sum)
}

##Chiu & Chao (2014)'s function for functional diversity based on Hill's numbers
#http://chao.stat.nthu.edu.tw/wordpress/paper/functional%20diversity.r
#Accessed 5 May 2017

Func2014=function (Dis,abun,q)  
  # input:
  # Dis: species pairwise functional distance matrix.
  # abun: one assemblage species abundance vector data or multi-assemblage
  #species by assemblages matrix data, where row is species and column is 
  # community or site(plot).
  # q: a non-negative value specifying diversity order.
  
  #output:
  # Q: Rao's quadratic entropy for each community.
  # FuncD= functional diversity of each community. 
  # Gamma= functional gamma diversity.
# Alpha= functional alpha diversity.
# Beta= functional beta diversity.
# FunCqN = functional overlap index CqN(similarity index) 
# FunUqN= functional overlap index UqN(similarity index)
{
  if(is.vector(abun)){
    Dis=as.matrix(Dis); 
    n=sum(abun);
    p=abun/n;I=which(p>0);p=p[I];
    Q=c(t(p)%*%Dis[I,I]%*%p);
    temp=p%*%t(p);
    if(q==1){FD= exp(sum(-Dis[I,I]*temp/Q*log(temp/Q)))}
    else{FD=(t(p^q)%*%Dis[I,I]%*%(p/Q)^q)^(1/(1-q));}
    #output=matrix(ncol=2,nrow=1);
    output=c(Q,FD)
    names(output)=c("Q","FuncD")
    return(output);
  }else{
    abun=as.matrix(abun);  
    N=ncol(abun);n=colSums(abun)
    FuncD=numeric(N);Q=numeric(N)
    Dis=as.matrix(Dis); 
    for(i in 1:N){            
      Q[i]=t(abun[,i]/n[i])%*%Dis%*%(abun[,i]/n[i]);
      temp=(abun[,i]/n[i])%*%t(abun[,i]/n[i]);
      I=which(temp>0);
      
      if(q==1){ FuncD[i]= exp(sum(-Dis[I]*temp[I]/Q[i]*log(temp[I]/Q[i]))); }
      else{ FuncD[i]=sum(Dis[I]*(temp[I]/Q[i])^q)^(1/(1-q)); }
    }
    
    gn=sum(abun);
    pop=abun/gn;  
    p=rowSums(pop);gI=which(p>0);p=p[gI];
    gQ=c(t(p)%*%Dis[gI,gI]%*%p);
    
    atemp=0;      
    for(i in 1:N){
      for(j in 1:N){
        pi=pop[,i];pj=pop[,j]                  
        pij=pi%*%t(pj);aI=which(pij>0);
        if(q==1){ atemp=atemp+ sum(-Dis[aI]*pij[aI]/gQ*log(pij[I]/gQ))}
        else{atemp=atemp+sum(Dis[aI]*(pij[aI]/gQ)^q) }
      }	
    }
    
    if(q==1){
      gFD= exp(sum(-Dis[gI,gI]*(p%*%t(p))/gQ*log( p%*%t(p)/gQ )))
      aFD= exp(atemp)/N^2
      bFD=gFD/aFD
      CqN= 1-log(bFD)/log(N^2);
      UqN=CqN;      
    }else{
      gFD= (t(p^q)%*%Dis[gI,gI]%*%(p/gQ)^q)^(1/(1-q));
      aFD= atemp^(1/(1-q))/N^2
      bFD=gFD/aFD;            
      CqN=1-(bFD^(1-q)-1)/(N^(2-2*q)-1);
      UqN=1-(bFD^(q-1)-1)/(N^(2*q-2)-1);
    }
    output=c(gQ,gFD,aFD,bFD,CqN,UqN)
    names(output)=c("Q","Gamma","Alpha","Beta","FunCqN","FunUqN")
    return( list(Q=Q,FuncD=FuncD,output))               
  }
}

##Monte Carlo permutation t-test to assess differences for each year

MCttest<-function(y, x, subset, nperm=999) {
	x<-x[subset]
	y<-y[subset]
	tobs<-t.test(y~x)$statistic
	tperm<-rep(NA, nperm)
	for (i in 1:nperm){
		xi<-sample(x, replace=FALSE)
		tperm[i]<-t.test(y~xi)$statistic
		}
	return(sum(tperm>abs(tobs))/(nperm+1)*2)
	}


####read data

###env data
#please check mandai_data_preparation to see how envi_combi_long.csv and dbh_long is generated

envi.long<-read.csv("envi_combi_long.csv")

envi.long$yt<-paste0(substr(envi.long$type,1,1), substr(envi.long$damage,1,1))
envi.long$yearf<-factor(envi.long$year)
envi.long$damage<-factor(envi.long$damage, levels=c("Unaffected", "Affected"))


#####################################
##Re-analyses with special conditions
  #SKIP IF ALL PLOTS RETAINED AS-IS

#if A27 and B10 were excluded because they were affected in 2012 and 13, respectively

envi<-read.csv("envi_5years.csv", row.names=1)
newly.affected<-envi$plot[!envi$damage_2011==envi$damage_2015]
unaffected.plots<-envi$plot[envi$damage_2011=="Unaffected"]
always.unaffected<-unaffected.plots[!(unaffected.plots %in% newly.affected)]

envi.long<-envi.long[!envi.long$plot %in% newly.affected,]

#if C2 was young instead of old

envi.long$type[envi.long$plot=="C2"]<-"Young"

#if C2 was excluded

envi.long<-envi.long[envi.long$plot!="C2",]

######################################


###community data
dbh<-read.csv("dbh_long.csv")
dbh$species<-sub("cf. ", "", dbh$species)
dbh$species<-sub("Aporosa frutescens ", "Aporosa frutescens", dbh$species)
dbh$y.plot<-paste0(dbh$plot,"-",substr(dbh$year, 3, 4))
dbh$presence<-ifelse(dbh$dbh>=1.0, 1, 0)

tree.matrix<-xtabs(presence~y.plot+species, data=dbh[dbh$plot %in% envi.long$plot,])
tree.matrix<-tree.matrix[,dimnames(tree.matrix)$species!="Unknown"]
tree.matrix<-matrix(tree.matrix, ncol=ncol(tree.matrix), nrow=nrow(tree.matrix), dimnames=dimnames(tree.matrix))

D0<-specnumber(tree.matrix)
D1<-exp(diversity(tree.matrix, index="shannon"))
D2<-diversity(tree.matrix, index="invsimpson")

###functional traits

trait <- read.csv("mandai_trait_compiled.csv", row.names=1)
trait <- trait[order(rownames(trait)),]
trait <- trait[rowSums(is.na(trait))<3,]

# dbh$has.trait <- dbh$species %in% rownames(trait.hastree)
# prop.table(with(dbh, tapply(UID, has.trait, function(x) length(unique(x)))))
# length(rownames(trait.hastree)) / length(unique(dbh$species))

#par(mfrow=c(3,2))
#for (i in 1:ncol(trait)) hist(trait[,i])

trait[,c(1,3,5)]<-log(trait[,c(1,3,5)])

tree.matrix.hastrait<-tree.matrix[,dimnames(tree.matrix)$species %in% rownames(trait)]
# write.csv(tree.matrix.hastrait / rowSums(tree.matrix.hastrait), "tree_matrix_hastrait.csv")

trait.hastree<-trait[rownames(trait) %in% dimnames(tree.matrix)$species,]
# write.csv(trait.hastree, "trait_has_tree.csv")

#FD<-dbFD(x=trait.hastree, a=tree.matrix.hastrait, corr="cailliez") #originally just for CWM
	#but buggy function, and I found another one to use specifically for CWM, functcomp()
	#used in the PRC section

trait.hastree.dist<-gowdis(trait.hastree)

FD0<-Func2014(trait.hastree.dist, t(tree.matrix.hastrait), 0)
FD1<-Func2014(trait.hastree.dist, t(tree.matrix.hastrait), 1)
FD2<-Func2014(trait.hastree.dist, t(tree.matrix.hastrait), 2)


####plot maps

wind<-read.csv("wind.csv", header=TRUE)
wind$date1+wind$time1->wind$time2
wind$time1*24->wind$time3

wind11<-wind[wind$date1==11,]

stations<-c("S100","S104","S122","S106","S103","S50","S109","S24","S44","S86","S111","S107","S116","S102","S60","S108")

par(mfrow=c(4,4), mar=c(2,2,0.5,0.5), mgp=c(2,1,0))
for (i in 1:length(stations)){
	plot(speed_max~time3, data=wind11, subset=station==stations[i], type="n", ylim=c(0,25), ylab="", xlab="")
	lines(speed_max~time3, data=wind11, subset=station==stations[i], col="red")
	lines(speed_avg~time3, data=wind11, subset=station==stations[i], col="gray")
	text(2.5, 22.5, LETTERS[i], font=1, cex=2)
	}

#19 Apr 2017 note from KY: excluded for now to concentrate on analysis



####Abiotic changes

Anova(N.lme<-lmer(N~type*damage*yearf+(1|plot), data=envi.long))
N.lme.refit<-lmer(N~yearf+(1|plot), data=envi.long)
N.glht<-summary(glht(N.lme.refit, linfct=mcp(yearf="Tukey")))

Anova(P.lme<-lmer(P~type*damage*yearf+(1|plot), data=envi.long))
P.lme.refit<-lmer(P~yearf+(1|plot), data=envi.long)
P.glht<-summary(glht(P.lme.refit, linfct=mcp(yearf="Tukey")))

Anova(K.lme<-lmer(K~type*damage*yearf+(1|plot), data=envi.long))
K.lme.refit<-lmer(K~yearf+(1|plot), data=envi.long)
K.glht<-summary(glht(K.lme.refit, linfct=mcp(yearf="Tukey")))

envi.long$CC.logit<-with(envi.long, log((CC/100)/(1-CC/100)))
Anova(CC.lme<-lmer(CC.logit~type*damage*yearf+(1|plot), data=envi.long))

damage.year.inter<-with(envi.long, interaction(damage, yearf))

CC.lme.refit<-lmer(CC.logit~damage.year.inter-1+(1|plot), data=envi.long)

damage.year.linfct<-glht(CC.lme.refit, linfct=mcp(damage.year.inter="Tukey"))$linfct[c(1,2,11,18,19,26,31,32,37,40,41,44,45),]

CC.glht<-summary(glht(CC.lme.refit, linfct=damage.year.linfct), test=adjusted("bonferroni"))

CC.glht.tab<-data.frame(diff=CC.glht$test$coef, se=CC.glht$test$sigma, z=CC.glht$test$tstat)
CC.glht.tab$pval<-replace(CC.glht$test$pvalues, CC.glht$test$pvalues<0.001, "<0.001")

write.csv(CC.glht.tab, "CCglhttab.csv", quote=FALSE, row.names=TRUE)

Anova(LL.lme<-lmer(LLavg~type*damage*yearf+(1|plot), data=envi.long))
LL.lme.refit<-lmer(LLavg~yearf+(1|plot), data=envi.long)
LL.glht<-summary(glht(LL.lme.refit, linfct=mcp(yearf="Tukey")))

envi.long$WD.logit<-with(envi.long, log((WD/100)/(1-WD/100)))
Anova(WD.lme<-lmer(WD.logit~type*damage*yearf+(1|plot), data=envi.long))
WD.lme.refit<-lmer(WD.logit~damage.year.inter-1+(1|plot), data=envi.long)
WD.glht<-summary(glht(WD.lme.refit, linfct=damage.year.linfct), test=adjusted("bonferroni"))
WD.glht.tab<-data.frame(diff=WD.glht$test$coef, se=WD.glht$test$sigma, z=WD.glht$test$tstat)
WD.glht.tab$pval<-replace(WD.glht$test$pvalues, WD.glht$test$pvalues<0.001, "<0.001")

write.csv(WD.glht.tab, "WDglhttab.csv", quote=FALSE, row.names=TRUE)

Anova(basal.lme<-lmer(basal~type*damage*yearf+(1|plot), data=envi.long))
basal.lme.refit<-lmer(basal~damage.year.inter-1+(1|plot), data=envi.long)
basal.glht<-summary(glht(basal.lme.refit, linfct=damage.year.linfct), test=adjusted("bonferroni"))
basal.glht.tab<-data.frame(diff=basal.glht$test$coef, se=basal.glht$test$sigma, z=basal.glht$test$tstat)
basal.glht.tab$pval<-replace(basal.glht$test$pvalues, basal.glht$test$pvalues<0.001, "<0.001")

write.csv(basal.glht.tab, "basalglhttab.csv", quote=FALSE, row.names=TRUE)

Anova(stem.d.lme<-glmer(stem.d~type*damage*yearf+(1|plot), data=envi.long, control=glmerControl(optimizer="bobyqa"), family=poisson))
(stem.d.chat<-sum(residuals(stem.d.lme)^2)/(nrow(envi.long)-length(fixef(stem.d.lme))-1))
#overdispersed, fit negative binomial instead
Anova(stem.d.lme<-glmer.nb(stem.d~type*damage*yearf+(1|plot), data=envi.long, control=glmerControl(optimizer="bobyqa")))
type.damage.year.inter<-with(envi.long, interaction(type, damage, yearf))
stem.d.lme.refit<-glmer.nb(stem.d~type.damage.year.inter-1+(1|plot), data=envi.long, control=glmerControl(optimizer="bobyqa"))

type.damage.year.linfct<-glht(stem.d.lme.refit, linfct=mcp(type.damage.year.inter="Tukey"))$linfct[c(2,4,21,23,41,58,72,74,87,89,103,116,126,128,137,139,149,158,164,166,171,173,179,184,186,189),]

stem.glht<-summary(glht(stem.d.lme.refit, linfct=type.damage.year.linfct), test=adjusted("bonferroni"))

stem.glht.tab<-data.frame(diff=stem.glht$test$coef, se=stem.glht$test$sigma, z=stem.glht$test$tstat)
stem.glht.tab$pval<-replace(stem.glht$test$pvalues, stem.glht$test$pvalues<0.001, "<0.001")

write.csv(stem.glht.tab, "stemglhttab.csv", quote=FALSE, row.names=TRUE)

ANOVA.results<-cbind(Variable=c(rep("N", 7),rep("P", 7),rep("K", 7),rep("Canopy cover", 7),rep("Leaf litter", 7),rep("Woody debris", 7),rep("Total basal area", 7),rep("Stem count", 7)),
	Term=rep(rownames(Anova(N.lme)),8),
	data.frame(rbind(
		Anova(N.lme), Anova(P.lme), Anova(K.lme), Anova(CC.lme), Anova(LL.lme), Anova(WD.lme), Anova(basal.lme), Anova(stem.d.lme)
		))
	)
rownames(ANOVA.results)<-NULL
ANOVA.results$Pr..Chisq.<-replace(ANOVA.results$Pr..Chisq.,ANOVA.results$Pr..Chisq.<0.001, "<0.001")

write.csv(ANOVA.results, "ANOVA_results.csv", quote=FALSE, row.names=FALSE)

NPKLL.pval<-c(N.glht$test$pvalues, P.glht$test$pvalues, K.glht$test$pvalues, LL.glht$test$pvalues)

write.csv(
	data.frame(var=c(rep("N", length(N.glht$test$coef)), rep("P", length(N.glht$test$coef)), rep("K", length(N.glht$test$coef)), rep("CC", length(N.glht$test$coef))),
		comp=names(N.glht$test$coef),
		diff=c(N.glht$test$coef, P.glht$test$coef, K.glht$test$coef, LL.glht$test$coef),
		se=c(N.glht$test$sigma, P.glht$test$sigma, K.glht$test$sigma, LL.glht$test$sigma),
		z=c(N.glht$test$tstat, P.glht$test$tstat, K.glht$test$tstat, LL.glht$test$tstat),
		pval=replace(NPKLL.pval, NPKLL.pval<0.001, "<0.001")
		), "NPKLLglhttab.csv", quote=FALSE, row.names=FALSE
	)

png("C:/Users/A0135555/Dropbox/Botany Lab/Mandai MS/1st MS/figs/1_env_dynamics.png", 180, 180, "mm", res=600)
par(mfrow=c(3,3), mar=c(3,3,1,1), mgp=c(2,0.8,0))
N.plot<-with(envi.long, draw.plot(var=N, year= year, yt=yt, 
	labely="Total N (mg/kg)"))
#text(c(2011.4, 2012.2, 2013.2, 2014, 2015), with(N.plot, tapply(mean, year, max))+c(-280,200,200,340,320),
#	labels=c("ac", "b", "c", "bc", "bc"), xpd=TRUE)
text(par()$usr[1]-0.18*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(a)", xpd=TRUE)

P.plot<-with(envi.long, draw.plot(var=P, year= year, yt=yt,
	labely="Extractable P (mg/kg)"))
#text(c(2011.4, 2012.2, 2013.3, 2013.9, 2015.1), with(P.plot, tapply(mean, year, max))+c(-2, 0, -2, -2.5, -3),
#	labels=c("a", "b", "a", "a", "a"), xpd=TRUE)
text(par()$usr[1]-0.18*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(b)", xpd=TRUE)

K.plot<-with(envi.long, draw.plot(var=K, year= year, yt=yt,
	labely="Extractable K (mg/kg)"))
#text(c(2011.4, 2012.3, 2013.2, 2013.9, 2015), with(K.plot, tapply(mean, year, max))+c(5, 5, 5, 8, 7),
#	labels=c("ac", "a", "b", "bc", "bc"), xpd=TRUE)
text(par()$usr[1]-0.18*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(c)", xpd=TRUE)

with(envi.long, draw.plot(var=CC, year= year, yt=yt,
	labely="Canopy cover (%)"))
text(par()$usr[1]-0.18*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(d)", xpd=TRUE)

with(envi.long, draw.plot(var=LLavg, year= year, yt=yt,
	labely="Leaf litter depth (cm)"))
text(par()$usr[1]-0.18*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(e)", xpd=TRUE)

with(envi.long, draw.plot(var=WD, year= year, yt=yt,
	labely="Woody debris cover (%)"))
text(par()$usr[1]-0.18*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(f)", xpd=TRUE)

with(envi.long, draw.plot(var=basal, year= year, yt=yt,
	labely="Total basal area (sq. cm)"))
text(par()$usr[1]-0.18*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(g)", xpd=TRUE)

with(envi.long, draw.plot(var=stem.d, year= year, yt=yt,
	labely="No. of stems"))
text(par()$usr[1]-0.18*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(h)", xpd=TRUE)

plot.new()
legend("topleft", legend=c("Young Unaffected", "Old Unaffected", "Young Affected", "Old Affected"),pch=draw.plot.pch, lty=draw.plot.lty, col=draw.plot.col, 
       lwd=1.4, cex=1, bty="n", y.intersp =1 )
dev.off()


####Biotic changes:
###Alpha diversity

Anova(D0.lme<-glmer(D0~type*damage*yearf+(1|plot), data=envi.long, family=poisson, control=glmerControl(optimizer="bobyqa")))
(D0.chat<-sum(residuals(D0.lme)^2)/(nrow(envi.long)-length(fixef(D0.lme))-1))
##underdispersed and tests will be over-conservative, but not considered a problem for now.
D0.lme.refit<-glmer(D0~type.damage.year.inter-1+(1|plot), data=envi.long, control=glmerControl(optimizer="bobyqa"), family=poisson)
D0.glht<-summary(glht(D0.lme.refit, linfct=type.damage.year.linfct), test=adjusted("bonferroni"))
D0.glht.tab<-data.frame(diff=D0.glht$test$coef, se=D0.glht$test$sigma, z=D0.glht$test$tstat)
D0.glht.tab$pval<-replace(D0.glht$test$pvalues, D0.glht$test$pvalues<0.001, "<0.001")

write.csv(D0.glht.tab, "D0glhttab.csv", quote=FALSE, row.names=TRUE)

Anova(D1.lme<-lmer(D1~type*damage*yearf+(1|plot), data=envi.long))
D1.lme.refit<-lmer(D1~type.damage.year.inter-1+(1|plot), data=envi.long)
D1.glht<-summary(glht(D1.lme.refit, linfct=type.damage.year.linfct), test=adjusted("bonferroni"))
D1.glht.tab<-data.frame(diff=D1.glht$test$coef, se=D1.glht$test$sigma, z=D1.glht$test$tstat)
D1.glht.tab$pval<-replace(D1.glht$test$pvalues, D1.glht$test$pvalues<0.001, "<0.001")

write.csv(D1.glht.tab, "D1glhttab.csv", quote=FALSE, row.names=TRUE)

Anova(D2.lme<-lmer(D2~type*damage*yearf+(1|plot), data=envi.long))
D2.lme.refit<-lmer(D2~type.damage.year.inter-1+(1|plot), data=envi.long)
D2.glht<-summary(glht(D2.lme.refit, linfct=type.damage.year.linfct), test=adjusted("bonferroni"))
D2.glht.tab<-data.frame(diff=D2.glht$test$coef, se=D2.glht$test$sigma, z=D2.glht$test$tstat)
D2.glht.tab$pval<-replace(D2.glht$test$pvalues, D2.glht$test$pvalues<0.001, "<0.001")

write.csv(D2.glht.tab, "D2glhttab.csv", quote=FALSE, row.names=TRUE)

Anova(FD0.lme<-lmer(FD0$FuncD~type*damage*yearf+(1|plot), data=envi.long))
FD0.lme.refit<-lmer(FD0$FuncD~type.damage.year.inter-1+(1|plot), data=envi.long)
FD0.glht<-summary(glht(FD0.lme.refit, linfct=type.damage.year.linfct), test=adjusted("bonferroni"))
FD0.glht.tab<-data.frame(diff=FD0.glht$test$coef, se=FD0.glht$test$sigma, z=FD0.glht$test$tstat)
FD0.glht.tab$pval<-replace(FD0.glht$test$pvalues, FD0.glht$test$pvalues<0.001, "<0.001")

write.csv(FD0.glht.tab, "FD0glhttab.csv", quote=FALSE, row.names=TRUE)

Anova(FD1.lme<-lmer(FD1$FuncD~type*damage*yearf+(1|plot), data=envi.long))
FD1.lme.refit<-lmer(FD1$FuncD~damage.year.inter-1+(1|plot), data=envi.long)
FD1.glht<-summary(glht(FD1.lme.refit, linfct=damage.year.linfct), test=adjusted("bonferroni"))
FD1.glht.tab<-data.frame(diff=FD1.glht$test$coef, se=FD1.glht$test$sigma, z=FD1.glht$test$tstat)
FD1.glht.tab$pval<-replace(FD1.glht$test$pvalues, FD1.glht$test$pvalues<0.001, "<0.001")

write.csv(FD1.glht.tab, "FD1glhttab.csv", quote=FALSE, row.names=TRUE)

Anova(FD2.lme<-lmer(FD2$FuncD~type*damage*yearf+(1|plot), data=envi.long))
FD2.lme.refit<-lmer(FD2$FuncD~type.damage.year.inter-1+(1|plot), data=envi.long)
FD2.glht<-summary(glht(FD2.lme.refit, linfct=type.damage.year.linfct), test=adjusted("bonferroni"))
FD2.glht.tab<-data.frame(diff=FD2.glht$test$coef, se=FD2.glht$test$sigma, z=FD2.glht$test$tstat)
FD2.glht.tab$pval<-replace(FD2.glht$test$pvalues, FD2.glht$test$pvalues<0.001, "<0.001")

write.csv(FD2.glht.tab, "FD2glhttab.csv", quote=FALSE, row.names=TRUE)

ANOVA.results.div<-cbind(Variable=c(rep("D0", 7),rep("D1", 7),rep("D2", 7),rep("FD0", 7),rep("FD1", 7),rep("FD2", 7)),
	Term=rep(rownames(Anova(D0.lme)),6),
	data.frame(rbind(
		Anova(D0.lme), Anova(D1.lme), Anova(D2.lme), Anova(FD0.lme), Anova(FD1.lme), Anova(FD2.lme)
		))
	)
rownames(ANOVA.results.div)<-NULL
ANOVA.results.div$Pr..Chisq.<-replace(ANOVA.results.div$Pr..Chisq.,ANOVA.results.div$Pr..Chisq.<0.001, "<0.001")

write.csv(ANOVA.results.div, "ANOVA_results_div.csv", quote=FALSE, row.names=FALSE)


png("2_SDFD_dynamics.png", 180, 120, "mm", res=600)
par(mfrow=c(2,3), mar=c(3,3,1,1), mgp=c(2,0.8,0), cex.axis=1, cex.lab=1.2)

with(envi.long, draw.plot(var=D0, year= year, yt=yt, ylim.drawplot=c(3,33),
	labely="SD0"))
text(par()$usr[1]-0.18*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(a)", xpd=TRUE)

with(envi.long, draw.plot(var=D1, year= year, yt=yt, ylim.drawplot=c(3,33),
	labely="SD1"))
text(par()$usr[1]-0.18*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(c)", xpd=TRUE)

with(envi.long, draw.plot(var=D2, year= year, yt=yt, ylim.drawplot=c(3,33),
	labely="SD2"))
text(par()$usr[1]-0.18*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(d)", xpd=TRUE)

with(envi.long, draw.plot(var=FD0$FuncD, year= year, yt=yt, ylim.drawplot=c(0,170),
	labely="FD0"))
text(par()$usr[1]-0.18*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(b)", xpd=TRUE)

with(envi.long, draw.plot(var=FD1$FuncD, year= year, yt=yt, ylim.drawplot=c(0,170),
	labely="FD1"))
text(par()$usr[1]-0.18*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(e)", xpd=TRUE)

with(envi.long, draw.plot(var=FD2$FuncD, year= year, yt=yt, ylim.drawplot=c(0,170),
	labely="FD2"))
text(par()$usr[1]-0.18*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(f)", xpd=TRUE)

legend("topright", legend=c("Young Unaffected", "Old Unaffected", "Young Affected", "Old Affected"),pch=draw.plot.pch, lty=draw.plot.lty, col=draw.plot.col, 
       lwd=1.4, cex=1, bty="n", y.intersp =1 )
dev.off()


####Community changes

damage.type<-factor(with(envi.long, paste(type, damage)), levels=c("Old Unaffected", "Young Unaffected", "Old Affected", "Young Affected"))

stems.affected.year<-tapply(rowSums(tree.matrix[envi.long$damage=="Affected",]), envi.long$yearf[envi.long$damage=="Affected"], sum)

imptval<-rep(NA, ncol(tree.matrix))

for (i in 1:ncol(tree.matrix)) imptval[i]<-max(tapply(tree.matrix[envi.long$damage=="Affected",i], envi.long$year[envi.long$damage=="Affected"], sum)[-1]/stems.affected.year[-1])

##Principal Response Curves (PRC) based on taxonomic identities
tax.prc<-prc(log(tree.matrix+1), damage.type, envi.long$yearf)

tax.prc.coef1<-summary(tax.prc, axis=1)$coef
tax.prc.coef1<-rbind(tax.prc.coef1, rep(0, ncol(tax.prc.coef1)))
rownames(tax.prc.coef1)[4]<-"Old Unaffected"
tax.prc.coef1<-tax.prc.coef1[c(1,4,3,2),]

tax.prc.sp1<-summary(tax.prc, axis=1)$sp

sp.select1<-order(imptval, decreasing=TRUE)[1:20]
#colSums(tree.matrix)>200

tax.prc.coef2<-summary(tax.prc, axis=2)$coef
tax.prc.coef2<-rbind(tax.prc.coef2, rep(0, ncol(tax.prc.coef2)))
rownames(tax.prc.coef2)[4]<-"Old Unaffected"
tax.prc.coef2<-tax.prc.coef2[c(1,4,3,2),]

tax.prc.sp2<-summary(tax.prc, axis=2)$sp

sp.select2<-order(imptval, decreasing=TRUE)[1:20]

tax.prc.coef3<-summary(tax.prc, axis=3)$coef
tax.prc.coef3<-rbind(tax.prc.coef3, rep(0, ncol(tax.prc.coef3)))
rownames(tax.prc.coef3)[4]<-"Old Unaffected"
tax.prc.coef3<-tax.prc.coef3[c(1,4,3,2),]

tax.prc.sp3<-summary(tax.prc, axis=3)$sp

sp.select3<-order(imptval, decreasing=TRUE)[1:20]

png("3_PRC_taxo.png", 180, 220, "mm", res=600)
par(mfrow=c(3,1), mar=c(3,3,1,15), mgp=c(2,0.8,0), cex=1)

matplot(colnames(tax.prc.coef1), t(tax.prc.coef1),
	col=draw.plot.col, pch=draw.plot.pch, 
	ylim=c(-0.4, 0.15), 
	xlim=c(2011, 2015.8),
	ylab="PRC 1", xlab="Year")
for (i in 1:nrow(tax.prc.coef1)) lines(colnames(tax.prc.coef1), tax.prc.coef1[i,], col=draw.plot.col[i], lty=draw.plot.lty[i])
matplot(colnames(tax.prc.coef1), t(tax.prc.coef1),
	col="white", pch=16, cex=2, add=TRUE)
matplot(colnames(tax.prc.coef1), t(tax.prc.coef1),
	col=draw.plot.col, pch=draw.plot.pch, add=TRUE)
linestack(tax.prc.sp1[sp.select1]/5, names(tax.prc.sp1)[sp.select1], add=TRUE, at=par()$usr[2], cex=0.7, font=3)
axis(2, at=axTicks(2), labels=axTicks(2)*5, pos=par()$usr[2])
text(par()$usr[1]-0.14*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(a)", xpd=TRUE)

matplot(colnames(tax.prc.coef2), t(tax.prc.coef2),
	col=draw.plot.col, pch=draw.plot.pch, 
	ylim=c(-0.15, 0.35), xlim=c(2011, 2015.8),
	ylab="PRC 2", xlab="Year")
for (i in 1:nrow(tax.prc.coef2)) lines(colnames(tax.prc.coef2), tax.prc.coef2[i,], col=draw.plot.col[i], lty=draw.plot.lty[i])
matplot(colnames(tax.prc.coef2), t(tax.prc.coef2),
	col="white", pch=16, cex=2, add=TRUE)
matplot(colnames(tax.prc.coef2), t(tax.prc.coef2),
	col=draw.plot.col, pch=draw.plot.pch, add=TRUE)
linestack(tax.prc.sp2[sp.select2]/5, names(tax.prc.sp2)[sp.select2], add=TRUE, at=par()$usr[2], cex=0.7, font=3)
axis(2, at=axTicks(2), labels=axTicks(2)*5, pos=par()$usr[2])
text(par()$usr[1]-0.14*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(b)", xpd=TRUE)

matplot(colnames(tax.prc.coef3), t(tax.prc.coef3),
	col=draw.plot.col, pch=draw.plot.pch, 
	ylim=c(-0.25, 0.35), xlim=c(2011, 2015.8),
	ylab="PRC 3", xlab="Year")
for (i in 1:nrow(tax.prc.coef3)) lines(colnames(tax.prc.coef3), tax.prc.coef3[i,], col=draw.plot.col[i], lty=draw.plot.lty[i])
matplot(colnames(tax.prc.coef3), t(tax.prc.coef3),
	col="white", pch=16, cex=2, add=TRUE)
matplot(colnames(tax.prc.coef3), t(tax.prc.coef3),
	col=draw.plot.col, pch=draw.plot.pch, add=TRUE)
linestack(tax.prc.sp3[sp.select3]/5, names(tax.prc.sp3)[sp.select3], add=TRUE, at=par()$usr[2], cex=0.7, font=3)
axis(2, at=axTicks(2), labels=axTicks(2)*5, pos=par()$usr[2])
text(par()$usr[1]-0.14*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(c)", xpd=TRUE)

legend("bottom", legend=c("Young Unaffected", "Old Unaffected", "Young Affected", "Old Affected"),pch=draw.plot.pch, lty=draw.plot.lty, col=draw.plot.col, 
       lwd=1.4, cex=0.7, bty="n")
dev.off()

##PRCs based on community-weighted mean functional trait values

CWM<-functcomp(x=trait.hastree, a=tree.matrix.hastrait, CWM.type="all")

CWM.scaled<-scale(CWM)

fun.prc<-prc(CWM.scaled, damage.type, envi.long$yearf)

fun.prc.coef1<-summary(fun.prc, axis=1)$coef
fun.prc.coef1<-rbind(fun.prc.coef1, rep(0, ncol(fun.prc.coef1)))
rownames(fun.prc.coef1)[4]<-"Old Unaffected"
fun.prc.coef1<-fun.prc.coef1[c(1,4,3,2),]

fun.prc.trait1<-summary(fun.prc, axis=1)$sp
names(fun.prc.trait1)[6] <- "Hmax"

fun.prc.coef2<-summary(fun.prc, axis=2)$coef
fun.prc.coef2<-rbind(fun.prc.coef2, rep(0, ncol(fun.prc.coef2)))
rownames(fun.prc.coef2)[4]<-"Old Unaffected"
fun.prc.coef2<-fun.prc.coef2[c(1,4,3,2),]

fun.prc.trait2<-summary(fun.prc, axis=2)$sp
names(fun.prc.trait1)[6] <- "Hmax"

png("4_PRC_cwm.png", 180, 150, "mm", res=600)
par(mfrow=c(2,1), mar=c(3,3,1,15), mgp=c(2,0.8,0))

matplot(colnames(fun.prc.coef1), t(fun.prc.coef1),
	col=draw.plot.col, pch=draw.plot.pch, 
	ylim=c(-0.3,0.5), xlim=c(2011,2015.8),
	ylab="PRC 1", xlab="Year")
for (i in 1:nrow(fun.prc.coef1)) lines(colnames(fun.prc.coef1), fun.prc.coef1[i,], col=draw.plot.col[i], lty=draw.plot.lty[i])
matplot(colnames(fun.prc.coef1), t(fun.prc.coef1),
	col="white", pch=16, cex=2, add=TRUE)
matplot(colnames(fun.prc.coef1), t(fun.prc.coef1),
	col=draw.plot.col, pch=draw.plot.pch, add=TRUE)
linestack(fun.prc.trait1/4, names(fun.prc.trait1), add=TRUE, at=par()$usr[2], cex=0.7)
axis(2, at=axTicks(2), labels=axTicks(2)*4, pos=par()$usr[2])
text(par()$usr[1]-0.14*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(a)", xpd=TRUE)
legend("top", legend=c("Young Unaffected", "Old Unaffected", "Young Affected", "Old Affected"),pch=draw.plot.pch, lty=draw.plot.lty, col=draw.plot.col, 
       lwd=1.4, cex=0.8, bty="n")

matplot(colnames(fun.prc.coef2), t(fun.prc.coef2),
	col=draw.plot.col, pch=draw.plot.pch,
	ylim=c(-0.74,0.5), xlim=c(2011, 2015.8),
	ylab="PRC 2", xlab="Year")
for (i in 1:nrow(fun.prc.coef2)) lines(colnames(fun.prc.coef2), fun.prc.coef2[i,], col=draw.plot.col[i], lty=draw.plot.lty[i])
matplot(colnames(fun.prc.coef2), t(fun.prc.coef2),
	col="white", pch=16, cex=2, add=TRUE)
matplot(colnames(fun.prc.coef2), t(fun.prc.coef2),
	col=draw.plot.col, pch=draw.plot.pch, add=TRUE)
linestack(fun.prc.trait2/2, names(fun.prc.trait2), add=TRUE, at=par()$usr[2], cex=0.7)
axis(2, at=axTicks(2), labels=axTicks(2)*2, pos=par()$usr[2])
text(par()$usr[1]-0.14*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(b)", xpd=TRUE)

dev.off()

##old-growth and young secondary colonizer indices

OG<-mean(tax.prc.coef1[4,2:5]-tax.prc.coef1[3,2:5])*tax.prc.sp1
YS<-mean(tax.prc.coef2[3,2:5]-tax.prc.coef2[4,2:5])*tax.prc.sp2

plot(OG, YS, col="gray", cex=100*(exp(imptval)-1), pch=16, xlab="Colonizers", ylab="Young secondary colonizers")
abline(h=0); abline(v=0)
abline(0,1, lty=2)
points(YS~OG, subset=order(imptval, decreasing=TRUE)[1:20], col="white", bg="gray40", cex=100*(exp(imptval)-1), pch=21, lwd=2)
identify(OG, YS, colnames(tree.matrix), xpd=TRUE, cex=0.7, font=3)

par(mfrow=c(2,2), mar=c(4,4,1,1), mgp=c(2.5, 1, 0))

plot(tax.prc.sp1, tax.prc.sp2, col="gray", cex=100*(exp(imptval)-1), pch=16, xlab="Colonizers", ylab="Young secondary colonizers")
abline(h=0); abline(v=0)
points(tax.prc.sp2~tax.prc.sp1, subset=order(imptval, decreasing=TRUE)[1:20], col="white", bg="gray40", cex=100*(exp(imptval)-1), pch=21, lwd=2)
identify(tax.prc.sp1, tax.prc.sp2, names(tax.prc.sp1), xpd=TRUE, cex=0.7, font=3)

plot(tax.prc.sp3, tax.prc.sp2, col="gray", cex=100*(exp(imptval)-1), pch=16, xlab="Old-growth Unaffected versus Young Unaffected", ylab="Young secondary colonizers")
abline(h=0); abline(v=0)
points(tax.prc.sp2~tax.prc.sp3, subset=order(imptval, decreasing=TRUE)[1:20], col="white", bg="gray40", cex=100*(exp(imptval)-1), pch=21, lwd=2)
identify(tax.prc.sp3, tax.prc.sp2, names(tax.prc.sp1), xpd=TRUE, cex=0.7, font=3)


plot(tax.prc.sp1, tax.prc.sp3, col="gray", cex=100*(exp(imptval)-1), pch=16, xlab="Colonizers", ylab="Old-growth Unaffected versus Young Unaffected")
abline(h=0); abline(v=0)
points(tax.prc.sp3~tax.prc.sp1, subset=order(imptval, decreasing=TRUE)[1:20], col="white", bg="gray40", cex=100*(exp(imptval)-1), pch=21, lwd=2)
identify(tax.prc.sp1, tax.prc.sp3, names(tax.prc.sp1), xpd=TRUE, cex=0.7, font=3)


##partial-RDA following PRC formulation to test significance of axes
tax.prc.rda<-rda(log(tree.matrix+1)~damage*type*yearf+Condition(yearf), data=envi.long)

tax.prc.rda.aov<-anova(tax.prc.rda, by="axis")
with(tax.prc.rda.aov, Variance[1:3]/sum(Variance))
  #inertia constrained

fun.prc.rda<-rda(CWM.scaled~damage*type*yearf+Condition(yearf), data=envi.long)

fun.prc.rda.aov<-anova(fun.prc.rda, by="axis")
with(fun.prc.rda.aov, Variance[1:3]/sum(Variance))

write.csv(rbind(data.frame(tax.prc.rda.aov), data.frame(fun.prc.rda.aov)), "ANOVA_results_comp.csv",
	quote=FALSE, row.names=TRUE)


#Monte Carlo pairwise comparisons for significant axes

tax.prc.scores<-scores(tax.prc, choices=1:3)$sites
fun.prc.scores<-scores(fun.prc, choices=1:3)$sites

MCttest(tax.prc.scores[,1], envi.long$damage, envi.long$type=="Young"&envi.long$year==2011)
MCttest(tax.prc.scores[,1], envi.long$damage, envi.long$type=="Young"&envi.long$year==2012)
MCttest(tax.prc.scores[,1], envi.long$damage, envi.long$type=="Young"&envi.long$year==2013)
MCttest(tax.prc.scores[,1], envi.long$damage, envi.long$type=="Young"&envi.long$year==2014)
MCttest(tax.prc.scores[,1], envi.long$damage, envi.long$type=="Young"&envi.long$year==2015)

MCttest(tax.prc.scores[,2], envi.long$damage, envi.long$type=="Young"&envi.long$year==2011)
MCttest(tax.prc.scores[,2], envi.long$damage, envi.long$type=="Young"&envi.long$year==2012)
MCttest(tax.prc.scores[,2], envi.long$damage, envi.long$type=="Young"&envi.long$year==2013)
MCttest(tax.prc.scores[,2], envi.long$damage, envi.long$type=="Young"&envi.long$year==2014)
MCttest(tax.prc.scores[,2], envi.long$damage, envi.long$type=="Young"&envi.long$year==2015)

MCttest(tax.prc.scores[,3], envi.long$damage, envi.long$type=="Young"&envi.long$year==2011)
MCttest(tax.prc.scores[,3], envi.long$damage, envi.long$type=="Young"&envi.long$year==2012)
MCttest(tax.prc.scores[,3], envi.long$damage, envi.long$type=="Young"&envi.long$year==2013)
MCttest(tax.prc.scores[,3], envi.long$damage, envi.long$type=="Young"&envi.long$year==2014)
MCttest(tax.prc.scores[,3], envi.long$damage, envi.long$type=="Young"&envi.long$year==2015)

MCttest(tax.prc.scores[,1], envi.long$damage, envi.long$type=="Old"&envi.long$year==2011)
MCttest(tax.prc.scores[,1], envi.long$damage, envi.long$type=="Old"&envi.long$year==2012)
MCttest(tax.prc.scores[,1], envi.long$damage, envi.long$type=="Old"&envi.long$year==2013)
MCttest(tax.prc.scores[,1], envi.long$damage, envi.long$type=="Old"&envi.long$year==2014)
MCttest(tax.prc.scores[,1], envi.long$damage, envi.long$type=="Old"&envi.long$year==2015)

MCttest(tax.prc.scores[,2], envi.long$damage, envi.long$type=="Old"&envi.long$year==2011)
MCttest(tax.prc.scores[,2], envi.long$damage, envi.long$type=="Old"&envi.long$year==2012)
MCttest(tax.prc.scores[,2], envi.long$damage, envi.long$type=="Old"&envi.long$year==2013)
MCttest(tax.prc.scores[,2], envi.long$damage, envi.long$type=="Old"&envi.long$year==2014)
MCttest(tax.prc.scores[,2], envi.long$damage, envi.long$type=="Old"&envi.long$year==2015)

MCttest(tax.prc.scores[,3], envi.long$damage, envi.long$type=="Old"&envi.long$year==2011)
MCttest(tax.prc.scores[,3], envi.long$damage, envi.long$type=="Old"&envi.long$year==2012)
MCttest(tax.prc.scores[,3], envi.long$damage, envi.long$type=="Old"&envi.long$year==2013)
MCttest(tax.prc.scores[,3], envi.long$damage, envi.long$type=="Old"&envi.long$year==2014)
MCttest(tax.prc.scores[,3], envi.long$damage, envi.long$type=="Old"&envi.long$year==2015)

MCttest(fun.prc.scores[,1], envi.long$damage, envi.long$type=="Young"&envi.long$year==2011)
MCttest(fun.prc.scores[,1], envi.long$damage, envi.long$type=="Young"&envi.long$year==2012)
MCttest(fun.prc.scores[,1], envi.long$damage, envi.long$type=="Young"&envi.long$year==2013)
MCttest(fun.prc.scores[,1], envi.long$damage, envi.long$type=="Young"&envi.long$year==2014)
MCttest(fun.prc.scores[,1], envi.long$damage, envi.long$type=="Young"&envi.long$year==2015)

MCttest(fun.prc.scores[,2], envi.long$damage, envi.long$type=="Young"&envi.long$year==2011)
MCttest(fun.prc.scores[,2], envi.long$damage, envi.long$type=="Young"&envi.long$year==2012)
MCttest(fun.prc.scores[,2], envi.long$damage, envi.long$type=="Young"&envi.long$year==2013)
MCttest(fun.prc.scores[,2], envi.long$damage, envi.long$type=="Young"&envi.long$year==2014)
MCttest(fun.prc.scores[,2], envi.long$damage, envi.long$type=="Young"&envi.long$year==2015)

MCttest(fun.prc.scores[,1], envi.long$damage, envi.long$type=="Old"&envi.long$year==2011)
MCttest(fun.prc.scores[,1], envi.long$damage, envi.long$type=="Old"&envi.long$year==2012)
MCttest(fun.prc.scores[,1], envi.long$damage, envi.long$type=="Old"&envi.long$year==2013)
MCttest(fun.prc.scores[,1], envi.long$damage, envi.long$type=="Old"&envi.long$year==2014)
MCttest(fun.prc.scores[,1], envi.long$damage, envi.long$type=="Old"&envi.long$year==2015)

MCttest(fun.prc.scores[,2], envi.long$damage, envi.long$type=="Old"&envi.long$year==2011)
MCttest(fun.prc.scores[,2], envi.long$damage, envi.long$type=="Old"&envi.long$year==2012)
MCttest(fun.prc.scores[,2], envi.long$damage, envi.long$type=="Old"&envi.long$year==2013)
MCttest(fun.prc.scores[,2], envi.long$damage, envi.long$type=="Old"&envi.long$year==2014)
MCttest(fun.prc.scores[,2], envi.long$damage, envi.long$type=="Old"&envi.long$year==2015)

#plot of partial RDA

damage.type.col<-with(envi.long,
ifelse(damage=="Unaffected",
	ifelse(type=="Young", "green", "darkgreen"),
	ifelse(type=="Young", "orange", "red"))
	)

damage.type.pch<-with(envi.long,
ifelse(damage=="Unaffected",
	ifelse(type=="Young", 2, 17),
	ifelse(type=="Young", 1, 16))
	)

envi.plot<-unique(envi.long$plot)


##Taxonomic composition

par(mfrow=c(2,3), mar=c(4,4,1,1), mgp=c(2.5, 1, 0))

boxplot(tax.prc.scores[,2]~damage.type+envi.long$yearf, axes=FALSE, col=draw.plot.col[c(2,1,4,3)])
axis(1, c(2.5,6.5,10.5,14.5,18.5), c(2011,2012,2013,2014,2015))
plot(tax.prc.scores[,1], tax.prc.scores[,2], type="n", xlab="RDA 1", ylab="RDA 2")
for (i in 1:length(envi.plot)) {
	index<-envi.long$plot==envi.plot[i]
	lines(tax.prc.scores[index,1], tax.prc.scores[index,2], col="gray")
	}
text(tax.prc, display="sites", choices=c(1,2),
	labels=as.numeric(envi.long$yearf)-1,
	col=damage.type.col,
	cex=0.7)

plot(tax.prc.scores[,3], tax.prc.scores[,2], type="n", xlab="RDA 3", ylab="RDA 2")
for (i in 1:length(envi.plot)) {
	index<-envi.long$plot==envi.plot[i]
	lines(tax.prc.scores[index,3], tax.prc.scores[index,2], col="gray")
	}
text(tax.prc, display="sites", choices=c(3,2),
	labels=as.numeric(envi.long$yearf)-1,
	col=damage.type.col,
	cex=0.7)

plot.new()

boxplot(tax.prc.scores[,1]~damage.type+envi.long$yearf, axes=FALSE, col=draw.plot.col[c(2,1,4,3)], horizontal=TRUE)
axis(2, c(2.5,6.5,10.5,14.5,18.5), c(2011,2012,2013,2014,2015))
boxplot(tax.prc.scores[,3]~damage.type+envi.long$yearf, axes=FALSE, col=draw.plot.col[c(2,1,4,3)], horizontal=TRUE)
axis(2, c(2.5,6.5,10.5,14.5,18.5), c(2011,2012,2013,2014,2015))


##Functional composition

par(mfrow=c(2,2), mar=c(4,4,1,1), mgp=c(2, 0.8, 0))

boxplot(fun.prc.scores[,2]~damage.type+envi.long$yearf, axes=FALSE, col=draw.plot.col[c(2,1,4,3)])
axis(1, c(2.5,6.5,10.5,14.5,18.5), c(2011,2012,2013,2014,2015))

plot(fun.prc.scores[,1], fun.prc.scores[,2], type="n", xlab="RDA 1", ylab="RDA 2")
for (i in 1:length(envi.plot)) {
	index<-envi.long$plot==envi.plot[i]
	lines(fun.prc.scores[index,1], fun.prc.scores[index,2], col="gray")
	}
text(fun.prc, display="sites",
	labels=as.numeric(envi.long$yearf)-1,
	col=damage.type.col,
	cex=0.7)

plot.new()

boxplot(fun.prc.scores[,1]~damage.type+envi.long$yearf, axes=FALSE, col=draw.plot.col[c(2,1,4,3)], horizontal=TRUE)
axis(2, c(2.5,6.5,10.5,14.5,18.5), c(2011,2012,2013,2014,2015))



basal.mat<-tapply(1/dbh$basal, list(dbh$y.plot, dbh$species), sum, na.rm=TRUE)

CWM.basal<-functcomp(trait.hastree, basal.mat[,
	colnames(basal.mat) %in% rownames(trait.hastree) & colnames(basal.mat)!="Unknown"],
	CWM.type="all")

CWM.basal.scaled<-scale(CWM.basal)

fun.prc<-prc(CWM.basal.scaled, damage.type, envi.long$yearf)

fun.prc.coef1<-summary(fun.prc, axis=1)$coef
fun.prc.coef1<-rbind(fun.prc.coef1, rep(0, ncol(fun.prc.coef1)))
rownames(fun.prc.coef1)[4]<-"Old Unaffected"
fun.prc.coef1<-fun.prc.coef1[c(1,4,3,2),]

fun.prc.trait1<-summary(fun.prc, axis=1)$sp
names(fun.prc.trait1)[6] <- "Hmax"

fun.prc.coef2<-summary(fun.prc, axis=2)$coef
fun.prc.coef2<-rbind(fun.prc.coef2, rep(0, ncol(fun.prc.coef2)))
rownames(fun.prc.coef2)[4]<-"Old Unaffected"
fun.prc.coef2<-fun.prc.coef2[c(1,4,3,2),]

fun.prc.trait2<-summary(fun.prc, axis=2)$sp
names(fun.prc.trait1)[6] <- "Hmax"

par(mfrow=c(2,1), mar=c(3,3,1,15), mgp=c(2,0.8,0))

matplot(colnames(fun.prc.coef1), t(fun.prc.coef1),
	col=draw.plot.col, pch=draw.plot.pch, 
	ylim=c(-0.4,0.5), xlim=c(2011,2015.8),
	ylab="PRC 1", xlab="Year")
for (i in 1:nrow(fun.prc.coef1)) lines(colnames(fun.prc.coef1), fun.prc.coef1[i,], col=draw.plot.col[i], lty=draw.plot.lty[i])
matplot(colnames(fun.prc.coef1), t(fun.prc.coef1),
	col="white", pch=16, cex=2, add=TRUE)
matplot(colnames(fun.prc.coef1), t(fun.prc.coef1),
	col=draw.plot.col, pch=draw.plot.pch, add=TRUE)
linestack(fun.prc.trait1/4, names(fun.prc.trait1), add=TRUE, at=par()$usr[2], cex=0.7)
axis(2, at=axTicks(2), labels=axTicks(2)*4, pos=par()$usr[2])
text(par()$usr[1]-0.14*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(a)", xpd=TRUE)

matplot(colnames(fun.prc.coef2), t(fun.prc.coef2),
	col=draw.plot.col, pch=draw.plot.pch,
	ylim=c(-0.7,0.85), xlim=c(2011, 2015.8),
	ylab="PRC 2", xlab="Year")
for (i in 1:nrow(fun.prc.coef2)) lines(colnames(fun.prc.coef2), fun.prc.coef2[i,], col=draw.plot.col[i], lty=draw.plot.lty[i])
matplot(colnames(fun.prc.coef2), t(fun.prc.coef2),
	col="white", pch=16, cex=2, add=TRUE)
matplot(colnames(fun.prc.coef2), t(fun.prc.coef2),
	col=draw.plot.col, pch=draw.plot.pch, add=TRUE)
linestack(fun.prc.trait2/2, names(fun.prc.trait2), add=TRUE, at=par()$usr[2], cex=0.7)
axis(2, at=axTicks(2), labels=axTicks(2)*2, pos=par()$usr[2])
text(par()$usr[1]-0.14*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(b)", xpd=TRUE)

legend("bottom", legend=c("Young Unaffected", "Old Unaffected", "Young Affected", "Old Affected"),pch=draw.plot.pch, lty=draw.plot.lty, col=draw.plot.col, 
       lwd=1.4, cex=0.8, bty="n")

fun.prc.rda<-rda(CWM.basal.scaled~damage*type*yearf+Condition(yearf), data=envi.long)
fun.prc.rda.aov<-anova(fun.prc.rda, by="axis")



##trash code

data.frame(
	row.names=c("type", "damage", "yearf", "type:damage", "type:yearf", "damage:yearf", "type:damage:yearf"),
	Chisq=c(
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=2)$result$chi2[1],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=3)$result$chi2[1],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=4:7)$result$chi2[1],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=8)$result$chi2[1],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=9:12)$result$chi2[1],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=13:16)$result$chi2[1],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=17:20)$result$chi2[1]
		),
	df=c(
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=2)$result$chi2[2],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=3)$result$chi2[2],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=4:7)$result$chi2[2],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=8)$result$chi2[2],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=9:12)$result$chi2[2],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=13:16)$result$chi2[2],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=17:20)$result$chi2[2]
		),
	pval=c(
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=2)$result$chi2[3],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=3)$result$chi2[3],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=4:7)$result$chi2[3],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=8)$result$chi2[3],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=9:12)$result$chi2[3],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=13:16)$result$chi2[3],
		wald.test(vcov(D0.lme), fixef(D0.lme), Terms=17:20)$result$chi2[3]
		)
	)

FD10<-FD1$FuncD/FD0$FuncD
FD20<-FD2$FuncD/FD0$FuncD
FD21<-FD2$FuncD/FD1$FuncD

aggregate(FD10~year+type+damage, data=envi.long, mean)
aggregate(FD20~year+type+damage, data=envi.long, mean)

par(mfrow=c(2,3), mar=c(3,3,1,1), mgp=c(2,0.8,0), cex.axis=1, cex.lab=1.2)
plot(FD0$FuncD~FD$FRic, col=damage.type.col, pch=damage.type.pch, ylab="FD0", xlab="FRic")
plot(FD1$FuncD~FD$FRic, col=damage.type.col, pch=damage.type.pch, ylab="FD1", xlab="FRic")
plot(FD2$FuncD~FD$FRic, col=damage.type.col, pch=damage.type.pch, ylab="FD2", xlab="FRic")
plot(FD10~FD$FEve, col=damage.type.col, pch=damage.type.pch, ylab="FD1 / FD0", xlab="FEve")
plot(FD20~FD$FEve, col=damage.type.col, pch=damage.type.pch, ylab="FD2 / FD0", xlab="FEve")
plot(FD21~FD$FEve, col=damage.type.col, pch=damage.type.pch, ylab="FD2 / FD1", xlab="FEve")


FuncRed2<-D2/FD2$FuncD

with(envi.long, draw.plot(var=FuncRed2, year=year, yt=yt, labely="Taxonomic Redundancy"))

plot(FuncRed0~FuncRed2, col=damage.type.col)

Anova(FRic.lme<-lmer(FD$FRic~type*damage*yearf+(1|plot), data=envi.long))
FRic.lme.refit<-lmer(FD$FRic~type.damage.year.inter-1+(1|plot), data=envi.long, REML=FALSE)
FRic.glht<-summary(glht(FRic.lme.refit, linfct=type.damage.year.linfct), test=adjusted("bonferroni"))
FRic.glht.tab<-data.frame(diff=FRic.glht$test$coef, se=FRic.glht$test$sigma, z=FRic.glht$test$tstat)
FRic.glht.tab$pval<-replace(FRic.glht$test$pvalues, FRic.glht$test$pvalues<0.001, "<0.001")

write.csv(FRic.glht.tab, "FRicglhttab.csv", quote=FALSE, row.names=TRUE)

Anova(FEve.lme<-lmer(FD$FEve~type*damage*yearf+(1|plot), data=envi.long))
FEve.lme.refit<-lmer(FD$FEve~damage.year.inter-1+(1|plot), data=envi.long)
FEve.glht<-summary(glht(FEve.lme.refit, linfct=damage.year.linfct), test=adjusted("bonferroni"))
FEve.glht.tab<-data.frame(diff=FEve.glht$test$coef, se=FEve.glht$test$sigma, z=FEve.glht$test$tstat)
FEve.glht.tab$pval<-replace(FEve.glht$test$pvalues, FEve.glht$test$pvalues<0.001, "<0.001")

write.csv(FEve.glht.tab, "FEveglhttab.csv", quote=FALSE, row.names=TRUE)

Anova(FDis.lme<-lmer(FD$FDis~type*damage*yearf+(1|plot), data=envi.long))
FDis.lme.refit<-lmer(FD$FDis~type.damage.year.inter-1+(1|plot), data=envi.long)
FDis.glht<-summary(glht(FDis.lme.refit, linfct=type.damage.year.linfct), test=adjusted("bonferroni"))
FDis.glht.tab<-data.frame(diff=FDis.glht$test$coef, se=FDis.glht$test$sigma, z=FDis.glht$test$tstat)
FDis.glht.tab$pval<-replace(FDis.glht$test$pvalues, FDis.glht$test$pvalues<0.001, "<0.001")

write.csv(FDis.glht.tab, "FDisglhttab.csv", quote=FALSE, row.names=TRUE)

with(envi.long, draw.plot(var=FD$FRic, year= year, yt=yt,
	labely="FRic"))
text(par()$usr[1]-0.18*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(d)", xpd=TRUE)

with(envi.long, draw.plot(var=FD$FEve, year= year, yt=yt,
	labely="FEve"))
text(par()$usr[1]-0.18*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(e)", xpd=TRUE)

with(envi.long, draw.plot(var=FD$FDis, year= year, yt=yt,
	labely="FDis"))
text(par()$usr[1]-0.18*(par()$usr[2]-par()$usr[1]), par()$usr[3]+1*(par()$usr[4]-par()$usr[3]), "(f)", xpd=TRUE)

library(AICcmodavg)

new.data<-with(envi.long, expand.grid(levels(type), levels(damage), levels(yearf)))
colnames(new.data)<-c("type", "damage", "yearf")

new.data.pred<-modavgPred(list(FDis.lme), newdata=new.data)

pred<-cbind(new.data, pred=new.data.pred$mod.avg.pred, se=new.data.pred$uncond.se, lower.CL=new.data.pred$lower.CL, upper.CL=new.data.pred$upper.CL)

aggregate(FD0$FuncD~type+damage+year, data=envi.long, mean)
aggregate(FD1$FuncD~type+damage+year, data=envi.long, mean)
aggregate(FD2$FuncD~type+damage+year, data=envi.long, mean)

offset=0.05

pred$year<-with(pred, as.numeric(as.character(yearf))+ifelse(
	type=="Young", ifelse(damage=="Affected", -3*offset, -1*offset),
		ifelse(damage=="Affected", 1*offset, 3*offset)))


plot(pred~year, data=pred,
	col=ifelse(damage=="Affected", 
		ifelse(type=="Young", "gold2", "red"),
		ifelse(type=="Young", "green", "darkgreen")),
	pch=ifelse(damage=="Affected",
		ifelse(type=="Young", 1, 16),
		ifelse(type=="Young", 2, 17)),
	ylim=c(min(pred-se), max(pred+se)))
with(pred, arrows(year, pred+se, year, pred-se, angle=90, code=3, length=0.01))

lines(pred~year, data=pred, subset=type=="Old"&damage=="Affected", col="red", lty=1)
lines(pred~year, data=pred, subset=type=="Young"&damage=="Affected", col="gold2", lty=1)
lines(pred~year, data=pred, subset=type=="Old"&damage=="Unaffected", col="darkgreen", lty=2)
lines(pred~year, data=pred, subset=type=="Young"&damage=="Unaffected", col="green", lty=2)


#PRC-inspired plot
rda1<-scores(tax.prc.rda)$sites[,1]
rda2<-scores(tax.prc.rda)$sites[,2]

sp.select1<-abs(scores(tax.prc.rda)$species[,1])>0.2&colSums(tree.matrix)>10
sp.select2<-abs(scores(tax.prc.rda)$species[,2])>0.2&colSums(tree.matrix)>10

par(mfrow=c(2,1), mar=c(3,3,1,13), mgp=c(2,0.8,0))
with(envi.long, draw.plot(var=rda1, year=year, yt=yt, 
	labely="RDA 1"))
abline(h=0, col="gray")
linestack(scores(tax.prc.rda)$species[sp.select1,1], rownames(scores(tax.prc.rda)$species)[sp.select1], add=TRUE, at=2015.35, cex=0.7)
with(envi.long, draw.plot(var=rda2, year=year, yt=yt, 
	labely="RDA 2"))
abline(h=0, col="gray")
linestack(scores(tax.prc.rda)$species[sp.select2,2], rownames(scores(tax.prc.rda)$species)[sp.select2], add=TRUE, at=2015.35, cex=0.7)

###pairwise functional similarity indices for CqN and UqN similarity indices
#input:
#abun:species by community matrix dataframe, where row is species and
#column is community or site.
# q: order q should be any non-negative value.
#output: pairwise functional similarity matrix 
# CqN: Sorensen-type overlap index (from a local view)(similarity index) 
# UqN: Jaccard-type overlap index (from a regional view)(similarity index)

pairFunc=function(Dis, abun,q){
  N=ncol(abun);
  CqN=matrix(1,ncol=N,nrow=N);UqN=CqN;
  for(i in 1:(N-1)){
    for(j in i:N){
      o=Func2014(Dis,abun[,c(i,j)],q);
      CqN[i,j]=o[[3]][5];CqN[j,i]=CqN[i,j];
      UqN[i,j]=o[[3]][6];UqN[j,i]=UqN[i,j];
    }
  }
  return(list(CqN=CqN,UqN=UqN));
}


CWM.scaled.dist<-dist(CWM.scaled)

CWM.betadisper<-betadisper(CWM.scaled.dist, with(envi.long, paste(damage, type, year)))

boxplot(~year, data=envi.long, subset=damage=="Affected"&type=="Young")

with(envi.long, draw.plot(var=CWM.betadisper$distances, year=year, yt=yt, labely="betadisper"))


with(envi.long, draw.plot(var=scores(fun.prc.rda)$sites[,1], year=year, yt=yt, labely="RDA1 (functional)"))


trait.hastree$rda1<-tax.prc.sp1[na.omit(match(names(tax.prc.sp1), rownames(trait.hastree)))]

par(mfrow=c(2,2), mar=c(4,4,1,1), mgp=c(2.5,1,0))
plot(rda1~Ht, data=trait.hastree, ylab="RDA 1", pch=16, col="gray", cex=2)
plot(rda1~LTh, data=trait.hastree, xlab="RDA 1", pch=16, col="gray", cex=2)
plot(rda1~SDM, data=trait.hastree, xlab="RDA 1", pch=16, col="gray", cex=2)
plot(rda1~SLA, data=trait.hastree, xlab="RDA 1", pch=16, col="gray", cex=2)

trait.hastree$rda2<-tax.prc.sp2[na.omit(match(names(tax.prc.sp2), rownames(trait.hastree)))]

par(mfrow=c(2,2), mar=c(4,4,1,1), mgp=c(2.5,1,0))
plot(rda2~WD, data=trait.hastree, ylab="RDA 2", pch=16, col="gray", cex=2)
plot(rda2~LTh, data=trait.hastree, ylab="RDA 2", pch=16, col="gray", cex=2)
plot(rda2~SDM, data=trait.hastree, ylab="RDA 2", pch=16, col="gray", cex=2)
plot(rda2~LDMC, data=trait.hastree, ylab="RDA 2", pch=16, col="gray", cex=2)

#typical PRC plots
plot(tax.prc, select=colSums(tree.matrix)>200, xlab="Year", legpos=NA, ylab="PRC1", type="b", pch=draw.plot.pch[2:4], col=draw.plot.col[2:4], axis=1)
plot(tax.prc, select=colSums(tree.matrix)>200, xlab="Year", legpos=NA, ylab="PRC2", type="b", pch=draw.plot.pch[2:4], col=draw.plot.col[2:4], axis=2)

