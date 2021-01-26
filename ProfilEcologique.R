#####################################################################
# frequence(x0,yobs)
# Dessin de la fréquence de l'espèce yobs le long du gradient x0
# Fonction par Emilien Kuhn
# janvier 2013
#####################################################################

frequence <- function(x0, yobs, xname, nclasses) {
	# Arguments x0 (explicatrice), yobs (à expliquer), xname (nom de la variable ecologique), nclasses (nombre de classes a realiser)
  	
	# histogramme de la varible explicatrice : creation des classes
	istx0=hist(x0,breaks=nclasses,plot=FALSE)
	
	# Preparation du resultat
	freqyobs=istx0$mids
	
	# CALCUL DES FREQUENCES
	for (i in 2:length(istx0$breaks)) {
  	tmp=yobs[which(x0>istx0$breaks[(i-1)] & x0<=istx0$breaks[i])]
  	freqyobs[(i-1)]=length(which(tmp==1))/length(tmp)
  }
  
	# Trace les frequences calculees
	par(mar=c(5.1,5.1,1.1,1.1))
  barplot(freqyobs,space=0,xlab=xname,ylab="Frequence des presences",cex.lab=2,col="chocolate2",ylim=c(0,1))
  axis(1,at=seq(0,length(freqyobs),by=1),labels=istx0$breaks)
}




#####################################################################
# freqgrad(x0,yobs)
# Dessin du profil écologique de l'espèce yobs le long du gradient x0
# Fonction écrite par Daniel Chessel, adaptée par Christophe Coudun, puis par Emilien Kuhn (R)
# avril 2012
#####################################################################

freqgrad <- function(x0, yobs, ncla = 10)
{
	# Arguments x0 (explicatrice), yobs (à expliquer), ncla (classes)
	# PREPARATION de VECTEURS sans VALEURS MANQUANTES
	yobs <- yobs[is.na(x0) == F]
	# On elimine les releves avec donnees manquantes
	x0 <- x0[is.na(x0) == F]
	# DISTRIBUTION de la VARIABLE QUANTITATIVE
	q0 <- quantile(x0, probs = seq(0, 1, 1/ncla), na.rm = F)
	# Quantiles de la distribution de x0 (0, 10, 20...100%)
	c.cla <- quantile(x0, probs = seq(0 + 1/2/ncla, 1 - 1/2/ncla, 1/ncla),
		na.rm = F)
	# Quantiles x0 (5, 15, 25...95%)
	xmin <- min(x0)
	# Valeurs extremes de x0
	xmax <- max(x0)
	# DISCRETISATION de la VARIABLE QUANTITATIVE et OCCURRENCES dans chaque CLASSE
	q1 <- cut(x0, q0, include.lowest = T)
	# Decoupage de la variable x0 en classes avec intervalles fixes par q0
	t0 <- table(q1, yobs)
	# Table de contingence des presences et absences dans chaque classe de x0
	t0[, 1] <- (t0[, 1] + t0[, 2])
	# Nombre de releves totaux dans chaque intervalle (absence + presence)
	freq <- t0[, 2]/t0[, 1]
	# Calcul de la frequence de l'espece dans chaque classe
	toto1 <- rep(0, ncla)
	# Creation de deux variables toto receptrices d'info
	toto2 <- rep(0, ncla)
	for(i in 1:ncla) {
		# Pour chaque intervalle de la variable explicatrice
		succes <- t0[i, 2]
		# Nombre de presences de l'espece dans l'intervalle
		essai <- t0[i, 1]
		# Nombre de releves dans l'intervalle
		if(essai > 10) {
			a0 <- prop.test(succes, essai)$conf.int
			# Tests de proportion et intervalle de confiance (par interv.)
			toto1[i] <- a0[1]
			# Borne inf. intervalle de confiance
			toto2[i] <- a0[2]
			# Borne sup. intervalle de confiance
			if(a0[1] > freq[i]) toto1[i] <- NA
			if(a0[2] < freq[i])
				toto2[i] <- NA
		}
		else {
			toto1[i] <- NA
			toto2[i] <- NA
		}
	}
	ymin <- min(toto1, na.rm = T)
	# Borne inf. des bornes inf.
	ymax <- max(toto2, na.rm = T)
	# Borne sup. des bornes sup.
	x11()
	plot(xmin, 0, ylim = c(0, 1), xlim = c(xmin, xmax), type = "n", xlab = "", ylab = "")
	# On peut prendre ymin et ymax pour borner y
	points(c.cla, freq, pch = 16)
	# Dessin des frequences calculees dans chaque intervalle centrees sur milieu
	for(i in 1:ncla) {
		if(!is.na(toto1[i])) {
			lines(x = c(c.cla[i], c.cla[i]), y = c(toto1[i], toto2[i]))
		}
		if(i > 1) {
			abline(v = q0[i], lty = 2)
		}
	}
}



######################################################################################################
# vi.s(coeff,grad)
# Extraction des valeurs indicatrices d'une courbe de réponse issue d'une régression logistique simple
# Extraction de l'optimum, de l'amplitude et de la probabilité maximale de la courbe
# Fonction écrite par Christophe Coudun et adaptée a R par Emilien Kuhn
# novembre 2012
######################################################################################################

vi.s <- function(coeff, grad, s.dist = 80, method = "gep", graph = F) {
	# Découpage du gradient écologique pour calculer la courbe
	x.cr <- seq(min(grad), max(grad), length = 1000)
	# Coefficients de régression
	b0 <- coeff[1]
	b1 <- coeff[2]
	if(is.na(b1)) {
		b1 <- 0
	}
	b2 <- coeff[3]
	if(is.na(b2)) {
		b2 <- 0
	}
	# Calcul des probabilités associées aux valeurs de x.cr
	y.cr <- 1/(1 + exp( - b0 - b1 * x.cr - b2 * x.cr^2))
	# Probabilité maximale de la courbe de réponse
	pmax <- max(y.cr)
	# Optimum (valeur de grad correspondant à pmax)
	opt <- x.cr[y.cr == max(y.cr)][1]
	# Amplitude
	# Selon Ter Braak et Looman (1986)
	##################################
	if(method == "tbl") {
		if(b2 < 0) {
			u <-  - b1/(2 * b2)
			# Optimum
			t <- 1/sqrt(-2 * b2)
			# Tolérance
			h <- 1/(1 + exp(b1^2/(4 * b2) - b0))
			# Probabilité maximale
			vi <- list(u = u, t = t, h = h)
		}
		else {
			vi <- list(u = NA, t = NA, h = NA)
		}
		if(graph) {
			plot(x.cr, y.cr, type = "l", xlab = "", ylab = "", ylim = c(0, 1))
			lines(x.cr[x.cr >= u - t & x.cr <= u + t], y.cr[x.cr >= u - t & x.cr <= u + t], lwd = 4, col = 2)
			points(opt, pmax, pch = 20, col = 1,cex=3)
		}
	}
	# Selon Heegaard (2002)
	#######################
	if(method == "heg") {
		pcentral <- pmax/exp(0.5)
		pouter <- pmax/exp(2)
		cbi <- min(x.cr[y.cr > pcentral])[1]
		# central border inferieur
		cbs <- max(x.cr[y.cr > pcentral])[1]
		# central border superieur
		tol <- cbs - cbi
		# largeur central border (tolerance)
		obi <- min(x.cr[y.cr > pouter])[1]
		# outer border inferieur
		obs <- max(x.cr[y.cr > pouter])[1]
		# outer border superieur
		ran <- obs - obi
		# largeur outer border (range)
		if(graph) {
			plot(x.cr, y.cr, type = "l", xlab = "", ylab = "", ylim = c(0, 1))
			lines(x.cr[x.cr >= obi & x.cr <= obs], y.cr[x.cr >= obi & x.cr <= obs], lwd = 4, col = 2)
			lines(x.cr[x.cr >= cbi & x.cr <= cbs], y.cr[x.cr >= cbi & x.cr <= cbs], lwd = 4, col = 3)
			points(opt, pmax, pch = 20, col = 1,cex=3)
		}
		vi <- list(cbi = cbi, cbs = cbs, tol = tol, obi = obi, obs = obs, ran = ran)
	}
	# Selon Gégout et Pierrat (1998)
	################################
	if(method == "gep") {
		s <- 0
		k <- 0
		inclus <- rep(0, 1000)
		# Tant qu'on n'a pas les 80 % de distribution a inclure
		while(s <= 0.01 * s.dist * sum(y.cr)) {
			k <- k + 1
			s <- s + y.cr[rev(order(y.cr))[k]]
			# order(y.cr)[1] renvoie la case de y.cr ou on observe la proba minimum de y.cr 
			# order(y.cr)[1000] renvoie la case de y.cr ou on observe la proba maximum de y.cr
			# On retourne le vecteur order(y.cr) avec la fonction "rev"
			# On lit alors les cases de y.cr avec la proba rangee en ordre decroissant
			# On incremente la somme en ajoutant les probas decroissantes
			inclus[k] <- x.cr[rev(order(y.cr))[k]]
		}
		# !!! Re-crire plus clairement cette explication !!!
		# Vecteur des valeurs de x conservees pour atteindre 80 % de distribution 
		inclus <- inclus[inclus != 0]
		min.amp <- min(inclus)
		max.amp <- max(inclus)
		amp.st <- (max(inclus) - min(inclus))/(0.01 * s.dist * (max(grad) - min(grad)))
		amp <- max.amp - min.amp
		vi <- list(opt = opt, amp.st = amp.st, amp = amp, pmax = pmax, min.amp = min.amp, max.amp = max.amp)
		if(b2 > 0) {
			vi <- list(opt = NA, amp.st = amp.st, amp = NA, pmax = NA, min.amp = NA, max.amp = NA)
		}
		if(b1 == 0 & b2 == 0) {
			vi <- list(opt = NA, amp.st = 1, amp = NA, pmax = pmax, min.amp = NA, max.amp = NA)
		}
		if(graph) {
			plot(x.cr, y.cr, type = "l", xlab = "", ylab = "", ylim = c(0, 1))
			lines(x.cr[x.cr >= min.amp & x.cr <= max.amp], y.cr[x.cr >= min.amp & x.cr <= max.amp], lwd = 4, col = 2)
			points(opt, pmax, pch = 20, col = 1,cex=3)
		}
	}
	return (vi)
}

############################################################################################
# roc(pred,obs)
# Evaluation de la qualité des prédictions d'un modèle par rapport aux données d'observation
# Evaluation des principales mesures de qualité dont l'aire sous la courbe ROC
# Fonction écrite par Christophe Coudun a R par Emilien Kuhn
# novembre 2012
############################################################################################

roc <- function(pred, obs, npoints = 100) {
	# Arguments : probas predites, valeurs observees et npoints
	# INITIALISATION des COORDONNEES
	xroc <- rep(0, npoints + 2)
	# Recoit les abscisses de roc (1-Sp)
	yroc <- rep(0, npoints + 2)
	# Recoit les ordonnees de roc (Sn)
	# VECTEURS des SEUILS testes SUCCESSIVEMENT
	seuil <- seq(0, 1, length = npoints + 1)
	seuil <- seuil * (max(pred)[1])
	# Seuil de proba d'acceptabilite de la presence
	for(k in 1:(npoints + 1)) {
		predPA <- pred
		# Transformation des donnees de probas en presence absence
		# MATRICES de CONFUSION SUCCESSIVES
		predPA[predPA < seuil[k]] <- 0
		predPA[predPA >= seuil[k]] <- 1
		conf <- table(predPA, obs)
		# Matrice de confusion associee a "seuil"
		if(min(dim(conf)[1], dim(conf)[2]) != 1) {
			# Si la matrice de confusion est une matrice carree 2*2
			a <- conf[2, 2]
			b <- conf[2, 1]
			c <- conf[1, 2]
			d <- conf[1, 1]
			Sn <- a/(a + c)
			# Sensitivity (pourcentage de presences bien predites)
			Sp <- d/(b + d)
		}
		else {
			Sp <- 1
			Sn <- 0
		}
		# COORDONNEES des POINTS de la COURBE ROC
		xroc[k] <- 1 - Sp
		yroc[k] <- Sn
	}
	# GRAPHE ROC
	xroc[npoints + 2] <- 1
	# Il faut ajouter le point (1,1) a la courbe ROC
	yroc[npoints + 2] <- 1
	par(mar=c(5.1,5.1,1.1,1.1))
	plot(sort(xroc), sort(yroc),type="l",lwd=2,xlab="1 - Specificite",ylab="Sensibilite",cex.lab=3,xlim=c(0,1),ylim=c(0,1))
	abline(a=0,b=1,lty="dotted",lwd=2)
	# Possibilite de tracer le graphe de la courbe ROC
	# AUC (AIRE sous la COURBE)
	x <- xroc[order(xroc)]
	# Ordination des valeurs xroc et yroc
	y <- yroc[order(xroc)]
	auc <- rep(0, npoints + 1)
	for(k in 1:(npoints + 1)) {
		auc[k] <- ((y[k] + y[k + 1]) * (x[k + 1] - x[k]))/2
	}
	auc <- sum(auc)

	# Aire (Area Under the Curve)
	# QUALITE de PREDICTION du MODELE
	if(auc <= 0.7) {
		# Mediocre
		qual <- 3
	}
	if(auc > 0.7 & auc < 0.9) {
		# Moyen
		qual <- 2
	}
	if(auc > 0.9) {
		# Bon
		qual <- 1
	}
	# SEUIL OPTIMAL (optimisation sensitivity et specificity)
	sopt <- seuil[(yroc - xroc) == max(yroc - xroc)[1]][1]
	# Extraction du seuil optimal
	# STATISTIQUES associees au SEUIL OPTIMAL
	predPA <- pred
	predPA[predPA < sopt] <- 0
	predPA[predPA >= sopt] <- 1
	conf <- table(predPA, obs)
	# Matrice de confusion associee a sopt
	a <- conf[2, 2]
	b <- conf[2, 1]
	c <- conf[1, 2]
	d <- conf[1, 1]
	n <- a + b + c + d
	S <- (a + d)/n
	# Succes global de la prediction
	Sn <- a/(a + c)
	# Sensitivity (pourcentage de presences bien predites)
	Sp <- d/(b + d)
	# Specificite (pourcentage d'absences bien predites)
	or <- (a * d)/(c * b)
	# Odds ratio (quotient bien predits/mal predits)
	npp <- d + (c + d)
	# Negative predictive power (pourcentage de absences predites alors que vraies presences)
	NMI <- ( - a * log(a) - b * log(b) - c * log(c) - d * log(d) + (a +	b) * log(a + b) + (c + d) * log(c + d))/(n * log(n) - ((a + c) * (a + b) + (b + d) * (c + d))/n)
	kappa <- ((a + d) - ((a + c) * (a + b) + (b + d) * (c + d))/n)/(n -	((a + c) * (a + b) + (b + d) * (c + d))/n)
	roc.all <- rbind(xroc, yroc, yroc - xroc)
	# Coordonnees des points de la courbe ROC
	# ATTRIBUTION des RESULTATS
	list(courbe.roc = roc.all, seuil = round(sopt[1], 3), seuil.fract = sopt[1]/max(pred), auc = round(auc, 3), succes = round(S * 100, 0), sensitivity = round(Sn * 100, 0), specificity = round( Sp * 100, 0), confusion = conf, OR = or, npp = npp, NMI = NMI, kappa = kappa, qual = qual)
}
