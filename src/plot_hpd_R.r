plot.hpd <- function(coda.object,
		wanted=varnames(coda.object),
		prob.wide=0.95,
		prob.narrow=0.50,
		xlab="HPD interval",
		names=NULL,
		type="p",
        centered=F,
        pos=0,  ## adjust the label's position
		...)
{
	which.wanted=ifow(is.integer(wanted), wanted, match(wanted, varnames(coda.object)))
	num.wanted=length(which.wanted)

	chain <- mcmc.stack(coda.object)
	mu    <- colMeans(chain[,which.wanted])
	med   <- apply(coda::HPDinterval(chain, prob=0.01)[which.wanted,],
			1, mean)
	hpd.wide    <- coda::HPDinterval(chain, prob=prob.wide)[which.wanted,]
	hpd.narrow  <- coda::HPDinterval(chain, prob=prob.narrow)[which.wanted,]

	if (is.null(names)) names <- varnames(chain)[which.wanted]
	else names <- rep(names, length.out=length(wanted))
	ypos <- plot.ci(med, hpd.narrow, hpd.wide, names=names, xlab=xlab, col.midvals="white", pch.midvals="|", type=type, pos=pos, ...)
	if ("p"==type)
	{
		points(mu, ypos, pch="|")
	}
	invisible(ypos)
}