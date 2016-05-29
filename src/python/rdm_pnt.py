def rdm_pnt(nbr_cluster):
	from numpy.random import uniform
	from numpy.random import multivariate_normal
	from numpy import savetxt

	import sys

	x = uniform(50, 950)
	y = uniform(50, 950)
	z = uniform(-500, 500) + 20.0
	cov = [[uniform(-5, 5),z], [z,uniform(-5, 5)]]
	s=multivariate_normal([x, y], cov, nbr_cluster[1])
	xo= s[:,0]
	yo= s[:,1]
	print s
	sys.stdout.flush()

	for i in range(nbr_cluster[0]-1):
		x = uniform(50, 950)
		y = uniform(50, 950)
		z = uniform(-500, 500) + 20.0
		cov = [[uniform(-5, 5),z], [z,uniform(-5, 5)]]
		s=multivariate_normal([x, y], cov, nbr_cluster[1])
		xo.extend(s[:,0])
		yo.extend(s[:,1])
		print xo
		sys.stdout.flush()

	savetxt('data.txt', s)

	return 1