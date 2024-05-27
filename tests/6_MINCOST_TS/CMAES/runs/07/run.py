# 07
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[9.204405941,	1.933291667,	0.051039286,	4.497695391,	0.351875752], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
