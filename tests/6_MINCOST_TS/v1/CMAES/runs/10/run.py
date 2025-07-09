# 10
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[7.418118812,	5.47825,	2.884857143,	0.166725251,	7.63993988], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
