# 30
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[1.497227723	,	9.0521875	,	4.6175	,	2.324408818	,	5.214789579], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
