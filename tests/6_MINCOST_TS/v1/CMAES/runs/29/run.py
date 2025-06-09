# 29
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[4.267871287	,	6.852770833	,	1.207017857	,	7.77739479	,	7.141442886], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
