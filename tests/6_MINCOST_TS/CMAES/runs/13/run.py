# 13
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[2.029257426	,	0.73543125	,	3.866714286	,	8.037815631	,	2.120741483], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
