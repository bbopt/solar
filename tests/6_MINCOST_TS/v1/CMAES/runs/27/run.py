# 27
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[8.112227723	,	1.263302083	,	0.804642857	,	2.451923848	,	4.767174349], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
