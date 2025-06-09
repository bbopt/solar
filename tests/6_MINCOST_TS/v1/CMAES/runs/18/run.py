# 18
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[5.367524752	,	3.297166667	,	7.907678571	,	4.000701403	,	6.233146293], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
