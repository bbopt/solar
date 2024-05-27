# 17
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[7.117425743	,	6.22775	,	6.760892857	,	8.713426854	,	2.432284569], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
