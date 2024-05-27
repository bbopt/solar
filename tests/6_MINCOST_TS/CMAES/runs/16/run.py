# 16
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[9.504306931	,	0.156110417	,	8.941607143	,	1.004296593	,	1.47751503], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
