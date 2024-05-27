# 24
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[5.088316832	,	8.725375	,	7.047535714	,	9.248857715	,	2.28002004], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
