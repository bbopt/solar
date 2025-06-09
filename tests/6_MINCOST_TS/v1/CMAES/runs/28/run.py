# 28
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[6.132623762	,	3.5643125	,	4.512607143	,	3.562645291	,	9.974088176], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
