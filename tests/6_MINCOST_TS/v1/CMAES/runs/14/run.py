# 14
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[6.967227723	,	4.280145833	,	5.747214286	,	0.502911824	,	5.896232465], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
