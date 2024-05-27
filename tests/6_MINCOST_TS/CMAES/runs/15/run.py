# 15
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[9.843712871	,	9.9573125	,	2.178139286	,	7.463406814	,	3.464849699], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
