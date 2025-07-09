# 21
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[2.654306931	,	0.3136	,	4.899464286	,	1.707823647	,	6.886873747], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
