# 22
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[2.740594059	,	2.448979167	,	3.195821429	,	6.338076152	,	9.476352705], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
