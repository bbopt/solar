# 12
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[2.140445545	,	8.4465625	,	5.928571429	,	9.567715431	,	8.795971944], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
