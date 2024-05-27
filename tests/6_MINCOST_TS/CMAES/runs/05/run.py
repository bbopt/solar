# 05
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[3.837128713,	3.781041667,	1.437739286,	2.867635271,	8.15250501], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
