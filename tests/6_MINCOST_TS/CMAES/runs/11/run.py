# 11
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[7.720148515,	2.819416667,	2.595828571,	3.30260521,	8.636773547], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
