# 03
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[3.202475248,	4.168395833,	3.985,	9.728256513,	0.991625251], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
