# 09
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[5.887128713,	4.6925,	9.435821429,	5.578376754,	6.655791583], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
