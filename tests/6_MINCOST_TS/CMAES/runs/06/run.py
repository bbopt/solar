# 06
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[4.017227723,	9.391375,	1.691835714,	9.072705411,	5.08759519], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
