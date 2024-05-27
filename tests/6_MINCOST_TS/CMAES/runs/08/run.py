# 08
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[4.596485149,	2.406166667,	6.291964286,	4.556753507,	3.901683367], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
