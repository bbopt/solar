# 20
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[8.502475248	,	5.190458333	,	9.87725	,	8.298436874	,	9.297014028], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
