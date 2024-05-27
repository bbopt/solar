# 25
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[8.420742574	,	9.451104167	,	8.763571429	,	5.107154309	,	2.818036072], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
