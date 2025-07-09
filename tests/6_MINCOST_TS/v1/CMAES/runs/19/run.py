# 19
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[5.604405941	,	7.296229167	,	3.590178571	,	6.673727455	,	0.112946693], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
