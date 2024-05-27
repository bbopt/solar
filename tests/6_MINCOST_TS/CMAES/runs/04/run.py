# 04
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[1.627475248,	1.648133333,	7.280107143,	6.558717435,	3.099659319], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
