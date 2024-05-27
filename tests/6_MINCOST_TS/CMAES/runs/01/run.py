import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[1.069158416,6.438916667,0.309107143,1.304803607,7.393046092], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
