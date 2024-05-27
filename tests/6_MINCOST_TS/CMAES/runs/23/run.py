# 23
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[8.947722772	,	7.904541667	,	9.182607143	,	6.044028056	,	1.722757515], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
