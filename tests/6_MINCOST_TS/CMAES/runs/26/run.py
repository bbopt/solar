# 26
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[6.600841584	,	7.637604167	,	8.257428571	,	5.402044088	,	4.25248497], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
