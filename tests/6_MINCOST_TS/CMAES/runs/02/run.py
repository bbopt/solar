# 02
import cma
import bb as bb

x, es = cma.fmin_con2(bb.f,x0=[3.465940594,	5.1396875,	6.630928571,	7.06492986,	8.259338677], sigma0=1, constraints=bb.c, options={'tolfun':0,'bounds':[0,10],'maxfevals':5000,'seed':1}, restarts=1000)
