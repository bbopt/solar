import os

n=5
m=6

bbe=0
last_x=n*[0.0]
bbx=n*[0.0]
obj=1e20
cstr=m*[0.0]

lb=[793.0,  2.0,  2.0, 0.01, 0.01]
ub=[995.0, 50.0, 30.0, 5.00, 5.00]

def c(x):
    global bbe, last_x, cstr, obj, bbx
    last_x=x
    bbe=bbe+1
    obj=1e20

    #print("Constraints: bbe=",bbe," x=",x)

    # unscaling: [0;10] --> [lb;ub]:
    for i in range(n):
        bbx[i]=lb[i]+(x[i]/10.0)*(ub[i]-lb[i])

    #print("Constraints: bbe=",bbe," bbx=",bbx)

        
    # blackbox evaluation:
    
    with open('x_tmp.txt', 'w') as x_file:
        for i in range(n):
            print(bbx[i]," ",file=x_file,end="")
        x_file.close()

    os.system("nice $SOLAR_HOME/bin/solar 6 x_tmp.txt > solar_output_tmp.txt")

    with open('solar_output_tmp.txt', 'r') as solar_output:
        s=solar_output.read()
        solar_output.close()

    s.strip("\n")
    tmp=s.split(" ")

    obj=float(tmp[0])
    
    for j in range(m):
        cstr[j]=float(tmp[j+1])
    
    return cstr

def f(x):

    #print("Objective  : bbe=",bbe," x=",x)

    if bbe==1:
        with open('out.txt', 'a') as out:
            print("bbe sol bbo obj",file=out)
            out.close()

    for i in range(n):
        if last_x[i] != x[i]:
            raise Exception("Error: not the right x in bb.f(): ",x," except of ",last_x)
        
    with open('out.txt', 'a') as out:
        print(bbe," ",bbx," ",cstr," ",obj,file=out)
        out.close()
        
    return obj

