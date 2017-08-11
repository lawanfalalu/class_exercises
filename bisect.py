#

def main():
    print(bisect())

def bisect():
    int1 = float(3)
    int2 = float(4.5)
    f1 = float(int1**2 - 3*int1 - 4)
    f2 = float(int2**2 -3*int2 - 4)
    prod = float(f1 * f2)
    err = 1e-16
    
    if prod < 0:
        while(abs(int1 - int2)>err):        
            midint = float( (int1+int2)/2)                   
            f3 = float(midint**2 -3*midint -4 )           
            int3 = f3 * f1        
            if  int3 < midint:
                int2 = midint
            else:
                int1 = midint
            print(midint)
    else:
        return 

def bisect2():
    int1 = float(3)
    int2 = float(4.5)
    f1 = float(int1**2 - 3*int1 - 4)
    f2 = float(int2**2 -3*int2 - 4)
    prod = float(f1 * f2)
    err = 1e-16
    

if __name__ == '__main__': main()