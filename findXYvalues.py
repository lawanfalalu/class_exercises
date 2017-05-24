import math
def findXYvalues(min,max):
    for x in range(min,max+1):
        f = math.pow(x,2)-3*x -4
        print(x,int(f))