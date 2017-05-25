import numpy as np
def patho(M,N,S):
    for m in range(M,N,S):
        print("m is %d "%(m))
        for j in range(20):
            print("j is %.4f "%((0.5-1/(m*np.abs(j-1.05)))))
patho(1,20,1)