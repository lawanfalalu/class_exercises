%prologcodes
factorial(0,1).
factorial(X,Y):- X1 is X - 1, Y1 is Y * factorial(X1,Y).