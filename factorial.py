#factorial
def fact(n):
	if n == 1
		return 1
	return fact(n-1)

	
def factClass(m):
	result = 1
	for i in range(1,m+1):
		result = result * i
	print("The result of %d factorial is %.4E "%(m,result))