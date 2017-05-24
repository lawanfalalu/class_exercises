#summation
def summation(n):
	if n == 0:
		return 0
	return n**2 + summation(n-1)