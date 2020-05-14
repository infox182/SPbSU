A = [[8, 1, 1, 1], [1, 10, 1, 1], [1, 1, 12, 1], [1, 1, 1, 14]]
b = [10, 12, 14, 16]
B = [[8, 1, 1, 1, 10], [1, 10, 1, 1, 12], [1, 1, 12, 1, 14], [1, 1, 1, 14, 16]]


def mult(ar, val):
	a = [0, 0, 0, 0, 0]
	for i in range(len(ar)):
		a[i] = ar[i]*val
	return a

def diff(ar1, ar2):
	a = [0, 0, 0, 0, 0]
	for i in range(len(ar1)):
		a[i] = ar1[i] - ar2[i]
	return a
	
for i in range(len(A)):
	for j in range(1, len(A)-i):
		B[j+i] = diff(B[j+i], mult(mult(B[i], 1/B[i][i]), B[j+i][i]))
		

x4 = B[3][4]/B[3][3]
x3 = (B[2][4] - x4*B[2][3])/B[2][2]
x2 = (B[1][4] - x4*B[1][3] - x3*B[1][2])/B[1][1]
x1 = (B[0][4] - x4*B[0][3] - x3*B[0][2] - x2*B[0][1])/B[0][0]

x = [x1, x2, x3, x4]

print(x)	
	
