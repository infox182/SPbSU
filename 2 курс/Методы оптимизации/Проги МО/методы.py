from math import exp, sqrt

#Метод Дихотомии
def dichotomy(f, a, b, e, d):
	e1 = (b - a) / 2

	if e1 <= e:
		x = (a + b) / 2
		ff = f(x, a, b)

		print('x* ≈ %f\nf* = %f' % (x, ff))
	else:
		x1 = (a + b - d) / 2
		x2 = (a + b + d) / 2

		if f(x1, a, b) <= f(x2, a, b):
			b = x2
		else:
			a = x1

		dichotomy(f, a, b, e, d)

#Метод золотого сечения
def golden_ratio(f, a, b, e, x1, x2):
	e1 = (b - a) / 2

	if e1 <= e:
		x = (a + b) / 2
		ff = f(x, a, b)

		print('x* ≈ %f\nf* = %f' % (x, ff))
	else:
		#x1 = a + ((3 - sqrt(5)) / 2) * (b - a)
		#x2 = a + ((sqrt(5) - 1) / 2) * (b - a)

		if f(x1, a, b) <= f(x2, a, b):
			b = x2 + 0
			x2 = x1 + 0
			x1 = a + ((3 - sqrt(5)) / 2) * (b - a)

		else:
			a = x1 + 0
			x1 = x2 + 0
			x2 = a + ((sqrt(5) - 1) / 2) * (b - a)

		golden_ratio(f, a, b, e, x1, x2)


#Данные
f = lambda x, a, b: a / exp(x) + b * x
a = 219 / 20
b = 3 / 2
e = 0.001
d = e * 1.5

print('Метод Дихотомии')
dichotomy(f, a, b, e, d)
print('Метод золотого сечения')
golden_ratio(f, a, b, e, a + ((3 - sqrt(5)) / 2) * (b - a), a + ((sqrt(5) - 1) / 2) * (b - a))