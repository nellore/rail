#!/usr/bin/env python

import random
import numpy as np


def interpolate(data):
	"""
		Returns a list of coefficients for a polynomial that approximates the values in data.

		data: Dictionary containing x and y values of data
		
		Return value: A list of coefficients describing a polynomial function that passes
			through all the data points
	"""

	num_pts = len(data)


	A = np.array([[0]*num_pts]*num_pts)
	b = np.array([0]*num_pts)
	
	for i in xrange(num_pts):
		x = data.keys()[i]
		y = data[x]

		b[i] = y

		for j in xrange(num_pts):
			A[i][num_pts-j-1] = pow(x, j)

	x = np.linalg.solve(A, b)
	return list(x)

def evaluate(coeffs, x):
	"""
		Return the value of the polynomial described by coeffs at x

		coeffs: Coefficients defining a polynomial. The equation for the polynomial is
			f(x) = coeffs[0] * x^n + coeffs[1] * x^(n-1) + ... + coeffs[n]
		x: Value at which to evaluate the function.

		Return value: The value of the polynomial at x.
	"""
	val = 0
	for j in xrange(len(coeffs)):
		val += coeffs[j] * pow(i, len(coeffs)-j-1)
	return val


num_pts = 21
xs = range(num_pts)
ys = [random.randint(0,20) for i in xrange(len(xs))]

data = dict()
for i in xrange(0, num_pts , 5):
	data[xs[i]] = ys[i]

coeffs = interpolate(data)
print coeffs
approx_ys = [0]*num_pts

for i in xrange(num_pts):
	#for j in xrange(len(coeffs)):
	#		approx_ys[i] += coeffs[j] * pow(i, len(coeffs)-j-1)
	approx_ys[i] = evaluate(coeffs, i)


for i in xrange(num_pts):
	print str(ys[i]) + '\t' + str(approx_ys[i])

