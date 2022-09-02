import Steenrod
import Brown_Gitler as BG


### Steenrod algebra

A = Steenrod.Steenrod_algebra(2) # only the mod 2 Steenrod algebra is currently supported
a = A.adem([[1, 2], [2, 1]])
b = A.adem([[2, 1], [4]])
print("Steenrod algebra calculations:")
print("({}) + ({}) = {}".format(a, b, a + b))
print("({}) * ({}) = {}".format(a, b, a * b))
print()


### Brown-Gitler modules

J = [BG.Brown_Gitler_module(n) for n in range(9)]

n = 6
print("Bases for J_{}:".format(n)); J[n].printBases(); print()
print("Actions of Steenrod squares on basis elements of J_{}:".format(n)); J[n].printActions(); print()


# Specify an element with a list of items [m, c], where c is a list of coefficients with respect to the given basis of J(n)^m. Example:
x = J[6].element([ [2,[1]], [3,[0,1]] ])
y = J[6].element([ [2,[1]], [6,[1]] ])
print("Examples with Brown-Gitler modules:")
print("({}) + ({}) = {}".format(x, y, x + y))
print("({}) * ({}) = {}".format(a, x, a * x))
print()


T = BG.Brown_Gitler_polynomial_algebra()
# Specify an element by providing a list consisting of a list for each monomial. The first entry of the list is the coefficient, and the rest are the exponents of x_0, x_1, x_2, etc. Example:
x = T.element([[1], [1, 2], [1, 0, 0, 3, 0, 1]])
y = T.element([[1, 2]])
print("Examples with the polynomial algebra F_2[x_0, x_1, x_2, ...]:")
print("({}) + ({}) = {}".format(x, y, x + y))
print("({}) * ({}) = {}".format(x, y, x * y))
print("({}) * ({}) = {}".format(a, x, a * x))
print()

n = 6
print("Basis for J_{} as a submodule of T: {}".format(n, ", ".join([str(x) for x in T.basis(n)])))
print("Actions of Steenrod squares on basis elements of J_{}:".format(n)); T.printActions(n); print()
