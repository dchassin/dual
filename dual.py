"""Dual number arithmetic"""

# References:
#
# [1] https://en.wikipedia.org/wiki/Dual_number
#
# [2] https://www.researchgate.net/profile/Farid-Messelmi/publication/264300733_Analysis_of_Dual_Functions/
#     links/53d7c11e0cf2a19eee7fcf1c/Analysis-of-Dual-Functions.pdf
#

__version__ = "0.0"
__author__ = "dchassin@stanford.edu"

from parse import parse
from math import exp, log, sin, cos, sinh, cosh

indeterminatechar = "w"
precision = 1e-8
parseformat = f"{{:g}}{{:+g}}{indeterminatechar}"
strformat = f"{{self.x}}{{self.y:+}}{indeterminatechar}"
reprformat = f"<dual:{{:.4g}}{{:+.4g}}{indeterminatechar}>"


class DualException(Exception):
    pass


class dual:

    def __init__(self,
            x = 0.0,
            y = 0.0,
            ):
        """Construct a dual number

        Arguments:
       
        x (dual|float|int|str) - dual to copy, real part, int of real part, or
        string representation
       
        y (float) - indeterminate part
        """
        if type(x) is str:
            r = parse(parseformat, x)
            self.x = r[0]
            self.y = r[1]
        elif type(x) is dual:
            self.x = x.re()
            self.y = x.im()
        elif type(x) in [int, float] and type(y) in [int, float]:
            self.x = float(x)
            self.y = float(y)

    def __str__(self):
        return strformat.format(self.x, self.y).replace(" ", "")

    def __repr__(self):
        return reprformat.format(self.x, self.y).replace(" ", "")

    def __neg__(self):
        return dual(-self.x, -self.y)

    def __pos__(self):
        return dual(abs(self.x), abs(self.y))

    def __invert__(self):
        return dual(1) / self

    def __add__(self, other):
        if not type(other) is dual:
            other = dual(other)
        return dual(self.x+other.x, self.y+other.y)

    def __sub__(self, other):
        if not type(other) is dual:
            other = dual(other)
        return dual(self.x-other.x, self.y-other.y)

    def __mul__(self, other):
        if not type(other) is dual:
            other = dual(other)
        return dual(self.x*other.x, self.x*other.y+self.y*other.x)

    def __truediv__(self, other):
        if not type(other) is dual:
            other = dual(other)
        return dual(self.x/other.x, (self.y*other.x-self.x*other.y)/(other.x*other.x))

    def __pow__(self, n):
        if not type(n) is int:
            raise DualException("n must be an integer")
        elif n > 1:
            return dual(self.x**n, n*self.x**(n-1)*self.y)
        elif n == 1:
            return dual(self.x, self.y)
        elif n == 0:
            return dual(1.0)
        else:
            raise DualException("n must be non-negative")

    def __iadd__(self, other):
        if not type(other) is dual:
            other = dual(other)
        self.x += other.x
        self.y += other.y
        return self

    def __isub__(self, other):
        if not type(other) is dual:
            other = dual(other)
        self.x -= other.x
        self.y -= other.y
        return self

    def __imul__(self, other):
        if not type(other) is dual:
            other = dual(other)
        self.y *= other.x
        self.y += self.x*other.y
        self.x *= other.x
        return self

    def __itruediv__(self, other):
        if not type(other) is dual:
            other = dual(other)
        self.y *= other.x
        self.y -= self.x*other.y
        self.y /= other.x*other.x
        self.x /= other.x
        return self

    def __ipow__(self, other):
        if not type(other) is dual:
            other = dual(other)
        z = dual(self.pow(other))
        self.x = z.x
        self.y = z.y
        return self

    def __eq__(self, other):
        if type(other) in [int, float, str]:
            other = dual(other)
        return abs(self.x-other.x) <= precision and abs(self.y-other.y) <= precision

    def __ne__(self, otherother):
        if type(other) in [int, float, str]:
            other = dual(other)
        return abs(self.x-other.x) > precision or abs(self.y-other.y) > precision

    def __lt__(self, other):
        if not type(other) is dual:
            other = dual(other)
        return self.x < other.x

    def __le__(self, other):
        if not type(other) is dual:
            other = dual(other)
        return self.x <= other.x

    def __gt__(self, other):
        if not type(other) is dual:
            other = dual(other)
        return self.x > other.x

    def __ge__(self, other):
        if not type(other) is dual:
            other = dual(other)
        return self.x >= other.x

    def re(self):
        """Get real part"""
        return self.x

    def im(self):
        """Get indeterminate part"""
        return self.y

    def abs(self):
        """Get the magnitude"""
        return self.x

    def arg(self):
        """Get the angle"""
        return self.y / self.x

    def conj(self):
        """Get the conjugate number"""
        return dual(self.x, -self.y)

    def exp(self):
        """Exponent"""
        ex = exp(self.x)
        return dual(ex, self.y*ex)

    def pow(self, n):
        """Power"""
        c = dual(n).re()
        if c > 0:
            a = self.x
            b = self.y
            if type(n) in [int]:
                return dual(a**n, n*a**(n-1)*b)
            elif type(n) in [float, dual, str]:
                d = dual(n).im()
                ex = exp(c*log(a))
                return dual(ex, ex*(c*b/a+d*log(a)))
            else:
                raise DualException("n must be either int, float, dual, or str")
        elif c == 0:
            return dual(1)
        else:
            raise DualException("n must be non-negative")

    def D(self):
        """Derivative"""
        return self.y

    def sin(self):
        """Sine"""
        return dual(sin(self.x), self.y*cos(self.x))

    def cos(self):
        """Cosine"""
        return dual(cos(self.x), -self.y*sin(self.x))

    def tan(self):
        """Tangent"""
        return self.sin() / self.cos()

    def sinh(self):
        """Hypebolic sine"""
        return dual(sinh(self.x), self.y*cosh(self.x))

    def cosh(self):
        """Hyperbolic cosine"""
        return dual(cosh(self.x), self.y*sinh(self.x))

    def tanh(self):
        """Hyperbolic tangent"""
        return self.sinh() / self.cosh()

    def log(self):
        """Logarithm"""
        return dual(log(self.x), self.y/self.x)


w = dual(0, 1)

if __name__ == "__main__":

    import unittest

    x = dual(1, 2)
    y = dual(3, 4)

    class TestDual(unittest.TestCase):

        def test_add(self):
            self.assertEqual(x+y, dual(4, 6))
            self.assertEqual(y+x, dual(4, 6))

        def test_sub(self):
            self.assertEqual(x-y, dual(-2, -2))
            self.assertEqual(y-x, dual(2, 2))

        def test_mul(self):
            self.assertEqual(x*y, dual(3, 10))
            self.assertEqual(y*x, dual(3, 10))

        def test_div(self):
            self.assertEqual(x/y, dual(1/3, 2/9))
            self.assertEqual(y/x, dual(3, -2))

        def test_pow(self):
            self.assertEqual(x**0, 1)
            self.assertEqual(x**1, x)
            self.assertEqual(x**2, x*x)
            self.assertEqual(x**3, x*x*x)
            self.assertEqual(x.pow(0), 1)
            self.assertEqual(x.pow(1), x)
            self.assertEqual(x.pow(2), x*x)
            self.assertEqual(x.pow(3), x*x*x)
            self.assertEqual(x.pow(y), dual(1, 6))
            self.assertEqual(y.pow(x), dual(3, 10.59167373200866))

        def test_exp(self):
            self.assertEqual(x.log().exp(), x)
            self.assertEqual(x.exp().log(), x)

        def test_trig(self):
            self.assertEqual((-x).sin(), -(x.sin()))
            self.assertEqual((-x).cos(), x.cos())
            self.assertEqual(x.tan(), x.sin()/x.cos())
            self.assertEqual((x+y).sin(), x.sin()*y.cos()+x.cos()*y.sin())
            self.assertEqual((x+y).cos(), x.cos()*y.cos()-x.sin()*y.sin())
            self.assertEqual(x.cos()**2+x.sin()**2, 1)
            self.assertEqual((w*x).cos(), 1)
            self.assertEqual((w*x).sin(), w*x)

        def test_hyp(self):
            self.assertEqual(x.sinh(), (x.exp()-(-x).exp())/2.0)
            self.assertEqual(x.cosh(), (x.exp()+(-x).exp())/2.0)
            self.assertEqual(x.tanh(), x.sinh()/x.cosh())
            self.assertEqual(x.cosh()**2-x.sinh()**2, 1)

    unittest.main()
