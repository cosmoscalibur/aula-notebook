# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 07:05:56 2016

@author: cosmoscalibur
"""
# Required for Jupyter Notebooks integration
from IPython.display import Latex, Markdown, display
# Required for building math expressions in text
from sympy import latex
# Required for Symbolic Algebra Computation
from sympy import symbols, Symbol, Function
from sympy.abc import x, f
from sympy import diff, integrate, simplify
# Algorithms
from interpolation import lagrange_polynomial, lagrange_interpolation


def tex_inline(sp_expr):
    return latex(sp_expr, mode='inline')
    
def tex_equation(sp_expr):
    return latex(sp_expr, mode='equation')
    
def tex_cat(list_expr):
    list_tex = [latex(expr) for expr in list_expr]
    return "$$" + str().join(list_tex) + "$$"
    
def dtx(expr):
    display(Latex(expr))
    
def dmd(expr):
    display(Markdown(expr))

def lagrange_interpolation_order_n(N):
    # n: Number of points to use in interpolation
    x_ = symbols('x_0:'+str(N))
    msg = "## Polinomio interpolante  \nComo primer paso se construyen los polinomios de Lagrange asociados" \
        + " a cada uno de los {} puntos $x_i$ solicitados para la interpolación. ".format(N) \
        + "Tenga en cuenta que los polinomios serán de grado {}.".format(N-1)
    dmd(msg)
    
    msg = "L_{i,n}(x) = \\prod_{j=0 \\ j\\neq i}^n \\frac{x-x_j}{x_i-x_j}"
    dtx("$$"+msg+"$$")
    cont = 0
    while cont < N:
        dtx(tex_cat(["L_{{" + str(cont) + "," + str(N) + "}}(x)=", \
        lagrange_polynomial(x_, cont, x)]))
        cont += 1
    msg = "Ahora, con los polinomios de Lagrange construidos y la evaluación de" \
        + " la función en los puntos, construimos el polinomio interpolante de" \
        + " grado {}, que suponemos como aproximación a la función.".format(N-1)
    dmd(msg)
    msg = "P_n(x) = \\sum_{i=0]^n f(x_i) L_{i,n}(x)"
    polynomial_approx = lagrange_interpolation(x_, f, x)
    dtx(tex_cat(["P_"+str(N)+"(x)=", polynomial_approx]))
    return polynomial_approx
    
def diff_methods(N):
    msg = "# Diferenciación numérica"
    dmd(msg)
    x_ = symbols('x_0:'+str(N))
    h = Symbol('h')
    polynomial_approx = lagrange_interpolation_order_n(N)
    msg = "## Derivada analítica \nUsando el polinomio interpolante, derivamos analíticamente respecto" \
        + " a la variable $x$, y obtenemos la siguiente aproximación.  "
    dmd(msg)
    diff_approx = diff(polynomial_approx, x)
    dtx(tex_cat(["\\frac{df(x)}{dx}\\approx\\frac{L_"+str(N)+"(x)}{dx}=", diff_approx]))
    msg = "Ahora, sustituimos la condición de puntos equidistantes con separación $h$ y factorizando."
    dmd(msg)
    msg_point = ", y posteriormente $x \\rightarrow x_0$."
    msg = "Diferencias hacia atras, $x_i \\rightarrow x_0 - ({}-i)h$".format(str(N-1)) + msg_point
    dmd(msg)
    point = {x:x_[0]}
    backward_points = [(x_[i], x_[0] - (N-1-i)*h) for i in range(N)]
    back_diff = simplify(diff_approx.subs(backward_points, simultaneous=True))
    back_diff = simplify(back_diff.subs(point))
    dtx(tex_cat(["\\left.\\frac{df(x)}{dx}\\right|_{x=x_0}\\approx", back_diff]))
    msg = "Diferencias hacia adelante, $x_i \\rightarrow x_0 + ih$" + msg_point
    dmd(msg)
    forward_points = [(x_[i], x_[0] + i*h) for i in range(N)]
    for_diff = simplify(diff_approx.subs(forward_points, simultaneous=True))
    for_diff = simplify(for_diff.subs(point))
    dtx(tex_cat(["\\left.\\frac{df(x)}{dx}\\right|_{x=x_0}\\approx", for_diff]))
    if N%2:
        msg = "Diferencias centradas, $x_i \\rightarrow x_0 + \\left(i - \\frac{{{} - 1}}{{2}}\\right)h$".format(N) + msg_point
        dmd(msg)
        center_points = [(x_[i], x_[0] + (i - (N - 1)//2)*h ) for i in range(N)]
        center_diff = simplify(diff_approx.subs(center_points, simultaneous=True))
        center_diff = simplify(center_diff.subs(point))
        dtx(tex_cat(["\\left.\\frac{df(x)}{dx}\\right|_{x=x_0}\\approx", center_diff]))
        
def integral_methods(N):
    msg = "# Integración numérica"
    dmd(msg)
    x_ = symbols('x_0:'+str(N))
    h = Symbol('h')
    polynomial_approx = lagrange_interpolation_order_n(N)
    msg = "## Integral analítica \nUsando el polinomio interpolante, integramos analíticamente respecto" \
        + " a la variable $x$ en el intervalo dado por $\\left[x_0, x_{}\\right]$ y obtenemos la siguiente aproximación.".format(N-1)
    dmd(msg)
    integral_approx = integrate(polynomial_approx, (x, x_[0], x_[-1]))
    dtx(tex_cat(["\\int_{{x_0}}^{{x_{0}}} f(x)dx\\approx \\int P_{0}(x)dx=".format(N-1), integral_approx]))
    msg = "Suponemos ahora puntos equidistantes de la forma $x_i \\rightarrow x_0 + ih$ y simplificamos," \
        + " factorizando $h$ con la fracción requerida para que el interior solo tenga números enteros."
    dmd(msg)
    x_points = [(x_[i], x_[0] + i*h) for i in range(N)]
    integral_method = simplify(integral_approx.subs(x_points))
    dtx(tex_cat(["\\int_{{x_0}}^{{x_{}}}f(x)dx\\approx", integral_method]))