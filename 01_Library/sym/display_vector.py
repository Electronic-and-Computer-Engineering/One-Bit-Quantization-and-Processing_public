from IPython.display import display, Math
import sympy as symp

def display_vector(name: str, vector_expr, symbol_only=False):
    """
    Display a vector (Matrix) in LaTeX, either symbolically or evaluated.
    """
    if not symbol_only:
        display(Math(fr"\mathbf{{{name}}} = " + symp.latex(vector_expr)))
    else:
        display(Math(fr"\mathbf{{{name}}} \in \mathbb{{R}}^{{{vector_expr.rows}}}"))