from IPython.display import display, Math
import sympy as symp

def display_vector(name: str, vector_expr, symbol_only=False, transpose=True):
    """
    Display a vector in LaTeX. If transpose=True, shows as row vector with ^T.
    """
    if symbol_only:
        dim = f"{vector_expr.rows}" if not transpose else f"{vector_expr.cols}"
        display(Math(fr"\mathbf{{{name}}} \in \mathbb{{R}}^{{{dim}}}"))
        return

    if transpose:
        # Convert column vector to row for display
        row_expr = vector_expr.T
        latex_vec = symp.latex(row_expr)
        display(Math(fr"\mathbf{{{name}}} = {latex_vec}^{{\mathsf{{T}}}}"))
    else:
        display(Math(fr"\mathbf{{{name}}} = " + symp.latex(vector_expr)))