from IPython.display import display, Math
import sympy as symp

def display_vector(name: str,
                   symbolic_vec,
                   values: list = None,
                   transpose=True,
                   precision=3,
                   show_symbolic=True,
                   show_values=True):
    """
    Display a symbolic vector in LaTeX, optionally with numerical values and transpose.

    Parameters:
    - name: variable name (e.g. 'v_x')
    - symbolic_vec: sympy.Matrix of symbols
    - values: list of numerical values (optional)
    - transpose: whether to display as a row vector with superscript T
    - precision: rounding precision for numeric display
    - show_symbolic: whether to show symbolic form
    - show_values: whether to show numeric (evaluated) form
    """
    assert show_symbolic or show_values, "At least one of show_symbolic or show_values must be True"
    
    latex_T = r"^{\mathsf{T}}" if transpose else ""
    vec_disp = symbolic_vec.T if transpose else symbolic_vec

    # Show symbolic version
    if show_symbolic:
        display(Math(fr"\mathbf{{{name}}} = {symp.latex(vec_disp)}{latex_T}"))

    # Show numeric version
    if show_values and values is not None:
        assert len(symbolic_vec) == len(values), "Length mismatch between symbolic vector and value list"
        subs_dict = {symbolic_vec[i]: values[i] for i in range(len(values))}
        evaluated = symbolic_vec.subs(subs_dict).evalf(precision)
        evaluated_disp = evaluated.T if transpose else evaluated
        display(Math(fr"\mathbf{{{name}}} = {symp.latex(evaluated_disp)}{latex_T}"))