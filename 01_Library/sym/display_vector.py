from IPython.display import display, Math
import sympy as symp

def display_vector(name: str, symbolic_vec, values: list = None, transpose=True, precision=3):
    """
    Display a symbolic vector in LaTeX. Optionally shows evaluated version with numerical values.

    Parameters:
    - name: Variable name (e.g., 'v_x')
    - symbolic_vec: sympy.Matrix of symbols
    - values: list of numerical values to substitute (optional)
    - transpose: if True, display as row vector with T superscript
    - precision: number of digits for numerical display
    """
    latex_transpose = r"^{\mathsf{T}}" if transpose else ""
    vec_disp = symbolic_vec.T if transpose else symbolic_vec

    # Display symbolic version
    display(Math(fr"\mathbf{{{name}}} = {symp.latex(vec_disp)}{latex_transpose}"))

    # Display numeric version if values are provided
    if values is not None:
        assert len(symbolic_vec) == len(values), "Length mismatch between vector and values"
        subs_dict = {symbolic_vec[i]: values[i] for i in range(len(values))}
        evaluated = symbolic_vec.subs(subs_dict).evalf(precision)
        evaluated_disp = evaluated.T if transpose else evaluated
        display(Math(fr"\mathbf{{{name}}} = {symp.latex(evaluated_disp)}{latex_transpose}"))