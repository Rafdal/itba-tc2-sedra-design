import sympy as sp

def S(y, x):
    return (sp.diff(y, x)*(x/y))

def SensTable(y, name='y', subs_dict = {}, latex=False, sympy=False):
    if not isinstance(name, str):
        name = sp.latex(name)

    # get all symbols of y
    y_symbols = list(y.free_symbols)
    table = []
    labels = []
    # Create a table with the sensitivities
    for x in y_symbols:
        Sy_x = S(y, x).simplify()
        if len(subs_dict) > 0:
            Sy_x = Sy_x.subs(subs_dict).simplify()
        expr = sp.Eq(sp.Symbol(f'S_{{{sp.latex(x)}}}^{{{name}}}'), Sy_x)
        if sympy:
            print(sp.Eq(sp.Symbol(f'S_{{{sp.latex(x)}}}^{{{name}}}'), Sy_x))
        elif not latex:
            display(expr)
        else:
            print(sp.latex(expr))
        table.append(Sy_x)
        labels.append(sp.latex(x))
    return table, labels