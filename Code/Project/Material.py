from sympy import *
import numpy as np
import sys
import matplotlib.pyplot as plt


class MaterialModel:
    def __init__(self, id, name, no_of_ranges, range_upper_limits, formulas_list):
        self.no_of_ranges = no_of_ranges
        self.id = id
        self.name = name
        self.range_upper_limits = np.array(range_upper_limits, dtype=Function)
        self.range_upper_limits = sympify(self.range_upper_limits)
        self.formulas = np.empty(no_of_ranges, dtype=Function)
        self.d_formulas = np.empty(no_of_ranges, dtype=Function)
        index = 0
        for formula in formulas_list:
            sympy_formula = sympify(formula)
            f_sympy_formula = lambdify(x, sympy_formula)
            self.formulas.put(index, f_sympy_formula)

            d_sympy_formula = sympy_formula.diff(Symbol('x'))
            f_d_sympy_formula = lambdify(x, d_sympy_formula)
            self.d_formulas.put(index, f_d_sympy_formula)

            index += 1
        self.plot()

    def get_stress(self, strain):
        for index in range(self.no_of_ranges - 1):
            # print("stress:",stress, " range_upper_limits[index]:",self.range_upper_limits[index])
            if strain < self.range_upper_limits[index]:
                # print('model id:',self.id,'range index', index)
                return self.formulas[index](strain)
        # print('model id:',self.id,'range index', 2)
        return self.formulas[-1](strain)

    def get_e(self, strain):
        for index in range(self.no_of_ranges - 1):
            if strain < self.range_upper_limits[index]:
                return self.d_formulas[index](strain)
        return self.d_formulas[-1](strain)

    def plot(self):

        plt.xlabel("strain")
        plt.ylabel("stress")
        plt.title(self.name)
        lower_limit = -0.005
        range_id = 0
        for formula in self.formulas:
            if range_id == self.no_of_ranges - 1:
                upper_limit = 0.005
            else:
                upper_limit = self.range_upper_limits[range_id]
                upper_limit = eval(str(upper_limit))
            range_id += 1
            x_values = np.linspace(lower_limit, upper_limit, 100)
            # print("lower:", lower_limit, "upper:", upper_limit)
            lower_limit = upper_limit
            y_values = formula(x_values)
            plt.plot(x_values, y_values)
        plt.show()


x = Symbol('x')
material_models = None


def load_material_models(js):
    global material_models
    no_of_material_models = js["no_of_material_models"]
    material_models = np.empty(no_of_material_models, dtype=MaterialModel)
    material_models_js = js["material_models"]
    for material_model_js in material_models_js:
        id = material_model_js["id"]
        name = material_model_js["name"]
        no_of_ranges = material_model_js["no_of_ranges"]
        range_upper_limits = material_model_js["range_upper_limits"]
        formulas = material_model_js["formulas"]
        material_model = MaterialModel(id, name, no_of_ranges, range_upper_limits, formulas)
        material_models.put(id, material_model)
