import json

from sympy import *
from sympy.parsing.sympy_parser import parse_expr

from prop_matrix import PropagatorCalculator
from shapes import ShapeFunction

import sys


class SolverInput:
    """
    Parses and encapsulates JSON input into an object with the following fields:
    `functions`, `shapes`, `ode`
    """

    def __init__(self, json_serialization):
        # define variables that are set as a part of the JSON deserialisation
        self.functions = []
        self.shapes = []
        self.ode = ""

        self.__dict__ = json.loads(json_serialization)


class SolverOutput:
    """
    The class is used to store computation results that are then easily converted into the json format
    """

    def __init__(self,
                 status,
                 solver,
                 propagator_elements,
                 ode_var_factor,
                 const_input,
                 ode_var_update_instructions):
        self.status = status
        self.solver = solver
        self.initial_values = []
        self.propagator_elements = propagator_elements
        self.ode_var_factor = ode_var_factor
        self.ode_var_update_instructions = ode_var_update_instructions
        self.shape_state_variables = []
        self.const_input = const_input
        self.updates_to_shape_state_variables = []
        self.shape_state_odes = []

    def decode_apostroph(self, ode):
        return

    def add_shape_state_variables(self, state_variables):
        self.shape_state_variables += state_variables

    def add_shape_state_odes(self, shape_state_odes):
        self.shape_state_odes += shape_state_odes

    def add_updates_to_shape_state_variables(self, updates_to_state_shape_variable):
        self.updates_to_shape_state_variables += updates_to_state_shape_variable

    def add_initial_values(self, initial_values):
        self.initial_values += initial_values


h = symbols("__h")


class OdeAnalyzer(object):
    """
    Orchestrates the execution of analysis activities which lead to a exact solution.
    """

    @staticmethod
    def is_linear_constant_coefficient_ode(ode_var, ode_rhs, shape_vars, shape_definitions, function_vars,
                                           function_definitions):
        for shape, shape_definition in zip(shape_vars, shape_definitions):
            exec ("{} = parse_expr(\"{}\")".format(shape, shape_definition))

        for function_var, function_definition in zip(function_vars, function_definitions):
            exec ("{} = parse_expr(\"{}\", local_dict=locals())".format(function_var, function_definition))

        # ode functions are defined in the local scope. it must be passed explicitly into the parse_exp,
        # otherwise, all functions in the ode would not be recognized as `complex` defined through the rhs symbols
        ode_var = parse_expr(ode_var)
        ode_rhs = parse_expr(ode_rhs, local_dict=locals())

        dvar = diff(ode_rhs, ode_var)
        dtdvar = diff(dvar, Symbol("t"))

        if simplify(dtdvar) == simplify(0):
            return True
        else:
            return False

    @staticmethod
    def compute_solution(input_json):
        input_ode_block = SolverInput(input_json)

        """
        The function computes a list with propagator matrices.
        :arguments A list starting with an ODE of the first order followed by shape definitions. An ODE is of the form 
        `var_name' = python_expression`. A shape function is of the form: `shape_name = python_expression`
        
        :returns JSON object containing all data necessary to compute an update step.
        """
        if len(input_ode_block.shapes) == 1:
            tmp = input_ode_block.shapes[0].split('=')
            shape_name = str(tmp[0])
            shape_expr = parse_expr(tmp[1].strip())
            if shape_expr.is_Function and str(shape_expr.func).startswith("delta"):
                # extract the name of the ODE and its defining expression
                ode_definition = input_ode_block.ode.split('=')  # it is now a list with 2 elements
                ode_var = ode_definition[0].replace("'", "").strip()
                ode_rhs = ode_definition[1].strip()

                function_vars = []  #
                function_definitions = []  # contains shape functions as ShapeFunction objects
                if "functions" in input_ode_block.__dict__:
                    for function_def in input_ode_block.functions:
                        tmp = function_def.split('=')
                        function_vars.append(tmp[0].strip())
                        function_definitions.append(tmp[1].strip())

                if OdeAnalyzer.is_linear_constant_coefficient_ode(ode_var,
                                                                  ode_rhs,
                                                                  [shape_name],
                                                                  [shape_expr],
                                                                  function_vars,
                                                                  function_definitions):
                    for function_var, function_definition in zip(function_vars, function_definitions):
                        exec ("{} = parse_expr(\"{}\", local_dict=locals())".format(function_var, function_definition))

                    # TODO discuss with Inga
                    ode_rhs_expr = parse_expr(ode_rhs, local_dict=locals())
                    const_input = simplify(1 / diff(ode_rhs_expr, parse_expr(shape_name)) * (
                        ode_rhs_expr - diff(ode_rhs_expr, parse_expr(ode_var)) * parse_expr(ode_var)) - (
                        parse_expr(shape_name)))

                    c1 = diff(ode_rhs_expr, parse_expr(ode_var))
                    # The symbol must be declared again. Otherwise, the right hand side will be used for the derivative
                    c2 = diff(ode_rhs_expr, parse_expr(shape_name))

                    tau_constant = shape_expr.args[1] # is is passed as the second argument of the delta function
                    ode_var_factor = exp(-h/tau_constant)
                    ode_var_update_instructions = [
                        ode_var + " = __ode_var_factor * " + str(ode_var),
                        ode_var + " += " + str(str(simplify(c2 / c1 * (exp(h * c1) - 1)))) + " * __const_input"]
                    result = SolverOutput(
                        "success",
                        "delta",
                        None,
                        {"__ode_var_factor": str(ode_var_factor)},
                        {"__const_input": str(const_input)},
                        ode_var_update_instructions)
                    return json.dumps(result.__dict__, indent=2)
                return None

        shape_functions = []  # contains shape functions as ShapeFunction objects
        for shape_function in input_ode_block.shapes:
            tmp = shape_function.split('=')
            lhs_var = tmp[0].strip()
            shape_expr = tmp[1].strip()
            shape_functions.append(ShapeFunction(lhs_var, shape_expr))

        if input_ode_block.ode is None:
            result = OdeAnalyzer.convert_shapes_to_odes(shape_functions)
            return json.dumps(result.__dict__, indent=2)

        # extract the name of the ODE and its defining expression
        ode_definition = input_ode_block.ode.split('=')  # it is now a list with 2 elements
        ode_var = ode_definition[0].replace("'", "").strip()
        ode_rhs = ode_definition[1].strip()

        function_vars = []  #
        function_definitions = []  # contains shape functions as ShapeFunction objects
        if "functions" in input_ode_block.__dict__:
            for function_def in input_ode_block.functions:
                tmp = function_def.split('=')
                function_vars.append(tmp[0].strip())
                function_definitions.append(tmp[1].strip())

        if OdeAnalyzer.is_linear_constant_coefficient_ode(ode_var,
                                                          ode_rhs,
                                                          [shape.name for shape in shape_functions],
                                                          [shape.shape_expr for shape in shape_functions],
                                                          function_vars,
                                                          function_definitions):

            return OdeAnalyzer.compute_exact_solution(function_definitions,
                                                      function_vars,
                                                      ode_rhs,
                                                      ode_var,
                                                      shape_functions)
        else:  # is_linear_constant_coefficient_ode evaluates to false
            result = OdeAnalyzer.convert_shapes_to_odes(shape_functions)
            return json.dumps(result.__dict__, indent=2)

    @staticmethod
    def compute_exact_solution(function_definitions, function_vars, ode_rhs, ode_var, shape_functions):
        calculator = PropagatorCalculator()
        prop_matrices, const_input, step_const = calculator.ode_to_prop_matrices(
            shape_functions,
            ode_var,
            ode_rhs,
            function_vars,
            function_definitions)
        propagator_elements, ode_var_factor, const_input, ode_var_update_instructions = \
            calculator.prop_matrix_to_prop_step(
                prop_matrices,
                const_input,
                step_const,
                shape_functions,
                ode_var)
        # build result JSON
        result = SolverOutput("success",
                              "exact",
                              propagator_elements,
                              ode_var_factor,
                              const_input,
                              ode_var_update_instructions)
        for shape in shape_functions:
            result.add_shape_state_variables(shape.additional_shape_state_variables())
            result.add_initial_values(shape.get_initial_values())
            result.add_updates_to_shape_state_variables(shape.get_updates_to_shape_state_variables())
        return json.dumps(result.__dict__, indent=2)

    @staticmethod
    def convert_shapes_to_odes(shape_functions):
        result = SolverOutput("success", "numeric", None, None, None, None)
        for shape in shape_functions:
            result.add_shape_state_odes(shape.nestml_ode_form)
            result.add_shape_state_variables(shape.additional_shape_state_variables())
            result.add_initial_values(shape.get_initial_values())
        return result


# MAIN ENTRY POINT ###
if __name__ == "__main__":
    result = OdeAnalyzer.compute_solution(sys.argv[1])
    f = open('result.tmp', 'w')
    f.write(result)
