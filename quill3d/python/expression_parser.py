import re
import numpy as np

number_pattern = '^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$'
variable_pattern = '^([x-z]|u[x-z]|v[x-z]|g|chi|t|xi|phi|theta)$'
operator_pattern = '^([+\-*/\^<>=]|neq|geq|leq|&&|\|\||and|or)$'
function_pattern = '^(sin|cos|exp|log|min|max|any|all)$'
split_by_operator_pattern = '([+\-*/\^<>()=]|neq|geq|leq|&&|\|\||and|or|)'
split_by_number_pattern = '([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)'


# Reverse Polish notation
def to_polish(expr_text, should_log=False):
    def operator_priority(op):
        if op in '()':
            return 0
        elif op in '&&||andor':
            return 1
        elif op in ('>', '<', '=', 'neq', 'leq', 'geq'):
            return 2
        elif op in '+-':
            return 3
        elif op in '*/':
            return 4
        elif op == '^':
            return 5
        elif re.match(function_pattern, op):
            return 6
    operators = []
    result = []

    if expr_text is None or expr_text == '':
        return []

    expr_text = expr_text.replace(' ', '').replace('==', '=').replace('**','^')\
        .replace('>=', 'geq').replace('<=', 'leq').replace('!=', 'neq').lower()

    # First find numbers (without leading plus/minus) and split by them
    tokens = re.split(split_by_number_pattern, expr_text)
    # print(tokens)
    tmp = [re.split(split_by_operator_pattern, t) if not re.match(number_pattern, t) else [t] for t in tokens]
    tokens = list(filter(lambda t: t != '', [t for t2 in tmp for t in t2]))
    # Insert zeros before unary negations
    minuses = [i for i, v in enumerate(tokens)
               if v == '-' and (i == 0 or i > 0 and (re.match(operator_pattern, tokens[i-1]) or tokens[i-1] == '('))]
    for i in range(len(minuses)):
        tokens.insert(minuses[i], '0')
        minuses = list(map(lambda x: x+1, minuses))
    # print(tokens)

    if should_log:
        print('Parsed: {0}'.format(tokens))

    for token in tokens:
        if re.match(number_pattern, token):
            result.append(token)
        elif re.match(variable_pattern, token):
            result.append(token)
        elif token == '(':
            operators.append(token)
        elif token == ')':
            if '(' not in operators:
                raise ValueError('Opening bracket is missing!')
            while len(operators) > 0:
                op = operators.pop()
                if op == '(':
                    break
                result.append(op)
        elif re.match(function_pattern, token):
            operators.append(token)
        elif re.match(operator_pattern, token):
            while len(operators) > 0:
                op = operators.pop()
                if operator_priority(token) > operator_priority(op):
                    operators.append(op)
                    break
                result.append(op)
            operators.append(token)
        else:
            raise ValueError('Unrecognized character: {0}'.format(token))
    while len(operators) > 0:
        op = operators.pop()
        if not re.match(operator_pattern, op) and not re.match(function_pattern, op):
            raise ValueError('Closing bracket is missing')
        result.append(op)

    if should_log:
        print('Rev. Polish notation: {0}'.format(result))
    return result


def get_vars(expr):
    return list(filter(lambda token: re.match(variable_pattern, token), expr))


def evaluate_var(name, var_evaluator=None):
    if var_evaluator is None:
        print('evaluate_var({0}): variable evaluator is None, returning zero'.format(name))
        return 0.0
    else:
        return var_evaluator(name)


def operand_to_float(operand, var_evaluator=None):
    if isinstance(operand, float) or isinstance(operand, np.ndarray):
        return operand
    elif re.match(variable_pattern, operand):
        return evaluate_var(operand, var_evaluator)
    elif re.match(number_pattern, operand):
        return float(operand)
    else:
        raise ValueError('Unable to parse {0}'.format(operand))


def evaluate_binary(operation, operand2, operand1, var_evaluator=None):
    operand1 = operand_to_float(operand1, var_evaluator)
    operand2 = operand_to_float(operand2, var_evaluator)

    if operation == '+':
        return operand1 + operand2
    elif operation == '-':
        return operand1 - operand2
    elif operation == '*':
        return operand1 * operand2
    elif operation == '/':
        return operand1 / operand2
    elif operation == '^':
        return operand1 ** operand2
    elif operation == '&&' or operation == 'and':
        return operand1 & operand2
    elif operation == '||' or operation == 'or':
        return operand1 | operand2
    elif operation == '=' or operation == '==':
        return operand1 == operand2
    elif operation == 'neq':
        return operand1 != operand2
    elif operation == '>':
        return operand1 > operand2
    elif operation == '<':
        return operand1 < operand2
    elif operation == 'geq':
        return operand1 >= operand2
    elif operation == 'leq':
        return operand1 <= operand2


def evaluate(expr, var_evaluator=None, should_log=False):
    if isinstance(expr, str):
        return evaluate(to_polish(expr, should_log))

    stack = []
    for token in expr:
        if re.match(number_pattern, token):
            stack.append(token)
        elif re.match(variable_pattern, token):
            stack.append(token)
        elif re.match(operator_pattern, token):
            if len(stack) < 2:
                raise ValueError('Unable to evaluate expression!')
            stack.append(evaluate_binary(token, stack.pop(), stack.pop(), var_evaluator))
        elif re.match(function_pattern, token):
            if token == 'sin':
                stack.append(np.sin(operand_to_float(stack.pop(), var_evaluator)))
            elif token == 'cos':
                stack.append(np.cos(operand_to_float(stack.pop(), var_evaluator)))
            elif token == 'exp':
                stack.append(np.exp(operand_to_float(stack.pop(), var_evaluator)))
            elif token == 'log':
                stack.append(np.log(operand_to_float(stack.pop(), var_evaluator)))
            elif token == 'min':
                val = operand_to_float(stack.pop(), var_evaluator)
                if isinstance(val, np.ndarray) and len(val) == 0:
                    stack.append(val)
                else:
                    stack.append(np.min(val))
            elif token == 'max':
                val = operand_to_float(stack.pop(), var_evaluator)
                if isinstance(val, np.ndarray) and len(val) == 0:
                    stack.append(val)
                else:
                    stack.append(np.max(val))
            elif token == 'any' or token == 'all':
                val = operand_to_float(stack.pop(), var_evaluator)
                if isinstance(val, np.ndarray):
                    stack.append(np.any(val) if token == 'any' else np.all(val))
                else:
                    stack.append(val)
        else:
            raise ValueError('Unrecognized function: {0}'.format(token))
    if len(stack) != 1:
        raise ValueError('Invalid expression')
    return stack.pop()


def test_parser():
    print(parser.evaluate('2.54E2 + 8/3 * (23 * (38.23*2/7)) - 2 * exp(2)', True))
    print(parser.evaluate('2 * exp(2)', True))
    print(parser.evaluate('(sin(1 - 2*(1+1)/2*cos(0)) + 1) / 3 + cos(0)', True))
    print(parser.evaluate('3-5.1e-6', True))
    print(parser.evaluate('sin(0) >= 0', True))
    print(parser.evaluate('(sin(1 - 3/2*exp(2)) + 1) / 2'))

