


def partial(f, args):

    if type(args) == int:
        def new_f(*f_args):
            return f(*f_args, args)
    else:
        def new_f(*f_args):
            return f(*f_args, *args)

    return new_f

def multiply(a, b):
    return a * b

def multiply_by_factors(a, b, c):
    return a * b * c


def main():
    b = 2
    multi = multiply
    multiply_by_two = partial(multi, args=(b))
    a=3
    print(multiply_by_two(a))

    multiply_by_six = partial(multiply_by_factors, args=(2, 3))
    print(multiply_by_six(a))

main()