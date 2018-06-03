"""
Takes a positive integer as an argument and prints all of its divisors.

example usage:
    python divisbles.py 100

""" 

import sys


def is_divisible(a, b):
    """Determines if integer a is divisible by integer b."""
    
    remainder = a % b
    # if there's no remainder, then a is divisible by b
    if not remainder:
        return True
    else:
        return False


def find_divisors(integer):
    """Find all divisors of an integer and return them as a list."""

    divisors = []
    # we know that an integer divides itself
    divisors.append(integer)
    # we also know that the biggest divisor other than the integer itself
    # must be at most half the value of the integer (think about it)
    divisor = integer / 2

    while divisor >= 0:
        if is_divisible(integer, divisor):
            divisors.append(divisor)
        divisor =- 1

    return divisors


if __name__ == '__main__':
    # do some checking of the user's input
    try:
        if len(sys.argv) == 2:
            # the following statement will raise a ValueError for
            # non-integer arguments
            test_integer = int(sys.argv[1])
            # check to make sure the integer was positive
            if test_integer <= 0:
                raise ValueError("integer must be positive")
        elif len(sys.argv) == 1:
            # print the usage if there are no arguments and exit
            print __doc__
            sys.exit(0)
        else:
            # alert the user they did not provide the correct number of
            # arguments
            raise ValueError("too many arguments")
    # catch the errors raised if sys.argv[1] is not a positive integer
    except ValueError, e:
        # alert the user to provide a valid argument
        print "Error: please provide one positive integer as an argument."
        sys.exit(2)

    divisors = find_divisors(test_integer)
    # print the results
    print "The divisors of %d are:" % test_integer
    for divisor in divisors:
        print divisor