import os
from test_norbs import test_norbs
from test_fermifab import test_fermistate, test_fermiop
from test_tensor_op import test_tensor_op

__all__ = ['test']

def test():
    test_norbs()
    test_fermistate()
    test_fermiop()
    test_tensor_op()
    print("All tests passed :)")


if __name__ == "__main__":
    test()