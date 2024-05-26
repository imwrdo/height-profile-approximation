import aproxMethods
import matplotlib.pyplot as plt

if __name__ == "__main__":
    for i in [10,30,70]:
        aproxMethods.LagrangeMethod(i,use_chebyshev=False)
        aproxMethods.LagrangeMethod(i,use_chebyshev=True)
        aproxMethods.splineMethod(i)