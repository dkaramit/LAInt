# LAInt
Local Adaptive 1D integration

Everything is in ```LAInt.hpp```. The class LocalAdaptive has two static methods: `partialIntegral` and `integrate`. The first calculates 
the integral in a cell, and subdivides it until two consecutive subdivisions give alsmost the same result. The second method, runs
the first in the entire integration region.


The signature of integrate is:

<code>
template< class LD > 

lai::LocalAdaptive< LD >::integrate(
                        funcType integrand, const std::vector<LD> &X, 
                        const LD & atol, const LD & rtol, const unsigned int & maxDivisions, 
                        ruleType rule)
</code>,

where `LD` is the numeric type,  `funcType = std::function<LD(LD x)>`, and `ruleType = std::function<LD(std::function<LD(LD)> integrand, const LD & x0, const LD & x1)>`. 

The arguments are:

1. integrand: the integrand.
1. X: integration cells. These are the ones that will be subdivided. 
1. atol,rtol:  absolute and relative tolerance. These are used to decide whether the integral in a cell has converged. 
1. maxDivisions: Maximum number of subdivisions of each cell. 

You can implement a new integration rule as 

<code>
template<class LD>
LD Rule(std::function<LD(LD)> integrand, const LD & x0, const LD & x1){
		return $\int_{x_0}^{x_1} \ dx$ integrand$(x)$;
}

</code>