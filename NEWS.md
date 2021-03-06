# igcop 0.2.0

* Added a `NEWS.md` file to track changes to the package.
* Completed the collection of distributional quantities to match that of the CopulaModel package.
* Computations of quantile functions are now more robust. 
* Changed the parameterization of the copula families to make them more numerically stable.
    * Small values of the parameter `alpha` (previously = `k - 1`) are now reliable.
* Placed the theta and alpha parameters of the IG copula as their own arguments. 

# igcop 0.1.0

* Computational tools for IG and IGL copulas are now available, using the original parameterization of these copulas as defined in Coia (2017) thesis.
