# Implied Willow Tree -- released June 2021

This is Matlab code for Implied Willow Tree, associated with the paper

Bing Dong and Wei Xu. Implied Willow Tree Method for pricing, hedging and 
the term structure of moments risk premia, 2021.

We propose an implied willow tree method to explore the implied 
stochastic process of asset prices from its option prices with various 
strike prices and maturities under the risk-neutral (Q) measure, and 
transform to physical (P) measure given a CRRA-utility function. 

-------------------------------------------------------------------------------
Note:

The folder sourceFiles contains matlab scripts and c code to reproduce
the results for the implied willow tree. Current package supports Matlabs under
Windows and MacOS under 64-bit. You have to recompile   “f_hhh.c” into mex file 
for Matlab if running under other operating systems.

-------------------------------------------------------------------------------
## Folder
1. [Demo] contains some scripts to show (i) construct implied willow tree under Q and P measure, 
            (ii) price European, American and Asian options and greeks under implied willow tree.

2. [ImpliedWT] contains some functions used to (i) construct implied willow tree under Q and P 
            measure, (ii) calculate impiled moments under Q and P measure.

3. [OptionsPricing] contains some functions used to calculate European, American, Asian and Digital 
            options prices and greeks by willow tree method and Black-Scholes formula.

4. [Utils] contains some utility functions when using Johnson Curve inverse transformation.


## demo

1. demo_construct_ImpWT_underQ.m is an example for construct implied willow tree under Q-measure. 
    In this demo, we generate options prices under GBM by Black-Scholes formular. 
    You can replace the options prices and strike prices. 

2. demo_construct_ImpWT_underP.m is an example for construct implied willow tree under P-measure. 

3. demo_ImpWT_OptionPricing_Greeks.m is an example for options pricing and greeks calculation 
    based on implied willow tree, including European, American and Asian options.
