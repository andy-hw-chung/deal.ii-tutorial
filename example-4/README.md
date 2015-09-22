# Example 4 - Run-time parameters

This example shows how to make use of the `deal.ii` `ParameterHandler` object to handle reading parameters at run-time. It does the same things as `example-3`, but now we can configure some parameters to change without having to compile each time we make a change.

In this example, the parameters are handled in the `Parameters` class. It has a `ParameterHandler` object which can be used to declare the parameters and then read them from the given file. A reference to `Parameter` object is passed to the `Application` which in turn is passed to the `LaplaceProblem`. 

The program looks for the name of the parameter file passed in through the command line. If none is given, the default parameter file name `"parameters.param"` is used. 
