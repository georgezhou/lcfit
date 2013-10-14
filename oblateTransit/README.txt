Instructions for this package

First thing first: Compile
This is a package need compile. 
make cc make the c part, which is only used for testing. 

make all make the python part.
So incase to use the program, 
>> make all

Then you will see a bunch of autometic compiled files.

To generate a simulated residual of a oblate planet. 

1) edit the paramters in test.py 
2) python test.py

To generate a simulated lightcurve, 
1) edit the parameters in lcgen.py 
function: trangen 
and function :main

2) python lcgen.py

----------------------------------
I don't got time to write a configure file for the parameters yet, which 
will be done in the near future.

