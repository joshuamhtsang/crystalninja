# crystalninja README.md


###################
### Description ###
###################

A tool for creating and modifying crystal structures relevant to materials science.
Outputs can be used for atomic simulation.


########################
### virtualenv setup ###
########################

This project uses Python virtualenv to simplify environment consistency.  Simply

> pip3 install -r requirements.txt

to install the needed dependencies.  If, during development, you introduce new 
package dependencies, then ensure you generate a new requirements.txt:

> pip3 freeze --user > requirements.txt


##################
### How to run ###
##################

To run, you need to create a *.py script that imports crystal ninja as well as numpy:

	import numpy as np
	import crystalninja as cn

Look at test.py for an example.  To execute, execute in terminal using Python 3:

> python3 test.py
