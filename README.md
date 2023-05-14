# econ_SEIR_example
Basic covid SEIR model

SEIR_code.R contains the code to create and SEIR model with set parameters. It uses a differential equation solver.

SEIR_pomp.R is code to calibrate an SEIR model. It uses the pomp package. First, the code writes the SEIR model. 
It then uses C snippets to speed up the calibration, and completes calibration.

InClass.R is the code for an in class exercise. For this, we add a simple vaccination scenario to the code in SEIR_code.R.
