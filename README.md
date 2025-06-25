# morphyd
MORPHology of HYDrides (morhpyd)

Authors: Teimur Aliev and Mikhail Kolesnik;
         Lebedev Physical Institute of the Russian Academy of Science, Moscow, Russia;
email: kolesnik.mikhail@gmail.com

The program morphyd simulates the nucleation and growth of hydrides in zirconium alloys during cooling of a specimen in 3D representation.
The result of the simulation allows estimating any metric of hydrides' morphology and thus predict the hydrogen embrittlement degree.

Call up TestScenario.m to start the simulation and postprocessing.
This file contains the history of temperature, geometry of the specimen, parameters of nucleation sites (density, spatial and defect energy distribution),
main control flow and the basic postprocessing.

NUCLGROWTH.m contains the governing equations of the mathematical model.

The code architecture is represented in code_architecture.pptx
