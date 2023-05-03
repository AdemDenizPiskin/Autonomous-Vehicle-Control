# Autonomous-Vehicle-Control
This project was part of my EE498 course. The object was lane keeping and obstacle avoidance using MPC. It is assumed that speed of the vehicle is constant. 

The implementation is a little bit simplified version of the algorithm in this article: Turri, V., Carvalho, A., Tseng, H.E., Johansson, K.H. and Borrelli, F., 2013, October. Linear model
predictive control for lane keeping and obstacle avoidance on low curvature roads. In 16th international
IEEE conference on intelligent transportation systems (ITSC 2013) (pp. 378-383). IEEE.

The difference being that the speed control is not implemented. The pdf file contains a detailed report of the simulation, relevent pages are 1-22. I originally wanted to
implement the full article but the project didn't required thee speed control part so I had to abonden it after I run out of time. There are still some functions and 
scripts about that part.

### main2.m 
The main code in the project that simulates lane keeping and obstecle avoidance.

### functions_file.m
A file that keeps most of the relevent functions for thee project such generating road and such.

### Fy_linearize.m
A script to calculate empirically some linearized parameters for the vehicle dynamics.

### semi_lin_VD.m
A function that models the semi-linear dynamics of the car.

### test_vehicle_dynamics.m
A simple test of the vehicle dynamics on a circular road

### Speed_control.m
A PI controller for the speed. Works well on it's own however has problem while integrating it with the MPC.

### main.m
The old original main code that would have implemented the speed control togetheer with the lane keeping. The controlleer unfortunately is not very stable tho.
