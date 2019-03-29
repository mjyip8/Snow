# CS348C Final Project: Snow
Used this as my final project for CS348c Animation and Simulation.

![alt text](https://raw.githubusercontent.com/mjyip8/Snow/master/artifacts/snow.gif)

## Implementation Details
Implementing the method described in 2013 paper, [A material point method for snow simulation](http://www.math.ucla.edu/~jteran/papers/SSCTS13.pdf). The simulation runs in 2D. 
- Used explicit time integration scheme to update the velocities. 
- Implemented PIC method and transferred particle data to grid at each time step.
- Applied forces from energy function in grid.
- Transferred grid velocities back to the particles.

## Notes
I did not really implement body collisions. Only handled collisions with the wall.
