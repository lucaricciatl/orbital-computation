# Mission to Io: Aerospace Robotics Unipi


## Overview

This project involves the use of two main files, "solar_system.m" and "generate_trj.m," to model the solar system and generate orbital trajectories for a mission to Io. The mission is focused on aerospace robotics, utilizing the capabilities of the Unipi platform.

## Trajectory Calculation

### Introduction
To calculate the trajectory of the probe, we adopted the approach of "patched conics" or the method of connected conics. The idea involves assembling the final trajectory using segments of conic curves (ellipses, hyperbolas, parabolas, or circles) by overlaying the final position of one with the initial position of the next.

In selecting the trajectory, we referenced the results of the master's thesis in aerospace engineering by Engineer Andrea Caruso on "Optimization of Interplanetary Trajectories" (supervisor Professor G. Mengali). Caruso, using the genetic algorithm (GA), found a trajectory to reach Jupiter with the minimum possible ∆𝑣 in a given time window, involving two fly-bys: the first at Venus and the second at Earth.

While we do not require an optimization study, we used Caruso's results as a starting point. Using NASA's Eyes software, we recreated a planetary configuration similar to that proposed in the thesis. Therefore, the trajectory we impose on the probe to reach Io, starting from an orbit 200 km above Earth's surface, involves two fly-bys: the first at Venus and the second at Earth.

### Trajectory Calculation Assumptions
In calculating the trajectory, we assumed that if the probe is within the Sphere of Influence (SOI) of a planet, it is influenced solely by the planet's gravitational force, neglecting interaction with the Sun. This simplification significantly streamlines calculations, reducing the problem to the study of the two-body motion, in accordance with the initial assumptions.

In the following sections, we differentiate between two situations: first, the part outside the Sphere of Influence, and later, we shift our focus to describe the motion of bodies within it.

## Trajectory Analysis
1. Interplanetary Transfer (Outside the SOI)

In this initial part, called the "heliocentric phase," we consider the Sun as the sole attractor for the spacecraft, positioned at the origin of the "frame sun" reference system. Additionally, due to the length of the interplanetary journey, we consider planets with negligible size, identified only by the position of their center of mass.

The fundamental problem in this part is to calculate, for each conic segment, the velocity required for the spacecraft to be at the predetermined point at the predetermined time. Once the velocity is obtained, we can consequently calculate the trajectory to go from the initial point to the final point in the specified time. This is nothing but the "Lambert's Problem."

The problem solved by Johann Heinrich Lambert involves finding the conic transfer trajectory that connects two arbitrary points in the solar system, and it is a key aspect of our trajectory calculation.

## computation 
Following the approach described in "Keplerian Elements for Approximate Position of Major Planet," we highlight the following steps:

Calculate the value of the mean anomaly (𝑀) as 𝑀 = 𝐿 − 𝜔, where 𝐿 and 𝜔 are obtained linearly:

𝐿 = 𝐿₀ + 𝐿̇ ∙ 𝑇 ; 𝜔 = 𝜔₀ + 𝜔̇ ∙ 𝑇

Using the "Kepler equation for elliptical orbits," derive the value of the eccentric anomaly 𝐸 through a recursive algorithm:

𝑀(𝑇) = 𝐸(𝑇) − 𝑒∗ ∙ 𝑠𝑖𝑛(𝐸(𝑇))

where 𝑒∗ is the eccentricity with a conversion factor applied: 𝑒∗ = 180/𝜋 ∙ 𝑒.

The numerical iterative algorithm involves calculating ∆𝑀, ∆𝐸, and updating 𝐸 until the exit condition |∆𝐸| < 10⁻⁶𝑑𝑒𝑔 is met.
Once the value of 𝐸 is calculated, the planet's position in the local reference system is determined using the relations:

{ 𝑥 = 𝑎(𝑐𝑜𝑠𝐸 − 𝑒)
𝑦 = 𝑎√(1 − 𝑒²) ∙ 𝑠𝑖𝑛𝐸
𝑧 = 0 }

To facilitate spatial position readings, transform the coordinates to the fixed inertial reference system "frame sun" using a Euler transformation "ZXZ" with rotation angles corresponding to the orbital parameter vector values.

The transformation is implemented in the MATLAB function 𝑇𝑧𝑥𝑧.m.

Assembling elementary rotations, we obtain the following ZXZ transformation matrix:

[𝑇(Ω,𝑖,𝜔)] = [𝑐(𝜔)𝑐(Ω)−𝑠(𝜔)𝑠(Ω)𝑐(𝑖) −𝑐(Ω)𝑠(𝜔)+𝑐(𝜔)𝑠(Ω)𝑐(𝑖) −𝑠(𝜔)𝑠(𝑖)
𝑐(𝜔)𝑠(Ω)+𝑠(𝜔)𝑐(Ω)𝑐(𝑖) −𝑠(𝜔)𝑠(Ω)+𝑐(𝜔)𝑐(Ω)𝑐(𝑖) 𝑠(𝑖)𝑐(𝜔)
𝑠(𝜔)𝑠(𝑖) −𝑠(𝑖)𝑐(𝜔) 𝑐(𝑖) ]

where 𝑐(∙) = 𝑐𝑜𝑠(∙) and 𝑠(∙) = 𝑠𝑖𝑛(∙).

In the project code, steps 3 and 4 are condensed into a single function called advance.m. Given the planet and the instant at which we want to calculate its position (𝑇𝑒𝑝ℎ), it outputs a vector containing Cartesian coordinates with respect to the "sun" reference system.

To implement the position of each planet over time, we placed this function in a for loop where the value 𝑇𝑒𝑝ℎ varies from the desired initial instant to the final one. The output of the loop is a matrix with 3 rows, equal to the number of spatial coordinates of the planet, for a number of columns corresponding to the value of the desired time interval (𝑇𝑒𝑝ℎ_𝑖𝑛𝑖𝑧𝑖𝑎𝑙𝑒 − 𝑇𝑒𝑝ℎ_𝑓𝑖𝑛𝑎𝑙𝑙𝑖𝑎𝑙𝑒). The columns of this matrix are the spatial coordinates, written in the "sun" reference system, that each planet occupies on different Julian days within the chosen time interval.

If we had applied Professor E. M. Standish's algorithm in its entirety, we would have had to consider that the orbits change over time, and depending on the length of the time interval, the deviation from J2000 would have taken significant values. To avoid this, we made the simplifying assumption commonly adopted, where the Keplerian parameters remain constant over the chosen time interval. 

## Files

1. solar_system.m
This file is responsible for creating a solar system model and defining orbits to extract the position and velocity of the planets. It serves as the foundation for generating accurate trajectories for the mission.

2. generate_trj.m
This file is used to generate orbital trajectories within the solar system and inside the Sphere of Influence (SOI) of specific celestial bodies. It plays a crucial role in planning the mission route and ensuring the spacecraft follows a predefined path.

## Instructions

Running the Simulation:
Execute the solar_system.m file to create the solar system model.
Run the generate_trj.m file to generate orbital trajectories for the mission.
Parameters and Configuration:
Ensure to review and modify any relevant parameters within the files to customize the mission according to specific requirements.
Check for comments within the code to understand the purpose and functionality of each section.
Output:
The output of the simulation may include position and velocity vectors, trajectory plots, or any other relevant information depending on the implementation.
Dependencies

Ensure that you have the necessary dependencies installed, such as MATLAB or Octave, to run the simulation successfully. Additionally, any specific libraries or toolboxes required for aerospace modeling should be installed and configured.

## Additional Notes

This project assumes a basic understanding of celestial mechanics and aerospace engineering principles.
Make sure to refer to the documentation of the Unipi platform for any specific guidelines or recommendations regarding robotics and control systems.


![alt text](https://github.com/ATLED-3301/Aerospace-mission-to-Io/blob/main/solar_system.jpg?raw=true))


![alt text](https://github.com/ATLED-3301/Aerospace-mission-to-Io/blob/main/spaceship_orbits.jpg?raw=true)


![alt text](https://github.com/ATLED-3301/Aerospace-mission-to-Io/blob/main/Io.png?raw=true)
