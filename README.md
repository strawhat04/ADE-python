# ADE-Python

## Introduction
**ADE-Python** is the finite volume method based Advection-Diffusion Equation Solver. Which can be very useful to simulate various transport phenomenon like flow of concentration,
temperature, enerygy or momemtum in a media.

<p align="center">
 <img src="https://latex.codecogs.com/png.latex?\dpi{80}&space;\fn_phv&space;\LARGE&space;\frac{\partial}{\partial&space;t}(\rho&space;\phi)&plus;\frac{\partial}{\partial&space;x_{j}}\left(\rho&space;u_{j}&space;\phi\right)=\frac{\partial}{\partial&space;x_{j}}\left(\Gamma&space;\frac{\partial&space;\phi}{\partial&space;x_{j}}\right)&plus;S" title="\LARGE \frac{\partial}{\partial t}(\rho \phi)+\frac{\partial}{\partial x_{j}}\left(\rho u_{j} \phi\right)=\frac{\partial}{\partial x_{j}}\left(\Gamma \frac{\partial \phi}{\partial x_{j}}\right)+S" />
</p>

The Python scipts are writted in very basic programming language, any newcomer can easily understand our code and contribute further.


We've implemented first order time implicit scheme because of its unconditionally stability, so you can use any length of time scale. 

To evaluate the value of fluxes at the interefernce of Control Volume, it uses power-scheme used by Suhas Patankar(1980). But you can use other below mentioned by chosing the 
A(P) fuction in the discretisation script

<a href="https://imgbb.com/"><img src="https://i.ibb.co/877zbtJ/Scheme.jpg" alt="Scheme" border="0"></a><br />


This project have been processed into two part i.e 2D unstructured grid and 3D structured grid. 

### 3D-Structured Grid
This part can only handly cuboidal geometery with cuboidal mesh elements. Simple numpy 3D array is used to create mesh of the geometry. We've used Line-by-Line solver which is an iterative solver, which is much slower than direct solvers but always give you realistic solution for any input conditions.

### 2D-Unstructured Grid
You can create any complex 2D geometry using the [Pygmsh](https://pypi.org/project/pygmsh/). Triangular Mesh elements are used and seperate class module is made to handle Mesh Topography information which increases complexity.
Direct Solver [numpy.linalg.solve](https://numpy.org/doc/stable/reference/generated/numpy.linalg.solve.html) is used to solve the linear equations for this part.

## Limitations
This solver gives significant False Diffusion for Peclet No greater than 10 due to upwind scheme of advection that'd been implemented. Although you can decrease this false diffusion
by using very fine meshing and compromising computational time

![gih](https://media.giphy.com/media/UhH90i5qBMcz0vdHRw/giphy.gif)

You can see the False Diffusion in advection direction. But don't worry this solver works good for simulations when your Peclet No. is less that 10.
