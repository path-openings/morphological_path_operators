morphological_path_operators
============================

ABSTRACT:

This repository contains code and documentation for morphological path operators, 
which are non-linear filters useful in image analysis and image processing

Created 2014-03-13

Path operators are filters from mathematical morphology suitable for the
filtering, enhancement, detection and segmentation of thin objects. Some 
examples are given in the repository.

1- 2D path operators
====================

This code was initially written in 2004 initially by Benjamin Appleton and Hugues Talbot, 
while at CSIRO Mathematical and Imaging Sciences. Ben was a PhD student at 
University of Queensland and Hugues was his advisor. 

It relies on the ImageMagick library for I/O.

Note: currently the master branch relies on ImageMagick-6 ; due to API changes this branch is not
compatible with ImageMagick-7. 

Work is ongoing for porting this code to ImageMagick-7 in the magick-7 branch. This is alpha code
as of this writing (August 2017).


2- 3D path operators
====================

3D-path operators are available with the RORPO code: https://github.com/path-openings/RORPO

However the RORPO code uses an approximate version of the robustness factor, and the
default code is slower by a factor of 20-30% approximately. With the approximate
robustness it is significant faster.

Please refer to the appropriate articles for an explanation on how robustness is handled
in this present code vs. the RORPO one.

3- Documentation
================

Books, theses, articles and tutorials, to be announced.

4- Misc
======

Python and Matlab bindings, tutorial code, applications. To be announced.
