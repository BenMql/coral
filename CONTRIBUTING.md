# Contributing to Coral

Greetings, potential contributor!

We are excited to find you here, with perhaps the intent to contribute to this open-source project. We summarily describe in this file how you can make a very much appreciated contribution to the Fortran solver itself, the documentation, data-processing routines, etc.



## How to contribute

A good starting point is perhaps [Git's guide to making open-source contributions](https://opensource.guide/how-to-contribute/).

At the time of writing, Coral is still in its relative infancy, by which we mean that there is a long list of planned modifications
that range from unifying coding style and source files templates, to implementing new features. Thus before you put some work into 
any form of contribution, please create an issue that motivates (found a bug; thought of a new feature; etc.) and describes your proposed contribution.
It is my hope that in a near future contributions can be initiated in a more independent way. However, for the time being, this is how 
conflicts with planned development would be avoided.

After making sure that the proposed contribution fits in the planned roadmap, it will be included after the pull request has been reviewed.

#### Guidelines
Disclaimer: some of the sources of coral do not follow the guidelines below. The reason why is that the development started a while ago, 
when I did not reckon that the project would be shared with a community. Once in a while, I clean one of these ancient sources so that it
follows a better coding hygiene, but there remains some files with poor syntax style. Nonetheless, please follow a few guidelines for new contributions:
+ Fortran code: Use modern fortran standard (f90 to f2008) and syntax (mostly lower case; long, meaningful variable names; line breaks and indents for legibility). Avoid deprecated constructions. 
 Comment and document (using any syntax, or Doxygen if you know it). Debug and put your code through a grinder (e.g. use `gfortran -pedantic -std=f2003`, use multiple compilers, etc.).
+ Python: please use exclusively `python3`. I suggest following [these guidelines](https://www.python.org/dev/peps/pep-0008/)



## Copyrights
While you retain copyrights on your contributions, these contributions have to be licensed under the GNU Public Licence v3, similarly to the rest of Coral. 

