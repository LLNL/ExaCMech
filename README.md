		ExaCMech
	
	 _______      ___    ___ ________  ________  _____ ______   _______   ________  ___  ___     
	|\  ___ \    |\  \  /  /|\   __  \|\   ____\|\   _ \  _   \|\  ___ \ |\   ____\|\  \|\  \    
	\ \   __/|   \ \  \/  / | \  \|\  \ \  \___|\ \  \\\__\ \  \ \   __/|\ \  \___|\ \  \\\  \   
	 \ \  \_|/__  \ \    / / \ \   __  \ \  \    \ \  \\|__| \  \ \  \_|/_\ \  \    \ \   __  \  
	  \ \  \_|\ \  /     \/   \ \  \ \  \ \  \____\ \  \    \ \  \ \  \_|\ \ \  \____\ \  \ \  \ 
	   \ \_______\/  /\   \    \ \__\ \__\ \_______\ \__\    \ \__\ \_______\ \_______\ \__\ \__\
	    \|_______/__/ /\ __\    \|__|\|__|\|_______|\|__|     \|__|\|_______|\|_______|\|__|\|__|
	             |__|/ \|__|                                                                     
	                                                                                             

BACKGROUND
======

ExaCMech is a GPU-friendly library of continuum mechanics (constitutive) models. Crystal-mechanics-based and porosity-mechanics-based models are a principal focus. 

Examples of work that has made use of substantially equivalent algorithms:
  * [Journal publication](http://dx.doi.org/10.1063/1.4971654) on use with a porosity mechanics type model
  * [Journal publication](http://dx.doi.org/10.1088/0965-0393/17/3/035003) on use with a crystal mechanics type model

BUILDING
======

The build system is cmake-based.

Dependencies:
* blt -- required
  - https://github.com/LLNL/blt
* snls -- required
  - https://github.com/LLNL/SNLS

TESTING
======

Run cmake with `-DENABLE_TESTS=ON` and do `make test`

DEVELOPMENT
======

The develop branch is the main development branch. Changes to develop are by pull request.

AUTHORS
======

The principal devleoper of SNLS is Nathan Barton, nrbarton@llnl.gov. 

LICENSE
======

NOTE : This software is not yet released. It is _anticipated_ that it will be released under the BSD-3-Clause license. 
