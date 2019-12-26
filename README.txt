/*	------------------------------------------------------------------------------------

                            * DBC454 *
	unbiased parallel density based clustering of large scale ITS data


    Copyright (C) SIB  - Swiss Institute of Bioinformatics,   2013-2019 Nicolas Guex
    Copyright (C) UNIL - University of Lausanne, Switzerland       2019 Nicolas Guex


    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.


	Code:       Nicolas Guex, 2013-2019
	Contact:    Nicolas.Guex@unil.ch
	Repository: https://github.com/sib-swiss/dbc454


	Article:    Density-based hierarchical clustering of pyro-sequences
	            on a large scale—the case of fungal ITS1

                   Bioinformatics. 2013 May 15; 29(10): 1268–1274.
                   https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtt149




	Machine :	Unix
	Language:	C
	Requires:	mpi, pthread

	Version information

	Version:	1.4.3	Jan.  2013	Article published in Bioinformatics (see above).
	Version:	1.4.4	Dec.  2019	Public release of code under GPL2+ license




	Compiling:   (you will need mpi on your system)
	
	mpicc -O3 -o dbc454 dbc454.c -Wall -lpthread -lm



	Testing:

	./unit_test.sh


	------------------------------------------------------------------------------------
*/

