/************************************************************************

	Copyright 2012-2013 Luciano Buglioni

	Contact: luciano.buglioni@gmail.com

	This file is a part of FluxSol

	FluxSol is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    FluxSol is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    For a copy of the GNU General Public License,
    see <http://www.gnu.org/licenses/>.

*************************************************************************/
#ifndef _SOLVER_H
#define _SOLVER_H

#include "./Type/Vec3d.h"
#include "./Field/Field.h"

#include "SistEcuac.h"

//Sparse libraries
//#include "laspack.h"
//#include "./Nastran/Varios.h"
#include "Utils.h"

#include <time.h>


namespace FluxSol{

template <typename number>
class Solver{


    protected:
    Vec3D vtol;     //Tolerancia de velocidades
    double ptol;    //Tolerancia de presiones
    double maxiter; //maxximo numero de iteraciones
    double actualiter;  //iteracion actual


	number rtol,abstol;
	const int matdim;		//CONST?
	int maxits;


    public:

    //ES VIRTUAL!!
    //Ojo que si no la defino, si no coloco nada y solo la declaro tengo problemas con la vtable en compilacion
    virtual void Resolver_Iteracion(){};

	Solver<number>():matdim(0)
	{
	    maxiter=100;
        ptol=1.0e-3;
        vtol=1.e-03;    //Todos los valores iguales
	}
	Solver<number>(const int &d):
	matdim(d)
	{}


};

	template <typename T>
    void Solve(EqnSystem <T> &);

}

#endif
