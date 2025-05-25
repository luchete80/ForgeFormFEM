Hi all,

after a long time I've finally succeeded in compiling a 64 bit parallel version of Petsc for Windows.
All tries with MPICH2 failed. It is known to have problems with newer 64 bit Windows (e.g. with UAC).
Therefore I gave Microsoft's MPI (which is also free of charge) a try.
Unfortunately there are some problems when using mingw compilers. (Some Fortran symbol
names are not compatible). So it does not work out of the box.
But the problems can easily be solved.

Here is the way it works:

Preparing Microsoft's MPI to work with x86_64-w64-mingw32-gfortran

0.) You need Cygwin. Please make sure that the Devel tools and Python are installed.

1.) Download Microsoft's "HPC Pack 2008 R2 MS-MPI Redistributable Package":
http://www.microsoft.com/en-us/download/details.aspx?id=14737

2.) Install it in a Directory without any spaces in the path/folder name
(LO INSTALE EN C:\Microsoft-HPCPack-2008R2\Lib\amd64)

3.) Insert the following line at the beginning of the file InstallationPath\Inc\mpi.h
#include <stdint.h>

4.) In the file InstallationPath\Inc\mpif.h
replace INT_PTR_KIND()
by 8

5.) Copy C:\Windows\System32\msmpi.dll to InstallationPath\Lib\amd64\msmpi.dll

6.) In a Cygwin terminal please enter the following:
cd /cygdrive/c/your/HPC/installation/path/Lib/amd64
gendef msmpi.dll
x86_64-w64-mingw32-dlltool -d msmpi.def -l libmsmpi.a -D msmpi.dll

This creates a library (libmsmpi.a) which now works with x86_64-w64-mingw32-gfortran.


Compiling Petsc

7.) Download Petsc v3.1-p8: petsc-3.1-p8.tar.gz<http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.1-p8.tar.gz>  (or http://www.mcs.anl.gov/petsc/download/index.html)

8.) Untar it (In a Cygwin terminal go to the folder where you've downloaded petsc-3.1-p8.tar.gz<http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.1-p8.tar.gz>
and type: tar -xzvf petsc-3.1-p8.tar.gz) and enter this directory (cd petsc-3.1-p8).

9.) Edit the file src/sys/memory/mal.c:
Change in line 39 (unsigned long) to (unsigned long long)
NO LO ENCONTRE

10.) Edit the file include/mpiuni/mpi.h:
Change in line 112
#define MPIUNI_INTPTR long
to:
#define MPIUNI_INTPTR long long

11.) Set the PATH variable in the Cygwin terminal:
export PATH=/usr/x86_64-w64-mingw32/bin:/usr/x86_64-w64-mingw32/sys-root/mingw/bin:/usr/local/bin:/usr/bin:/bin: / cygdrive/c/your/HPC/installation/path/Bin:$PATH

12.) Configure PETSc with the command given below in the Cygwin shell. Please adapt your mpi path and choose with-scalar-type=complex
or with-scalar-type=real according your needs.
Please note: GetDP can solve problems with real or complex DOFs regardless the scalar-type of Petsc. In this context the scalar-type is just an optimization in terms of speed and memory consumption.

PERO ANTES
0) cd PEtSC DIRECTory
1) export PETSC_DIR=$PWD
2) export PETSC_ARCH=NOMBRE-CONFIGURACION

Ejemplo NOMBRE-CONFIG: arch-mswin-c-mpi-hpc-release-64 
O 
–with-cc=
3-
export PATH=/usr/x86_64-w64-mingw32/bin:/usr/x86_64-w64-mingw32/sys-root/mingw/bin:/usr/local/bin:/usr/bin:/bin:/cygdrive/c/Microsoft-HPCPack-2008R2/Bin:$PATH

PC DEL TRABAJO
export PATH=/usr/x86_64-w64-mingw32/bin:/usr/x86_64-w64-mingw32/sys-root/mingw/bin:/usr/local/bin:/usr/bin:/bin:/cygdrive/d/Luciano/Programas/Microsoft-HPCPack-2008R2/Bin:$PATH

MUESTRA LO SIGUIENTE:
*******************************************************************************
         UNABLE to CONFIGURE with GIVEN OPTIONS    (see configure.log for details):
-------------------------------------------------------------------------------
Did not find package METIS needed by parmetis.
Enable the package using --with-metis or --download-metis
*******************************************************************************

PARA ESTO SE AGREGA EL –-WITH-AR

./configure --CC=x86_64-w64-mingw32-gcc.exe --CXX=x86_64-w64-mingw32-g++.exe --FC="x86_64-w64-mingw32-gfortran.exe -fno-range-check" --CPP=x86_64-w64-mingw32-cpp.exe --with-debugging=0 --with-clanguage=cxx --with-shared=0 --with-x=0 --useThreads=0 --download-f-blas-lapac=1 --download-mumps=1 --download-parmetis=1 --download-scalapack=1 --download-blacs=1 --with-scalar-type=real --with-mpi-include=/cygdrive/c/Microsoft-HPCPack-2008R2/Inc --with-mpi-lib=/cygdrive/c/Microsoft-HPCPack-2008R2/Lib/amd64/libmsmpi.a --with-ar=/usr/bin/ar


CON BLAS Y LAPACK DE FORTRAN PARA BAJAR, SE CORRIGIO –download-f-blas-lapack=1 agregandose la k. por otra parte, SE BAJA PARMETIS Y SE COLOCA EN EXTERNAL PACKAGES y se coloca --download-parmetis=0, SI SE COLOCA BAJAR PARMETIS PIDE A METIS, AUNQUE UNO LO ASIGNE EN UN DIRECTORIO NO FUNCIONA

ESTO ES CON RELEASE
./configure --CC=x86_64-w64-mingw32-gcc.exe --CXX=x86_64-w64-mingw32-g++.exe --FC="x86_64-w64-mingw32-gfortran.exe -fno-range-check" --CPP=x86_64-w64-mingw32-cpp.exe --with-debugging=0 --with-clanguage=cxx --with-shared-libraries=0 --with-x=0 --useThreads=0 --download-f-blas-lapack=0 --download-mumps=1 --download-parmetis=0 –-with-metis-include=/cygdrive/d/Luciano/Numerico/Libs/metis-5.1.0/include –-with-metis-lib=/cygdrive/d/Luciano/Numerico/Libs/metis-5.1.0/bin-x86_64-mingw/libmetis.a --download-scalapack=1 --download-blacs=1 --with-scalar-type=real --with-mpi-include=/cygdrive/d/Luciano/Programas/Microsoft-HPCPack-2008R2/Inc --with-mpi-lib=/cygdrive/d/Luciano/Programas/Microsoft-HPCPack-2008R2/Lib/amd64/libmsmpi.a --with-ar=/usr/bin/ar


===============================================================================
             Configuring PETSc to compile on your system
===============================================================================
===============================================================================                                                                                                                                                                                                      ***** WARNING: CC (set to x86_64-w64-mingw32-gcc.exe) found in environment variables - ignoring                                                                                                                                                                                 use ./configure CC=$CC if you really want to use that value ******                                                                                                                                                                                                      ===============================================================================                                                                                                                                                                                                ===============================================================================                                                                                                                                                                                                      ***** WARNING: CXX (set to x86_64-w64-mingw32-g++.exe) found in environment variables - ignoring                                                                                                                                                                                use ./configure CXX=$CXX if you really want to use that value ******                                                                                                                                                                                                    ===============================================================================                                                                                                                                                                                                ===============================================================================                                                                                                                                                                                                                WARNING! Compiling PETSc with no debugging, this should                                                                                                                                                                                                                              only be done for timing and production runs. All development should                                                                                                                                                                                                            be done when configured using --with-debugging=1                                                                                                                                                                                                         ===============================================================================                                                                                                                                                                                                ===============================================================================                                                                                                                                                                                                      Compiling FBLASLAPACK; this may take several minutes                                                                                                                                                                                                                     ===============================================================================                                                                                                                                                                                                ===============================================================================                                                                                                                                                                                                      Trying to download http://www.netlib.org/scalapack/scalapack-2.0.2.tgz for SCALAPACK                                                                                                                                                                                     ===============================================================================  

EL –AR LO AGREGUE YO, TB AGREGUE EL DOWNLOAD METIS, DEBE ESTAR ANTES QUE EL PARMETIS, YA HABIENDO BAJADO TODO

./configure --CC=x86_64-w64-mingw32-gcc.exe --CXX=x86_64-w64-mingw32-g++.exe --FC="x86_64-w64-mingw32-gfortran.exe -fno-range-check" --CPP=x86_64-w64-mingw32-cpp.exe --with-debugging=0 --with-clanguage=cxx --with-shared-libraries=0 --with-x=0 --useThreads=0 --download-f-blas-lapack=1 --download-mumps=0 --download-parmetis=0 –-with-metis-include=/cygdrive/d/Luciano/Numerico/Libs/metis-5.1.0/include –-with-metis-lib=/cygdrive/d/Luciano/Numerico/Libs/metis-5.1.0/bin-x86_64-mingw/libmetis.a --download-blacs=1 --with-scalar-type=real --with-mpi-include=/cygdrive/d/Luciano/Programas/Microsoft-HPCPack-2008R2/Inc --with-mpi-lib=/cygdrive/d/Luciano/Programas/Microsoft-HPCPack-2008R2/Lib/amd64/libmsmpi.a --with-ar=/usr/bin/ar

SI NO QUIERO BAJAR NADA (ANDA, PERO NO ME AGREGA LAS LIBRERIAS)
./configure --CC=x86_64-w64-mingw32-gcc.exe --CXX=x86_64-w64-mingw32-g++.exe --FC="x86_64-w64-mingw32-gfortran.exe -fno-range-check" --CPP=x86_64-w64-mingw32-cpp.exe --with-debugging=0 --with-clanguage=cxx --with-shared-libraries=0 --with-x=0 --useThreads=0 --download-f-blas-lapack=1 --download-mumps=0 --download-parmetis=0 –-with-metis-include=/cygdrive/d/Luciano/Numerico/Libs/metis-5.1.0/include –-with-metis-lib=/cygdrive/d/Luciano/Numerico/Libs/metis-5.1.0/bin-x86_64-mingw/libmetis.a --download-blacs=1 --with-scalar-type=real --with-mpi-include=/cygdrive/d/Luciano/Programas/Microsoft-HPCPack-2008R2/Inc --with-mpi-lib=/cygdrive/d/Luciano/Programas/Microsoft-HPCPack-2008R2/Lib/amd64/libmsmpi.a --with-ar=/usr/bin/ar

SI QUIERO USAR EN DIRECTORIOS YA BAJADOS (ESTA ES LA QUE UTILICÉ EN LA NB)

export EXT_LIBS_DIR=/cygdrive/c/EPSol-Libs


./configure --CC=x86_64-w64-mingw32-gcc.exe --CXX=x86_64-w64-mingw32-g++.exe --FC="x86_64-w64-mingw32-gfortran.exe -fno-range-check" --CPP=x86_64-w64-mingw32-cpp.exe --with-debugging=0 --with-clanguage=cxx --with-shared-libraries=0 --with-x=0 --useThreads=0 --download-f-blas-lapack=1 --download-mumps=0 --download-parmetis=0 –-with-metis-include=$EXT_LIBS_DIR/metis-5.1.0/include –-with-metis-lib=$EXT_LIBS_DIR/metis-5.1.0/bin-x86_64-mingw/libmetis.a –with-parmetis-dir=$EXT_LIBS_DIR/parmetis-4.0.3 –with-mumps-dir=$EXT_LIBS_DIR/MUMPS_4.10.0-p3 –with-scalapack-dir=$EXT_LIBS_DIR/scalapack-2.0.2 --download-blacs=1 --with-scalar-type=real --with-mpi-include=/cygdrive/d/Luciano/Programas/Microsoft-HPCPack-2008R2/Inc --with-mpi-lib=/cygdrive/d/Luciano/Programas/Microsoft-HPCPack-2008R2/Lib/amd64/libmsmpi.a --with-ar=/usr/bin/ar

MEJOR:
export EXT_LIBS_DIR=/cygdrive/c/EPSol-Libs


13.) Before Petsc is compiled the file ./cygwin-cxx-opt/include/petscconf.h  has to be edited.
Replace ( BEWARE: check the number of underscores, don't touch to PETSC_HAVE__SLEEP )

#ifndef PETSC_HAVE_SLEEP
#define PETSC_HAVE_SLEEP 1
#endif

with

#undef PETSC_HAVE_SLEEP

and replace

#ifndef PETSC_HAVE_GETPAGESIZE
#define PETSC_HAVE_GETPAGESIZE 1
#endif

with

#undef PETSC_HAVE_GETPAGESIZE


14.) Now let's compile it.
make PETSC_DIR=/your/path/to/petsc-3.1-p8 PETSC_ARCH=cygwin-cxx-opt all

15.) Drink some coffee and pray that everything works ;-)

TAMBIEN HAY QUE IR AL CODIGO DE:
petsc-3.4.4\src\sys\objects

Y AGREGAR EN options.c
//LUCIANO
#include <stdio.h>
#include <string.h>


APARTE, ANTES QUE NADA (UNA VEZ BAJADA LA LIBRERIA), se hacen estos dos cambios:
1)	TB MODIFIQUE /petscdir/src/dm/impls/plex.c agregandole el include petscsys.h
2)	Por ultimo a /petscdir/petscsys.h le DESCOMENTÉ un inline en línea 2339:
PETSC_STATIC_INLINE PetscErrorCode PetscSegBufferGetInts


------------------------ LUEGO DE COMPILAR (SI NO SE QUIERE WINDOWS) ---------------
lueggo se cambio dps de configurar el archivo petscconf.h (dentro de debug)
se saco el def de 
//MODIFIED BY LUCIANO, SET TO 0
//#ifndef PETSC_USE_WINDOWS_GRAPHICS
//#define PETSC_USE_WINDOWS_GRAPHICS 0
//#endif
Si no se hace esto no se encuentra PetscDrawCreate_Win32 no se encuentra

SE PUEDE USAR PERO SE TIENE QUE INCLUIR EL ARCHIVO DE WINDOWS

That's all.

Hopefully this is helpful for some of you.


Have a nice day!
Michael




Michael Asam
Infineon Technologies AG
ATV BP PD1 M1, Munich

Infineon Technologies AG
Vorsitzender des Aufsichtsrats: Wolfgang Mayrhuber
Vorstand: Peter Bauer (Vorsitzender), Dominik Asam, Arunjai Mittal, Dr. Reinhard Ploss
Sitz der Gesellschaft: Neubiberg
Registergericht: München HRB 126492




EL SUMARIO DE COMPILACION ES ESTE
===============================================================================                                                                                                                                                                                                TESTING: alternateConfigureLibrary from PETSc.packages.mpi4py(config/PETSc/packages/mpi4py.py:49)                                                                                                                                                                              Compilers:
  C Compiler:         x86_64-w64-mingw32-gcc.exe  -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -O
  C++ Compiler:       x86_64-w64-mingw32-g++.exe  -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -O
  Fortran Compiler:   x86_64-w64-mingw32-gfortran.exe -fno-range-check  -Wall -Wno-unused-variable -Wno-unused-dummy-argument -O
Linkers:
  Static linker:   /usr/bin/ar cr
BLAS/LAPACK: -Wl,-rpath,/petsc-3.4.4/arch-mswin-mpi1-release/lib -L/petsc-3.4.4/arch-mswin-mpi1-release/lib -lflapack -Wl,-rpath,/petsc-3.4.4/arch-mswin-mpi1-release/lib -L/petsc-3.4.4/arch-mswin-mpi1-release/lib -lfblas
MPI:
  Includes: -I/cygdrive/d/Luciano/Programas/Microsoft-HPCPack-2008R2/Inc
  Library:  -Wl,-rpath,/cygdrive/d/Luciano/Programas/Microsoft-HPCPack-2008R2/Lib/amd64 -L/cygdrive/d/Luciano/Programas/Microsoft-HPCPack-2008R2/Lib/amd64 -lmsmpi
cmake:
  Arch:
pthread:
  Library:  -lpthread
PETSc:
  PETSC_ARCH: arch-mswin-mpi1-release
  PETSC_DIR: /petsc-3.4.4
  Clanguage: Cxx
  Scalar type: real
  Precision: double
  Memory alignment: 16
  shared libraries: disabled
  dynamic loading: disabled
xxx=========================================================================xxx
 Configure stage complete. Now build PETSc libraries with (cmake build):
   make PETSC_DIR=/petsc-3.4.4 PETSC_ARCH=arch-mswin-mpi1-release all
 or (experimental with python):
   PETSC_DIR=/petsc-3.4.4 PETSC_ARCH=arch-mswin-mpi1-release ./config/builder.py
xxx=========================================================================xxx



EN EL TRABAJO TENGO ESTA SALIDA

===============================================================================
             Configuring PETSc to compile on your system
===============================================================================
===============================================================================                                                                                                                                 WARNING! Compiling PETSc with no debugging, this should                                                                                                                                               only be done for timing and production runs. All development should                                                                                                                             be done when configured using --with-debugging=1                                                                                                                          ===============================================================================                                                                                                                       Compiling FBLASLAPACK; this may take several minutes                                                                                                                                      ===============================================================================                                                                                                                 TESTING: alternateConfigureLibrary from PETSc.packages.mpi4py(config/PETSc/packages/mpi4py.py:49)                                                                                               Compilers:
  C Compiler:         x86_64-w64-mingw32-gcc.exe  -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -O
  C++ Compiler:       x86_64-w64-mingw32-g++.exe  -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -O
  Fortran Compiler:   x86_64-w64-mingw32-gfortran.exe -fno-range-check  -Wall -Wno-unused-variable -Wno-unused-dummy-argument -O
Linkers:
  Static linker:   /usr/bin/ar cr
BLAS/LAPACK: -Wl,-rpath,/cygdrive/c/cygwin64/petsc-3.4.4/arch-mswin-mpi1-release/lib -L/cygdrive/c/cygwin64/petsc-3.4.4/arch-mswin-mpi1-release/lib -lflapack -Wl,-rpath,/cygdrive/c/cygwin64/petsc-3.4.4/arch-mswin-mpi1-release/lib -L/cygdrive/c/cygwin64/petsc-3.4.4/arch-mswin-mpi1-release/lib -lfblas
MPI:
  Includes: -I/cygdrive/d/Luciano/Programas/Microsoft-HPCPack-2008R2/Inc
  Library:  -Wl,-rpath,/cygdrive/d/Luciano/Programas/Microsoft-HPCPack-2008R2/Lib/amd64 -L/cygdrive/d/Luciano/Programas/Microsoft-HPCPack-2008R2/Lib/amd64 -lmsmpi
cmake:
  Arch:
pthread:
  Library:  -lpthread
PETSc:
  PETSC_ARCH: arch-mswin-mpi1-release
  PETSC_DIR: /cygdrive/c/cygwin64/petsc-3.4.4
  Clanguage: Cxx
  Scalar type: real
  Precision: double
  Memory alignment: 16
  shared libraries: disabled
  dynamic loading: disabled
xxx=========================================================================xxx
 Configure stage complete. Now build PETSc libraries with (cmake build):
   make PETSC_DIR=/cygdrive/c/cygwin64/petsc-3.4.4 PETSC_ARCH=arch-mswin-mpi1-release all
 or (experimental with python):
   PETSC_DIR=/cygdrive/c/cygwin64/petsc-3.4.4 PETSC_ARCH=arch-mswin-mpi1-release ./config/builder.py
xxx=========================================================================xxx

