/*=========================================================================
                           ANEUSURFISOLATION
  =========================================================================
Module    : AneuSurfIsolation 
File      : BioFunctions.cpp
Copyright : (C)opyright 2021++
Authors   : D. Millan, A. Rosolen
Modified  : N. Muzi
Purpose   : Functions from BioMechanics library, for standalone 
            implementation of bioAneuIsolation application.
            This functions were used with permission from their authors. 
Date      :
Version   :
Changes   :

    This software is distributed WITHOUT ANY WARRANTY; without even 
    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
    PURPOSE.  See the above copyright notices for more information.
=========================================================================*/

#ifndef __aneuFunctions_h
#define __aneuFunctions_h

// Standard C++ header files
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <ios>
#include <fstream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <cassert>

// VTK headers
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkIdList.h>

/** \brief ABS: computes the absolute value of a number*/
#define ABS(x) ( ((x) < 0) ? -(x) : (x) )

/** \brief MAX: computes the maximum value between two numbers*/
#define MAX(x,y) ( (x) > (y) ? (x) : (y) )

/** \brief MIN: computes the minimum value between two numbers */
#define MIN(x,y) ( (x) < (y) ? (x) : (y) )

/** @name Auxiliary functions non related with the Isolation itself.
 *        These functions were taken from BioMechanics library, and 
 *        are used to perform very specific tasks.                 */

//@{
/** HeapSort algorithm to order two arrays in place up to element N
*  \param iwa the index work array whose elements are to be sorted (index's array)
*  \param wa  the associated work array
*  \param N   index of the last element to be sorted
*/
void HeapSort( int N, double *wa, int *iwa );

/** Fills an array wa of 1xN, with a given value
*  \param N    Number of elements in array
*  \param val  value used to fill the array
*  \param iwa  array to be filled with N values.
*/
void Fill ( int N, int val, int *wa );

/** Copy the elements of an array wa1 in wa2, bot of 1xN
*  \param N    Number of elements in array
*  \param wa1  origin array
*  \param wa2  destiny array
*/
void Copy ( int N, double *wa1, double *wa2 );

/** Merges duplicates or points that are closer than a given tolerance.
*   Requires VTK 
* 
*  \param [in]  curvePdt    polydata containing the curve, can be the output of a filter
*  \param [in]  tol         tolerance (desired minimum distance between nodes)
*  \param [out] curvePts    curve points
*  \param [out] curveIds    sorted Ids for the curve points
*/
void AneuCleanCurve ( vtkSmartPointer<vtkPolyData> curvePdt, double tol,
                      vtkSmartPointer<vtkPoints> curvePts, vtkSmartPointer<vtkIdList> curveIds);

/** This function checks the extension of a filename and return an error value.
*   If an error appears then it prints a message describing what happened.
* 
*  \param [in] input_name		The filename to be checked
*  \param [in] input_extension	The expected extension for this file
*  \return ierr		0: all is OK \n
* 			        1: No extension in the input filename \n
*        			2: The extension '%s' of the input filename is wrong
*/ 
int CheckFilename(const char *input_name, const char *input_extension);


/** Gets the input name in changes the extension to 'extension'.
* 
* \warning the ouput variable is set locally, then its value could be modified
* by other program, to avoid this copy the result in a global variable.
* 
*
* @param input_name ...
* @param input_extension ...
* @param output_extension ...
* @param output_name ...
* @return int flag indicating the error type (0 success)
**/
int GetOutputName( const char *input_name,
                   const char *input_extension,
                   const char *output_extension,
                   char *output_name );

//@}
#endif
