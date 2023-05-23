/*=========================================================================
  
Module    : Aneurysms - Biomechanics Applications
File      : AuxAneuIsolation.h
Copyright : (C)opyright 2022++
            See COPYRIGHT statement in top level directory.
Authors   : N. Muzi, D. Millan
Purpose   : Auxiliary file with the declaration of all functions used in
            BioAneuIsolation APP (see bioAneuIsolation.cpp and 
            auxAneuIsolation.cpp files). Each function is accompanied 
            by a brief description of its purpose, its inputs and outputs. 
Date      : August 2022
Version   : 1
Changes   :

    This software is distributed WITHOUT ANY WARRANTY; without even 
    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
    PURPOSE.  See the above copyright notices for more information.
=========================================================================*/

#ifndef __AuxAneuIsolation_h
#define __AuxAneuIsolation_h

// C standard lib
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//C++ STL
#include <iostream>
#include <iterator>
#include <set>

//VTK lib
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPlane.h>
#include <vtkIdList.h>
#include <vtkCleanPolyData.h>
#include <vtkCutter.h>
#include <vtkPolyDataConnectivityFilter.h>

/** @name Auxiliary functions used by the BioAneuIsolation program, which 
*         isolates an aneurysm model and the adyacent parent vessel from
*         the entire blood vessel tree.*/

//@{
/** This function computes the neck plane, using the neck polygon from the database.
 *   
 *  @param[in]   neck           neck polyData 
 *  @param[out]  neckPlane      neck vtkPlane 
 */    
void GetNeckPlane(vtkSmartPointer<vtkPolyData> neck, vtkSmartPointer<vtkPlane> neckPlane);

/** This function classifies all centerlines in one of three kinds:
 *  1 - The centerline goes under the neck, and continues through the main blood vessels
 *  2 - The centerline pass through the neck, and continues through a small vessel from it.
 *  3 - The centerline branches out before reaching the neck
 * 
 *  @param[in]   verbose        if true, print information on screen
 *  @param[in]   line2          boolean variable, indicates if there is a "type 2" centerline
 *  @param[in]   line3          boolean variable, indicates if there is a "type 3" centerline
 *  @param[in]   data_name      name of the file where clines data will be stored
 *  @param[in]   clines         centerlines PolyData
 *  @param[in]   neckPlane      neck vtkPlane
 *  @param[out]  clineTags      tags int array 
 */
void CenterLinesLocator(bool verbose, bool *line2, bool *line3, vtkSmartPointer<vtkPolyData> clines, vtkSmartPointer<vtkPlane> neckPlane, int clineTags[]);

/** This function finds all "type 3" centerlines and picks a point to clip them, as
 *  a pretreatment for the initial model mesh. 
 *   
 *  @param[in]   verbose        if true, print information on screen
 *  @param[in]   clines         centerlines PolyData
 *  @param[in]   clineTags      tags int array
 *  @param[out]  clipPoints     list of all points where the program will clip the model
 *  @param[out]  clipDirections list of directions, where the i position in this list correlates to the j position on clipPoints 
 *  @param[in]   clipFactor     distance from the neck where the clip will be performed
 */
void FindBranchClipPoints( bool verbose, vtkSmartPointer<vtkPolyData> clines, int clineTags[], 
                           vtkSmartPointer<vtkIdList> clipPoints, vtkSmartPointer<vtkIdList> clipDirections, double clipFactor);

/** This function calculates and returns the "clip radius" for a given point.
 *    
 *  @param[in]   x_p            i_point where the radius will be measured
 *  @param[in]   x_p1c          starting point from the neck  
 *  @param[in]   model          aneurysm model  
 */
double GetClipRadius(double x_p[3],double x_p1c[3], vtkSmartPointer<vtkPolyData> model);

/** This function picks, for lateral aneurysms only, the clipping points where the AneurysmClip function
 *  will perform the clip. It also assigns a direction for the clip.
 *
 *  @param[in]   verbose        if true, print information on screen
 *  @param[in]   linetype2      bool variable, its value depends if a "type 2" centerline exists
 *  @param[in]   model          aneurysm model
 *  @param[in]   clines         centerlines PolyData
 *  @param[in]   neck           neck polyData 
 *  @param[in]   nplane         neck vtkPlane
 *  @param[in]   clineTags      tags (int array)
 *  @param[in]   clipFactor     distance from the neck where the clip will be performed
 *  @param[out]  clipPoints     list of all points where the program will clip the model
 *  @param[out]  clipDirections list of directions, where the i position in this list correlates to the j position on clipPoints
 */
void FindLateralClipPoints(   bool verbose, bool linetype2, vtkSmartPointer<vtkPolyData> model, vtkSmartPointer<vtkPolyData> clines, 
                              vtkSmartPointer<vtkPolyData> neck, vtkSmartPointer<vtkPlane> nplane, int clineTags[], double clipFactor, 
                              vtkSmartPointer<vtkIdList> clipPoints, vtkSmartPointer<vtkIdList> clipDirections);

/** This function picks, for terminal aneurysms only, the clipping points where the AneurysmClip function
 *  will perform the clip. It also assigns a direction for the clip.
 *
 *  @param[in]   verbose        if true, print information on screen
 *  @param[in]   linetype2      bool variable, its value depends if a "type 2" centerline exists
 *  @param[in]   clines         centerlines PolyData
 *  @param[in]   neck           neck polyData 
 *  @param[in]   nplane         neck vtkPlane
 *  @param[in]   clineTags      tags (int array)
 *  @param[in]   clipFactor   distance from the neck where the clip will be performed
 *  @param[out]  clipPoints     list of all points where the program will clip the model
 *  @param[out]  clipDirections list of directions, where the i position in this list correlates to the j position on clipPoints
 */
void FindTerminalClipPoints(   bool verbose, bool linetype2, vtkSmartPointer<vtkPolyData> model, vtkSmartPointer<vtkPolyData> clines, 
                               vtkSmartPointer<vtkPolyData> neck, vtkSmartPointer<vtkPlane> nplane, int clineTags[], double clipFactor, 
                               vtkSmartPointer<vtkIdList> clipPoints, vtkSmartPointer<vtkIdList> clipDirections);

/** This function picks, for terminal aneurysms only, the point in the main inlet vessel where the AneurysmClip function
 *  will perform the clip. It also assigns a direction for the clip.
 *
 *  @param[in]   verbose        if true, print information on screen
 *  @param[in]   clines         centerlines PolyData
 *  @param[in]   neck           neck polyData 
 *  @param[in]   nplane         neck vtkPlane
 *  @param[in]   clineTags      tags (int array)
 *  @param[in]   clipFactor   distance from the neck where the clip will be performed
 *  @param[out]  clipPoints     list of all points where the program will clip the model
 *  @param[out]  clipDirections list of directions, where the i position in this list correlates to the j position on clipPoints
 */
void FindTerminalBifurcation( bool verbose, vtkSmartPointer<vtkPolyData> clines, vtkSmartPointer<vtkPolyData> neck,
                              vtkSmartPointer<vtkPlane> nplane, int clineTags[], double clipFactor,
                              vtkSmartPointer<vtkIdList> clipPoints, vtkSmartPointer<vtkIdList> clipDirections);

/** This function clips the model polydata using previously selected points and a direction for them.
 *
 *  @param[in]   clines         centerlines PolyData
 *  @param[in]   model          model polyData 
 *  @param[in]   clipPoints     list of all points where the program will clip the model
 *  @param[in]   clipDirections list of directions, where the i position in this list correlates to the j position on clipPoints
 *  @param[out]  clean          clean polydata after the clip
 */
void AneurysmClip(vtkSmartPointer<vtkPolyData> clines, vtkSmartPointer<vtkPolyData> model, vtkSmartPointer<vtkIdList> clipPoints,   
                  vtkSmartPointer<vtkIdList> clipDirections, vtkSmartPointer<vtkCleanPolyData> clean);

/** This function extracts the closest connected component to the neck origin, discarding the 
 *  rest of the blood vessel tree.
 * 
 *  @param[in]   clean          clean polydata after the clip
 *  @param[in]   nplane         neck vtkPlane
 *  @param[out]  outputPdt      isolated model polydata
 */
void ClipConnecExtractor (  vtkSmartPointer<vtkCleanPolyData> clean, vtkSmartPointer<vtkPlane> neckPlane, 
                            vtkSmartPointer<vtkPolyData> outputPdt);

/** This function renders and displays the isolated model when executing the program, giving
 *  a quick visual output of the result
 * 
 *  @param[in]   output         isolated model polydata
 *  @param[in]   clines         centerlines PolyData
 *  @param[in]   neck           neck polyData 
 *  @param[in]   nplane         neck vtkPlane
 */
void vtkAneuRender(vtkSmartPointer<vtkPolyData> output, vtkSmartPointer<vtkPolyData> clines, vtkSmartPointer<vtkPolyData> neck,
                   vtkSmartPointer<vtkPlane> neckPlane);

/** This function computes the ratio between the area of the circumbscribed circle 
 * around the clip "section", and the area of the inscribed circle in said section
 *
 *  @param[in]   slice          slice polydata
 *  @param[in]   center         coordinates of the center of the slice (point in the centerlines)
 *  @result      area_ratio     returned cicumscribed_area/inscribed_area value
 */
double CircumscribedInscribedClipRatio(vtkSmartPointer<vtkPolyData> slice, double center[]);

//@}
#endif
