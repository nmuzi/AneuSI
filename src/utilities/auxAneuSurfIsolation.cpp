/*=========================================================================
  
Module    : Aneurysms - Biomechanics Applications
File      : AuxAneuIsolation.cpp
Copyright : (C)opyright 2022++
            See COPYRIGHT statement in top level directory.
Authors   : N. Muzi, D. Millan
Purpose   : Auxiliary file with the definition of all functions used in
            bioAneuIsolation APP (see bioAneuIsolation.cpp and 
            auxAneuIsolation.h files). 
Date      : August 2022
Version   : 1
Changes   :

    This software is distributed WITHOUT ANY WARRANTY; without even 
    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
    PURPOSE.  See the above copyright notices for more information.
=========================================================================*/

#include "auxAneuSurfIsolation.h"

#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkErrorCode.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkIdList.h>
#include <vtkCenterOfMass.h>
#include <vtkPlane.h>
#include <vtkMath.h>
#include <vtkPointLocator.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkTriangleFilter.h>

// Clipping
#include <vtkPolyDataNormals.h>
#include <vtkPolygon.h>
#include <vtkClipPolyData.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkLoopSubdivisionFilter.h> 
#include <vtkCutter.h>
#include <vtkPlaneCutter.h>
#include <vtkGenericCell.h>

// RENDERING
#include <vtkLookupTable.h>
#include <vtkPlaneSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkCamera.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>

//Standard
#include <iostream>     // std::cout
#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <fstream>      // Stream class to both read and write from/to files.
#include <set>

// BioFunctions
#include "aneuConfigFile.h"
#include "aneuFunctions.h"

#define MAX(x,y) ( (x) > (y) ? (x) : (y) )

void GetNeckPlane(vtkSmartPointer<vtkPolyData> neck, vtkSmartPointer<vtkPlane> neckPlane)
{
    // number of points describing the neck curve/plane
    int nNPts = neck->GetNumberOfPoints();
    
    double alpha;         //absolute value of the angle between vectors on the neck line whose origin is the center of the neck
    double neckPoint[3];  //3D coordinates of the center of mass of the neck curve
    double neckNormal[3]; //normal to the cut plane defining the neck
    double point[3];
    double vect0[3];
    double vect1[3];
    
    // Compute center of mass of neck points set
    vtkSmartPointer<vtkCenterOfMass> neckCenter = vtkSmartPointer<vtkCenterOfMass>::New();
    neckCenter->SetInputData(neck);
    neckCenter->SetUseScalarsAsWeights(false);
    neckCenter->Update();
        
    neckCenter->GetCenter(neckPoint); // Assign values to central point components.
    
    neck->GetPoint(0, point); 
    vect0[0] = point[0] - neckPoint[0];
    vect0[1] = point[1] - neckPoint[1];
    vect0[2] = point[2] - neckPoint[2];
    vtkMath::Normalize(vect0);

    for(int PointId = 1; PointId < nNPts;	PointId++ ) 
    {
        neck->GetPoint(PointId, point);
        vect1[0] = point[0] - neckPoint[0];
        vect1[1] = point[1] - neckPoint[1];
        vect1[2] = point[2] - neckPoint[2];
        vtkMath::Normalize(vect1);

        alpha = ABS(vtkMath::DegreesFromRadians(vtkMath::AngleBetweenVectors(vect0, vect1)));
        
        if ( alpha > 30 && alpha < 135 )
        {
            vtkMath::Cross(vect0, vect1, neckNormal);
            PointId = nNPts; // end the loop
        }
    }
    
    vtkMath::Normalize(neckNormal);    
    
    //set main ouput
    neckPlane->SetOrigin(neckPoint);
    neckPlane->SetNormal(neckNormal);
        
    return;
}

void CenterLinesLocator(bool verbose, bool *line2, bool *line3, vtkSmartPointer<vtkPolyData> clines, vtkSmartPointer<vtkPlane> neckPlane, int clineTags[])
{
    
    vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
    
    vtkDataArray *radius;
    vtkIdList*    pointIds;
    vtkIdType     idLocal, idGlobal;
    
    double pointDist;       // Distance between the selected point and neck center
    double pointRadius;     // Vessel inner radius at the selected point
    double alpha;           // Angle between the cline tangent at the selected point and the neck normal
    double closestPoint[3]; // selected point coordinates
    double nextPoint[3];    // following point coordinates
    double pointVector[3];  // centerline approximate tangent at the selected point
    double neckNormal[3];    
    
    radius = clines->GetPointData()->GetArray("MaximumInscribedSphereRadius");
    
    // Print on screen
    if (verbose == true) 
    {
        cout << "###################################" << endl;
        cout << "Classifying centerlines...\n" << endl;
    }
    
    for ( int i=0; i<clines->GetNumberOfLines(); i++ )
    {
        // "Store" clines points in polydata
        vtkSmartPointer<vtkPolyData> lineData =  vtkSmartPointer<vtkPolyData>::New();
        lineData->SetPoints(clines->GetCell(i)->GetPoints());
        pointIds = clines->GetCell(i)->GetPointIds();
        
        // Build the locator for the centerline        
        pointLocator->SetDataSet(lineData);
        pointLocator->BuildLocator();
        
        // Find the closest points of clines to the neck's center of mass
        idLocal  = pointLocator->FindClosestPoint(neckPlane->GetOrigin());  // Local ID of the closest point
        idGlobal = pointIds->GetId(idLocal);                                // Global ID of the closest point
        
        // Get the coordinates of the closest point
        lineData->GetPoint(idLocal, closestPoint);        
        
        // Calculate the distance between the selected point and the origin
        pointDist = sqrt(vtkMath::Distance2BetweenPoints(closestPoint, neckPlane->GetOrigin()));
        
        pointRadius = radius->GetTuple1(idGlobal);
            
        if ( pointDist < 2.0*pointRadius )  // Discard "type 3" centerlines from the analysis
        {
            lineData->GetPoint(idLocal,closestPoint);       // Get closest points coordinates
            lineData->GetPoint(idLocal+1,nextPoint);        // Get next point coordinates
            
            vtkMath::Subtract(closestPoint, nextPoint, pointVector);
            vtkMath::Normalize(pointVector);
                        
            neckPlane->GetNormal(neckNormal);               // get normal
            
            // Find angle between the selected point and 
            alpha = ABS(vtkMath::DegreesFromRadians(vtkMath::AngleBetweenVectors(neckNormal, pointVector)));
            
            if ( alpha > 45 && alpha < 135 )
            {
                clineTags[i] = 1;                
            }
            else
            {
                clineTags[i] = 2;
                *line2 = true;
            }            
        }
        else
        {
            clineTags[i] = 3;      
            *line3 = true;
        }
        printf("Centerline = %d  |  Tag = %d\n", i, clineTags[i] );
    
    }
    
    cout << "###################################" << endl;
    return;
}

void FindBranchClipPoints( bool verbose, vtkSmartPointer<vtkPolyData> clines,  int clineTags[],
                           vtkSmartPointer<vtkIdList> clipPoints, vtkSmartPointer<vtkIdList> clipDirections, double clipFactor)
{
    // Variables
    
    int i;            // iterator
    int iend, skip;   // finish while loop, and skip for iteration
    int step;         // change with each iteration
    
    int    refLineId;   // Id of the reference centerline    
    int    pcId, maxId;
    double x_pL1[3];    // Coordinates of point in type 1 centerline
    double x_pL3[3];    // Coordinates of point in type 3 centerline
    double x_p1[3];
    double x_p2[3];
    double dist;        
    double pRad;    

    vtkSmartPointer<vtkIdList>       refLineIds   = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList>       selLineIds   = vtkSmartPointer<vtkIdList>::New();
    
    // Initialize variables
    
    vtkDataArray *radius = clines->GetPointData()->GetArray("MaximumInscribedSphereRadius");    
    iend    = 0;
    step    = 0;
    
    // Print on screen
    if (verbose)
        printf("Clip 'Type 3' lines\n");
    
    // Select a type 1 line, used as reference    
    while ( iend == 0 )
    {
        if (clineTags[step] == 1)
        {
            iend = 1;
            refLineId = step;
        }
        
        step++;
    }
    
    refLineIds = clines->GetCell(refLineId)->GetPointIds();        
    
    // --------------------------------- 
    // Get all type 3 centerlines
    
    skip = 0;
    
    for ( i=0; i<clines->GetNumberOfLines(); i++ )
    {   
        if ( clineTags[i] == 3 )
        {
            step = 0;
            iend = 0;
            
            while (iend == 0)
            {
                clines->GetCell(refLineId)->GetPoints()->GetPoint(step,x_pL1);
                clines->GetCell(i)->GetPoints()->GetPoint(step,x_pL3);
            
                pRad = radius->GetTuple1(refLineIds->GetId(step));
            
                dist = sqrt(vtkMath::Distance2BetweenPoints(x_pL1,x_pL3));
                
                if (dist > pRad/6)
                {
                    iend = 1;                    
                }
                else
                {   
                    step++;
                }                
            }
            
            // From the bifurcation, move 1 diameter along type 3 centerline
            
            iend = 0;
            dist = 0.0;
            pRad = radius->GetTuple1(refLineIds->GetId(step));  // Use the reference line, as it is practically the same point
            
            while (iend == 0)
            {
                clines->GetCell(i)->GetPoints()->GetPoint(step,x_p1);
                clines->GetCell(i)->GetPoints()->GetPoint(step+1,x_p2);
                
                dist += sqrt(vtkMath::Distance2BetweenPoints(x_p1,x_p2));
                
                if (dist >= 2.0*pRad || step+1 == clines->GetCell(i)->GetNumberOfPoints()-1 )                    
                {
                    if (step+1 == clines->GetCell(i)->GetNumberOfPoints()-1)
                    {
                        printf("Reached the end of the 'Type 3' artery, skipping bifurcation\n");
                        iend = 1;
                        skip = 1;
                    }
                    else 
                    {
                        pcId  = clines->GetCell(i)->GetPointIds()->GetId(step);
                        maxId = clines->GetCell(i)->GetPointIds()->GetId(clines->GetCell(i)->GetNumberOfPoints()-1);
                        iend  = 1;
                    }
                }
                else 
                {
                    step++;
                }                
            }
            
            if (skip == 0)
            {
                // Move a clipFactor*radius distance from the point
                iend = 0;
                dist = 0.0;
                
                pRad = radius->GetTuple1(pcId);                
                
                while (iend == 0)
                {
                    clines->GetPoint(pcId,x_p1);
                    clines->GetPoint(pcId+1,x_p2);
                    
                    dist += sqrt(vtkMath::Distance2BetweenPoints(x_p1,x_p2));
                    
                    if (dist >= pRad*clipFactor || pcId+1 == maxId )                    
                    {
                        if (pcId+1 == maxId)
                        {
                            printf("WARNING: Reached the end of the 'Type 3' artery, skipping bifurcation\n");
                            iend = 1;                            
                        }
                        else 
                        {
                            clipPoints->InsertNextId(pcId);
                            clipDirections->InsertNextId(1);
                            iend = 1;                            
                        }
                    }
                    else 
                    {
                        pcId++;
                    }
                }
            }            
        }
    }
    
    if (verbose == true)
    {
        cout << "###################################" << endl;
    }
    
    return;
}


double GetClipRadius(double x_p[3],double x_p1c[3], vtkSmartPointer<vtkPolyData> model)
{   
    vtkSmartPointer<vtkPolyDataConnectivityFilter> connec  = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    vtkSmartPointer<vtkCleanPolyData>              clean   = vtkSmartPointer<vtkCleanPolyData>::New();
    vtkSmartPointer<vtkCutter>                     cutter  = vtkSmartPointer<vtkCutter>::New();
    vtkSmartPointer<vtkPlane>                      plane   = vtkSmartPointer<vtkPlane>::New();
    
    double  pnormal[3]; 
    double  nPt[3];
    
    // Get clip radius
    // Define slice plane
    vtkMath::Subtract(x_p,x_p1c,pnormal);  
    vtkMath::Normalize(pnormal);
    
    pnormal[0] = -pnormal[0];
    pnormal[1] = -pnormal[1];
    pnormal[2] = -pnormal[2];
    
    plane->SetOrigin(x_p1c);
    plane->SetNormal(pnormal);
    
    // We will use the slices to calculate the radius of the "cylinder"
    /// --------------------------------------------------------------------------------
    //cut with and infinite plate
    cutter->SetInputData(model);
    cutter->SetCutFunction(plane);    
    cutter->GenerateValues(1, 0, 0);
    cutter->GenerateTrianglesOff();
    cutter->Update();
            
    //clean the cut, basically remove duplicated points
    double reltol = 0.0001;
    clean->SetInputData(cutter->GetOutput());
    clean->SetTolerance(reltol);
    clean->ConvertStripsToPolysOn();
    clean->ConvertPolysToLinesOn();
    clean->ConvertLinesToPointsOn();
    clean->Update();
    
    //from the cut commonly there are many curves, choose the nearest one to the 
    //selected centerline center
    connec->SetInputData(clean->GetOutput());
    connec->SetExtractionModeToClosestPointRegion();
    connec->SetClosestPoint(x_p1c);
    connec->Update();
    
    // Get furthest point distance
    double fDist = 0.0;
    double tempDist = 0.0;
    
    for ( int j=0; j<connec->GetOutput()->GetPoints()->GetNumberOfPoints(); j++ )
    {
        connec->GetOutput()->GetPoints()->GetPoint(j,nPt);
        tempDist = sqrt(vtkMath::Distance2BetweenPoints(x_p1c,nPt));
        if (tempDist > fDist)
        {
            fDist= tempDist;
        }                
    }                
    
    return fDist;    
}

void FindLateralClipPoints(   bool verbose, bool linetype2, vtkSmartPointer<vtkPolyData> model ,vtkSmartPointer<vtkPolyData> clines, 
                              vtkSmartPointer<vtkPolyData> neck, vtkSmartPointer<vtkPlane> nplane, int clineTags[], double clipFactor, 
                              vtkSmartPointer<vtkIdList> clipPoints, vtkSmartPointer<vtkIdList> clipDirections)
{
    // Basic information from input data (needed for array definition)       
    int  nPts;                          // Number of points in neck
    nPts = neck->GetNumberOfPoints();       
    
    // -----------------------------------------------------
    // Variables 
    int     i, j;
    int     iend, step;             // iterators
    double  dist, sDist;            // Distance between points, and selected/reference distance. Used in multiple loops
    int     line1Id;                // Selected "type 1" centerline        
    double  x_n[3];                 // Line point coordinates and neck point coordinates    
    double  x_p1n[3], x_p2n[3];     // neck "ends" coordinates, and auxiliary point coordinates
    int     id1, id2;               // auxiliary Ids
    int     p1cId, p2cId;           // Id of clipping "starting points"
    double  x_p1c[3], x_p2c[3];     // Coordinates of clipping "starting points"
    double  x_f[3],   x_p[3];       // "following" and "previous" point
    double  distI;                  // Distance between F and P points
    double  p1cRad, p2cRad;         // p1c and p2c radius     
    double  diam;                   // vessel diameter
    int     p1cPts, p2cPts;         // Number of p1c and p2c points
    int     p1cTag, p2cTag;         // Tags for p1c and p2c
    double  cRad, pRad;             // "potential" clip radius
    double  area_ratio;
    int     nearId;                 // Possible near point Id (group search)
    double  nearPt[3];              // Possible near point coordinates (group search)
    int     sel1Pts, sel2Pts;       // Number of selected points (p1cSelection and p2cSelection)
    
    vtkDataArray *radius;
    
    vtkSmartPointer<vtkPolyData>     lineData    = vtkSmartPointer<vtkPolyData>::New();        // Selected type 1 cline Pdt
    vtkSmartPointer<vtkIdList>       lineIds     = vtkSmartPointer<vtkIdList>::New();          // Global Ids of the type 1 cline
    vtkSmartPointer<vtkPointLocator> lineLocator = vtkSmartPointer<vtkPointLocator>::New();    // Locator, used to find the closest point
    
    vtkSmartPointer<vtkIdList> p1cList      = vtkSmartPointer<vtkIdList>::New();     // Ids of all p1c points
    vtkSmartPointer<vtkIdList> p2cList      = vtkSmartPointer<vtkIdList>::New();     // Ids of all p2c points
    vtkSmartPointer<vtkIdList> p1clineMinId = vtkSmartPointer<vtkIdList>::New();     // Max global Id of line for all p1c points
    vtkSmartPointer<vtkIdList> p2clineMaxId = vtkSmartPointer<vtkIdList>::New();     // Max global Id of line for all p1c points
    vtkSmartPointer<vtkIdList> p1cSelMinId  = vtkSmartPointer<vtkIdList>::New();     // Max global Id on selected p1c points
    vtkSmartPointer<vtkIdList> p2cSelMaxId  = vtkSmartPointer<vtkIdList>::New();     // Max global Id on selected p1c points
    vtkSmartPointer<vtkIdList> p1cSelection = vtkSmartPointer<vtkIdList>::New();     // Ids of p1c points for clipping
    vtkSmartPointer<vtkIdList> p2cSelection = vtkSmartPointer<vtkIdList>::New();     // Ids of p2c points for clipping

    // -----------------------------------------------------
    radius = clines->GetPointData()->GetArray("MaximumInscribedSphereRadius"); // Extract radius 
    
    // Print info on screen    
    if (verbose == true)
    {
        cout << "Clip 'Type 1' lines                " << endl;        
    }   
    
    // Pick the first "Type 1" centerline. Some of the first centerlines could be "type 3"
    iend = 0;
    step = 0;
    
    while ( iend == 0 ) 
    {
        if ( clineTags[step] == 1 )
        {
            line1Id = step;     // Store the number of cline
            iend    = 1;        // End the while loop  
        }
        else step++;
    }
        
    lineIds = clines->GetCell(line1Id)->GetPointIds();    
    
    lineData->SetPoints(clines->GetCell(line1Id)->GetPoints()); // Store points in new polydata
    lineData->BuildLinks();
    
    lineLocator->SetDataSet(lineData);
    lineLocator->BuildLocator();
    
    int tempId;
    
    id1 = lineData->GetNumberOfPoints()-1;
    id2 = 0;
    
    // Get closest point to neck on the centerline¸    
    for ( i=0; i<nPts; i++ )
    {
        neck->GetPoint(i,x_n);
        
        tempId = lineLocator->FindClosestPoint(x_n);    // Get closest point id in centerline
        
        id2 = MAX(tempId,id2);
        id1 = MIN(tempId,id1);        

    }
    
    double x_p1cInit[3], x_p2cInit[3];
    
    lineData->GetPoints()->GetPoint(id1,x_p1cInit);
    lineData->GetPoints()->GetPoint(id2,x_p2cInit);    
    
    // ---------------------------------------------------
    // Find p1c and p2c for each centerline
    // Move a clipFactor and store potential clip points    
    
    double ratio = 2.5;    
    
    for ( i=0; i<clines->GetNumberOfLines(); i++ )
    {
        if (clineTags[i] == 1)
        {
            // ----------------------------------------------
            // Find closest point to each neck end
            
            lineData->SetPoints(clines->GetCell(i)->GetPoints());
            
            lineIds = clines->GetCell(i)->GetPointIds(); 
            
            lineLocator->SetDataSet(lineData);
            lineLocator->BuildLocator();
            
            p1cId = lineLocator->FindClosestPoint(x_p1cInit); // Get p1c possible id
            p2cId = lineLocator->FindClosestPoint(x_p2cInit); // Get p2c possible id
            
            // -----------------------------------------
            // Move a "clipFactor" from each point and store the point in list            
            // Initialize variables
            iend  = 0;
            step  = p1cId-1;
            distI = 0.0;
            
            p1cRad = radius->GetTuple1(lineIds->GetId(p1cId));   // Get p1c radius (to calculate diameter)
            lineData->GetPoint(p1cId,x_p1c);                     // Get p1c coordinates
            
            diam = 2.0*p1cRad;
            
            // For p1c move "backwards" along the line
            while (iend == 0)
            {   
                lineData->GetPoint(step,x_p);
                distI += sqrt(vtkMath::Distance2BetweenPoints(x_p1c, x_p));
                
                if (distI > clipFactor*diam*0.8)
                {
                    pRad = radius->GetTuple1(lineIds->GetId(step));  // point radius 
                
                    if ( step < 1 )                        
                    {
                        printf("WARNING: Reached the beggining of the artery tree. Skipping clip.\n");
                        iend = 1;                        
                    }
                    else
                    {
                        cRad = GetClipRadius(x_p, x_p1c, model); // clip radius                        
                        
                        area_ratio = (cRad*cRad)/(pRad*pRad);
                    
                        if (distI < clipFactor*diam || area_ratio > ratio) // The final distance is a multiple of the diameter
                        {
                            x_p1c[0] = x_p[0];
                            x_p1c[1] = x_p[1];
                            x_p1c[2] = x_p[2];
                            step--;
                        }
                        else if (distI >= clipFactor*diam)
                        {
                            if ( step < 1 )
                            {
                                printf("WARNING: Bifurcation reached, got to the beggining of the artery tree. Skipping clip.\n");
                                iend = 1;                
                            }
                            else
                            {
                                p1cList->InsertNextId(lineIds->GetId(step)); // Keep the global id
                                p1clineMinId->InsertNextId(lineIds->GetId(0));
                                iend = 1;
                            }
                        }               
                    }
                }
                else                    
                {            
                    x_p1c[0] = x_p[0];
                    x_p1c[1] = x_p[1];
                    x_p1c[2] = x_p[2];
                    step--;            
                }
            }
            
            // move "forward" along the line
            // P2C (point "after" the neck)
            // Initialize variables
            iend   = 0;
            step   = p2cId+1;
            distI  = 0.0;
            p2cRad = radius->GetTuple1(lineIds->GetId(p2cId));   // Get p2c radius (to calculate diameter)
            lineData->GetPoint(p2cId,x_p2c);                     // Get p2c coordinates´
            
            diam = p2cRad*2;
            
            // move "forward" along the line
            while (iend == 0)
            {
                
                lineData->GetPoint(step,x_f);
                distI += sqrt(vtkMath::Distance2BetweenPoints(x_p2c, x_f));
                
                if (distI > clipFactor*diam*0.8)
                {
                
                    pRad = radius->GetTuple1(lineIds->GetId(step));
                    
                    if (step >= lineData->GetNumberOfPoints() - 1)
                    {
                        printf("WARNING: Reached the end of the artery tree in line %d. Skipping clip.\n",i);
                        iend = 1;                           
                    }
                    else
                    {
                        cRad = GetClipRadius(x_f, x_p2c, model);
                        
                        area_ratio = (cRad*cRad)/(pRad*pRad);
                                       
                        if (distI < clipFactor*diam || area_ratio > ratio) // The final distance is a multiple of the diameter
                        {
                            x_p2c[0] = x_f[0];
                            x_p2c[1] = x_f[1];
                            x_p2c[2] = x_f[2];
                            step++;
                        }
                        else if (distI > clipFactor*diam)
                        {
                            if ( step > lineData->GetNumberOfPoints() - 1 )
                            {
                                printf("WARNING: Bifurcation reached in line %d, got to the end of the artery tree. Skipping clip.\n",i);
                                iend=1;                                
                            }
                            else
                            {
                                p2cList->InsertNextId(lineIds->GetId(step)); // Keep the global id
                                p2clineMaxId->InsertNextId(lineIds->GetId(lineData->GetNumberOfPoints()-1));
                                
                                iend=1;                
                            }
                        }
                    }
                }
                else 
                {
                    x_p2c[0] = x_f[0];
                    x_p2c[1] = x_f[1];
                    x_p2c[2] = x_f[2];
                    step++;
                }
            }
        }        
    }
    
    // -----------------------------------------------------
    // Discard all p1c and p2c points outside a small radius        
    // Initialize variables
    p1cPts = p1cList->GetNumberOfIds();
    p2cPts = p2cList->GetNumberOfIds();
    
    int p1cGroups[p1cPts];
    int p2cGroups[p2cPts];
    
    Fill(p1cPts, 0, p1cGroups);
    Fill(p2cPts, 0, p2cGroups);
    
    p1cTag  = 1;
    p2cTag  = 1;        
    // Clasify p1c in groups
    for ( i=0; i<p1cPts; i++ )
    {
        if ( p1cGroups[i] == 0 )
        {
            p1cId  = p1cList->GetId(i);
            p1cRad = radius->GetTuple1(p1cId);    // We use global Ids
            clines->GetPoint(p1cId,x_p1c);        // Get Point coordinates
         
            // Searching near points and assigning tags 
            for ( j=0; j<p1cPts; j++ )
            {
                if ( p1cGroups[j] == 0 )  // Check if the point has no tag
                {
                    nearId = p1cList->GetId(j);         // Get other "critic" point
                    clines->GetPoint(nearId,nearPt);    // Get coordinates
                    dist   = sqrt(vtkMath::Distance2BetweenPoints(x_p1c,nearPt));
                    
                    if (dist < p1cRad/3)
                    {   
                        p1cGroups[j] = p1cTag;  // Assign tag 
                    }                    
                }                
            }
            
            p1cTag++;
        }        
    }
    
    // Classify p2c in groups
    for ( i=0; i<p2cPts; i++ )
    {
        if ( p2cGroups[i] == 0 )
        {
            
            p2cId = p2cList->GetId(i);
            p2cRad = radius->GetTuple1(p2cId);    // We use global Ids
            clines->GetPoint(p2cId,x_p2c);        // Get Point coordinates
            
            // Searching near points and assigning tags 
            for ( j=0; j<p2cPts; j++ )
            {
                if ( p2cGroups[j] == 0 )
                {
                    nearId = p2cList->GetId(j);
                    clines->GetPoint(nearId,nearPt);
                    dist = sqrt(vtkMath::Distance2BetweenPoints(x_p2c,nearPt));
                    
                    if ( dist < p2cRad/3 )
                    {
                        p2cGroups[j] = p2cTag;
                    }                    
                }                
            }
        
            p2cTag++;
        }
    }
    
    // Now select only one for each group.     
    // For p1c
    // Initialize variables
    int iter;
    int minId, maxId;
    
    minId = 0;
    maxId = 0;
    
    iter = 1;
    id1  = 0;
    
    while ( iter < p1cTag )  // Iterate through all p1c points, until we reach the max tag
    {
        sDist = 0.0;
        dist  = 0.0;
        
        for ( i=0; i<p1cList->GetNumberOfIds(); i++ )
        {
            if ( p1cGroups[i] == iter )
            {
                p1cId = p1cList->GetId(i);
                clines->GetPoint(p1cId,x_p1c);
                                
                dist = vtkMath::Distance2BetweenPoints(x_p1c,x_p1n); // Get distance to neck "end"
                
                if ( dist > sDist )
                {
                    sDist = dist;               // Keep the furthest point
                    id1 = p1cList->GetId(i);
                    minId = p1clineMinId->GetId(i);
                }
            }            
        }
    
        p1cSelection->InsertNextId(id1);
        p1cSelMinId->InsertNextId(minId);
        iter++;
    }
    
    // For p2c
    // initialize variables
    iter = 1;
    id2  = 0;
    
    while ( iter < p2cTag )
    {
        sDist = 0.0;
        dist  = 0.0;
        
        for ( i=0; i<p2cList->GetNumberOfIds(); i++ ) 
        {
            if ( p2cGroups[i] == iter )
            {
                p2cId = p2cList->GetId(i);
                clines->GetPoint(p2cId,x_p2c);
                
                dist = vtkMath::Distance2BetweenPoints(x_p2c,x_p2n);
                
                if ( dist > sDist )
                {
                    sDist = dist;
                    id2 = p2cList->GetId(i);
                    maxId = p2clineMaxId->GetId(i);
                }
            }
        }
        
        p2cSelection->InsertNextId(id2);
        p2cSelMaxId->InsertNextId(maxId);
        iter++;
    }
    
    // #############################################
    double distance;    
    double x_p1[3], x_p2[3];    
    
    // Number of points on each group
    sel1Pts = p1cSelection->GetNumberOfIds();
    sel2Pts = p2cSelection->GetNumberOfIds();
    
    // Skip "sharp" ends - move backwards
    int move_backwards[sel1Pts];
    Fill(sel1Pts,-1,move_backwards);
        
    for (i=0; i<sel1Pts; i++)
    {
        if (move_backwards[i] == -1)
        {
            id1  = p1cSelection->GetId(i);
            pRad = radius->GetTuple1(id1);
            clines->GetPoint(id1, x_p1);
            
            for (j=i+1; j<sel1Pts; j++)
            {
                
                id2 = p2cSelection->GetId(j);
                clines->GetPoint(id2, x_p2);
                pRad = MAX(pRad, radius->GetTuple1(id2));
                
                distance = sqrt(vtkMath::Distance2BetweenPoints(x_p1,x_p2));                
              
                if (distance < 4.0*pRad)
                {
                    move_backwards[i] = 1;
                    move_backwards[j] = 1;
                }
            }
        }
    }    
    
    // Skip "sharp" ends - move forward
    int move_forward[sel2Pts];
    Fill(sel2Pts,-1,move_forward);
    
    for (i=0; i<sel2Pts; i++)
    {
        if (move_forward[i] == -1)
        {
            move_forward[i] = 0;
            
            id1  = p2cSelection->GetId(i);
            pRad = radius->GetTuple1(id1);
            clines->GetPoint(id1, x_p1);
            
            for (j=i+1; j<sel2Pts; j++)
            {    
                id2 = p2cSelection->GetId(j);
                clines->GetPoint(id2, x_p2);

                pRad = MAX(radius->GetTuple1(id1), radius->GetTuple1(id2));

                distance = sqrt(vtkMath::Distance2BetweenPoints(x_p1,x_p2));                
                
                if (distance < 4.0*pRad)
                {
                    move_forward[i] = 1;
                    move_forward[j] = 1;
                }
            }
        }
    }
    
    // Set Isolation points
    for ( i=0; i<sel1Pts; i++ )
    {
        id1   = p1cSelection->GetId(i);
        minId = p1cSelMinId->GetId(i);
        
        if (move_backwards[i] == 0)
        {
            clipPoints->InsertNextId(id1);
            clipDirections->InsertNextId(0);            
        }
        else
        {
            iend = 0;
            distance = 0.0;
            pRad = radius->GetTuple1(id1);
            
            while (iend == 0)
            {
                id2 = id1 - 1;
                clines->GetPoint(id1, x_p1);                
                clines->GetPoint(id2, x_p2);
                
                distance += sqrt(vtkMath::Distance2BetweenPoints(x_p1,x_p2));
                
                if ( id2 == minId || distance >= 1.25*pRad )
                {
                    clipPoints->InsertNextId(id1+1);
                    clipDirections->InsertNextId(0);
                    iend = 1;
                }
                else
                {
                   id1--; 
                }
            }   
        }
    }
    
    for ( i=0; i<sel2Pts; i++ )
    {
        id1 = p2cSelection->GetId(i);
        maxId = p2cSelMaxId->GetId(i);
        
        if (move_forward[i] == 0)
        {
            clipPoints->InsertNextId(id1);
            clipDirections->InsertNextId(1);            
        }
        else
        {
            printf("Bifurcation found, advance along the centerline\n");
            
            iend     = 0;
            distance = 0.0;
            
            pRad = radius->GetTuple1(id1);
            clines->GetPoint(id1, x_p1);                
            
            while (iend == 0)
            {
                id2 = id1 + 1;
                
                clines->GetPoint(id1, x_p1);                
                clines->GetPoint(id2, x_p2);
                distance += sqrt(vtkMath::Distance2BetweenPoints(x_p1,x_p2));
                
                if ( id2 == maxId || distance >= 1.25*pRad )
                {
                    clipPoints->InsertNextId(id1-1);
                    clipDirections->InsertNextId(1);
                    iend = 1;
                }
                else
                {
                   id1++; 
                }    
            }
        }
    }
    
    if (verbose == true)
    {        
        cout << "###################################" << endl;
    }
    
    /// --------------------------------------------------------
    // CLIP CENTERLINES THAT PASS THROUGH THE DOME
    if (linetype2 == true)
    {
        if (verbose == true)
        {
        cout << "Clip 'Type 2' lines                " << endl;        
        }       
        
        // Variables (create only if there is a "type 2" centerline
        int     i;
        int     line2Id; 
        int     sPtId;              // Starting point Id (closest to the neck origin)
        int     p1Id;
        double  x_sPt[3], x_nPt[3]; // Starting point and neck point, respectively.
        double  x_p1[3], x_p2[3];   // Points from the centerline (we use 2 for moving along the cline)
        double  dist, refDist;      
        double  p1rad;              // Inscripted sphere radius for point 1
        int     p1GlobalId;
        
        // Search all type 2 centerlines and get the closest point to the neck center. Store the clines number        
        vtkSmartPointer<vtkIdList> type2clines = vtkSmartPointer<vtkIdList>::New();
        
        for (int i=0; i<clines->GetNumberOfLines(); i++)
        {
            if (clineTags[i] == 2)
            {
                type2clines->InsertNextId(i);
            }
            else continue;
        }
        
        for (int j=0; j<type2clines->GetNumberOfIds(); j++)
        {
            
            line2Id = type2clines->GetId(j); // Get line id
            
            vtkSmartPointer<vtkPolyData> line2data = vtkSmartPointer<vtkPolyData>::New();
            line2data->SetPoints(clines->GetCell(line2Id)->GetPoints());
            
            vtkSmartPointer<vtkIdList> line2Ids = vtkSmartPointer<vtkIdList>::New(); 
            line2Ids = clines->GetCell(line2Id)->GetPointIds(); // Store global ids
            
            vtkSmartPointer<vtkPointLocator> line2Locator = vtkSmartPointer<vtkPointLocator>::New();
            line2Locator->SetDataSet(line2data);
            line2Locator->BuildLocator();
            
            refDist = 0.0;
            dist    = 0.0;
            
            // Get closest point to the neck's center
            sPtId = line2Locator->FindClosestPoint(nplane->GetOrigin()); // Starting point Id
            line2data->GetPoint(sPtId,x_sPt);
            
            for ( i=0; i<nPts; i++ )
            {
                neck->GetPoint(i,x_nPt);
                dist = sqrt(vtkMath::Distance2BetweenPoints(x_sPt,x_nPt));
                
                if (dist > refDist)
                {
                    refDist = dist;
                }                
            }
            
            p1Id  = sPtId;
            p1rad = dist;
            
            while (p1rad > 0.4*refDist)
            {
                p1Id++;
                p1GlobalId = line2Ids->GetId(p1Id);
                p1rad = radius->GetTuple1(p1GlobalId);
            }            
            
            // After choosing the point, move 2 diameter forward and store that point            
            dist = 0.0;            
            iend = 0;
            while ( iend == 0 )
            {
                if ( p1Id < line2data->GetNumberOfPoints()-2 )   
                {
                    if ( dist < 4*p1rad )
                    {
                        line2data->GetPoint(p1Id,x_p1);
                        line2data->GetPoint(p1Id+1,x_p2);
                        dist += sqrt(vtkMath::Distance2BetweenPoints(x_p1,x_p2));                
                        p1Id++;
                    }
                    else 
                    {
                        clipPoints->InsertNextId(line2Ids->GetId(p1Id));    
                        clipDirections->InsertNextId(1);
                        iend = 1;
                    }
                    
                }
                else
                {
                    printf("WARNING: Reached the end of 'type 2' artery in line %d. Skipping clip.\n",line2Id);
                    iend = 1;
                }                
            }            
        }
        
        if (verbose == true)
        {        
            cout << "###################################" << endl;
        }        
    }
    
    return;    
}

void FindTerminalClipPoints(  bool verbose, bool linetype2, vtkSmartPointer<vtkPolyData> model, vtkSmartPointer<vtkPolyData> clines, 
                              vtkSmartPointer<vtkPolyData> neck, vtkSmartPointer<vtkPlane> nplane, int clineTags[], double clipFactor,
                              vtkSmartPointer<vtkIdList> clipPoints, vtkSmartPointer<vtkIdList> clipDirections)
{
    /* In this particular case, we only need one "critical" point, as for terminal 
     * aneurysms the "starting" point of all "type 1" and "type 2" centerlines will 
     * be at the bifurcation where the aneurysm was formed.                         
     */        
    
    int    step;
    int    pcId, id;  
    int    nPts;
    int    iend;
    double area_ratio;
    double distI, pRadius;
    double pt[3], nextpt[3];
    double npt[3], cpt[3];
    double dist, diam;                // vessel diameter
    double cRad, pRad;
    
    double ratio = 2.5;
    pcId = 0;
    
    vtkSmartPointer<vtkIdList>       lineIds     = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkPolyData>     lineData    = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPointLocator> lineLocator = vtkSmartPointer<vtkPointLocator>::New();

    vtkSmartPointer<vtkIdList> pcList       = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> pclineMaxId  = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> pcSelMaxId   = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> pcSelection  = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> critPts      = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> clineIds     = vtkSmartPointer<vtkIdList>::New();
 
    nPts = neck->GetNumberOfPoints();
    
    vtkDataArray *radius = clines->GetPointData()->GetArray("MaximumInscribedSphereRadius"); // Keep the radius
    
     // Print in screen
    if (verbose == true)
    {   
       cout << "Clip type 1 centerlines after bifurcation" << endl;
    }
    
    for (int i=0; i<clines->GetNumberOfLines(); i++)
    {
        if (clineTags[i] == 1)
        {
            lineIds = clines->GetCell(i)->GetPointIds();
            
            lineData->SetPoints(clines->GetCell(i)->GetPoints());
            lineData->BuildLinks();            
            
            lineLocator->SetDataSet(lineData);
            lineLocator->BuildLocator();
            
            // Find closest neck point to cline (x_p1n)                        
            distI = 100.0;
            id    = 0;
            
            for (int j=0; j<neck->GetNumberOfPoints(); j++)
            {
                neck->GetPoint(j,npt);                      // Get neck point j coordinates
                id = lineLocator->FindClosestPoint(npt);    // Get closest point on centerline i
                lineData->GetPoint(id,cpt);                 // Get closest point coordinates
                
                dist = sqrt(vtkMath::Distance2BetweenPoints(npt,cpt)); // Calculate distance between points
                
                if (dist < distI)
                {
                    distI = dist;                    
                    pcId  = id;
                }                
            }
            
            // Move from that point the desired distance, and store point in list
            // Initialize variables
            iend   = 0;
            step   = pcId+1;
            distI  = 0.0;
            
            pRadius = radius->GetTuple1(lineIds->GetId(pcId));    // Get point radius (to calculate diameter)
            diam    = pRadius*2.0;
            
            lineData->GetPoint(step,pt);                          // Get point coordinates
            
            // move "forward" along the line
            while (iend == 0)
            {
                
                lineData->GetPoint(step,nextpt);
                distI += sqrt(vtkMath::Distance2BetweenPoints(pt, nextpt));
                
                if (distI > clipFactor*diam*0.8)
                {                
                    pRad = radius->GetTuple1(lineIds->GetId(step));
                    
                    if ( step >= lineData->GetNumberOfPoints() - 1)                        
                    {
                        printf("WARNING: Reached the end of the artery tree in line %d. Skipping clip.\n",i);
                        iend = 1;                        
                    }
                    else
                    {
                        cRad = GetClipRadius(nextpt, pt, model);
                        
                        area_ratio = (cRad*cRad)/(pRad*pRad);
                    
                        if (distI < clipFactor*diam || area_ratio > ratio) // The final distance is a multiple of the diameter
                        {
                            pt[0] = nextpt[0];
                            pt[1] = nextpt[1];
                            pt[2] = nextpt[2];
                            step++;
                        }
                        else if (distI >= clipFactor*diam)
                        {
                            if (step >= lineData->GetNumberOfPoints() - 1)
                            {
                                printf("WARNING: Bifurcation reached in line %d, got to the end of the artery tree. Skipping clip.\n",i);
                                iend=1;      
                            }
                            else
                            {
                                pcList->InsertNextId(lineIds->GetId(step)); // Keep the global id
                                pclineMaxId->InsertNextId(lineIds->GetId(lineData->GetNumberOfPoints()-1));
                                
                                iend=1;
                            }
                        }
                    }
                }  
                else 
                {
                    pt[0] = nextpt[0];
                    pt[1] = nextpt[1];
                    pt[2] = nextpt[2];
                    step++;
                }
            }            
        }   
    }
        
    // -----------------------------------------------------
    // Discard all pc outside a small radius        
    // Initialize variables
    
    int i,j;
    int id1, id2;
    int pcPts, pcTag;
    int nearId;
    int selPts;
    
    double pcRad;
    double sDist;
    double x_pc[3], x_pn[3];
    double nearPt[3];
    double distance;    
    double x_p1[3], x_p2[3];    
    
    neck->GetCenter(x_pn);
    
    pcPts = pcList->GetNumberOfIds();
    
    int pcGroups[pcPts];    
    
    Fill(pcPts, 0, pcGroups);    
    
    pcTag  = 1;
            
    // Clasify p1c in groups
    for ( i=0; i<pcPts; i++ )
    {
        if ( pcGroups[i] == 0 )
        {
            pcId  = pcList->GetId(i);
            pcRad = radius->GetTuple1(pcId);    // We use global Ids
            clines->GetPoint(pcId,x_pc);        // Get Point coordinates
         
            // Searching near points and assigning tags 
            for ( j=0; j<pcPts; j++ )
            {
                if ( pcGroups[j] == 0 )  // Check if the point has no tag
                {
                    nearId = pcList->GetId(j);         // Get other "critic" point
                    clines->GetPoint(nearId,nearPt);    // Get coordinates
                    dist   = sqrt(vtkMath::Distance2BetweenPoints(x_pc,nearPt));
                    
                    if (dist < pcRad/3)
                    {   
                        pcGroups[j] = pcTag;  // Assign tag 
                    }                    
                }                
            }
            
            pcTag++;
        }        
    }
    
    // Now select only one for each group.     
    // For p1c
    // Initialize variables

    int maxId;
        
    id1   = 0;
    step  = 1;
    maxId = 0;
    
    while ( step < pcTag )  // Iterate through all p1c points, until we reach the max tag
    {
        sDist = 0.0;
        dist  = 0.0;
        
        for ( i=0; i<pcList->GetNumberOfIds(); i++ )
        {
            if ( pcGroups[i] == step )
            {
                pcId = pcList->GetId(i);
                clines->GetPoint(pcId,x_pc);
                                
                dist = vtkMath::Distance2BetweenPoints(x_pc,x_pn); // Get distance to neck "center"
                
                if ( dist > sDist )
                {
                    sDist = dist;               // Keep the furthest point
                    id1   = pcList->GetId(i);
                    maxId = pclineMaxId->GetId(i);
                }
            }            
        }
    
        pcSelection->InsertNextId(id1);
        pcSelMaxId->InsertNextId(maxId);
        step++;
    }
    
    // Number of points on each group
    selPts = pcSelection->GetNumberOfIds();
    
    // Skip "sharp" ends - move forward
    int move_forward[selPts];
    Fill(selPts,-1,move_forward);
    
    for (i=0; i<selPts; i++)
    {
        if (move_forward[i] == -1)
        {
            move_forward[i] = 0;
            
            id1  = pcSelection->GetId(i);
            pRad = radius->GetTuple1(id1);
            clines->GetPoint(id1, x_p1);
            
            for (j=i+1; j<selPts; j++)
            {    
                id2 = pcSelection->GetId(j);
                clines->GetPoint(id2, x_p2);

                pRad = MAX(radius->GetTuple1(id1), radius->GetTuple1(id2));

                distance = sqrt(vtkMath::Distance2BetweenPoints(x_p1,x_p2));                
                
                if (distance < 4.0*pRad)
                {
                    move_forward[i] = 1;
                    move_forward[j] = 1;
                }
            }
        }
    }
    
    // Set Isolation points
    for ( i=0; i<selPts; i++ )
    {
        id1 = pcSelection->GetId(i);
        maxId = pcSelMaxId->GetId(i);
        
        if (move_forward[i] == 0)
        {
            clipPoints->InsertNextId(id1);
            clipDirections->InsertNextId(1);     
        }
        else
        {
            iend     = 0;
            distance = 0.0;
            
            pRad = radius->GetTuple1(id1);
            clines->GetPoint(id1, x_p1);                
            
            while (iend == 0)
            {
                id2 = id1 + 1;
                
                clines->GetPoint(id1, x_p1);                
                clines->GetPoint(id2, x_p2);
                distance += sqrt(vtkMath::Distance2BetweenPoints(x_p1,x_p2));
                
                if ( id2 == maxId || distance >= 1.25*pRad )
                {
                    if ( id2 == maxId )
                    {
                        printf("Reached the end of the artery. Skipping clip.\n");
                    }
                    
                    clipPoints->InsertNextId(id1);
                    clipDirections->InsertNextId(1);
                    iend = 1;
                }
                else
                {
                    id1++;                     
                }    
            }
        }
    }
    
    if (verbose == true)
    {        
        cout << "###################################" << endl;
    }    
    
    /// --------------------------------------------------------
    // CLIP CENTERLINES THAT PASS THROUGH THE DOME
    
    if (linetype2 == true)
    {
        if (verbose == true)
        {
        cout << "Clip 'Type 2' lines                " << endl;        
        }       
        
        // Variables (create only if there is a "type 2" centerline
        int     i;
        int     line2Id; 
        int     sPtId;              // Starting point Id (closest to the neck origin)
        int     p1Id;
        double  x_sPt[3], x_nPt[3]; // Starting point and neck point, respectively.
        double  x_p1[3], x_p2[3];   // Points from the centerline (we use 2 for moving along the cline)
        double  dist, refDist;      
        double  p1rad;              // Inscripted sphere radius for point 1
        int     p1GlobalId;
        
        // Search all type 2 centerlines and get the closest point to the neck center. Store the clines number        
        vtkSmartPointer<vtkIdList> type2clines = vtkSmartPointer<vtkIdList>::New();
        
        for (int i=0; i<clines->GetNumberOfLines(); i++)
        {
            if (clineTags[i] == 2)
            {
                type2clines->InsertNextId(i);
            }
            else continue;
        }
        
        for (int j=0; j<type2clines->GetNumberOfIds(); j++)
        {
            
            line2Id = type2clines->GetId(j); // Get line id
            
            vtkSmartPointer<vtkPolyData> line2data = vtkSmartPointer<vtkPolyData>::New();
            line2data->SetPoints(clines->GetCell(line2Id)->GetPoints());
            
            vtkSmartPointer<vtkIdList> line2Ids = vtkSmartPointer<vtkIdList>::New(); 
            line2Ids = clines->GetCell(line2Id)->GetPointIds(); // Store global ids
            
            vtkSmartPointer<vtkPointLocator> line2Locator = vtkSmartPointer<vtkPointLocator>::New();
            line2Locator->SetDataSet(line2data);
            line2Locator->BuildLocator();
            
            refDist = 0.0;
            dist    = 0.0;
            
            // Get closest point to the neck's center
            sPtId = line2Locator->FindClosestPoint(nplane->GetOrigin()); // Starting point Id
            line2data->GetPoint(sPtId,x_sPt);
            
            for ( i=0; i<nPts; i++ )
            {
                neck->GetPoint(i,x_nPt);
                dist = sqrt(vtkMath::Distance2BetweenPoints(x_sPt,x_nPt));
                
                if (dist > refDist)
                {
                    refDist = dist;
                }                
            }
            
            p1Id  = sPtId;
            p1rad = dist;
            
            while (p1rad > 0.4*refDist)
            {
                p1Id++;
                p1GlobalId = line2Ids->GetId(p1Id);
                p1rad = radius->GetTuple1(p1GlobalId);
            }            
            
            // After choosing the point, move 2 diameter forward and store that point            
            dist = 0.0;            
            iend = 0;
            while ( iend == 0 )
            {
                if ( p1Id < line2data->GetNumberOfPoints()-2 )   
                {
                    if ( dist < 4*p1rad )
                    {
                        line2data->GetPoint(p1Id,x_p1);
                        line2data->GetPoint(p1Id+1,x_p2);
                        dist += sqrt(vtkMath::Distance2BetweenPoints(x_p1,x_p2));                
                        p1Id++;
                    }
                    else 
                    {
                        clipPoints->InsertNextId(line2Ids->GetId(p1Id));    
                        clipDirections->InsertNextId(1);
                        iend = 1;
                    }
                    
                }
                else
                {
                    printf("WARNING: Reached the end of 'type 2' artery in line %d. Skipping clip.\n",line2Id);
                    iend = 1;
                }                
            }            
        }
        
        if (verbose == true)
        {        
            cout << "###################################" << endl;
        }
    }
    
    return;
    
}

void FindTerminalBifurcation(  bool verbose, vtkSmartPointer<vtkPolyData> clines, vtkSmartPointer<vtkPolyData> neck, 
                               vtkSmartPointer<vtkPlane> nplane, int clineTags[], double clipFactor,
                               vtkSmartPointer<vtkIdList> clipPoints, vtkSmartPointer<vtkIdList> clipDirections)
{
    /* In this particular case, we only need one "critical" point, as for terminal 
     * aneurysms the "starting" point of all "type 1" and "type 2" centerlines will 
     * be at the bifurcation where the aneurysm was formed.                         
     */        
    
    int    step;
    int    pcId, id;  
    int    iend;
    double distI, pRad;
    double npt[3], cpt[3];
    double dist, diam;                // vessel diameter    
    
    pcId = 0;
    
    vtkSmartPointer<vtkIdList>       lineIds     = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkPolyData>     lineData    = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPointLocator> lineLocator = vtkSmartPointer<vtkPointLocator>::New();

    vtkSmartPointer<vtkIdList> pcList       = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> clineIds     = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> pcSelLineId  = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> pcSelection  = vtkSmartPointer<vtkIdList>::New();
    
    vtkDataArray *radius = clines->GetPointData()->GetArray("MaximumInscribedSphereRadius"); // Keep the radius
    
     // Print in screen
    if (verbose == true)
    {   
       cout << "Clip type 1 centerlines before the bifurcation" << endl;
    }
    
    for (int i=0; i<clines->GetNumberOfLines(); i++)
    {
        if (clineTags[i] == 1)
        {
            lineIds = clines->GetCell(i)->GetPointIds();
            
            lineData->SetPoints(clines->GetCell(i)->GetPoints());
            lineData->BuildLinks();            
            
            lineLocator->SetDataSet(lineData);
            lineLocator->BuildLocator();
            
            // Find closest neck point to cline (x_p1n)                        
            distI = 100.0;
            id    = 0;
            
            for (int j=0; j<neck->GetNumberOfPoints(); j++)
            {
                neck->GetPoint(j,npt);                      // Get neck point j coordinates
                id = lineLocator->FindClosestPoint(npt);    // Get closest point on centerline i
                lineData->GetPoint(id,cpt);                 // Get closest point coordinates
                
                dist = sqrt(vtkMath::Distance2BetweenPoints(npt,cpt)); // Calculate distance between points
                
                if (dist < distI)
                {
                    distI = dist;                    
                    pcId  = id;
                }                
            }
            
            pcList->InsertNextId(lineIds->GetId(pcId)); // Keep the global id of the point
            clineIds->InsertNextId(i);                  // Keep the line number
        }
    }
    
    // Keep only one for each "branch". We want to keep only 2 centerlines
    // Initialize variables
    
    int i,j;
    int pId, lId;
    int pcPts, pcTag;
    int nearId;
    
    double pcRad;
    double sDist;
    double x_pc[3], x_pn[3];
    double nearPt[3];
    
    neck->GetCenter(x_pn);
    pcPts = pcList->GetNumberOfIds();
    
    int pcGroups[pcPts];    
    
    Fill(pcPts, 0, pcGroups);    
    
    pcTag  = 1;
            
    // Clasify p1c in groups
    for ( i=0; i<pcPts; i++ )
    {
        if ( pcGroups[i] == 0 )
        {
            pcId  = pcList->GetId(i);
            pcRad = radius->GetTuple1(pcId);    // We use global Ids
            clines->GetPoint(pcId,x_pc);        // Get Point coordinates
         
            // Searching near points and assigning tags 
            for ( j=0; j<pcPts; j++ )
            {
                if ( pcGroups[j] == 0 )  // Check if the point has no tag
                {
                    nearId = pcList->GetId(j);         // Get other "critic" point
                    clines->GetPoint(nearId,nearPt);    // Get coordinates
                    dist   = sqrt(vtkMath::Distance2BetweenPoints(x_pc,nearPt));
                    
                    if (dist < pcRad/3)
                    {   
                        pcGroups[j] = pcTag;  // Assign tag 
                    }                    
                }                
            }
            
            pcTag++;
        }        
    }
    
    step = 1;
    pId  = 0;
    lId  = 0;
    
    while ( step < pcTag )  // Iterate through all pc points, until we reach the max tag
    {
        sDist = 0.0;
        dist  = 0.0;
        
        for ( i=0; i<pcList->GetNumberOfIds(); i++ )
        {
            if ( pcGroups[i] == step )
            {
                pcId = pcList->GetId(i);
                clines->GetPoint(pcId,x_pc);
                                
                dist = vtkMath::Distance2BetweenPoints(x_pc,x_pn); // Get distance to neck "end"
                
                if ( dist > sDist )
                {
                    sDist = dist;               // Keep the furthest point
                    lId = clineIds->GetId(i);
                    pId = pcList->GetId(i);
                }
            }            
        }
    
        pcSelLineId->InsertNextId(lId);
        pcSelection->InsertNextId(pId);     // Saving just in case we have only 1 "type 1" centerline
        step++;
    }
    
    //     
    int minId;
    double x_p1[3], x_p2[3];
    
    iend = 0;    
    step = 0;
    dist = 0;    
    
    if (pcSelLineId->GetNumberOfIds() == 1)
    {
        int pId;
        
        pId = pcSelection->GetId(0);
        
        pRad = radius->GetTuple1(pId);
            
        while (iend == 0)
        {
            clines->GetPoints()->GetPoint(pId,x_p1);
            clines->GetPoints()->GetPoint(pId-1,x_p2);
            
            dist += sqrt(vtkMath::Distance2BetweenPoints(x_p1,x_p2));

            if (dist >= pRad*3)
            {
                pcId = pId;
                iend = 1;
            }
            else
            {
                pId--;
            }
        }
    }
    else
    {
        int lid1, lid2;
        int counter = 0;
        
        // Select 2 "type 1" centerlines that start at the same point
        // This is needed because not always happen (1 case in AneuriskWeb database)
        while (iend == 0)
        {
            lineIds = clines->GetCell(pcSelLineId->GetId(counter))->GetPointIds(); // Only use global ids from the first one
            lid1 = pcSelLineId->GetId(counter);   // line 1
            lid2 = pcSelLineId->GetId(counter+1);
            
            clines->GetCell(lid1)->GetPoints()->GetPoint(0,x_p1);
            clines->GetCell(lid2)->GetPoints()->GetPoint(0,x_p2);
            
            pRad = radius->GetTuple1(lineIds->GetId(0));            
            dist = sqrt(vtkMath::Distance2BetweenPoints(x_p1,x_p2));
            
            if (dist < pRad/6)
            {
                iend = 1;
            }
            else
            {
                counter++;
            }            
        }
        
        iend = 0;
        lineIds = clines->GetCell(pcSelLineId->GetId(counter))->GetPointIds(); // Use global ids
        
        while (iend == 0)
        {
            clines->GetCell(lid1)->GetPoints()->GetPoint(step,x_p1);
            clines->GetCell(lid2)->GetPoints()->GetPoint(step,x_p2);
            
            pRad = radius->GetTuple1(lineIds->GetId(step));
            
            dist = sqrt(vtkMath::Distance2BetweenPoints(x_p1,x_p2));
            
            if (dist > pRad/6)
            {
                iend = 1;
            }
            else
            {   
                step++;
            }
        }
        
        pcId  = lineIds->GetId(step); // Global id of the bifurcation point    
    }
    
    minId = lineIds->GetId(0);
    
    pRad = radius->GetTuple1(pcId);    
    diam = 2.0*pRad;

    iend = 0;
    dist = 0.0;    
    
    while (iend == 0)
    {
        clines->GetPoint(pcId, x_p1);                
        clines->GetPoint(pcId-1, x_p2);
        
        dist += sqrt(vtkMath::Distance2BetweenPoints(x_p1,x_p2));
        
        if ( pcId-1 == minId || dist >= clipFactor*diam )
        {
            if ( pcId-1 == minId )
            {
                printf("Reached the beginning of the artery. Skipping clip.\n");
            }
                    
            clipPoints->InsertNextId(pcId+1);
            clipDirections->InsertNextId(0);
            iend = 1;
        }
        else
        {
            pcId--; 
        }
    }       
    
    return;

}

void AneurysmClip(vtkSmartPointer<vtkPolyData> clines, vtkSmartPointer<vtkPolyData> model, vtkSmartPointer<vtkIdList> clipPoints,  
                   vtkSmartPointer<vtkIdList> clipDirections, vtkSmartPointer<vtkCleanPolyData> clean)
{
    int    p1Id;
    int    dir;
    int    numPts;
    double p1a[3], p1b[3];          // Cylinder base and normal points
    double p1tang[3], p2tang[3];    // Tangents   
    double reltol;                  // Tolerance for cleaning
    int    iend;    
    double area_ratio;
    
    // Make clipping 
    vtkSmartPointer<vtkCleanPolyData> recurClean = vtkSmartPointer<vtkCleanPolyData>::New();
    recurClean->SetInputData(model);
    recurClean->Update();
    
    for (int i=0; i<clipPoints->GetNumberOfIds(); i++)
    {
        iend = 0;
        dir = clipDirections->GetId(i);
                
        vtkSmartPointer<vtkPolyData>      slice  = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPlane>         top    = vtkSmartPointer<vtkPlane>::New();
        vtkSmartPointer<vtkClipPolyData>  clip   = vtkSmartPointer<vtkClipPolyData>::New();
        vtkSmartPointer<vtkCutter>        cutter = vtkSmartPointer<vtkCutter>::New();
        vtkSmartPointer<vtkCleanPolyData> clean  = vtkSmartPointer<vtkCleanPolyData>::New();
        
        vtkSmartPointer<vtkPolyDataConnectivityFilter> connec = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
        
        // Get Points    
        p1Id = clipPoints->GetId(i);    // p1a id            
        
        int iteration = 0;
        
        while (iend == 0)
        {
            clines->GetPoint(p1Id,p1a);     // p1a coordinates 
            // Find "tangent" to centerline in p1a, it will be the cutting plane normal
            if (dir == 0) 
            {        
                clines->GetPoint(p1Id-1,p1b);   // p1b coordinates
            }
            else 
            {
                clines->GetPoint(p1Id+1,p1b);   // p1b coordinates
            }
                
            vtkMath::Subtract(p1b,p1a,p1tang);  // Get p1tang
            vtkMath::Normalize(p1tang);
            
            p2tang[0] = -p1tang[0];
            p2tang[1] = -p1tang[1];
            p2tang[2] = -p1tang[2];
            
            top->SetOrigin(p1a);
            top->SetNormal(p2tang);     
            
            // We will use the slices to calculate the radius of the "cylinder"
            /// --------------------------------------------------------------------------------
            //cut with and infinite plate
            cutter->SetInputData(model);
            cutter->SetCutFunction(top);    
            cutter->GenerateValues(1, 0, 0);
            cutter->GenerateTrianglesOff();
            cutter->Update();
            
            //from the cut commonly there are many curves, choose the nearest one to the 
            //input point
            connec->SetInputData(cutter->GetOutput());
            connec->SetExtractionModeToClosestPointRegion();
            connec->SetClosestPoint(p1a);
            connec->Update();
            
            connec->GetOutput()->BuildLinks();
            connec->GetOutput()->BuildCells();
            
            //clean the cut, remove points that are "too close"
            reltol = 0.001;
            clean->SetInputData(connec->GetOutput());
            clean->SetTolerance(reltol);
            clean->ConvertStripsToPolysOff();
            clean->ConvertPolysToLinesOff();
            clean->ConvertLinesToPointsOff();
            clean->PointMergingOn();        
            clean->Update();
            
            // Measure if there is a bifurcation on the clipping point
            // If so, move a little bit further along the centerline
            area_ratio = CircumscribedInscribedClipRatio(clean->GetOutput(),p1a);
            
            if (area_ratio < 2.5)
            {
                iend = 1;
            }
            else
            {
                iteration++;
                
                if (dir == 0)
                {
                    p1Id = p1Id + 5;
                }
                else
                {
                    p1Id = p1Id - 5;
                }
            }                        
        }
        
        // Create connectivity of top, bottom and "cylinder"        
        clean->GetOutput()->BuildLinks(); // needed to sort
        clean->GetOutput()->BuildCells();
        
        vtkSmartPointer<vtkPoints> aneuboxPts = vtkSmartPointer<vtkPoints>::New();      // Points of box extracted from the model 
        vtkSmartPointer<vtkPoints> clipboxPts = vtkSmartPointer<vtkPoints>::New();      // Points of box to perform the clip        
        
        vtkSmartPointer<vtkPoints> slicePts   = vtkSmartPointer<vtkPoints>::New();      // Points extracted from slice
        vtkSmartPointer<vtkIdList> topPts     = vtkSmartPointer<vtkIdList>::New();      // Ids of points in the "top lid" (upstream)
        vtkSmartPointer<vtkIdList> botPts     = vtkSmartPointer<vtkIdList>::New();      // Ids of points in the "bottom lid" (downstream)
        vtkSmartPointer<vtkIdList> topCells   = vtkSmartPointer<vtkIdList>::New();        
        
        numPts = clean->GetOutput()->GetPoints()->GetNumberOfPoints();        
        
        topPts->InsertNextId(0);
        topCells->InsertNextId(0);
        
        //now the edges and points are sorted
        int ptId=0, cellId=0;
        for( int id=0; id < numPts; id++ )
        {           
            vtkSmartPointer<vtkIdList> ptIds   = vtkSmartPointer<vtkIdList>::New();
            vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();

            clean->GetOutput()->GetPointCells(topPts->GetId(id),cellIds);
            if(cellIds->GetId(0) == topCells->GetId(id))
                cellId = cellIds->GetId(1);
            else
                cellId = cellIds->GetId(0);

            //By each cell are stored the edge cells in a sorted way
            topCells->InsertNextId(cellId);

            clean->GetOutput()->GetCellPoints(topCells->GetId(id+1),ptIds);
            if(ptIds->GetId(0) == topPts->GetId(id))
                ptId = ptIds->GetId(1);
            else
                ptId = ptIds->GetId(0);

            // Fill IDs Lists
            topPts->InsertNextId(ptId);                                                     
        }
        
        for (int i=0; i<topPts->GetNumberOfIds(); i++)
        {
            ptId = topPts->GetId(i);
            botPts->InsertNextId(ptId + numPts);
        }
         
        
        // Create polygons
        vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
        polys->InsertNextCell(topPts);        
        polys->InsertNextCell(botPts);

        // Create connectivity for the "cylinder"
        vtkSmartPointer<vtkIdList> t1 = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> t2 = vtkSmartPointer<vtkIdList>::New();
        
        for ( int j=0; j<numPts; j++ )
        {
            // First triangle of the rectangle
            t1->InsertNextId(topPts->GetId(j));     
            t1->InsertNextId(topPts->GetId(j+1));
            t1->InsertNextId(botPts->GetId(j));
        
            // Second triangle of the rectangle
            t2->InsertNextId(topPts->GetId(j+1));
            t2->InsertNextId(botPts->GetId(j+1));
            t2->InsertNextId(botPts->GetId(j));
        
            // Create polygons
            polys->InsertNextCell(t1);
            polys->InsertNextCell(t2);
            
            t1->Reset();
            t2->Reset();
        }
        
        // Displace points radially, avoid using system generated normals
        double refPt[3], dispPt[3], sCenter[3];
        double ray[3];
        
        vtkSmartPointer<vtkCenterOfMass> sliceCenter = vtkSmartPointer<vtkCenterOfMass>::New();
        sliceCenter->SetInputData(clean->GetOutput());
        sliceCenter->SetUseScalarsAsWeights(false);
        sliceCenter->Update();
        
        sliceCenter->GetCenter(sCenter);
        
        // Iterate trough the rest of points
        for ( int i=0; i<numPts; i++ ) //All except the first element. We can use the last one as the first element repeats at the end of the list
        {

            clean->GetOutput()->GetPoints()->GetPoint(i,refPt);
            vtkMath::Subtract(refPt,sCenter,ray);
            
            dispPt[0] = refPt[0] + 0.1*ray[0];
            dispPt[1] = refPt[1] + 0.1*ray[1];
            dispPt[2] = refPt[2] + 0.1*ray[2];
            
            aneuboxPts->InsertNextPoint(dispPt);
        }
        
        // Add "the other lid" points        
        for ( int j=0; j<numPts; j++ )             
        {
            aneuboxPts->GetPoint(j,refPt); 
            
            dispPt[0] = refPt[0] + 0.5*p1tang[0];
            dispPt[1] = refPt[1] + 0.5*p1tang[1];
            dispPt[2] = refPt[2] + 0.5*p1tang[2];
            
            aneuboxPts->InsertNextPoint(dispPt);                      
        }
        
        // Use normals to displace points and make the "clipBox"
        vtkSmartPointer<vtkPolyData> clipboxPdt = vtkSmartPointer<vtkPolyData>::New();
        clipboxPdt->SetPoints(aneuboxPts);
        clipboxPdt->SetPolys(polys);
        clipboxPdt->BuildLinks();
        clipboxPdt->BuildCells();
        
        // Compute normals
        vtkSmartPointer<vtkPolyDataNormals> boxNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
        boxNormals->SetInputData(clipboxPdt);
        boxNormals->SplittingOff();
        boxNormals->ComputePointNormalsOn();
        boxNormals->ConsistencyOn();
        boxNormals->FlipNormalsOn();        
        boxNormals->NonManifoldTraversalOff();   
        boxNormals->SetOutputPointsPrecision(1);
        boxNormals->Update();        
        
        // For testing purposes only, uncomment if needed
        vtkSmartPointer<vtkPolyDataWriter> boxWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
        boxWriter->SetInputData(boxNormals->GetOutput());
        boxWriter->SetFileName("clipBox.vtk");
        boxWriter->Write();
        boxWriter->Update();  
        
        // Implicit function that will be used to slice the mesh
        vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
        implicitPolyDataDistance->SetInput(boxNormals->GetOutput());
        
        clip->SetClipFunction(implicitPolyDataDistance);
        clip->SetInputData(recurClean->GetOutput());
        clip->SetOutputPointsPrecision(1);
        clip->SetInsideOut(false);
        clip->SetValue(0.1);
        clip->GenerateClippedOutputOn();
        clip->Update();
        
        recurClean->SetInputData(clip->GetOutput());
        recurClean->SetOutputPointsPrecision(1);
        recurClean->ConvertStripsToPolysOn();
        recurClean->ConvertPolysToLinesOn();
        recurClean->ConvertLinesToPointsOn();
        recurClean->Update();   
        
    }
    
    reltol = 0.0001;
    
    clean->SetInputData(recurClean->GetOutput());
    clean->SetTolerance(reltol);
    clean->ConvertStripsToPolysOn();
    clean->ConvertPolysToLinesOn();
    clean->ConvertLinesToPointsOn();
    clean->Update();
    
    return;
}

void ClipConnecExtractor ( vtkSmartPointer<vtkCleanPolyData> clean, vtkSmartPointer<vtkPlane> neckPlane, 
                            vtkSmartPointer<vtkPolyData> outputPdt)
{
    
    
    vtkSmartPointer<vtkPolyDataConnectivityFilter> connec = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    connec->SetInputData(clean->GetOutput());
    connec->SetExtractionModeToClosestPointRegion();
    connec->SetClosestPoint(neckPlane->GetOrigin());
    connec->Update();
    
    vtkSmartPointer<vtkCleanPolyData> cleanConnec = vtkSmartPointer<vtkCleanPolyData>::New();
    cleanConnec->SetInputData(connec->GetOutput());
    cleanConnec->Update();
    
    outputPdt->SetPoints(cleanConnec->GetOutput()->GetPoints());
    outputPdt->SetPolys(cleanConnec->GetOutput()->GetPolys());
    outputPdt->BuildLinks();
    outputPdt->BuildCells();
}

void vtkAneuRender(vtkSmartPointer<vtkPolyData> output, vtkSmartPointer<vtkPolyData> clines, vtkSmartPointer<vtkPolyData> neck, 
                   vtkSmartPointer<vtkPlane> neckPlane)
{
    vtkSmartPointer<vtkPlaneSource> VTKplane = vtkSmartPointer<vtkPlaneSource>::New();
    VTKplane->SetCenter(neckPlane->GetOrigin());
    VTKplane->SetNormal(neckPlane->GetNormal());
    VTKplane->Update();
    
    vtkSmartPointer<vtkPolyDataMapper>  mMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mMapper->SetInputData(output);
    mMapper->ScalarVisibilityOff();

    vtkSmartPointer<vtkActor> mActor = vtkSmartPointer<vtkActor>::New();
    mActor->SetMapper(mMapper);
    mActor->GetProperty()->SetColor(1.0000, 0.3882, 0.2784);
    mActor->GetProperty()->SetOpacity(0.6);

    //-------------------------------------------------
    vtkSmartPointer<vtkPolyDataMapper>  cMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    cMapper->SetInputData(clines);
    cMapper->ScalarVisibilityOff();

    vtkSmartPointer<vtkActor> cActor = vtkSmartPointer<vtkActor>::New();
    cActor->SetMapper(cMapper);
    cActor->GetProperty()->SetColor(0.1, 0.1, 1.0);
    cActor->GetProperty()->SetOpacity(1.0);
    cActor->GetProperty()->SetLineWidth(4.0);

    //-------------------------------------------------
    vtkSmartPointer<vtkPolyDataMapper> pMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pMapper->SetInputData(VTKplane->GetOutput());

    vtkSmartPointer<vtkActor> pActor = vtkSmartPointer<vtkActor>::New();
    pActor->SetMapper(pMapper);
    pActor->PickableOn();
    pActor->SetOrigin(neckPlane->GetOrigin());
    pActor->SetScale(5.0);
    pActor->GetProperty()->SetColor(0.3, 1.0, 0.2);
    pActor->GetProperty()->SetOpacity(0.6);

    //-------------------------------------------------
    vtkSmartPointer<vtkPolyDataMapper> nMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    nMapper->SetInputData(neck);
    
    vtkSmartPointer<vtkActor> nActor = vtkSmartPointer<vtkActor>::New();
    nActor->SetMapper(nMapper);
    nActor->PickableOn();
    nActor->GetProperty()->SetColor(0.3, 0.3, 1.0);
    nActor->GetProperty()->SetOpacity(1.0);
    nActor->GetProperty()->SetEdgeColor(0,0,0);
    nActor->GetProperty()->SetRenderLinesAsTubes(true);
    nActor->GetProperty()->SetLineWidth(4.0);
            
    //-------------------------------------------------
    vtkSmartPointer<vtkRenderer>               ren    = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow>           renWin = vtkSmartPointer<vtkRenderWindow>::New();
    vtkSmartPointer<vtkRenderWindowInteractor> iren   = vtkSmartPointer<vtkRenderWindowInteractor>::New();

    renWin->AddRenderer(ren);
    iren->SetRenderWindow(renWin);
    
    //Add the objects to visualize
    ren->AddActor(mActor);
    ren->AddActor(cActor);
    ren->AddActor(pActor);
    ren->AddActor(nActor);

    ren->SetBackground(1,1,1);
    
    renWin->SetSize(700, 700);

    // We'll zoom in a little by accessing the camera and invoking a "Zoom"
    // method on it.
    renWin->Render(); 

    // This starts the event loop and as a side effect causes an initial render.
    iren->Start();
}

double CircumscribedInscribedClipRatio(vtkSmartPointer<vtkPolyData> slice, double center[])
{
    int    nPts;
    double min_radius, max_radius;
    double dist;
    double area_ratio;
    double x_p[3];
    
    max_radius = 0.0;
    min_radius = 100.0;
    
    nPts = slice->GetNumberOfPoints();

    
    
    for (int i=0; i<nPts; i++)
    {
        slice->GetPoints()->GetPoint(i,x_p);
        
        dist = sqrt(vtkMath::Distance2BetweenPoints(x_p,center));
        
        if (dist >= max_radius)
        {
            max_radius = dist;            
        }
        if (dist <= min_radius)
        {
            min_radius = dist;
        }
    }
        
    area_ratio = (max_radius*max_radius)/(min_radius*min_radius);
    
    return area_ratio;
}

