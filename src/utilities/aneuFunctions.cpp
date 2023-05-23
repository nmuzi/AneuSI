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

#include "aneuFunctions.h"

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
#include <vtkCellArray.h>

void HeapSort( int N, double *wa2, int *iwa2 )
{
    //	iwa2[]	- the work array whose elements are to be sorted (index's array)
    // 	wa2[]	- the associated work array to the iwa2 array
    // 	N		- index of the last element to be sorted 
    
	int n= N, i = n/2, parent, child, aux;
	double t;

	for (;;)					// Loops until wa is sorted 
	{ 
		if (i > 0)				// First stage - Sorting the heap 
		{
			i--;				// Save its index to i 
			t = wa2[i];			// Save parent value to t 
			aux=iwa2[i];
		} 
		else					// Second stage - Extracting elements in-place 
		{
			n--;				// Make the new heap bioaller 
			if (n == 0)			// When the heap is empty, we are done 
				return; 
			t = wa2[n];			// Save last value (it will be overwritten) 
			aux=iwa2[n];
			wa2[n] = wa2[0];	// Save largest value at the end of wa 
			iwa2[n] = iwa2[0];
		}

		parent = i;				// We will start pushing down t from parent 
		child = i*2 + 1;		// parent's left child 

		// Shift operation - pushing the value of t down the heap 
		while (child < n)
		{
			if (child + 1 < n  &&  wa2[child + 1] > wa2[child]) 
			{
				child++;					// Choose the largest child 
			}
			if (wa2[child] > t)				// If any child is bigger than the parent 
			{ 
				wa2[parent] = wa2[child];	// Move the largest child up 
				iwa2[parent] = iwa2[child];
				parent = child;				// Move parent pointer to this child 
				child = parent*2 + 1;		// Find the next child 
			} 
			else 
			{
				break;						// t's place is found 
			}
		}
		wa2[parent] = t;					// We save t in the heap 
		iwa2[parent] = aux;
	}
}

void Fill ( int N, int val, int *wa )
        {
            int i;
            for ( i=0; i < N; i++ )
                wa[i] = val;
            return;
        }
        
void Copy ( int N, double *wa1, double *wa2 )
        {
            int i;
            for ( i=0; i < N; i++ )
                wa2[i] = wa1[i];
            return;
        }

int CheckFilename( const char *input_name, const char *input_extension )
{
    const char * ch_ptr;

    ch_ptr = std::strrchr (input_name, '.');
    if (ch_ptr == NULL)
    {
        printf( "ERROR::No extension in input filename : %s\n", input_name );
        return 1;
    }

    if ( input_extension != NULL )
    {
        if ( std::strncmp(input_extension, ch_ptr, ch_ptr - input_extension) != 0 )
        {
            printf( "ERROR::The extension '%s' of the input filename is wrong\n", ch_ptr );
            printf( "     ::It should be '%s'\n", input_extension );            
            return 2;
        }
    }
    return 0;
}

void AneuCleanCurve ( vtkSmartPointer<vtkPolyData> curvePdt, double tol,
                      vtkSmartPointer<vtkPoints> curvePts, vtkSmartPointer<vtkIdList> curveIds)
{
     /*  This script merges duplicated nodes or nodes closer to the desired tolerance. Also it
        sorts the curve points. */
    
    int numPts;
    
    // Pre process needed to sort
    curvePdt->BuildLinks(); 
    curvePdt->BuildCells();    
    numPts = curvePdt->GetPoints()->GetNumberOfPoints();        
    
    // Numbers
    double dist;
    int    ptId, cellId;
    int    ii;
    int    id1, id2;   // n stands for "new"    
    
    // Arrays
    double points[numPts*3];    // array where all points will be stored
    double x1[3], x2[3], xn[3]; // point coordinates
    double pt[3];
    int    ids[numPts+1];    
        
    // VTK
    vtkSmartPointer<vtkIdList> curveCells = vtkSmartPointer<vtkIdList>::New();  // Curve elements    
    vtkSmartPointer<vtkIdList> auxIds     = vtkSmartPointer<vtkIdList>::New();    
    
    auxIds->InsertNextId(0);
    curveCells->InsertNextId(0);
    
    //now the edges and points are sorted
    ptId=0, cellId=0;
    for( int id=0; id < numPts; id++ )
    {           
        vtkSmartPointer<vtkIdList> ptIds   = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();

        curvePdt->GetPointCells(auxIds->GetId(id),cellIds);
        if(cellIds->GetId(0) == curveCells->GetId(id))
            cellId = cellIds->GetId(1);
        else
            cellId = cellIds->GetId(0);

        //By each cell are stored the edge cells in a sorted way
        curveCells->InsertNextId(cellId);

        curvePdt->GetCellPoints(curveCells->GetId(id+1),ptIds);
        if(ptIds->GetId(0) == auxIds->GetId(id))
            ptId = ptIds->GetId(1);
        else
            ptId = ptIds->GetId(0);

        // Fill IDs Lists
        auxIds->InsertNextId(ptId);                    
    }
    
    // Put 0 at the beggining. This way we avoid repeating the 0. 
    for (int i=0; i<numPts+1; i++)
    {   
        ids[i] = auxIds->GetId(i);        
    }
    
    // Save point coordinates in our own list
    for (int i=0; i<numPts; i++)
    {
        curvePdt->GetPoints()->GetPoint(i,pt);
        ii = i*3;
        
        points[ii]   = pt[0];
        points[ii+1] = pt[1];
        points[ii+2] = pt[2];
    }
    
    int id, del;
    int iend = 0;
    
    while (iend == 0)
    {
        del = 0;
        iend = 1;
        
        for (int i=0; i<numPts; i++)
        {
            if (ids[i] != -1)
            {
                id1 = ids[i];
                id2 = ids[i+1];
                
                x1[0] = points[id1*3];
                x1[1] = points[id1*3+1];
                x1[2] = points[id1*3+2];
                
                x2[0] = points[id2*3];
                x2[1] = points[id2*3+1];
                x2[2] = points[id2*3+2];            
                
                dist = sqrt(vtkMath::Distance2BetweenPoints(x1,x2));            
                
                if (dist <= tol)
                {
                    vtkMath::Add(x1,x2,xn);
                    vtkMath::MultiplyScalar(xn,0.5);
                    
                    points[id1*3]   = xn[0];
                    points[id1*3+1] = xn[1];
                    points[id1*3+2] = xn[2];
                    
                    id = ids[i+1];  // Id to replace
                    
                    // Delete one id
                    for (int j=i+1; j<numPts; j++)
                    {
                        ids[j] = ids[j+1];
                    }            
                    // Reduce all ids accordingly
                    for (int j=0; j<numPts+1; j++)
                    {
                        if (ids[j]>id)
                        {
                            ids[j] = ids[j]-1;                    
                        }                
                    }
                    // Move all points on the list to consider the change in the Ids
                    for (int j=id; j<numPts; j++)
                    {   
                        points[j*3]   = points[(j+1)*3];
                        points[j*3+1] = points[(j+1)*3+1];
                        points[j*3+2] = points[(j+1)*3+2];                
                    }
                    del++;
                    
                    ids[numPts+1-del] = -1;
                    iend = 0;
                }
            }        
        }            
    }    
    
    int nPts=-1; // new number of points, take repeated 0 into account
    for (int i=0; i<numPts+1; i++)
    {
        if (ids[i]>-1)
        {
            curveIds->InsertNextId(ids[i]);
            nPts++;
        }
    }
    for (int j=0; j<nPts; j++)
    {
        pt[0] = points[3*j];
        pt[1] = points[3*j+1];
        pt[2] = points[3*j+2];
        
        curvePts->InsertNextPoint(pt);
    }
}

int GetOutputName(  const char *input_name,
                    const char *input_extension,
                    const char *output_extension,
                    char *output_name )
{
    //char  buffer[MAX_LINE_LENGTH];
    
    const char * ch_ptr;

    ch_ptr = std::strrchr (input_name, '.');

    std::strncpy (output_name, input_name, ch_ptr - input_name);

    output_name[ch_ptr - input_name] = '\0';

    if (input_extension == NULL )
    {
        printf( "ERROR::The extension for the input filename is NULL\n" );
        return 1;
    }
    
    if (output_extension == NULL )
    {
        printf( "ERROR::The extension for the output filename is NULL\n" );
        return 2;
    }

    std::strcat(output_name, output_extension);
    return 0;
}

