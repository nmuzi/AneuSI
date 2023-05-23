/*=========================================================================
                           ANEUSURFISOLATION
  =========================================================================
File      : AneuSurfIsolation.cpp
Copyright : (C)opyright 2021++            
Authors   : N. Muzi, D. Millan
Purpose   : This application isolates an aneurysm and the adyacent parent 
            vessel from a blood vessel tree. It receives a .cfg file (see
            app source) with the inputs/outputs:
            
            Inputs:
                - Polydata meshes of the aneurysm model, the ostium and 
                the centerlines. 
                - Clip distance (as a multiple of the inner radius of the
                vessel).
                - Aneurysm Type: 1 for lateral, 2 for terminal.
                - Options: verbose (display information on screen),
                           bRender (displays the isolated model after the 
                           clip). 
            Outputs: 
                - Polydata mesh of the isolated model. 
                
Date      : February 2023
Version   : 1.0
Changes   : 

    This software is distributed WITHOUT ANY WARRANTY; without even 
    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
    PURPOSE.  See the above copyright notices for more information.
=========================================================================*/

// C standard lib
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>     // std::cout
#include <fstream>      // Stream class to both read and write from/to files.
#include <iterator>
#include <set>

//VTK lib
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkErrorCode.h>
#include <vtkPlane.h>
#include <vtkIdList.h>
#include <vtkCleanPolyData.h>

// BioMec
#include "auxAneuSurfIsolation.h"    // Contains the main functions.
#include "aneuConfigFile.h"          // Class for reading .cfg configuration files
#include "aneuFunctions.h"           // Auxiliary functions from MOCCAI Group's libraries

int Usage();

int ReadInputData( const char *config_filename, bool *verbose, bool *bRender, int *aneuType, std::string *output_name, double *clipFactor,
                   vtkSmartPointer<vtkPolyData> model, vtkSmartPointer<vtkPolyData> clines, vtkSmartPointer<vtkPolyData> neck);

//------------------------------------START MAIN-----------------------------------------
int main(int argc, char *argv[])
{
    if (argc < 2)
        return Usage();
    
    // #################################################################
    //                          .CFG INPUTS 
    // #################################################################
    
    //CONFIGURATION FILE containing all the info
    const char *config_filename  = NULL;
    
    //names of inputs/output files
    std::string output_name;    // Clipped aneurysm polydata output filename
    
    bool    verbose  = false;   // Display information explaining what is being done
    bool    bRender  = false;   // Display an image of the output file
    
    int     aneuType = 1;       // Aneurysm type. 1:lateral 2:terminal
    double  clipFactor = 1;       // Set the clip distance from the neck borders
    
    vtkSmartPointer<vtkPolyData> model  = vtkSmartPointer<vtkPolyData>::New();      // Model polydata
    vtkSmartPointer<vtkPolyData> clines = vtkSmartPointer<vtkPolyData>::New();      // Centerlines polydata 
    vtkSmartPointer<vtkPolyData> neck   = vtkSmartPointer<vtkPolyData>::New();      // Automatedneck polydata
    vtkSmartPointer<vtkPolyData> output = vtkSmartPointer<vtkPolyData>::New();      // Output polydata
    
    // READ INPUTS 
    // Input parameters are read from .cfg file
    config_filename  = argv[1];
    argc--;
    argv++;
    
    if ( ReadInputData(config_filename, &verbose, &bRender, &aneuType, &output_name, &clipFactor, model, clines, neck) )
    {
        cerr << "Input reading fails. Check it." << endl;        
        return Usage();
    }
    
    printf("ClipFactor = %6.4f\n",clipFactor);
    
    // #################################################################
    //                        MESH PROCESSING 
    // #################################################################
    
    // Variable declaration (general, for every case)  
    vtkSmartPointer<vtkPlane>          nplane        = vtkSmartPointer<vtkPlane>::New();            // Neck plane    
    vtkSmartPointer<vtkIdList>         clipPoints    = vtkSmartPointer<vtkIdList>::New();           // Store clipping points global ids
    vtkSmartPointer<vtkIdList>         clipDirection = vtkSmartPointer<vtkIdList>::New();           // Direction of the clipping, 0:Backwards, 1:Forward    
    vtkSmartPointer<vtkCleanPolyData>  clean         = vtkSmartPointer<vtkCleanPolyData>::New();    // Clean polydata after clipping.
    vtkSmartPointer<vtkPolyDataWriter> aneWriter     = vtkSmartPointer<vtkPolyDataWriter>::New();   // Write polydata in disk
    
    int clineTags[clines->GetNumberOfLines()];      // "Tag" for each centerline    
    
    bool line2 = false;      // Used to check if exits a "type 2" centerline (enters into the dome) 
    bool line3 = false;      // Used to check if exits a "type 3" centerline (ends before passing behind the dome)    
    
    GetNeckPlane(neck, nplane);     // Get neck plane
    
    CenterLinesLocator(verbose, &line2, &line3, clines, nplane, clineTags);  // Classify centerlines.
    
    if (verbose)
    {
        if (aneuType == 1)
        {
            cout << "Aneurysm Type: Lateral" << endl;
        }
        else
        {
            cout << "Aneurysm Type: Terminal" << endl;        
        }
    }
    
    if (line3)
    {
        // Variable declaration (specific for "type 3" centerlines, if they exist)  
        vtkSmartPointer<vtkCleanPolyData> cleanPreModel      = vtkSmartPointer<vtkCleanPolyData>::New();  // Clean Polydata without bifurcations
        vtkSmartPointer<vtkIdList>        clipLine3Points    = vtkSmartPointer<vtkIdList>::New();         // Store bifurcations clipping points global ids
        vtkSmartPointer<vtkIdList>        clipLine3Direction = vtkSmartPointer<vtkIdList>::New();         // Direction of the clipping, 0:Backwards, 1:Forward
        
        FindBranchClipPoints(verbose, clines, clineTags, clipLine3Points, clipLine3Direction,clipFactor);
        AneurysmClip(clines, model, clipLine3Points, clipLine3Direction, cleanPreModel);
        ClipConnecExtractor(cleanPreModel, nplane, model);
        
    }
    
    if (aneuType == 1)
    {
        FindLateralClipPoints(verbose, line2, model, clines, neck, nplane, clineTags, clipFactor, clipPoints, clipDirection);
    }
    else
    {
        FindTerminalClipPoints(verbose, line2, model, clines, neck, nplane, clineTags, clipFactor, clipPoints, clipDirection);
        FindTerminalBifurcation(verbose, clines, neck, nplane, clineTags, clipFactor, clipPoints, clipDirection);
    }
    
    //---------------------------------------------------------
    AneurysmClip(clines, model, clipPoints, clipDirection, clean);
    
    // Keep only the closest connected component to the neck center    
    ClipConnecExtractor(clean, nplane, output);
    
    // Save output
    aneWriter->SetInputData(output);
    aneWriter->SetFileName(output_name.c_str());
    aneWriter->Write();
    aneWriter->Update();

    if(bRender) 
    {
        vtkAneuRender(output,clines,neck,nplane);
    }
    
    return 0;
}


// ------------------------------- FUNCTIONS ------------------------------------------

int Usage()
{
    cerr << "Usage: BioAneuIsolation config_file.cfg                                    " <<endl;
    cerr << "                                                                           " <<endl;      
    cerr << "This application isolates an aneurysm model from the vessel tree.          " <<endl;
    cerr << "The input is a cfg file with the following inputs/outputs:                 " <<endl;
    cerr << "Inputs:                                                                    " <<endl;
    cerr << "- Polydata meshes of the aneurysm model, the ostium and the centerlines.   " <<endl;                                            
    cerr << "- Clip distance (as a multiple of the inner radius of the vessel).         " <<endl;
    cerr << "- Aneurysm Type: 1 for lateral, 2 for terminal.                            " <<endl;
    cerr << "- Options: verbose (display information on screen),                        " <<endl;
    cerr << "  bRender (displays the isolated model after the clip).                    " <<endl;
    cerr << "Outputs:                                                                   " <<endl;
    cerr << "- Polydata mesh of the isolated model.                                     " <<endl;    
    return 1; 
}

int ReadInputData( const char *config_filename, bool *verbose, bool *bRender, int *aneuType, std::string *output_name, double *clipFactor,
                   vtkSmartPointer<vtkPolyData> model, vtkSmartPointer<vtkPolyData> clines, vtkSmartPointer<vtkPolyData> neck)
{
    std::string model_name;     // Artery and aneurysm mesh filename
    std::string clines_name;    // Centerlines polydata filename
    std::string neck_name;      // Aneurysm neck polydata filename
    
    if (  CheckFilename( config_filename, ".cfg" ) )
    {
        printf("ERROR::Configuration filename is not properly defined (should be .cfg)\n");
        return Usage();
    }
    aneuConfigFile config( config_filename );
    
    // Read input parameters from the .cfg file
    config.readInto( model_name, "model_name" );
    if ( model_name.empty()==true || CheckFilename( model_name.c_str(), ".vtk" ) )
    {
        printf("ERROR::model filename either is not set or doesn't have '.vtk' extension\n");
        return Usage();
    }
    
    config.readInto( clines_name, "clines_name" );
    if ( clines_name.empty()==true || CheckFilename( clines_name.c_str(), ".vtk" ) )
    {
        printf("ERROR::centerlines filename either is not set or doesn't have '.vtk' extension\n");
        return Usage();
    }
     
     config.readInto( neck_name, "neck_name" );
    if ( neck_name.empty()==true || CheckFilename( neck_name.c_str(), ".vtk" ) )
    {
        printf("ERROR::neck filename either is not set or doesn't have '.vtk' extension\n");
        return Usage();
    }
     
     config.readInto( *output_name, "output_name" );
    if ( output_name->empty()==true || CheckFilename( output_name->c_str(), ".vtk" ) )
    {
        printf("ERROR::Output filename either is not set or doesn't have '.vtk' extension\n");
        return Usage();
    }
    
    printf("inputs: %s  %s  %s\n", model_name.c_str(), clines_name.c_str(), neck_name.c_str());
    
    printf("outputs: %s¨\n", output_name->c_str());
    
    // Read other parameters
    config.readInto(    *clipFactor,    "clipFactor" );
    config.readInto(    *aneuType,    "aneuType" );
    
    
    // Read optional parameters
    config.readInto(    *verbose,    "verbose" );
    config.readInto(    *bRender,    "bRender" );    
    
    // Read model polydata from input file
    vtkSmartPointer<vtkPolyDataReader> mreader = vtkSmartPointer<vtkPolyDataReader>::New();
    mreader->SetFileName(model_name.c_str());
    mreader->Update();
    if(mreader->GetErrorCode() != 0)
    {
        cout << "ERROR CODE :: " 
             << vtkErrorCode::GetStringFromErrorCode(mreader->GetErrorCode()) << endl;
        return 1;
    }

    //define MODEL polydata
    model->ShallowCopy(mreader->GetOutput());

    
    if ( *verbose == true )
    {
        cout<< "VTK polygonal data surface has:\n\t" 
            << model->GetNumberOfPoints()  << " Points\n\t"
            << model->GetNumberOfCells ()  << " Cells"
            << "\n\t\tverts : " << model->GetNumberOfVerts ()
            << "\n\t\tlines : " << model->GetNumberOfLines ()
            << "\n\t\tpolys : " << model->GetNumberOfPolys ()
            << "\n\t\tstrips: " << model->GetNumberOfStrips() << endl;
    }
   
    // Read centerlines polydata from input file
    vtkSmartPointer<vtkPolyDataReader> creader = vtkSmartPointer<vtkPolyDataReader>::New();
    creader->SetFileName(clines_name.c_str());
    creader->Update();
    
    if(creader->GetErrorCode() != 0)
    {
        cout << "ERROR CODE :: " 
             << vtkErrorCode::GetStringFromErrorCode(creader->GetErrorCode()) << endl;
        return 1;
    }

    // define centerlines polydata
    clines->ShallowCopy(creader->GetOutput());
    clines->BuildLinks();           // Create upward links from points to cells that use each point.  
    clines->BuildCells();           // Create data structure that allows random access of cells
    
    if ( *verbose == true )
    {
        cout<< "VTK polygonal data centerlines has:\n\t" 
            << clines->GetNumberOfPoints()  << " Points\n\t"
            << clines->GetNumberOfCells ()  << " Cells"
            << "\n\t\tverts : " << clines->GetNumberOfVerts ()
            << "\n\t\tlines : " << clines->GetNumberOfLines ()
            << "\n\t\tpolys : " << clines->GetNumberOfPolys ()
            << "\n\t\tstrips: " << clines->GetNumberOfStrips() << endl;
    }    
        
    // Read neck polydata from input file
    vtkSmartPointer<vtkPolyDataReader> nreader = vtkSmartPointer<vtkPolyDataReader>::New();
    nreader->SetFileName(neck_name.c_str());
    nreader->Update();
        
    if(nreader->GetErrorCode() != 0)
    {
        cout << "ERROR CODE :: " 
             << vtkErrorCode::GetStringFromErrorCode(nreader->GetErrorCode()) << endl;
        return 1;
    }
    
    // define AutoNeck polydata
    neck->ShallowCopy(nreader->GetOutput());
    
    if ( *verbose == true )
    {
        cout<< "VTK polygonal data plane neck has:\n\t" 
            << neck->GetNumberOfPoints()  << " Points\n\t"
            << neck->GetNumberOfCells ()  << " Cells"
            << "\n\t\tverts : " << neck->GetNumberOfVerts ()
            << "\n\t\tlines : " << neck->GetNumberOfLines ()
            << "\n\t\tpolys : " << neck->GetNumberOfPolys ()
            << "\n\t\tstrips: " << neck->GetNumberOfStrips() << endl;
    }
    
    return 0;    
}



