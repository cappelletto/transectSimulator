// Extracted from http://geomalgorithms.com/a03-_inclusion.html

// Copyright 2000 softSurfer, 2012 Dan Sunday
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.
 
///Basic C and C++ libraries
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdlib.h>

#include <stdio.h>
#include <math.h>

#include "winding_test.h"

using namespace std;

#define VERSION 1

int main(int argc, char *argv[]) {

//*********************************************************************************
/*	PARSER section */
/*  Uses built-in OpenCV parsing method cv::CommandLineParser. It requires a string containing the arguments to be parsed from
	the command line. Further details can be obtained from opencv webpage
*/
    String keys =
            "{@input |<none>  | Input file containing shape polygon}"    // input image is the first argument (positional)
                    "{@output |<none> | Output file name, where transect coordinates will be stored}" // output prefix is the second argument (positional)
                    "{n      |10  | Desired number of transects}"
                    "{l      |30  | Desired transect length}"
                    "{l      |1  | Desired transect width}"
                    "{help h usage ?  |      | show this help message}";      // optional, show help optional

    CommandLineParser cvParser(argc, argv, keys);
    cvParser.about("random transect geeration module");	//adds "about" information to the parser method

	//if the number of arguments is lower than 3, or contains "help" keyword, then we show the help
	if (argc < 3 || cvParser.has("help")) {
        cout << "Random transect generator, delimited into a polygon describing the container shape" << endl;
        cvParser.printMessage();
        cout << endl << "\tExample:" << endl;
        cout << "\t$ transect_gen -n=100 -l=30 -w=5 input.csv transect_xy.txt" << endl;
        cout <<
        "\tThis will find the polygon decribed in input.csv, generate 100 random transects of size 30x5 (units)" << endl << endl;
        return 0;
    }

    String InputFile = cvParser.get<cv::String>(0);		//String containing the input file path+name from cvParser function
    String OutputFile = cvParser.get<cv::String>(1);	//String containing the output file template from cvParser function
    ostringstream OutputFileName;						// output string that will contain the desired output file name

	int   N = cvParser.get<int>("n"); 	// gets argument -n=X, where X is the number of transects to generate
    float L = cvParser.get<float>("l");		// gets argument -l=S, where S is the length of the transect
    float W = cvParser.get<float>("w");		// gets argument -w=S, where S is the width of the transect

	// Check if occurred any error during parsing process
    if (! cvParser.check()) {
        cvParser.printErrors();
        return -1;
    }

    cout << "***************************************" << endl;
    cout << "Input: " << InputFile << endl;

	printf ("Transect generator - version: %d\n", VERSION);
	return 0;
}


