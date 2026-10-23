/* ===============================================================================
     File: @(#)contour_interpolate.cpp
     Author: Christoph Hormann <chris_hormann@gmx.de>
     Creation Date: $ 24 May 2006, 12:23:02 $
     Last modified: $ 10 Jul 2007, 21:33:31 $
   -------------------------------------------------------------------------------
     Description: contour line to DEM interpolation 

		 Version history:

		   0.3: initial public version

   -------------------------------------------------------------------------------

     Copyright (C) 2006-2007 Christoph Hormann

     This program is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with this program.  If not, see <http://www.gnu.org/licenses/>.

   ===============================================================================
 */

#include "DemProc.h"

int main(int argc,char **argv)
{
	std::fprintf(stderr,"%s\n", PROGRAM_TITLE);
	std::fprintf(stderr,"-------------------------------------------------------\n");

	cimg_usage("Usage : contour_interpolate [options]");

	// --- Read command line parameters ---

	const char *out_opts = cimg_option("-o",(char*)NULL,"Output File (with corner coordinates and resolution)");
	const bool Disp = cimg_option("-d",false,"Display");
	const float SeaLev = cimg_option("-s",0.0,"Sea Level");
	const float SeaLevOut = cimg_option("-so",0.0,"Sea Level (Output)");
	const float MaxErr = cimg_option("-e",1000.0,"Error Bound");
	const bool ClipZero = cimg_option("-z",false,"Clip heights below sea level");
	const int SeaContWater = cimg_option("-sz",0,"Radius to set Sea Contours to water (0, 1 or 2)");
	const bool NoFill = cimg_option("-n",false,"Test run (no interpolation)");
	const int GrStep = cimg_option("-g",4,"Gradient Calculation stepping");
	char file_o[256];
	double x1;
	double x2;
	double y1;
	double y2;
	int SzX;
	int SzY;
	// first parameter is the file name, 
	// second to fifth are the area bounding box
	// and the last two the x and y resolution
	std::sscanf(out_opts,"%250[^:]:%lf:%lf:%lf:%lf:%d:%d", file_o, &x1, &x2, &y1, &y2, &SzX, &SzY);

	const bool helpflag = cimg_option("-h",false,"Display this help");
	if (helpflag) std::exit(0);

	if (file_o == NULL)
	{
		std::fprintf(stderr,"You must specify an output image (try '%s -help').\n\n",argv[0]);
		std::exit(0); 
	}

	DemProc ContourProcessor(file_o, x1, x2, y1, y2, SzX, SzY, SeaLev, SeaLevOut, MaxErr, Disp);

	// several input files are allowed, each consisting of contour file, 
	// sea water mask and undefined height contour file.  They are drawn 
	// one over the other, expected to be the same resolution but can be 
	// different size although they all have to start at the upper right corner
  for (int k=1; k<argc; k++)
		if (!cimg::strcasecmp(argv[k],"-i"))
		{
			for (int j=k+1; j<argc; j++)
			{
				char fnm[256];
				char fnm_w[256];
				char fnm_sw[256];
				std::sscanf(argv[j],"%250[^:]:%250[^:]:%250[^:]", fnm, fnm_w, fnm_sw);
				ContourProcessor.ReadContoursWhole(fnm, fnm_w, fnm_sw);
			}
		}

	if (!NoFill)
	{
		ContourProcessor.DEMPreProc();
		ContourProcessor.DEMMinMax();
		ContourProcessor.DEMGradient(GrStep);
		ContourProcessor.DEMPostProc();
	}

	ContourProcessor.SaveDEM(ClipZero, SeaContWater); 	
}
