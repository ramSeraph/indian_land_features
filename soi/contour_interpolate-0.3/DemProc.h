/* ===========================================================
     File: @(#)DemProc.h
     Author: chris_hormann@gmx.de
     Creation Date: $ 26 May 2006, 21:01:51 $
     Last modified: $ 10 Jul 2007, 21:43:22 $
   -----------------------------------------------------------
     Description: DEM processing class for contour line -> DEM conversion

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
   ===========================================================
 */

#ifndef _DEMPROC_H
#define _DEMPROC_H 1

#define cimg_use_magick 1

#include "CImg.h"

const char PROGRAM_TITLE[] = "contour line interpolation version 0.3";

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

using namespace cimg_library;

struct GradientLinePoint {
	float x,y;
	int xi,yi;
	float t;
};

struct TerrainPoint {
	int x,y;
	int type;
};	

struct HeightPoint {
	int x,y;
	float h;
};

struct Pnt2Dint	{
	int x,y;
	Pnt2Dint(int nx, int ny) : x(nx), y(ny) { }
	Pnt2Dint() : x(0), y(0) { }
};

/**
	 This class wraps all the contour interpolation
*/
class DemProc
{
 public:
	DemProc(const char *file_o, const double x1, const double x2, const double y1, const double y2, const int SizeX, const int SizeY, const float SeaLev, const float SeaLevOut, const float MaxErr, const bool Disp = true);
	~DemProc();

	/// read contours and draw them into the work image
	int ReadContoursWhole(const char *fnm, const char *fnm_w, const char *fnm_sw);
	/// does the initial interpolation forming the basis of the further steps
	int DEMPreProc();
	/// find and mark local minima/maxima
	int DEMMinMax();
	/// calculate and draw gradient lines
	int DEMGradient(int grid = 1);
	/// new interpolation based on the contours plus gradient lines etc.
	int DEMPostProc();
	/// write the results
	int SaveDEM(const bool ClipZero = false, const int SeaContWater = false);

 private:
	static const float OutsideHeightVal;
	static const float EmptyHeightValLoad;
	static const float EmptyHeightVal;
	static const float SWHeightVal;
	static const float SeaHeightVal;

	static const Pnt2Dint	GridDirections[];
	static const Pnt2Dint	GridDirections4[];
	static const Pnt2Dint	ContourDirections[];
	static const Pnt2Dint	ContourDirectionsLeft[];
	static const Pnt2Dint	ContourDirectionsRight[];

	/// checks if a certain coordinate is no more than Rim from the image borders
	bool IsInImage(int x, int y, int Rim = 0);
	/// fill image from a certain point until the next contours with the value of those contours
	void FillRec(CImg<float> *Img, int x, int y, int dx, int dy, float diff, float MaxRad = 1.0e6);
	/// directions for gradient line calculations
	bool MaxDir(CImg<float> *Img, int x, int y, int &Nx, int &Ny, float &HMax);
	bool MinDir(CImg<float> *Img, int x, int y, int &Nx, int &Ny, float &HMin);
	/// determine gradint of the image at a certain point
	float ImgGradient(CImg<float> *Img, float pxd, float pyd, float &Gradx, float &Grady, int &MinMax);
	/// do interpolation using Laplace function approximation
	void ImgLaplace2(CImg<float> *Mod, CImg<float> *Fixed, float &Radf, const int Cnt = 100, const int Rad = 5);
	/// do interpolation using direction search
	void DirSearchFill(CImg<float> *Mod, CImg<float> *Fixed, const int Rad = 5);

	std::string m_OutputFile;

	double m_MinX;
	double m_MaxX;
	double m_MinY;
	double m_MaxY;
	bool m_ContourLoaded;

	double m_DEMResX;
	double m_DEMResY;
	int m_DEMWidth;
	int m_DEMHeight;
	int m_MaxDist;
	int m_MinMax_Step;
	float m_MaxErr;
	float m_SeaLevel;
	float m_SeaLevelOut;

	float m_MaxDiff;

	OGRDataSource *m_poDS;
	GDALDataset *m_poDstDS; 
	GDALDataset *m_poDstDS_test; 

  CImg<float> *m_VectorImage;
  CImg<float> *m_DEMImage;
  CImg<float> *m_DEMImageS;
  CImg<float> *m_MaxImage;
  CImg<unsigned char> *m_MaxLineImage;
  CImg<unsigned char> *m_SeaMaskImage;
  CImg<float> *m_GradImage;
  CImg<float> *m_ErrImage;
  CImg<int> *m_SmallWaterImage;
  CImg<unsigned char> *m_WaterImage;
	int m_SWcount;

	CImgDisplay *m_Display1;
	CImgDisplay *m_Display2;
	bool m_Disp;
};

#endif /* _DEMPROC_H */

