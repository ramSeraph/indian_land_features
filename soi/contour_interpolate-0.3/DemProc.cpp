/* ===========================================================
     File: @(#)DemProc.cpp
     Author: chris_hormann@gmx.de
     Creation Date: $ 26 May 2006, 21:01:09 $
     Last modified: $ 10 Jul 2007, 21:37:13 $
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

#include "DemProc.h"

const float DemProc::OutsideHeightVal = -150000.0;
const float DemProc::EmptyHeightValLoad = -100000.0;
const float DemProc::EmptyHeightVal = 0.0;
const float DemProc::SeaHeightVal = -2.0;
const float DemProc::SWHeightVal = -120000.0;

const Pnt2Dint DemProc::GridDirections4[]=
{
	Pnt2Dint(1, 0),
	Pnt2Dint(0, -1),
	Pnt2Dint(-1, 0),
	Pnt2Dint(0, 1)
};

DemProc::DemProc(const char *file_o, const double x1, const double x2, const double y1, const double y2, const int SizeX, const int SizeY, const float SeaLev, const float SeaLevOut, const float MaxErr, const bool Disp)
{
	m_OutputFile = std::string(file_o);

	m_Display1 = NULL;
	m_Display2 = NULL;
	m_Disp = Disp;

  m_VectorImage = NULL;
  m_SeaMaskImage = NULL;
  m_DEMImage = NULL;
  m_DEMImageS = NULL;
  m_MaxImage = NULL;
  m_GradImage = NULL;
	m_MaxLineImage = NULL;
  m_ErrImage = NULL;

	m_poDS = NULL;
	m_poDstDS = NULL;

	m_MinMax_Step = 8;
	m_MaxDist = 120;
	m_MaxErr = MaxErr;

	m_SWcount = 0;
	m_MinX = x1;
	m_MaxX = x2;
	m_MinY = y1;
	m_MaxY = y2;
	m_DEMResX = (x2-x1)/(SizeX-1);
	m_DEMResY = (y2-y1)/(SizeY-1);
	m_DEMWidth = SizeX;
	m_DEMHeight = SizeY;
	m_SeaLevel = SeaLev;
	m_SeaLevelOut = SeaLevOut;

	std::fprintf(stderr,"---Initializing Interpolator...\n");
	std::fprintf(stderr,"-------------------------------------------------------\n");
	std::fprintf(stderr,"coordinates: x: %f -> %f\n", m_MinX, m_MaxX);
	std::fprintf(stderr,"             y: %f -> %f\n", m_MinY, m_MaxY);
	std::fprintf(stderr,"             size: %d x %d\n", m_DEMWidth, m_DEMHeight);
	std::fprintf(stderr,"             resolution: %f x %f\n", m_DEMResX, m_DEMResY);

  m_VectorImage = new CImg<float>(m_DEMWidth, m_DEMHeight, 1, 1, DemProc::EmptyHeightVal);
  m_SeaMaskImage = new CImg<unsigned char>(m_DEMWidth, m_DEMHeight, 1, 1, 0);
  m_DEMImage = new CImg<float>(m_DEMWidth, m_DEMHeight, 1, 1, DemProc::EmptyHeightVal);
  m_DEMImageS = new CImg<float>(m_DEMWidth, m_DEMHeight, 1, 1, DemProc::EmptyHeightVal);
	m_SmallWaterImage = new CImg<int>(m_DEMWidth, m_DEMHeight, 1, 1, 0);
  m_WaterImage = new CImg<unsigned char>(m_DEMWidth, m_DEMHeight, 1, 1, 0);

	m_ContourLoaded = false;

	m_MaxDist = 120;
	m_MinMax_Step = 40;

	GDALAllRegister();

	const char *pszFormat = "GTiff";
	GDALDriver *poDriver;

	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);

	if(poDriver == NULL)
	{
		std::fprintf(stderr,"GDAL initialisation failed.\n\n");
	}

	char **papszOptions = NULL;

	m_poDstDS = poDriver->Create(m_OutputFile.c_str(), m_DEMWidth, m_DEMHeight, 1, GDT_Float32, papszOptions);
	m_poDstDS_test = poDriver->Create("test_cont.tif", m_DEMWidth, m_DEMHeight, 1, GDT_Float32, papszOptions);
}

DemProc::~DemProc()
{
	if (m_Display1 != NULL)
		delete m_Display1;
	if (m_Display2 != NULL)
		delete m_Display2;

	if (m_VectorImage != NULL)
		delete m_VectorImage;
	if (m_SeaMaskImage != NULL)
		delete m_SeaMaskImage;
	if (m_DEMImage != NULL)
		delete m_DEMImage;
	if (m_DEMImageS != NULL)
		delete m_DEMImageS;
	if (m_MaxImage != NULL)
		delete m_MaxImage;
	if (m_GradImage != NULL)
		delete m_GradImage;
	if (m_MaxLineImage != NULL)
		delete m_MaxLineImage;

	if (m_poDS != NULL)
		OGRDataSource::DestroyDataSource( m_poDS );
	if (m_poDstDS != NULL)
		delete m_poDstDS;
	if (m_poDstDS_test != NULL)
		delete m_poDstDS_test;
}


bool DemProc::IsInImage(int x, int y, int Rim)
{
	return ((x < m_DEMWidth-Rim) && (y < m_DEMHeight-Rim) && (x >= Rim) && (y >= Rim));
}

void DemProc::FillRec(CImg<float> *Img, int x, int y, int dx, int dy, float diff, float MaxRad)
{
	(*Img)(x+dx,y+dy,0,0) = (*m_DEMImageS)(x+dx,y+dy,0,0) - diff;
	for (int i = 0; i < 4; i++)
	{
		int dxn = dx + GridDirections4[i].x;
		int dyn = dy + GridDirections4[i].y;
		if (IsInImage(x+dxn, y+dyn))
			if (((*m_VectorImage)(x+dxn, y+dyn,0,0) == DemProc::EmptyHeightVal) && ((*Img)(x+dxn, y+dyn,0,0) == DemProc::EmptyHeightVal))
			{
				float Dist = std::sqrt(dxn*dxn + dyn*dyn);
				if (Dist < MaxRad-1.1)
					FillRec(Img, x, y, dxn, dyn, diff, MaxRad);
			}
	}
}

bool DemProc::MaxDir(CImg<float> *Img, int x, int y, int &Nx, int &Ny, float &HMax)
{
	float H;
	Nx = x;
	Ny = y;
	bool Found = false;

	H = (*Img)(x+1,y,0,0);
	if (H-0.01 > HMax)
	{
		HMax = H;
		Nx = x+1;
		Ny = y;
		Found = true;
	}
	H = (*Img)(x-1,y,0,0);
	if (H-0.01 > HMax)
	{
		HMax = H;
		Nx = x-1;
		Ny = y;
		Found = true;
	}
	H = (*Img)(x,y+1,0,0);
	if (H-0.01 > HMax)
	{
		HMax = H;
		Nx = x;
		Ny = y+1;
		Found = true;
	}
	H = (*Img)(x,y-1,0,0);
	if (H-0.01 > HMax)
	{
		HMax = H;
		Nx = x;
		Ny = y-1;
		Found = true;
	}
	H = (*Img)(x-1,y-1,0,0);
	if (H-0.01 > HMax)
	{
		HMax = H;
		Nx = x-1;
		Ny = y-1;
		Found = true;
	}
	H = (*Img)(x+1,y+1,0,0);
	if (H-0.01 > HMax)
	{
		HMax = H;
		Nx = x+1;
		Ny = y+1;
		Found = true;
	}
	H = (*Img)(x+1,y-1,0,0);
	if (H-0.01 > HMax)
	{
		HMax = H;
		Nx = x+1;
		Ny = y-1;
		Found = true;
	}
	H = (*Img)(x-1,y+1,0,0);
	if (H-0.01 > HMax)
	{
		HMax = H;
		Nx = x-1;
		Ny = y+1;
		Found = true;
	}
	return Found;
}

bool DemProc::MinDir(CImg<float> *Img, int x, int y, int &Nx, int &Ny, float &HMin)
{
	float H;
	Nx = x;
	Ny = y;
	bool Found = false;

	H = (*Img)(x+1,y,0,0);
	if (H+0.01 < HMin)
	{
		HMin = H;
		Nx = x+1;
		Ny = y;
		Found = true;
	}
	H = (*Img)(x-1,y,0,0);
	if (H+0.01 < HMin)
	{
		HMin = H;
		Nx = x-1;
		Ny = y;
		Found = true;
	}
	H = (*Img)(x,y+1,0,0);
	if (H+0.01 < HMin)
	{
		HMin = H;
		Nx = x;
		Ny = y+1;
		Found = true;
	}
	H = (*Img)(x,y-1,0,0);
	if (H+0.01 < HMin)
	{
		HMin = H;
		Nx = x;
		Ny = y-1;
		Found = true;
	}
	H = (*Img)(x-1,y-1,0,0);
	if (H+0.01 < HMin)
	{
		HMin = H;
		Nx = x-1;
		Ny = y-1;
		Found = true;
	}
	H = (*Img)(x+1,y+1,0,0);
	if (H+0.01 < HMin)
	{
		HMin = H;
		Nx = x+1;
		Ny = y+1;
		Found = true;
	}
	H = (*Img)(x+1,y-1,0,0);
	if (H+0.01 < HMin)
	{
		HMin = H;
		Nx = x+1;
		Ny = y-1;
		Found = true;
	}
	H = (*Img)(x-1,y+1,0,0);
	if (H+0.01 < HMin)
	{
		HMin = H;
		Nx = x-1;
		Ny = y+1;
		Found = true;
	}
	return Found;
}


float DemProc::ImgGradient(CImg<float> *Img, float pxd, float pyd, float &Gradx, float &Grady, int &MinMax)
{
	float Hgt = Img->linear_pix2d(pxd,pyd);
	float Hgt_xm = Img->linear_pix2d(pxd-0.1,pyd);
	float Hgt_xp = Img->linear_pix2d(pxd+0.1,pyd);
	float Hgt_ym = Img->linear_pix2d(pxd,pyd-0.1);
	float Hgt_yp = Img->linear_pix2d(pxd,pyd+0.1);

	if ((Hgt_xp > Hgt) && (Hgt_xm > Hgt) && (Hgt_yp > Hgt) && (Hgt_ym > Hgt))
	{
		Gradx = 0;
		Grady = 0;
		MinMax = -1;
		return Hgt;
	}

	if ((Hgt_xp < Hgt) && (Hgt_xm < Hgt) && (Hgt_yp < Hgt) && (Hgt_ym < Hgt))
	{
		Gradx = 0.0;
		Grady = 0.0;
		MinMax = 1;
		return Hgt;
	}

	Gradx = (Hgt_xp-Hgt_xm);
	Grady = (Hgt_yp-Hgt_ym);
	float GradLen = sqrt(Gradx*Gradx + Grady*Grady);
	if (GradLen > 1.0e-6)
	{
		Gradx /= GradLen;
		Grady /= GradLen;
	}
	else
	{
		Gradx = 0.0;
		Grady = 0.0;
		MinMax = 2;
		return Hgt;
	}

	MinMax = 0;
	return Hgt;
}

void DemProc::ImgLaplace2(CImg<float> *Mod, CImg<float> *Fixed, float &Radf, const int Cnt, const int Rad)
{
	// prepare the weighting factors for the initial spread filter kernel
	CImg<float> WgtArr(Rad*2+1, Rad*2+1);

	for (int dy = -Rad; dy<=Rad; dy++)
	{
		for (int dx = -Rad; dx<=Rad; dx++)
		{
			if ((dx == 0) && (dy == 0))
				WgtArr(dx+Rad,dy+Rad) = 5.0;
			else
				WgtArr(dx+Rad,dy+Rad) = 1.0/sqrt(dx*dx+dy*dy);
		}
	}

	if (m_Disp)
	{
		if (m_Display1 == NULL)
			m_Display1 = new CImgDisplay(*Mod,PROGRAM_TITLE,0,3);
		if (m_Display2 == NULL)
			m_Display2 = new CImgDisplay(Fixed->get_normalize(0.0, 255.0),PROGRAM_TITLE,0,3);
		else
			Fixed->get_normalize(0.0, 255.0).display(*m_Display2);
	}

	std::fprintf(stderr,"Laplace initial spread...\n");

	CImg<float> Orig(*Fixed);

	// do the spread step a few times (increase this or use a larger filter kernel)
	for (int i = 0; i < 5; i++)
	{
		cimg_forXY((*Mod),qx,qy)
		{
			if ((*Fixed)(qx,qy,0,0) == DemProc::EmptyHeightVal)
			{
				float WSum = 0.0;
				(*Mod)(qx,qy,0,0) = 0.0;				
				for (int py = (qy-Rad); py<=(qy+Rad); py++)
				{
					if ((py >= 0) && (py < Orig.dimy()))
					{
						int dy = (py-qy);
						for (int px = (qx-Rad); px<=(qx+Rad); px++)
						{
							if ((px >= 0) && (px < Orig.dimx()))
							{
								if (Orig(px,py,0,0) != DemProc::EmptyHeightVal)
								{
									int dx = (px-qx);
									(*Mod)(qx,qy,0,0) += WgtArr(dx+Rad,dy+Rad)*Orig(px,py,0,0);
									WSum += WgtArr(dx+Rad,dy+Rad);
								}
							}
						}
					}
				}
				if (WSum > 0.0)
				{
					(*Mod)(qx,qy,0,0) /= WSum;
				}
			}
			else
				(*Mod)(qx,qy,0,0) = (*Fixed)(qx,qy,0,0);
		}

		if (m_Disp)
			Mod->get_normalize(0.0, 255.0).display(*m_Display1);

		cimg_forXY(Orig,x,y) { Orig(x,y,0,0) = (*Mod)(x,y,0,0); }

	}

	std::fprintf(stderr,"Laplace spread...\n");

	// do the actual Laplace approximation by alternating 
	// smoothing and setting the countour pixels to their values
	// it might well be that using the same approach as above
	// (not changing the contour pixels themseleves) is more efficient
	// but this was fast to code using existing CImg::blur()
	// and is reasonably fast.
	for (int i = 0; i < Cnt; i++)
	{
		cimg_forXY(*Mod,x,y)
		{
			if ((*Fixed)(x,y,0,0) != DemProc::EmptyHeightVal)
				(*Mod)(x,y,0,0) = (*Fixed)(x,y,0,0);
		}
		Mod->blur(Radf);

		if (m_Disp)
			Mod->get_normalize(0.0, 255.0).display(*m_Display1);

		std::fprintf(stderr,"\r  Iteration %d of %d...  ", i, Cnt);
	}

	float MaxDiff = 0.0;

	// calculate the error after the last blur.
	cimg_forXY(*Mod,x,y)
	{
		if ((*Fixed)(x,y,0,0) != DemProc::EmptyHeightVal)
		{
			MaxDiff = cimg::max(MaxDiff, cimg::abs((*Mod)(x,y,0,0)-(*Fixed)(x,y,0,0)));
			(*Mod)(x,y,0,0) = (*Fixed)(x,y,0,0);
		}
	}

	Radf = cimg::max(1.0,0.75*MaxDiff/m_MaxErr);

	std::fprintf(stderr,"\n  Finished (max error: %f -> %f)!\n", MaxDiff, Radf);
}


void DemProc::DirSearchFill(CImg<float> *Mod, CImg<float> *Fixed, const int Rad)
{
	CImg<float> WgtArr(Rad*2+1, Rad*2+1);

	for (int dy = -Rad; dy<=Rad; dy++)
	{
		for (int dx = -Rad; dx<=Rad; dx++)
		{
			if ((dx == 0) && (dy == 0))
				WgtArr(dx+Rad,dy+Rad) = 5.0;
			else
				WgtArr(dx+Rad,dy+Rad) = 1.0/sqrt(dx*dx+dy*dy);
		}
	}

	if (m_Disp)
	{
		if (m_Display1 == NULL)
			m_Display1 = new CImgDisplay(*Mod,PROGRAM_TITLE,0,3);
		if (m_Display2 == NULL)
			m_Display2 = new CImgDisplay(Fixed->get_normalize(0.0, 255.0),PROGRAM_TITLE,0,3);
		else
			Fixed->get_normalize(0.0, 255.0).display(*m_Display2);
	}

	std::fprintf(stderr,"Direction search fill...\n");

	CImg<float> Orig(*Fixed);

	// interpolate by sampling in 32 directions and
	// doing weighted average of the found heights
	cimg_forY((*Mod),qy)
	{
		cimg_forX((*Mod),qx)
		{
			if ((*Fixed)(qx,qy,0,0) == DemProc::EmptyHeightVal)
			{
				float WSum = 0.0;
				(*Mod)(qx,qy,0,0) = 0.0;	
				for (int i = 0; i < 32; i++)
				{
					for (int r = 1; r <= Rad; r++)
					{
						int px = qx + int(sin(2.0*cimg::PI*float(i)/32.0)*r);
						int py = qy + int(cos(2.0*cimg::PI*float(i)/32.0)*r);
						if (IsInImage(px, py))
							if (Orig(px,py,0,0) != DemProc::EmptyHeightVal)
							{
								int dy = (py-qy);
								int dx = (px-qx);
								(*Mod)(qx,qy,0,0) += WgtArr(dx+Rad,dy+Rad)*Orig(px,py,0,0);
								WSum += WgtArr(dx+Rad,dy+Rad);
								break;
							}
					}
				}
				if (WSum > 0.0)
					(*Mod)(qx,qy,0,0) /= WSum;
			}
			else
				(*Mod)(qx,qy,0,0) = (*Fixed)(qx,qy,0,0);
		}
		std::fprintf(stderr,"\r  line %d of %d...  ", qy, Mod->dimy());
	}

	if (m_Disp)
		Mod->get_normalize(0.0, 255.0).display(*m_Display1);

	std::fprintf(stderr,"\n  Finished!\n");

}

int DemProc::ReadContoursWhole(const char *fnm, const char *fnm_w, const char *fnm_sw)
{
	// special magick for SRTM based input: no water and undefined height files
	if (std::strncmp(fnm_w, "srtm", 4) == 0)
	{
		std::fprintf(stderr,"---Drawing contour file %s (whole area, srtm)...\n", fnm);
		CImg<float> Contours = CImg<float>().get_load(fnm);

		cimg_forXY(Contours,px,py)
		{
			if (IsInImage(px,py))
			{
				if ((*m_VectorImage)(px,py) == 0)
				{
					if ((Contours(px,py) < -1000.0) || (Contours(px,py) > 10000.0))
						(*m_VectorImage)(px,py) = DemProc::EmptyHeightVal;
					else
						(*m_VectorImage)(px,py) = Contours(px,py);
				}
			}
		}
	}
	else
	{
		std::fprintf(stderr,"---Drawing contour file %s (whole area)...\n", fnm);
		CImg<float> Contours = CImg<float>().get_load_cimg(fnm);

		if (fnm_w != NULL)
		{
			std::fprintf(stderr,"---Drawing water mask file %s...\n", fnm_w);
			CImg<unsigned char> Contours_W = CImg<unsigned char>().get_load(fnm_w);

			cimg_forXY(Contours_W,px,py)
			{
				if (IsInImage(px,py))
				{
					if (Contours_W(px,py) > 10)
					{
						if ((*m_WaterImage)(px,py) == 0)
							(*m_WaterImage)(px,py) = Contours_W(px,py);
					}
				}
			}
		}

		cimg_forXY(Contours,px,py)
		{
			if (IsInImage(px,py))
			{
				if (((*m_VectorImage)(px,py) == 0.0) || (((*m_WaterImage)(px,py) > 110) && (Contours(px,py) != DemProc::EmptyHeightValLoad)))
				{
					if ((Contours(px,py) == DemProc::EmptyHeightVal) && (DemProc::EmptyHeightVal != EmptyHeightValLoad))
					{
						(*m_VectorImage)(px,py) = DemProc::EmptyHeightVal+0.001;
					}
					else if (Contours(px,py) == DemProc::SWHeightVal)
					{
						(*m_VectorImage)(px,py) = DemProc::EmptyHeightVal;
					}
					else if (Contours(px,py) == DemProc::EmptyHeightValLoad)
						(*m_VectorImage)(px,py) = DemProc::EmptyHeightVal;
					else if (Contours(px,py) != DemProc::OutsideHeightVal)
					{
						if (Contours(px,py) == DemProc::SeaHeightVal)
						{
							(*m_VectorImage)(px,py) = m_SeaLevel+0.001;
							(*m_SeaMaskImage)(px,py) = 1;
						}
						else
							(*m_VectorImage)(px,py) = Contours(px,py);
					}
				}
				else
				{
					if (Contours(px,py) == DemProc::SeaHeightVal)
					{
						(*m_VectorImage)(px,py) = m_SeaLevel+0.001;
						(*m_SeaMaskImage)(px,py) = 1;
					}
				}
			}
		}

		if (fnm_sw != NULL)
		{
			std::fprintf(stderr,"---Drawing undefined height contour file %s...\n", fnm_sw);
			CImg<int> Contours_SW = CImg<int>().get_load_cimg(fnm_sw);

			int SWc = m_SWcount;

			cimg_forXY(Contours_SW,px,py)
			{
				if (IsInImage(px,py))
				{
					if (Contours_SW(px,py) != 0)
					{
						(*m_SmallWaterImage)(px,py) = Contours_SW(px,py)+SWc;
						if (Contours_SW(px,py)+SWc > m_SWcount) m_SWcount = Contours_SW(px,py)+SWc;
					}
				}
			}
		}
	}
	return 0;
}

int DemProc::DEMPreProc()
{
	std::fprintf(stderr,"---Preprocessing DEM...\n");

	DirSearchFill(m_DEMImage, m_VectorImage, 120);

	if (m_SWcount > 0)
	{
		std::fprintf(stderr,"Drawing %d small lakes...\n", m_SWcount);
		for (int i = 1; i <= m_SWcount; i++)
		{
			std::fprintf(stderr,"\r  %d of %d...", i, m_SWcount);
			float HgtAvg = 0.0;
			int HgtCnt = 0;
			cimg_forXY(*m_VectorImage,px,py)
			{
				if ((*m_SmallWaterImage)(px,py) == i)
				{
					HgtAvg += (*m_DEMImage)(px,py);
					HgtCnt++;
				}
			}
			std::fprintf(stderr,"...  ");
			HgtAvg /= HgtCnt;
			cimg_forXY(*m_VectorImage,px,py)
			{
				if ((*m_SmallWaterImage)(px,py) == i)
				{
					(*m_VectorImage)(px,py) = (*m_DEMImage)(px,py) = HgtAvg;
				}
			}
		}
		std::fprintf(stderr,"\nFinished!\n");
	}

	m_DEMImage->blur(4.8);

	CImg<float> img_err_o(*m_VectorImage);
	CImg<float> img_err(*m_VectorImage);

	cimg_forXY(img_err_o,px,py)
	{
		if ((*m_VectorImage)(px,py,0,0) != DemProc::EmptyHeightVal)
			img_err_o(px,py,0,0) = (*m_VectorImage)(px,py,0,0)-(*m_DEMImage)(px,py,0,0);
	}

	float Lev = 2.2;
	float Lev2 = Lev;

	ImgLaplace2(&img_err, &img_err_o, Lev2, 100, 5);

	img_err.blur(Lev);
	Lev2 = 1.0;

	(*m_DEMImage) += img_err*0.85;

	(*m_DEMImageS) = (*m_DEMImage);

	for (int j = 0; j < 2; j++)
	{
		float Lev = (1.5 - 0.8*j);

		CImg<float> img_err2_o(*m_VectorImage);
		CImg<float> img_err2(*m_VectorImage);

		cimg_forXY(img_err2_o,px,py)
		{
			if ((*m_VectorImage)(px,py,0,0) != DemProc::EmptyHeightVal)
				img_err2_o(px,py,0,0) = (*m_VectorImage)(px,py,0,0)-(*m_DEMImage)(px,py,0,0);
		}

		Lev2 = Lev;

		ImgLaplace2(&img_err2, &img_err2_o, Lev2, 100, 15);

		img_err2.blur(Lev/Lev2);

		(*m_DEMImage) += img_err2*0.85;

		float MaxDiff = 0.0;

		cimg_forXY(*m_DEMImage,x,y)
		{
			if ((*m_VectorImage)(x,y,0,0) != DemProc::EmptyHeightVal)
			{
				MaxDiff = cimg::max(MaxDiff, cimg::abs((*m_VectorImage)(x,y,0,0)-(*m_DEMImage)(x,y,0,0)));
			}
		}

		std::fprintf(stderr,"(max error: %f)!\n", MaxDiff);
	}

	m_MaxDiff = 0.0;
	cimg_forXY((*m_DEMImage),x,y)
	{
		if ((*m_VectorImage)(x,y,0,0) != DemProc::EmptyHeightVal)
		{
			m_MaxDiff = cimg::max(m_MaxDiff, cimg::abs((*m_VectorImage)(x,y,0,0)-(*m_DEMImage)(x,y,0,0)));
		}
	}

	std::fprintf(stderr,"---Finished (max error: %f)!\n", m_MaxDiff);

	return 0;
}

int DemProc::DEMMinMax()
{
  m_MaxImage = new CImg<float>(*m_VectorImage);
  m_MaxLineImage = new CImg<unsigned char>(*m_VectorImage, false);
	cimg_forXY((*m_MaxLineImage),x,y)
	{
		(*m_MaxLineImage)(x,y) = 0;
	}

	if (m_Disp)
	{
		if (m_Display1 == NULL)
			m_Display1 = new CImgDisplay(*m_MaxLineImage,PROGRAM_TITLE,0,3);
		else
			m_MaxLineImage->display(*m_Display2);
		
		if (m_Display2 == NULL)
			m_Display2 = new CImgDisplay(m_VectorImage->get_normalize(0.0, 255.0),PROGRAM_TITLE,0,3);
		else
			m_VectorImage->get_normalize(0.0, 255.0).display(*m_Display2);
	}

	std::fprintf(stderr,"---Calculating Minima/Maxima (step=%d)...\n", m_MinMax_Step);

	int MinCnt = 0;
	int MaxCnt = 0;

	std::vector<TerrainPoint> MinMaxPnts;

	TerrainPoint Pnt;

	for (int y = 1; y < m_DEMHeight-1; y += m_MinMax_Step)
	{
		for (int x = 1; x < m_DEMWidth-1; x += m_MinMax_Step)
		{
			if ((*m_MaxLineImage)(x,y,0,0) == 0)
			{
				float Hgt = (*m_DEMImage)(x,y,0,0);
				int MX = x;
				int MY = y;
				int MX2;
				int MY2;
				bool Found = true;
				while (MinDir(m_DEMImage, MX, MY, MX2, MY2, Hgt))
				{
					(*m_MaxLineImage)(MX,MY,0,0) += 2;
					MX = MX2;
					MY = MY2;
					if (!IsInImage(MX, MY, 1))
					{
						Found = false;
						break;
					}
					if ((*m_MaxLineImage)(MX,MY,0,0) > 0)
					{
						Found = false;
						break;
					}
				}
				if (Found)
				{
					if ((*m_VectorImage)(x,y,0,0) == DemProc::EmptyHeightVal)
					{
						(*m_MaxLineImage)(MX,MY,0,0) = 1;
						Pnt.x = MX;
						Pnt.y = MY;
						Pnt.type = -1;
						MinMaxPnts.push_back(Pnt);						
						(*m_MaxImage)(MX,MY,0,0) = (*m_DEMImage)(MX,MY,0,0) - m_MaxDiff*2.5;
						MinCnt++;
						std::fprintf(stderr,"found Min: (%d/%d) (count: %d)...\n", MX, MY, MinCnt);
						std::fprintf(stdout, "Minimum(%d,%d,%f)\n", MX, MY, (*m_DEMImage)(MX,MY,0,0));
					}
				}

				Hgt = (*m_DEMImage)(x,y,0,0);
				MX = x;
				MY = y;
				Found = true;
				while (MaxDir(m_DEMImage, MX, MY, MX2, MY2, Hgt))
				{
					(*m_MaxLineImage)(MX,MY,0,0) += 2;
					MX = MX2;
					MY = MY2;
					if (!IsInImage(MX, MY, 1))
					{
						Found = false;
						break;
					}
					if ((*m_MaxLineImage)(MX,MY,0,0) > 0)
					{
						Found = false;
						break;
					}
				}
				if (Found)
				{
					if ((*m_VectorImage)(x,y,0,0) == DemProc::EmptyHeightVal)
					{
						(*m_MaxLineImage)(MX,MY,0,0) = 1;
						Pnt.x = MX;
						Pnt.y = MY;
						Pnt.type = 1;
						MinMaxPnts.push_back(Pnt);	
						(*m_MaxImage)(MX,MY,0,0) = (*m_DEMImage)(MX,MY,0,0) + m_MaxDiff*2.5;
						MaxCnt++;
						std::fprintf(stderr,"found Max: (%d/%d) (count: %d)...\n", MX, MY, MaxCnt);
						std::fprintf(stdout, "Maximum(%d,%d,%f)\n", MX, MY, (*m_DEMImage)(MX,MY,0,0));
					}
				}
				
			}
		}
		if (m_Disp)
			m_MaxLineImage->get_normalize(0, 255).display(*m_Display1);
	}

	std::fprintf(stderr,"---Finished (count: %d/%d)!\n", MinCnt, MaxCnt);

	if ((MinCnt > 0) && (MaxCnt > 0))
	{
		std::fprintf(stderr,"---Post processing Minima/Maxima...\n");

		for(std::vector<TerrainPoint>::iterator iter = MinMaxPnts.begin(); iter != MinMaxPnts.end(); iter++)
		{
			std::fprintf(stderr,"doing MinMax (%d, %d)...\n", (*iter).x, (*iter).y);
			int dist[] = { 1, 1, 1, 1 };
			int MinDist = -1;
			while (true)
			{
				int Ends = 0;
				for (int i = 0; i < 4; i++)
				{
					if (dist[i] != 0)
					{
						int px = (*iter).x + dist[i]*GridDirections4[i].x;
						int py = (*iter).y + dist[i]*GridDirections4[i].y;
						if (IsInImage(px, py))
						{
							if ((*m_VectorImage)(px,py,0,0) == DemProc::EmptyHeightVal)
								dist[i]++;
							else 
							{
								if (MinDist < 0) MinDist = dist[i];
								Ends++;
							}
						}
						else
						{
							dist[i] = 0;
							Ends++;
						}
					}
					else Ends++;
				}
				if (Ends == 4) break;
			}
			float Difference = 0;
			int DiffCnt = 0;
			for (int i = 0; i < 4; i++)
			{
				if (dist[i] != 0)
				{
					int px = (*iter).x + dist[i]*GridDirections4[i].x;
					int py = (*iter).y + dist[i]*GridDirections4[i].y;
					if (IsInImage(px, py))
					{
						Difference += (*m_DEMImageS)(px,py,0,0) - (*m_VectorImage)(px,py,0,0);
						DiffCnt++;
					}
				}
			}
			Difference /= DiffCnt;

			std::fprintf(stderr,"Fill (type: %d, Difference: %f, MinDist: %d)\n", (*iter).type, Difference, MinDist);

			if (MinDist > 5) MinDist = 5;

			if (std::abs(Difference) > 1.0e-8)
				FillRec(m_MaxImage, (*iter).x, (*iter).y, 0, 0, Difference, (float)MinDist);

			std::fprintf(stderr,"Fill End\n");
		}

		std::fprintf(stderr,"---Finished!\n");
	}
	return 0;
}


int DemProc::DEMGradient(int grid)
{
  m_GradImage = new CImg<float>(*m_MaxImage);

	CImg<unsigned char> img_grad_d(m_DEMWidth, m_DEMHeight, 1, 1, 0);

	if (m_Disp)
	{
		if (m_Display1 == NULL)
			m_Display1 = new CImgDisplay(img_grad_d,PROGRAM_TITLE,0,3);
		if (m_Display2 == NULL)
			m_Display2 = new CImgDisplay(m_MaxImage->get_normalize(0.0, 255.0),PROGRAM_TITLE,0,3);
		else
			m_MaxImage->get_normalize(0.0, 255.0).display(*m_Display2);
	}

	CImg<int> img_grad_cnt(m_DEMWidth, m_DEMHeight, 1, 1, 0);
	img_grad_cnt = 0;

	std::fprintf(stderr,"---Calculating gradient lines...\n");

	for (int y = 0; y < m_DEMHeight; y += grid)
	{
		for (int x = 0; x < m_DEMWidth; x += grid)
		{
			if (((*m_MaxImage)(x,y,0,0) != DemProc::EmptyHeightVal))// && (y < 20))
			{
				{
					int pxi = x;
					int pyi = y;
					int pxi2;
					int pyi2;
					float px = x;
					float py = y;
					float pt = 0.0;

					std::vector<GradientLinePoint> Line;
					GradientLinePoint LinePoint;

					bool EndLine = false;
					float StartHgt = (*m_DEMImage)(x,y,0,0);
					float Hgt = StartHgt, HgtN;

					do
					{
						int MinMax;
						float Gx, Gy;
						HgtN = ImgGradient(m_DEMImage, px, py, Gx, Gy, MinMax);
						if (HgtN > Hgt) 
						{
							Line.clear();
							EndLine  = true;
							std::fprintf(stderr,"Gradient inconsistency (Min).\n");
							break;
						}
						Hgt = HgtN;
						if (MinMax != 0) EndLine  = true;
						px -= Gx*1.0;
						py -= Gy*1.0;
						pxi2 = (int)floor(px+0.5);
						pyi2 = (int)floor(py+0.5);

						if ((pxi2 != pxi) || (pyi2 != pyi))
						{
							if (!IsInImage(pxi2, pyi2, 1))
							{
								EndLine  = true;
								break;
							}
							if ((*m_MaxImage)(pxi2,pyi2,0,0) != DemProc::EmptyHeightVal) { EndLine  = true; break; }
							if ((pxi2 != pxi) && (pyi2 != pyi))
							{
								if ((*m_MaxImage)(pxi,pyi2,0,0) != DemProc::EmptyHeightVal) { EndLine  = true; break; }
								if ((*m_MaxImage)(pxi2,pyi,0,0) != DemProc::EmptyHeightVal) { EndLine  = true; break; }
							}
							pxi = pxi2;
							pyi = pyi2;
						}

						pt += 1.0;

						LinePoint.x = px;
						LinePoint.y = py;
						LinePoint.xi = pxi;
						LinePoint.yi = pyi;
						LinePoint.t = pt;
						Line.push_back(LinePoint);

						if (pt > 1000) { EndLine  = true; break; }
					}
					while (!EndLine);

					if (Line.size() > 1)
					{
						float EndHgt = (*m_DEMImage)(pxi,pyi,0,0);
						float xa = 0.0;
						float ya = 0.0;
						for(std::vector<GradientLinePoint>::iterator iter = Line.begin(); iter != Line.end(); iter++)
						{
							if ((*m_MaxImage)((*iter).xi,(*iter).yi,0,0) == DemProc::EmptyHeightVal)
							{
								(*m_GradImage)((*iter).xi,(*iter).yi,0,0) += StartHgt + ((*iter).t/pt) * (EndHgt-StartHgt);
								img_grad_cnt((*iter).xi,(*iter).yi,0,0) += 1;
								img_grad_d((*iter).xi,(*iter).yi,0,0) = 128;
								if ((xa != 0.0) || (ya != 0.0))
									std::fprintf(stdout, "GradientLineDown(%f,%f,%f,%f)\n", xa, ya, (*iter).x, (*iter).y);
								xa = (*iter).x;
								ya = (*iter).y;
							}
						}
					}
				}
				{
					int pxi = x;
					int pyi = y;
					int pxi2;
					int pyi2;
					float px = x;
					float py = y;
					float pt = 0.0;

					std::vector<GradientLinePoint> Line;
					GradientLinePoint LinePoint;

					bool EndLine = false;
					float StartHgt = (*m_DEMImage)(x,y,0,0);
					float Hgt = StartHgt, HgtN;

					do
					{
						int MinMax;
						float Gx, Gy;
						HgtN = ImgGradient(m_DEMImage, px, py, Gx, Gy, MinMax);
						if (HgtN < Hgt) 
						{
							Line.clear();
							EndLine  = true;
							std::fprintf(stderr,"Gradient inconsistency (Max).\n");
							break;
						}
						Hgt = HgtN;
						if (MinMax != 0) EndLine  = true;
						px += Gx*1.0;
						py += Gy*1.0;
						pxi2 = (int)floor(px+0.5);
						pyi2 = (int)floor(py+0.5);

						if ((pxi2 != pxi) || (pyi2 != pyi))
						{
							if (!IsInImage(pxi2, pyi2, 1))
							{
								EndLine  = true;
								break;
							}
							if ((*m_MaxImage)(pxi2,pyi2,0,0) != DemProc::EmptyHeightVal) { EndLine  = true; break; }
							if ((pxi2 != pxi) && (pyi2 != pyi))
							{
								if ((*m_MaxImage)(pxi,pyi2,0,0) != DemProc::EmptyHeightVal) { EndLine  = true; break; }
								if ((*m_MaxImage)(pxi2,pyi,0,0) != DemProc::EmptyHeightVal) { EndLine  = true; break; }
							}
							pxi = pxi2;
							pyi = pyi2;
						}

						pt += 1.0;

						LinePoint.x = px;
						LinePoint.y = py;
						LinePoint.xi = pxi;
						LinePoint.yi = pyi;
						LinePoint.t = pt;
						Line.push_back(LinePoint);

						if (pt > 1000) { EndLine  = true; break; }
					}
					while (!EndLine);

					if (Line.size() > 1)
				 	{
						float EndHgt = (*m_DEMImage)(pxi,pyi,0,0);
						float xa = 0.0;
						float ya = 0.0;
						for(std::vector<GradientLinePoint>::iterator iter = Line.begin(); iter != Line.end(); iter++)
						{
							if ((*m_MaxImage)((*iter).xi,(*iter).yi,0,0) == DemProc::EmptyHeightVal)
							{
								(*m_GradImage)((*iter).xi,(*iter).yi,0,0) += StartHgt + ((*iter).t/pt) * (EndHgt-StartHgt);
								img_grad_cnt((*iter).xi,(*iter).yi,0,0) += 1;
								img_grad_d((*iter).xi,(*iter).yi,0,0) = 255;
								if ((xa != 0.0) || (ya != 0.0))
									std::fprintf(stdout, "GradientLineUp(%f,%f,%f,%f)\n", xa, ya, (*iter).x, (*iter).y);
								xa = (*iter).x;
								ya = (*iter).y;
							}
						}
					}
				}
			}
		}
		if (m_Disp)
			img_grad_d.display(*m_Display1);
	}

	cimg_forXY(*m_GradImage,x,y)
	{
		if (img_grad_cnt(x,y,0,0) != DemProc::EmptyHeightVal)
		{
			(*m_GradImage)(x,y,0,0) /= img_grad_cnt(x,y,0,0);
		}				
	}

	std::fprintf(stderr,"\n---Finished!\n");

	return 0;
}


int DemProc::DEMPostProc()
{
	std::fprintf(stderr,"---Postprocessing DEM...\n");

	DirSearchFill(m_DEMImage, m_GradImage, 120);

	m_DEMImage->blur(4.8);

	CImg<float> img_err_o(*m_GradImage);
	CImg<float> img_err(*m_GradImage);

	cimg_forXY(img_err_o,px,py)
	{
		if ((*m_GradImage)(px,py,0,0) != DemProc::EmptyHeightVal)
			img_err_o(px,py,0,0) = (*m_GradImage)(px,py,0,0)-(*m_DEMImage)(px,py,0,0);
	}

	float Lev = 4.2;
	float Lev2 = 1.0;
	ImgLaplace2(&img_err, &img_err_o, Lev, 100, 4);

	img_err.blur(4.2);

	(*m_DEMImage) += img_err*0.95;
	Lev2 = 1.0;

	for (int j = 0; j < 10; j++)
	{
		float Lev = ((3.8 - 3.2*(float)j/10.0)/2.0)/Lev2;

		CImg<float> img_err2_o(*m_GradImage);
		CImg<float> img_err2(*m_GradImage);

		cimg_forXY(img_err2_o,px,py)
		{
			if ((*m_GradImage)(px,py,0,0) != DemProc::EmptyHeightVal)
				img_err2_o(px,py,0,0) = (*m_GradImage)(px,py,0,0)-(*m_DEMImage)(px,py,0,0);
		}

		Lev2 = Lev;
		ImgLaplace2(&img_err2, &img_err2_o, Lev2, 100, 4);

		img_err2.blur(Lev/Lev2);

		(*m_DEMImage) += img_err2*0.85;
	}

	float MaxDiff = 0.0;
	float MaxDiff2 = 0.0;

  m_ErrImage = new CImg<float>(*m_VectorImage, false);

	cimg_forXY(*m_DEMImage,x,y)
	{
		if ((*m_GradImage)(x,y,0,0) != DemProc::EmptyHeightVal)
		{
			float Err = cimg::abs((*m_VectorImage)(x,y,0,0)-(*m_DEMImage)(x,y,0,0));
			if ((*m_VectorImage)(x,y,0,0) > DemProc::EmptyHeightVal+10.0)
			{
				(*m_ErrImage)(x,y,0,0) = Err;
				MaxDiff = cimg::max(MaxDiff, Err);
			}
			else
				(*m_ErrImage)(x,y,0,0) = 0.0;

			MaxDiff2 = cimg::max(MaxDiff2, cimg::abs((*m_GradImage)(x,y,0,0)-(*m_DEMImage)(x,y,0,0)));
		}
		else
			(*m_ErrImage)(x,y,0,0) = 0.0;
	}

	std::fprintf(stderr,"---Finished (max error: %f/%f)!\n", MaxDiff, MaxDiff2);

	return 0;
}


int DemProc::SaveDEM(const bool ClipZero, const int SeaContWater)
{
	std::fprintf(stderr,"Writing results...\n");

	cimg_forXY(*m_SeaMaskImage,x,y)
	{
		if (IsInImage(x, y, 2))
		{
			// all empty pixels next to a sea masked pixel 
			if (((*m_VectorImage)(x,y) <= DemProc::EmptyHeightVal) || ((std::abs((*m_VectorImage)(x,y,0,0)) < 0.05) && (SeaContWater > 0)) )
			if ((*m_SeaMaskImage)(x,y) == 0)
			{
				if ((*m_SeaMaskImage)(x+1,y) == 1) (*m_SeaMaskImage)(x,y) = 2;
				else if ((*m_SeaMaskImage)(x-1,y) == 1) (*m_SeaMaskImage)(x,y) = 2;
				else if ((*m_SeaMaskImage)(x,y+1) == 1) (*m_SeaMaskImage)(x,y) = 2;
				else if ((*m_SeaMaskImage)(x,y-1) == 1) (*m_SeaMaskImage)(x,y) = 2;
			}
			if (SeaContWater > 1)
			{
				if (std::abs((*m_VectorImage)(x,y,0,0)) < 0.005)
				{
					if ((*m_SeaMaskImage)(x+1,y) == 1) (*m_SeaMaskImage)(x,y) = 2;
					else if ((*m_SeaMaskImage)(x+2,y) == 1) (*m_SeaMaskImage)(x,y) = 2;
					else if ((*m_SeaMaskImage)(x-1,y) == 1) (*m_SeaMaskImage)(x,y) = 2;
					else if ((*m_SeaMaskImage)(x-2,y) == 1) (*m_SeaMaskImage)(x,y) = 2;
					else if ((*m_SeaMaskImage)(x,y+1) == 1) (*m_SeaMaskImage)(x,y) = 2;
					else if ((*m_SeaMaskImage)(x,y+2) == 1) (*m_SeaMaskImage)(x,y) = 2;
					else if ((*m_SeaMaskImage)(x,y-1) == 1) (*m_SeaMaskImage)(x,y) = 2;
					else if ((*m_SeaMaskImage)(x,y-2) == 1) (*m_SeaMaskImage)(x,y) = 2;
				}
			}
		}
	}

	cimg_forXY(*m_VectorImage,x,y)
	{
		if ((*m_VectorImage)(x,y,0,0) <= DemProc::EmptyHeightVal) // set empty points in contour map to -1000 (for vis)
		{
			(*m_VectorImage)(x,y,0,0) = -1000.0;
		}
		if ((*m_SeaMaskImage)(x,y) == 1) // mask sea
		{
			(*m_VectorImage)(x,y,0,0) = m_SeaLevelOut;
			(*m_DEMImage)(x,y,0,0) = m_SeaLevelOut;
		}
		else if ((*m_SeaMaskImage)(x,y) == 2) // mask sea (enlarged) with difference threshold
		{
			if ((*m_DEMImage)(x,y) < m_SeaLevelOut+0.05)
			{
				(*m_VectorImage)(x,y,0,0) = m_SeaLevelOut;
				(*m_DEMImage)(x,y,0,0) = m_SeaLevelOut;
			}
		}
		else if ((*m_DEMImage)(x,y) < m_SeaLevelOut-0.25) // rule without mask
		{
			(*m_VectorImage)(x,y,0,0) = m_SeaLevelOut;
			(*m_DEMImage)(x,y,0,0) = m_SeaLevelOut;
		}
		else if ((*m_DEMImage)(x,y) == m_SeaLevelOut)
		{
			(*m_DEMImage)(x,y,0,0) = m_SeaLevelOut+0.001;
		}

		if (ClipZero)
		{
			if ((*m_DEMImage)(x,y) < m_SeaLevelOut)
			{
				(*m_DEMImage)(x,y,0,0) = m_SeaLevelOut+0.001;
			}
		}
	}

	if (m_VectorImage != NULL)
		m_VectorImage->get_normalize(0.0, 255.0).save("test_cont.png");
	if (m_SeaMaskImage != NULL)
		m_SeaMaskImage->get_normalize(0, 255).save("test_sea.png");
	if (m_DEMImage != NULL)
		m_DEMImage->get_normalize(0.0, 255.0).save("test_dem.png");
	if (m_DEMImageS != NULL)
		m_DEMImageS->get_normalize(0.0, 255.0).save("test_dems.png");
	if (m_MaxImage != NULL)
		m_MaxImage->get_normalize(0.0, 255.0).save("test_max.png");
	if (m_MaxLineImage != NULL)
		m_MaxLineImage->save("test_maxd.png");
	if (m_GradImage != NULL)
		m_GradImage->get_normalize(0.0, 255.0).save("test_grad.png");
	if (m_SmallWaterImage != NULL)
		m_SmallWaterImage->get_normalize(0, 255).save("test_sw.png");
	if (m_WaterImage != NULL)
		m_WaterImage->get_normalize(0, 255).save("test_w.png");
	if (m_ErrImage != NULL)
		m_ErrImage->get_normalize(0.0, 255.0).save("test_err.png");

	GDALRasterBand *poBand;
	float *pRaster;

	// if no processing is done
	if (m_ErrImage == NULL)
		pRaster = m_VectorImage->ptr();
	else
		pRaster = m_DEMImage->ptr();

	poBand = m_poDstDS->GetRasterBand(1);
	poBand->RasterIO( GF_Write, 0, 0, m_DEMWidth, m_DEMHeight, 
										pRaster, m_DEMWidth, m_DEMHeight, GDT_Float32, 0, 0 ); 

	GDALRasterBand *poBand2;
	float *pRaster2 = m_VectorImage->ptr();

	poBand2 = m_poDstDS_test->GetRasterBand(1);
	poBand2->RasterIO( GF_Write, 0, 0, m_DEMWidth, m_DEMHeight, 
										pRaster2, m_DEMWidth, m_DEMHeight, GDT_Float32, 0, 0 ); 

	std::fprintf(stderr,"Results written to file %s\n", m_OutputFile.c_str());

	return 0;
}
