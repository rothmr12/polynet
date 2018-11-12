/*

	PolyNet: A mesoscale model of a polymer network.
	This is based on the original EPnet algorithm of David E. Hanson, with some
	modifications.
													John L. Barber,
													jlbarber@lanl.gov
													505-664-0605

*/

#include "PolyNet.h"
#include "ReadInputDeck.h"

/*

	Appends current time, lambda, stress tensor, and fraction of chains ruptured to StrainStressFile. 

*/

bool Network::WriteStrainStressData(bool FirstTimeFlag){

	if(FirstTimeFlag){
		StrainStressFileStream << "<r2 / N b_m> = " << r2overNbm << endl;
		double DomainVolume = (L[0][1] - L[0][0]) * (L[1][1] - L[1][0]) * (L[2][1] - L[2][0]);
		StrainStressFileStream << "Chain number density = " << NumChains / DomainVolume << endl;
		StrainStressFileStream << endl;
		StrainStressFileStream << "Time\t Strain\t sigma_{t,xx}\t sigma_{t,xy}\t sigma_{t,xz}\t sigma_{t,yx}\t sigma_{t,yy}\t sigma_{t,yz}\t sigma_{t,zx}\t sigma_{t,zy}\t sigma_{t,zz}\t Frac. Ruptured" << endl;
	}

	ComputeStressTensor();

	StrainStressFileStream << Time << "\t";
	StrainStressFileStream << Lambda[2] - 1 << "\t";
	for(int mu = 0; mu < 3; ++mu) for(int nu = 0; nu < 3; ++nu)
		StrainStressFileStream << Stress[mu][nu] << "\t";
	StrainStressFileStream << (1.*NumRupturedChains) / NumChains << endl;

	return true;
}





/*

	Writes an image of the network, including nodes and chains.

*/

bool Network::WriteNetworkImage(string ImageFileName){

	clock_t c_begin = clock();

	cout << "Writing network image " << ImageFileName << "." << ImageFormat << " ";
	cout.flush();

	int n = NumberOfPixels;
	if(n <= 0){
		cout << "Error in WriteNetworkImage: n <= 0." << endl;
		cout.flush();
		return false;
	}


	// Calculate the camera coordinate axes xhat, yhat, and zhat
	zhat = -ViewDirection;

		// Make Up perpendicular to ViewDirection, make Up a unit vector, and set yhat = Up
	Up -= Dot(Up, ViewDirection)*ViewDirection;
	if(Up.Norm2() <= 0){
		cout << "Error in WriteNetworkImage: " << Up << " is not a valid vertical direction." << endl;
		cout.flush();
		return false;
	}
	Up.Normalize();
	yhat = Up;

		// Set xhat = yhat cross zhat
	xhat = Cross3D(yhat, zhat);


	// Determine actual pixel array dimensions, as well as angular resolution in both directions
	int nx, ny;
	if(AngularSpanX >= AngularSpanY){
		nx = n;
		DTan = 2*tan(AngularSpanX/2) / nx;
		ny = int(ceil(nx*tan(AngularSpanY/2)/tan(AngularSpanX/2)));
	}
	else{
		ny = n;
		DTan = 2*tan(AngularSpanY/2) / ny;
		nx = int(ceil(ny*tan(AngularSpanX/2)/tan(AngularSpanY/2)));
	}


	// Allocate space for image array and depth array
	unsigned char* ImageData = new(nothrow) unsigned char[3*nx*ny];
	if(!ImageData){
		cout << "Error in WriteNetworkImage: Unable to allocate space to store image data." << endl;
		cout.flush();
		return false;
	}

	double* Depth = new(nothrow) double[nx*ny];
	if(!Depth){
		cout << "Error in WriteNetworkImage: Unable to allocate space to store the depth array." << endl;
		cout.flush();
		return false;
	}


	// Add background color and initialize depth array

	int i, j, k;
	for(i = 0; i < nx; ++i) for(j = 0; j < ny; ++j){
		ImageData[3*(i + nx*(ny-1-j)) + 0] = BackgroundColor[0];
		ImageData[3*(i + nx*(ny-1-j)) + 1] = BackgroundColor[1];
		ImageData[3*(i + nx*(ny-1-j)) + 2] = BackgroundColor[2];

		Depth[j + ny*i] = +DBL_MAX;
	}


	// Add chains
	double ExtensionDevs = 5.0;
	double MaxExtension, Extension, AverageExtension;
	double ChainLength, NumKuhnLengths, ColorScale;
	unsigned char ChainColor[3];

	int jmin, kmin;
	double s[3][2], smin, W[] = {L[0][1]-L[0][0], L[1][1]-L[1][0], L[2][1]-L[2][0]};
	bool BoundaryCrossFlag, StillSomeLeft;
	Point r0, r1, r2, r;

	if(ChainRadius > 0) for(auto & ElementC : ChainList){
		if(ElementC.tag != 0) continue;
		r0 = NodeList[ElementC.Node0].r;
		r1 = NodeList[ElementC.Node1].r;
		r = r1 - r0;

		BoundaryCrossFlag = false;
		if(PeriodicFlag) for(j = 0; j < 3; ++j){
			s[j][0] = s[j][1] = DBL_MAX;
			if(r.r[j] >= +0.5*W[j]){
				BoundaryCrossFlag = true;
				r.r[j] -= W[j];
				s[j][0] = (L[j][0] - r0.r[j])/r.r[j];
			}
			else if(r.r[j] < -0.5*W[j]){
				BoundaryCrossFlag = 1;
				r.r[j] += W[j];
				s[j][1] = (L[j][1] - r0.r[j])/r.r[j];
			}
			r1.r[j] = r0.r[j] + r.r[j];
		}

		ChainLength = r.Norm();
		NumKuhnLengths = MonomerLength*ChainList[i].N / KuhnLength;
		AverageExtension = sqrt( 8 / (3*PI*NumKuhnLengths) );
		MaxExtension = 1.0;
		Extension = ChainLength / (ChainList[i].N*MonomerLength);
		ColorScale = (Extension - AverageExtension) / (MaxExtension - AverageExtension);
		if(ColorScale <= 0) ColorScale = 0;
		if(ColorScale >= 1) ColorScale = 1;
		ChainColor[0] = (unsigned char)( LooseChainColor[0] + ColorScale*(TightChainColor[0] - LooseChainColor[0]) );
		ChainColor[1] = (unsigned char)( LooseChainColor[1] + ColorScale*(TightChainColor[1] - LooseChainColor[1]) );
		ChainColor[2] = (unsigned char)( LooseChainColor[2] + ColorScale*(TightChainColor[2] - LooseChainColor[2]) );

		if(!PeriodicFlag || !BoundaryCrossFlag){
			if(!RenderCylinder(ImageData, Depth, nx, ny, r0, r1, ChainRadius, ChainColor, PeriodicFlag))
				return false;
		}
		else{
			while(1){ // Make cylinder segments between r0 and all successive boundary crossings
				StillSomeLeft = false;
				smin = DBL_MAX;
				for(j = 0; j < 3; ++j) for(k = 0; k < 2; ++k){ // Find remaining boundary crossing nearest to r0
					if(s[j][k] <= 0 || s[j][k] >= 1) continue;
					StillSomeLeft = true;
					if(s[j][k] >= smin) continue;
					smin = s[j][k];
					jmin = j;
					kmin = k;
				}

				if(!StillSomeLeft) break;

				r2 = r0 + smin*r; 			// Make r2 = point of bounday crossing
				if(!RenderCylinder(ImageData, Depth, nx, ny, r0, r2, ChainRadius, ChainColor, PeriodicFlag)) // Make cylinder between r0 and r2
					return false;

				r0 = r2;					// Reset r0 to be = to r2
				if(kmin == 0){			// Subtract off domain width in dimension corresponding to this boundary crossing
					r0.r[jmin] += W[jmin];
					r1.r[jmin] += W[jmin];
				}
				else{
					r0.r[jmin] -= W[jmin];
					r1.r[jmin] -= W[jmin];
				}

				for(j = 0; j < 3; ++j) for(k = 0; k < 2; ++k) // Subtract smin from all boundary crossing distances
					if(s[j][k] < 1) s[j][k] -= smin;
			}

			if(!RenderCylinder(ImageData, Depth, nx, ny, r0, r1, ChainRadius, ChainColor, PeriodicFlag)) // Make final cylinder between r0 (last boundary crossing) and r1
				return false;
		}
	}


	// Add nodes

	unsigned char NodeColor[3];
	if(NodeRadius > 0) for(auto & ElementN : NodeList){
		if(ElementN.tag == 0){
			NodeColor[0] = Node0Color[0];
			NodeColor[1] = Node0Color[1];
			NodeColor[2] = Node0Color[2];
		}
		else{
			NodeColor[0] = Node1Color[0];
			NodeColor[1] = Node1Color[1];
			NodeColor[2] = Node1Color[2];
		}

//		if(tag != 0) continue;

		if(!RenderSphere(ImageData, Depth, nx, ny, ElementN.r, NodeRadius, NodeColor, PeriodicFlag))
			return false;
	}


	// Add border
	if(BorderRadius > 0) for(i = -1; i <= +1; ++i) for(j = -1; j <= +1; ++j) for(k = -1; k <= +1; ++k){
		if(int(fabs(i) + fabs(j) + fabs(k)) != 2) continue;

		r0.r[0] = ( i == 0 ? L[0][0] - BorderRadius : (L[0][1]+L[0][0] + i*W[0])/2 );
		r0.r[1] = ( j == 0 ? L[1][0] - BorderRadius : (L[1][1]+L[1][0] + j*W[1])/2 );
		r0.r[2] = ( k == 0 ? L[2][0] - BorderRadius : (L[2][1]+L[2][0] + k*W[2])/2 );
		r1.r[0] = ( i == 0 ? L[0][1] - BorderRadius : (L[0][1]+L[0][0] + i*W[0])/2 );
		r1.r[1] = ( j == 0 ? L[1][1] - BorderRadius : (L[1][1]+L[1][0] + j*W[1])/2 );
		r1.r[2] = ( k == 0 ? L[2][1] - BorderRadius : (L[2][1]+L[2][0] + k*W[2])/2 );

		if(!RenderCylinder(ImageData, Depth, nx, ny, r0, r1, BorderRadius, BorderColor, false))
			return false;
	}

	delete [] Depth;


	// Create PPM file
	string ImageFileNameTemp = ImageFileName + ".ppm";
	FILE *Out;
	Out = fopen(ImageFileNameTemp.c_str(), "w");
	if(!Out){
		cout << "Error in WriteNetworkImage: Cannot write image to " << ImageFileNameTemp << endl;
		cout.flush();
		return false;
	}

	fprintf(Out, "P6\n");
	fprintf(Out, "%d %d\n", nx, ny);
	fprintf(Out, "%d\n", 255);    
	fwrite(ImageData, sizeof(unsigned char), 3*nx*ny, Out);
	fprintf(Out, "\n");
	fclose(Out);

	delete [] ImageData;

	// If called for, use system commands to convert ppm file to desired
	// image format, and then erase ppm file
	if(!CapitalCompare(ImageFormat, "ppm"))
		Converter(ImageFileName, "ppm", ImageFormat);

	cout << ((double)(clock() - c_begin)) / CLOCKS_PER_SEC << " sec" << endl;

	return true;
}





bool Network::RenderSphere(unsigned char* ImageArray, double* Depth, int nx, int ny, const Point& r0, double R, unsigned char* Color, bool PeriodicFlag){

	// Check arguments
	if(nx <= 0 || ny <= 0){
		cout << "Error in RenderSphere: " << nx << " x " << ny << " is not a valid set of image dimensions" << endl;
		cout.flush();
		return false;
	}

	if(R <= 0) return true; // Sphere has no radius. Don't draw.


	// Determine rough position of sphere w.r.t. camera,
	// distance from camera to nearest part of sphere.
	Point r = r0 - ViewPoint;		// Vector from ViewPoint to center of sphere
	double rmag = r.Norm();
	double rDotVD = Dot(r, ViewDirection);
	double NearestDistance;

	if(rmag >= R){ 		// Viewpoint is outside of sphere.
		if(rDotVD < 0)	// Viewpoint is in front of sphere: Sphere is not seen.
			return true;
		NearestDistance = rmag - R;
	}
	else{				// Viewpoint is inside of sphere
		if(rDotVD >= 0)	// Viewpoint is behind center of sphere
			NearestDistance = R + rmag;
		else			// Viewpoint is in front of center of sphere
			NearestDistance = R - rmag;
	}


	// Build coordinate system for sphere surface
	Point zhat1 = -Normalize(r);
	Point yhat1 = Normalize( yhat - Dot(yhat, zhat1)*zhat1 );
	Point xhat1 = Cross3D(yhat1, zhat1);


	// Render points on surface of sphere
	unsigned char FinalColor[3];
	double Dr, rad, DPhi, Phi, x, y, z, TanX, TanY, ThisPointsDepth;
	int i, j, k, l, i1, j1, NPhi, Nr;
	double TanXMax = +tan(AngularSpanX/2);
	double TanXMin = -tan(AngularSpanX/2);
	double TanYMax = +tan(AngularSpanY/2);
	double TanYMin = -tan(AngularSpanY/2);
	Point Normal;

	Dr = 0.5*DTan*NearestDistance;
	Nr = int(ceil(R/Dr));
	Dr = R/Nr;

	for(i = 0; i < Nr; ++i){
		rad = (i+0.5)*Dr;
		DPhi = Dr/rad;
		NPhi = int(ceil(2*PI/DPhi));
		DPhi = 2*PI/NPhi;
		for(j = 0; j < NPhi; ++j){
			Phi = j*DPhi;
			x = rad*cos(Phi);
			y = rad*sin(Phi);
			for(k = -1; k <= +1; k += 2){	// -1 = back of sphere, +1 = front of sphere.
				z = k*sqrt(R*R - rad*rad);
				for(l = 0; l < 3; ++l)
					r.r[l] = r0.r[l] + x*xhat1.r[l] + y*yhat1.r[l] + z*zhat1.r[l];

//				if(PeriodicFlag){
//					if(r.r[0] < L[0][0] || r[0] > L[0][1]) continue;
//					if(r.r[1] < L[1][0] || r[1] > L[1][1]) continue;
//					if(r.r[2] < L[2][0] || r[2] > L[2][1]) continue;
//				}

				if(!FindImageTangents(r, TanX, TanY)) continue;
				if(TanX < TanXMin || TanX >= TanXMax) continue;
				if(TanY < TanYMin || TanY >= TanYMax) continue;
				i1 = int( (TanX - TanXMin)/DTan );
				j1 = int( (TanY - TanYMin)/DTan );

				r = r - ViewPoint;

				ThisPointsDepth = r.Norm();
				if(ThisPointsDepth >= Depth[j1 + ny*i1]) continue;
				Depth[j1 + ny*i1] = ThisPointsDepth;

				for(l = 0; l < 3; ++l)
					Normal.r[l] = (x*xhat1.r[l] + y*yhat1.r[l] + z*zhat1.r[l]) / R;
				IlluminationModel(Normal, Color, FinalColor);

				for(l = 0; l < 3; ++l)
					ImageArray[3*(i1 + nx*(ny-1-j1)) + l] = FinalColor[l];
			}
		}
	}

	return true;
}





bool Network::RenderCylinder(unsigned char* ImageArray, double* Depth, int nx, int ny, const Point& r0, const Point& r1, double R, unsigned char* Color, bool PeriodicFlag){

	// Check arguments
	if(nx <= 0 || ny <= 0){
		cout << "Error in RenderCylinder: " << nx << " x " << ny << " is not a valid set of image dimensions" << endl;
		cout.flush();
		return false;
	}

	if(R <= 0)		// Cylinder has no radius. Don't draw.
		return true;


	// Check if this cylinder is within the "viewing frustum"
	double TanXMax = +tan(AngularSpanX/2);
	double TanXMin = -tan(AngularSpanX/2);
	double TanYMax = +tan(AngularSpanY/2);
	double TanYMin = -tan(AngularSpanY/2);
	double TanX, TanY, TanX1, TanY1;

	bool BehindFlag0 = !FindImageTangents(r0, TanX, TanY);
	bool BehindFlag1 = !FindImageTangents(r1, TanX1, TanY1);

	if(BehindFlag0 && BehindFlag1)
		return true;

	else if(!BehindFlag0 && !BehindFlag1){
		if(TanX < TanXMin && TanX1 < TanXMin)
			return true;

		if(TanX > TanXMax && TanX1 > TanXMax)
			return true;

		if(TanY < TanYMin && TanY1 < TanYMin)
			return true;

		if(TanY > TanYMax && TanY1 > TanYMax)
			return true;
	}


	// Build coordinate system for cylinder surface
	double Length = Norm(r0 - r1);
	if(Length <= 0)			// Cylinder has zero length: Don't draw
		return true;

	double Z, r;
	Point xhat1, yhat1, zhat1, c;

	zhat1 = (r0 - r1) / Length;
	c = 0.5*(r1 + r0) - ViewPoint;	// Center of cylinder wrt ViewPoint
	Z = -Dot(c, zhat1); 			// "z" of ViewPoint in cylinder's cylindrical coords
	r = Norm(c + Z*zhat1);			// "r" of ViewPoint in cylinder's cylindrical coords

	if(r == 0)	// ViewPoint is on axis
		yhat1 = Up - Dot(Up, zhat1)*zhat1;
	else		// ViewPoint is off axis
		yhat1 = c - Dot(c, zhat1)*zhat1;
	yhat1.Normalize();
	xhat1 = Cross3D(yhat1, zhat1);


	// Determine distance from camera to nearest part of cylinder,
	// and z and phi resolutions on cylinder surface
	double rmin, zmin, NearestDistance, Dz, DPhi;
	int Nz, NPhi;
	int LargestDimension = (nx >= ny ? nx : ny);

	rmin = (r <= R ? R - r : r - R);
	zmin = (fabs(Z) <= Length/2 ? 0 : fabs(Z) - Length/2);
	NearestDistance = sqrt(rmin*rmin + zmin*zmin);

	Dz = 0.5*DTan*NearestDistance;
	Nz = int(ceil(Length/Dz));
	Dz = Length/Nz;

	DPhi = 0.5*DTan*NearestDistance/R;
	NPhi = int(ceil(2*PI/DPhi));
	DPhi = 2*PI/NPhi;


	// Render points on surface of cylinder
	unsigned char FinalColor[3];
	double Phi, x, y, z, ThisPointsDepth;
	int i, j, k, i1, j1;
	Point Position, Normal;

	for(i = 0; i < Nz; ++i){
		z = -Length/2 + (i+0.5)*Dz;
		for(j = 0; j < NPhi; ++j){
			Phi = j*DPhi;
			x = R*cos(Phi);
			y = R*sin(Phi);

			for(k = 0; k < 3; ++k)
				Position.r[k] = 0.5*(r0.r[k] + r1.r[k]) + x*xhat1.r[k] + y*yhat1.r[k] + z*zhat1.r[k];

//			if(PeriodicFlag){
//				if(Position.r[0] < L[0][0] || Position[0] > L[0][1]) continue;
//				if(Position.r[1] < L[1][0] || Position[1] > L[1][1]) continue;
//				if(Position.r[2] < L[2][0] || Position[2] > L[2][1]) continue;
//			}

			if(!FindImageTangents(Position, TanX, TanY)) continue;
			if(TanX < TanXMin || TanX >= TanXMax) continue;
			if(TanY < TanYMin || TanY >= TanYMax) continue;
			i1 = int( (TanX - TanXMin)/DTan );
			j1 = int( (TanY - TanYMin)/DTan );

			Position -= ViewPoint;

			ThisPointsDepth = Position.Norm();
			if(ThisPointsDepth >= Depth[j1 + ny*i1]) continue;
			Depth[j1 + ny*i1] = ThisPointsDepth;

			Normal = (x/R)*xhat1 + (y/R)*yhat1;
			IlluminationModel(Normal, Color, FinalColor);

			for(k = 0; k < 3; ++k)
				ImageArray[3*(i1 + nx*(ny-1-j1)) + k] = FinalColor[k];
		}
	}

	return true;
}





/*

	Given a normal vector (Normal) and color (Color) of a given surface, uses the known
	LightDirection via an illumination model to generate an apparent color (FinalColor) for
	the surface at that point.

*/

void Network::IlluminationModel(const Point& Normal, unsigned char* Color, unsigned char* FinalColor){

	int MaxInt = 255;
	double CosAngle = -Dot(LightDirection, Normal) / sqrt(LightDirection.Norm2() * Norm2(Normal));
	double f, C0, C1, C2; 

	C0 = 0.5*Color[0];
	C1 = 0.5*Color[1];
	C2 = 0.5*Color[2];

	if(CosAngle > 0 && CosAngle <= 0.5){
		f = 2*CosAngle;
		C0 += 0.5*f*Color[0];
		C1 += 0.5*f*Color[1];
		C2 += 0.5*f*Color[2];
	}
	else if(CosAngle > 0.5){
		f = 2*(CosAngle- 0.5);
		C0 += 0.5*( f*MaxInt + (1-f)*Color[0] );
		C1 += 0.5*( f*MaxInt + (1-f)*Color[1] );
		C2 += 0.5*( f*MaxInt + (1-f)*Color[2] );
	}

	if(C0 < 0)      C0 = 0;
	if(C0 > MaxInt) C0 = MaxInt;
	if(C1 < 0)      C1 = 0;
	if(C1 > MaxInt) C1 = MaxInt;
	if(C2 < 0)      C2 = 0;
	if(C2 > MaxInt) C2 = MaxInt;

	FinalColor[0] = (unsigned char)(C0);
	FinalColor[1] = (unsigned char)(C1);
	FinalColor[2] = (unsigned char)(C2);

	return;
}





/*

	Finds the tangents of the horizontal (TanX) and vertical (TanY) of the angular
	position of the point SurfacePoint with respect to the imaging coordinate system.

*/

bool Network::FindImageTangents(const Point& SurfacePoint, double& TanX, double& TanY){

	Point RelPoint = SurfacePoint - ViewPoint;
	double z = Dot(RelPoint, zhat);
	if(z >= 0)	// Then point is behind "camera", return 1
		return false;

	TanX = Dot(RelPoint, xhat)/(-z);
	TanY = Dot(RelPoint, yhat)/(-z);

	return true;
}
