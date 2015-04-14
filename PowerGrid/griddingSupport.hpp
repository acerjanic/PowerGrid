//
//  griddingSupport.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/9/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//


#ifndef PowerGrid_griddingSupport_hpp
#define PowerGrid_griddingSupport_hpp

//using namespace std::literals; //where to put?


//__global__ static void
//Deinterleave_data2d_kernel(
//                           cufftComplex *src, float *outR_d, float *outI_d,
//                           int imageX, int imageY)
//{
//    int Y = blockIdx.x;
//    int X = blockIdx.y;
//
//    int lIndex = X + Y*imageX;
//    outR_d[lIndex] = src[lIndex] REAL;
//    outI_d[lIndex] = src[lIndex] IMAG;
//}


//void
//deinterleave_data2d(
//                    cufftComplex *src, float *outR_d, float *outI_d,
//                    int imageX, int imageY)
//{
//    dim3 threads(1,1);
//    dim3 blocks(imageY,imageX);
//
//    Deinterleave_data2d_kernel<<<blocks,threads>>>
//    (src,outR_d,outI_d,imageX,imageY);
//}

template <typename T1>
void
deinterleave_data2d( ///NAIVE
                    std::complex<T1> *src, T1 *outR_d, T1 *outI_d,
                    int imageX, int imageY)
{
    int lIndex;
    for(Y=0; Y<imageY; Y++)
    {
    	for(X=0; X<imageX; X++)
    	{
    		lIndex = X+Y*imageX;
    		outR_d[lIndex] = src[lIndex].real();
    		outI_d[lIndex] = src[lIndex].imag();
    	}

    }
    
}

//__global__ static void
//Deinterleave_data3d_kernel(
//                           cufftComplex *src, float *outR_d, float *outI_d,
//                           int imageX, int imageY, int imageZ)
//{
//    int Z;
//    int Y;
//    int X;
//
//    int lIndex = Z + X*imageZ + Y*imageZ*imageX;
//    outR_d[lIndex] = src[lIndex] REAL;
//    outI_d[lIndex] = src[lIndex] IMAG;
//}


template <typename T1>
void
deinterleave_data3d(
					std::complex<T1> *src, T1 *outR_d, T1 *outI_d,
                    int imageX, int imageY, int imageZ)
{
    int lIndex,X,Y,Z;
    for(Y=0; Y<imageY; Y++)
    {
    	for(X=0; X<imageX; X++)
		{
			for(Z=0; Z<imageZ; Z++)
			{
				lIndex = Z + X*imageZ + Y*imageZ*imageX;
				outR_d[lIndex] = src[lIndex].real();
				outI_d[lIndex] = src[lIndex].imag();
			}

		}
    }
}

//
//__global__ static void
//Deapodization2d_kernel(
//                       cufftComplex *dst, cufftComplex *src,
//                       int imageX, int imageY,
//                       float kernelWidth, float beta, float gridOS)
//{
//    /*Justin's gridding code:
//     [kernelX kernelY kernelZ] =meshgrid([-Nx/2:Nx/2-1]/Nx,
//     [-Ny/2:Ny/2-1]/Ny,
//     [-Nz/2:Nz/2-1]/Nz);
//     gridKernel = (sin(sqrt(pi^2*kernelWidth^2*kernelX.^2 - beta^2))./ ...
//     sqrt(pi^2*kernelWidth^2*kernelX.^2 - beta^2)).*
//     (sin(sqrt(pi^2*kernelWidth^2*kernelY.^2 - beta^2))./
//     sqrt(pi^2*kernelWidth^2*kernelY.^2 - beta^2)).*
//     (sin(sqrt(pi^2*kernelWidth^2*kernelZ.^2 - beta^2))./
//     sqrt(pi^2*kernelWidth^2*kernelZ.^2 - beta^2));
//     */
//    int imageNumElems = imageX * imageY;
//
//    int Y = blockIdx.x;
//    int X = blockIdx.y;
//
//    float gridKernelY = float(Y - (imageY/2)) / (float)imageY;
//    float gridKernelX = float(X - (imageX/2)) / (float)imageX;
//
//    float common_exprX = (PI*PI*kernelWidth*kernelWidth*gridKernelX*gridKernelX - beta*beta);
//    float common_exprY = (PI*PI*kernelWidth*kernelWidth*gridKernelY*gridKernelY - beta*beta);
//
//    float common_exprX1;
//    float common_exprY1;
//
//    if(common_exprX>=0)
//        common_exprX1 = (sin(sqrt(common_exprX))/sqrt(common_exprX));
//    else
//        common_exprX1 = (sinh(sqrt(-1.0f*common_exprX))/sqrt(-1.0f*common_exprX));
//
//    if(common_exprY>=0)
//        common_exprY1 = (sin(sqrt(common_exprY))/sqrt(common_exprY));
//    else
//        common_exprY1 = (sinh(sqrt(-1.0f*common_exprY))/sqrt(-1.0f*common_exprY));
//
//    float gridKernel =  common_exprX1 * common_exprY1;
//
//    if(gridKernel==gridKernel)
//    {
//        int common_index = X + Y*imageX;
//        float gridOS2 = gridOS * gridOS;
//        //dst[common_index] REAL = ((float)imageNumElems * src[common_index]REAL) / gridKernel;
//        //dst[common_index] IMAG = ((float)imageNumElems * src[common_index]IMAG) / gridKernel;
//        dst[common_index] REAL = (src[common_index]REAL) / gridKernel * (1.0f / gridOS2);
//        dst[common_index] IMAG = (src[common_index]IMAG) / gridKernel * (1.0f / gridOS2);
//    }
//}

template <typename T1>
void
deapodization2d(
				std::complex<T1> *dst,std::complex<T1> *src,
                int imageX, int imageY,
                float kernelWidth, float beta, float gridOS)
{
//    assert( (!(imageX%2) && !(imageY%2)) ); //<< CHECK EVEN
    
    int imageNumElems = imageX * imageY * imageZ;

	int Y;
	int X;

	float gridKernelY;
	float gridKernelX;

	float common_exprX;
	float common_exprY;

	float common_exprX1;
	float common_exprY1;

	float gridKernel;
	int common_index;
	float gridOS2;

	for(Y=0; Y<imageY; Y++)
	{
		for(X=0; X<imageX; X++)
		{
				gridKernelY = float(Y - (imageY/2)) / (float)imageY;
			    gridKernelX = float(X - (imageX/2)) / (float)imageX;

			    common_exprX = (MRI_PI*MRI_PI*kernelWidth*kernelWidth*gridKernelX*gridKernelX - beta*beta);
			    common_exprY = (MRI_PI*MRI_PI*kernelWidth*kernelWidth*gridKernelY*gridKernelY - beta*beta);

			    if(common_exprX>=0)
			        common_exprX1 = (sin(sqrt(common_exprX))/sqrt(common_exprX));
			    else
			        common_exprX1 = (sinh(sqrt(-1.0f*common_exprX))/sqrt(-1.0f*common_exprX));

			    if(common_exprY>=0)
			        common_exprY1 = (sin(sqrt(common_exprY))/sqrt(common_exprY));
			    else
			        common_exprY1 = (sinh(sqrt(-1.0f*common_exprY))/sqrt(-1.0f*common_exprY));

			    gridKernel =  common_exprX1 * common_exprY1;

			    if(gridKernel==gridKernel) //???
			    {
			        common_index = X + Y*imageX;
			        gridOS2 = gridOS * gridOS;
			        dst[common_index] = (src[common_index]) / gridKernel * (1.0f / gridOS2) //???
//			        dst[common_index] = (src[common_index].real()) / gridKernel * (1.0f / gridOS2);
//			        dst[common_index] += ((src[common_index].imag()) / gridKernel * (1.0f / gridOS2))*1(std::complex_literals::if);
			    }
		}
	}
}

//__global__ static void
//Deapodization3d_kernel(
//                       cufftComplex *dst, cufftComplex *src,
//                       int imageX, int imageY, int imageZ,
//                       float kernelWidth, float beta, float gridOS)
//{
//    /*Justin's gridding code:
//     [kernelX kernelY kernelZ] =meshgrid([-Nx/2:Nx/2-1]/Nx,
//     [-Ny/2:Ny/2-1]/Ny,
//     [-Nz/2:Nz/2-1]/Nz);
//     gridKernel = (sin(sqrt(pi^2*kernelWidth^2*kernelX.^2 - beta^2))./ ...
//     sqrt(pi^2*kernelWidth^2*kernelX.^2 - beta^2)).*
//     (sin(sqrt(pi^2*kernelWidth^2*kernelY.^2 - beta^2))./
//     sqrt(pi^2*kernelWidth^2*kernelY.^2 - beta^2)).*
//     (sin(sqrt(pi^2*kernelWidth^2*kernelZ.^2 - beta^2))./
//     sqrt(pi^2*kernelWidth^2*kernelZ.^2 - beta^2));
//     */
//    int imageNumElems = imageX * imageY * imageZ;
//
//    int Z = threadIdx.x;
//    int Y = blockIdx.x;
//    int X = blockIdx.y;
//
//    float gridKernelZ = (float(Z) - ((float)imageZ/2.0f)) / (float)imageZ;
//    float gridKernelY = (float(Y) - ((float)imageY/2.0f)) / (float)imageY;
//    float gridKernelX = (float(X) - ((float)imageX/2.0f)) / (float)imageX;
//
//    float common_exprX = (PI*PI*kernelWidth*kernelWidth*gridKernelX*gridKernelX - beta*beta);
//    float common_exprY = (PI*PI*kernelWidth*kernelWidth*gridKernelY*gridKernelY - beta*beta);
//    float common_exprZ = (PI*PI*kernelWidth*kernelWidth*gridKernelZ*gridKernelZ - beta*beta);
//
//    float common_exprX1;
//    float common_exprY1;
//    float common_exprZ1;
//
//    if(common_exprX>=0)
//        common_exprX1 = (sin(sqrt(common_exprX))/sqrt(common_exprX));
//    else
//        common_exprX1 = (sinh(sqrt(-1.0f*common_exprX))/sqrt(-1.0f*common_exprX));
//
//    if(common_exprY>=0)
//        common_exprY1 = (sin(sqrt(common_exprY))/sqrt(common_exprY));
//    else
//        common_exprY1 = (sinh(sqrt(-1.0f*common_exprY))/sqrt(-1.0f*common_exprY));
//
//    if(common_exprZ>=0)
//        common_exprZ1 = (sin(sqrt(common_exprZ))/sqrt(common_exprZ));
//    else
//        common_exprZ1 = (sinh(sqrt(-1.0f*common_exprZ))/sqrt(-1.0f*common_exprZ));
//
//    float gridKernel =  common_exprX1 * common_exprY1 * common_exprZ1;
//
//    if(gridKernel==gridKernel)
//    {
//        int common_index = Z + X*imageZ + Y*imageZ*imageX;
//        float gridOS3 = gridOS * gridOS * gridOS;
//        //dst[common_index] REAL = ((float)imageNumElems*src[common_index]REAL) / (gridKernel);
//        //dst[common_index] IMAG = ((float)imageNumElems*src[common_index]IMAG) / (gridKernel);
//        dst[common_index] REAL = (src[common_index]REAL) / gridKernel * (1.0f / gridOS3);
//        dst[common_index] IMAG = (src[common_index]IMAG) / gridKernel * (1.0f / gridOS3);
//    }
//}

template <typename T1>
void
deapodization3d(
				std::complex<T1> *dst,std::complex<T1> *src,
                int imageX, int imageY, int imageZ,
                float kernelWidth, float beta, float gridOS)
{
    
    int Z;
	int Y;
	int X;

    int imageNumElems = imageX * imageY * imageZ;

    float gridKernelZ;
	float gridKernelX;
	float gridKernelY;

	float common_exprZ;
	float common_exprY;
	float common_exprX;

	float common_exprX1;
	float common_exprY1;
	float common_exprZ1;

	float gridOS3;
	int common_index;
    for (X=0;X<imageX;X++)
    {
    	for (Y=0;Y<imageY;Y++)
    	{
    		for (Z=0;Z<imageZ;Z++)
    		{
    			gridKernelZ = (float(Z) - ((float)imageZ/2.0f)) / (float)imageZ;
    		    gridKernelY = (float(Y) - ((float)imageY/2.0f)) / (float)imageY;
    		    gridKernelX = (float(X) - ((float)imageX/2.0f)) / (float)imageX;

    		    common_exprX = (MRI_PI*MRI_PI*kernelWidth*kernelWidth*gridKernelX*gridKernelX - beta*beta);
    		    common_exprY = (MRI_PI*MRI_PI*kernelWidth*kernelWidth*gridKernelY*gridKernelY - beta*beta);
    		    common_exprZ = (MRI_PI*MRI_PI*kernelWidth*kernelWidth*gridKernelZ*gridKernelZ - beta*beta);

    		    if(common_exprX>=0)
    		        common_exprX1 = (sin(sqrt(common_exprX))/sqrt(common_exprX));
    		    else
    		        common_exprX1 = (sinh(sqrt(-1.0f*common_exprX))/sqrt(-1.0f*common_exprX));

    		    if(common_exprY>=0)
    		        common_exprY1 = (sin(sqrt(common_exprY))/sqrt(common_exprY));
    		    else
    		        common_exprY1 = (sinh(sqrt(-1.0f*common_exprY))/sqrt(-1.0f*common_exprY));

    		    if(common_exprZ>=0)
    		        common_exprZ1 = (sin(sqrt(common_exprZ))/sqrt(common_exprZ));
    		    else
    		        common_exprZ1 = (sinh(sqrt(-1.0f*common_exprZ))/sqrt(-1.0f*common_exprZ));

    		    float gridKernel =  common_exprX1 * common_exprY1 * common_exprZ1;

    		    if(gridKernel==gridKernel) //?
    		    {
    		        int common_index = Z + X*imageZ + Y*imageZ*imageX;
    		        gridOS3 = gridOS * gridOS * gridOS;
    		        dst[common_index] = (src[common_index]) / gridKernel * (1.0f / gridOS3);
//    		        dst[common_index] = (src[common_index].real()) / gridKernel * (1.0f / gridOS3);
//    		        dst[common_index] += (src[common_index].imag()) / gridKernel * (1.0f / gridOS3)*1(std::complex_literals::if);
    		    }
    		}
    	}
    }
}

//__global__ static void
//CropCenterRegion2d_kernel(
//                          cufftComplex *dst, cufftComplex *src,
//                          int imageSizeX, int imageSizeY,
//                          int gridSizeX, int gridSizeY)
//{
//    int dY_dst = blockIdx.x;
//    int dX_dst = blockIdx.y;
//
//    int offsetY = (int)(((float)gridSizeY / 2.0f) - ((float)imageSizeY / 2.0f));
//    int offsetX = (int)(((float)gridSizeX / 2.0f) - ((float)imageSizeX / 2.0f));
//
//    int dY_src = dY_dst + offsetY;
//    int dX_src = dX_dst + offsetX;
//
//    int common_index_dst = dY_dst*imageSizeX + dX_dst;
//    int common_index_src = dY_src*gridSizeX  + dX_src;
//
//    dst[common_index_dst] REAL = src[common_index_src] REAL;
//    dst[common_index_dst] IMAG = src[common_index_src] IMAG;
//}
template <typename T1>
void
crop_center_region2d(
					std::complex<T1> *dst,std::complex<T1> *src,
                     int imageSizeX, int imageSizeY,
                     int gridSizeX, int gridSizeY)
{
    /* (gridSizeX,gridSizeY) is the size of 'src' */
//    assert( (!(gridSizeX%2) && !(gridSizeY%2) ) );
//    assert( (!(imageSizeX%2) && !(imageSizeY%2) ) );

    int offsetY;
    int offsetX;
    int dY_src;
    int dX_src;
    int common_index_dst;
    int common_index_src;
    for (int dY_dst=0;dY_dst<ImageSizeY;dY_dst++)
    {
    	for (int dX_dst=0;dX_dst<ImageSizeX;dX_dst++)
    	{
    	    offsetY = (int)(((float)gridSizeY / 2.0f) - ((float)imageSizeY / 2.0f));
    	    offsetX = (int)(((float)gridSizeX / 2.0f) - ((float)imageSizeX / 2.0f));

    	    dY_src = dY_dst + offsetY;
    	    dX_src = dX_dst + offsetX;

    	    common_index_dst = dY_dst*imageSizeX + dX_dst;
    	    common_index_src = dY_src*gridSizeX  + dX_src;

    	    dst[common_index_dst] = src[common_index_src];
    	}
    }

}
//
//__global__ static void
//CropCenterRegion3d_kernel(
//                          cufftComplex *dst, cufftComplex *src,
//                          int imageSizeX, int imageSizeY, int imageSizeZ,
//                          int gridSizeX, int gridSizeY, int gridSizeZ)
//{
//    int dY_dst = blockIdx.x;
//    int dX_dst = blockIdx.y;
//    int dZ_dst = threadIdx.x;
//
//    int offsetY = (int)(((float)gridSizeY / 2.0f) - ((float)imageSizeY / 2.0f));
//    int offsetX = (int)(((float)gridSizeX / 2.0f) - ((float)imageSizeX / 2.0f));
//    int offsetZ = (int)(((float)gridSizeZ / 2.0f) - ((float)imageSizeZ / 2.0f));
//
//    int dY_src = dY_dst + offsetY;
//    int dX_src = dX_dst + offsetX;
//    int dZ_src = dZ_dst + offsetZ;
//
//    int common_index_dst = dY_dst*imageSizeX*imageSizeZ + dX_dst*imageSizeZ + dZ_dst;
//    int common_index_src = dY_src*gridSizeX*gridSizeZ   + dX_src*gridSizeZ  + dZ_src;
//
//    dst[common_index_dst] REAL = src[common_index_src] REAL;
//    dst[common_index_dst] IMAG = src[common_index_src] IMAG;
//}

template <typename T1>
void
crop_center_region3d(
					 std::complex<T1> *dst,std::complex<T1> *src,
                     int imageSizeX, int imageSizeY, int imageSizeZ,
                     int gridSizeX, int gridSizeY, int gridSizeZ)
{
    /* (gridSizeX,gridSizeY,gridSizeZ) is the size of 'src' */
//    assert( (!(gridSizeX%2) && !(gridSizeY%2) && !(gridSizeZ%2)) );
//    assert( (!(imageSizeX%2) && !(imageSizeY%2) && !(imageSizeZ%2)) );
    
    int offsetY;
    int offsetX;
    int offsetZ;

	int dY_src;
	int dX_src;
	int dZ_src;
	int common_index_dst;
	int common_index_src;

    for (int dY_dst = 0; dY_dst < imageSizeY; dY_dst++)
    {
    	for (int dX_dst = 0; dX_dst < imageSizeX; dX_dst++)
    	{
    		for (int dZ_dst = 0; dZ_dst < imageSizeZ; dZ_dst++)
    		{
    		    offsetY = (int)(((float)gridSizeY / 2.0f) - ((float)imageSizeY / 2.0f));
    		    offsetX = (int)(((float)gridSizeX / 2.0f) - ((float)imageSizeX / 2.0f));
    		    offsetZ = (int)(((float)gridSizeZ / 2.0f) - ((float)imageSizeZ / 2.0f));

    		    dY_src = dY_dst + offsetY;
    		    dX_src = dX_dst + offsetX;
    		    dZ_src = dZ_dst + offsetZ;

    		    common_index_dst = dY_dst*imageSizeX*imageSizeZ + dX_dst*imageSizeZ + dZ_dst;
    		    common_index_src = dY_src*gridSizeX*gridSizeZ   + dX_src*gridSizeZ  + dZ_src;

    		    dst[common_index_dst]= src[common_index_src];
    		}
    	}
    }
}

void
cuda_fft2shift_grid(
                    cufftComplex *src, cufftComplex *dst,
                    int dimY, int dimX, int inverse)
{
    //(dimX,dimY) is the size of 'src'
    
    int pivotY = 0;
    int pivotX = 0;
    if(inverse) {
        pivotY = (int)floor(float(dimY / 2));
        pivotX = (int)floor(float(dimX / 2));
    } else {
        pivotY = (int)ceil(float(dimY / 2));
        pivotX = (int)ceil(float(dimX / 2));
    }
    
    dim3 threads(FFTSHIFT_TILE_SIZE_X, FFTSHIFT_TILE_SIZE_Y);
    dim3 blocks(pivotX / FFTSHIFT_TILE_SIZE_X, pivotX / FFTSHIFT_TILE_SIZE_Y);
    
    CudaFFTShift<<<blocks, threads>>>
    (src, pivotY, pivotX, dimX, dst);
}

void
cuda_fft3shift_grid(
                    cufftComplex *src, cufftComplex *dst,
                    int dimY, int dimX, int dimZ, int inverse)
{  
    //(dimX,dimY,dimZ) is the size of 'src'
    
    int pivotY = 0;
    int pivotX = 0;
    int pivotZ = 0;
    if(inverse) {
        pivotY = (int)floor(float(dimY / 2));
        pivotX = (int)floor(float(dimX / 2));
        pivotZ = (int)floor(float(dimZ / 2));
    } else {
        pivotY = (int)ceil(float(dimY / 2));
        pivotX = (int)ceil(float(dimX / 2));
        pivotZ = (int)ceil(float(dimZ / 2));
    }
    
    dim3 threads(pivotZ, 1);
    dim3 blocks(pivotY, pivotX);
    
    CudaFFT3Shift<<<blocks, threads>>> 
    (src, dst, pivotY, pivotX, pivotZ, dimX, dimZ);
}


#endif
