/*********************************************************************************
 *
 * Compute the PSNR of two yuv sequences
 *
 * Jie Liang
 * School of Engineering Science
 * Simon Fraser University
 *
 * Feb. 2005
 *
 *********************************************************************************/

#include <iostream>
#include <fstream>
using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

void usage() {
    cout << "Usage: yuvpsnr seq1.yuv seq2.yuv w h frames" << endl;
    cout << "w: width" << endl;
    cout << "h: height" << endl;
    cout << "frames: number of frames to be measured." << endl;

}

double GetPSNR(
    unsigned char *pcImgBuf1, 
    unsigned char *pcImgBuf2,
    int iImageArea)
{
    //computer MSE and PSNR
    double mse = 0;
    double psnr = 0;

   //computer MSE
    for (int i = 0; i < iImageArea; i++) {
        mse += (double) (pcImgBuf1[i] - pcImgBuf2[i]) * (pcImgBuf1[i] - pcImgBuf2[i]);
    }
    mse /= iImageArea;

    //PSNR
    return (10 * log10(255.0 * 255.0 / mse));
}

//--------------------------------------------------------
//
// main routine
//
//--------------------------------------------------------
int main(int argc,char **argv) {
 
    ifstream ifs1, ifs2;
    int iWidth, iHeight, i, iFrames, iImageArea;
    unsigned char *pcImgBuf1, *pcImgBuf2;
    double psnr[3], tmp[3];

    if(argc < 6) {
        usage();
        return -1;
    }

    //open input file
    ifs1.open(argv[1], ios::in|ios::binary);
    ifs2.open(argv[2], ios::in|ios::binary);
    if (!ifs1 || !ifs2)
    {
        cout << "Can't open input file " << endl;
        return -2;
    }

    iWidth = atoi(argv[3]);    
    iHeight = atoi(argv[4]);    
    iFrames = atoi(argv[5]);    

    //allocate buffers
    iImageArea = iWidth * iHeight;
    pcImgBuf1 = new unsigned char[iImageArea];
    pcImgBuf2 = new unsigned char[iImageArea];
    if (!pcImgBuf1 || !pcImgBuf2) {
        cout << "Fail to create image buffer." << endl;
        return -7;
    }

    //Initialize YUV PSNRs
    psnr[0] = psnr[1] = psnr[2] = 0;
    tmp[0] = tmp[1] = tmp[2] = 0;

    cout << "PSNR (in dB) of Y U V:" << endl;

    for (i = 0; i < iFrames; i++) {
        //------------------------------------
        //HW5: 
        // Read the YUV components of the two input videos.
        // Example: ifs1.read((char *)pcImgBuf1, iImageArea);
        // Compute the PSNRs of the three components in each frame, save them in tmp[0] to tmp[2].
        // Calculate the average PSNR of each component over all frames.
        //------------------------------------

        cout << i << ": " << tmp[0] << "  " << tmp[1] << "  " << tmp[2] << endl;
    }

    cout << "Average: " << psnr[0] / iFrames << "  " << psnr[1] / iFrames << "  " << psnr[2] / iFrames << endl;

 
    ifs1.close();
    ifs2.close();

    delete pcImgBuf1;
    delete pcImgBuf2;

    return 0;
}