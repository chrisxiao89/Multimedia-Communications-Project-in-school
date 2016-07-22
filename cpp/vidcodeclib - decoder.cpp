#include <iostream>
#include <fstream>
using namespace std;

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "codeclib.h"
#include "arith.h"
#include "bits.h"
#include "xform.h"
#include "quant.h"

//-----------------------------------
// ICodec class
//-----------------------------------
ICodec::ICodec(unsigned w, unsigned h)
  :_w(w), _h(h), _w2(w / 2), _h2(h / 2),
  m_iMBNumW (w / MBSIZE), m_iMBNumH (h / MBSIZE)
{
  unsigned i;

  //Init Image buffer: the (y, x)-th pixel of the image can be accessed  by _ydata[y][x].
  m_pfFrameBuf[0] = new float[_w * _h * 3 / 2];
  m_pfFrameBuf[1] = new float[_w * _h * 3 / 2];

  _ydata = new float *[_h];
  _udata = new float *[_h2];
  _vdata = new float *[_h2];

  _ydataRef = new float *[_h];
  _udataRef = new float *[_h2];
  _vdataRef = new float *[_h2];

  SetFrameBufPointer(m_pfFrameBuf[0], m_pfFrameBuf[1]);

  //MV buffer
  m_iMVy = new int *[m_iMBNumH];
  m_iMVy[0] = new int[m_iMBNumW * m_iMBNumH];
  for(i = 1; i < m_iMBNumH; i++) {
    m_iMVy[i] = m_iMVy[i-1] + m_iMBNumW;
  }

  m_iMVx = new int *[m_iMBNumH];
  m_iMVx[0] = new int[m_iMBNumW * m_iMBNumH];
  for(i = 1; i < m_iMBNumH; i++) {
    m_iMVx[i] = m_iMVx[i-1] + m_iMBNumW;
  }

  //context prob models
  for (i = 0; i < MAX_CONTEXT_NUM; i++) {
      _context[i] = new binModel(1, 1, _w); 
  }

}

ICodec::~ICodec(void)
{
  int i;

  delete m_pfFrameBuf[0];
  delete m_pfFrameBuf[1];

  delete _ydata;
  delete _udata;
  delete _vdata;

  delete _ydataRef;
  delete _udataRef;
  delete _vdataRef;
 
  delete m_iMVy[0];
  delete m_iMVy;
  delete m_iMVx[0];
  delete m_iMVx;

  for (i = 0; i < MAX_CONTEXT_NUM; i++) {
      delete _context[i];
  }
}

//Two set of buffers are needed for video coding:
//_ydata, _udata, _vdata hold YUV of the current frame
//_ydataRef, _udataRef, _vdataRef hold the reference frame
void ICodec::SetFrameBufPointer(
    float *pfBuf1,
    float *pfBuf2
)
{
  unsigned y;

  //_d, _udata, _vdata hold the current frame
  _ydata[0] = pfBuf1;
  for(y = 1; y < _h; y++) {
    _ydata[y] = _ydata[y-1] + _w;
  }

  _udata[0] = _ydata[0] + _w * _h;
  for(y = 1; y < _h2; y++) {
    _udata[y] = _udata[y-1] + _w2;
  }

  _vdata[0] = _udata[0] + _w2 * _h2;
  for(y = 1; y < _h2; y++) {
    _vdata[y] = _vdata[y-1] + _w2;
  }

  //_ydataRef, _udataRef, _vdataRef hold the reference frame
  _ydataRef[0] = pfBuf2;
  for(y = 1; y < _h; y++) {
    _ydataRef[y] = _ydataRef[y-1] + _w;
  }

  _udataRef[0] = _ydataRef[0] + _w * _h;
  for(y = 1; y < _h2; y++) {
    _udataRef[y] = _udataRef[y-1] + _w2;
  }

  _vdataRef[0] = _udataRef[0] + _w2 * _h2;
  for(y = 1; y < _h2; y++) {
    _vdataRef[y] = _vdataRef[y-1] + _w2;
  }

}

//Create a reconstruced frame of the current frame,
//so that the next frame can use it as reference frame.
//Current frame already contains the reconstructed prediction error.
void ICodec::GetReconstructedFrame()
{
    unsigned y, x;

    for(y = 0; y < _h; y += BLOCKSIZE) {
        for(x = 0; x < _w; x += BLOCKSIZE) {
            GetBlockRecon(_ydata, _ydataRef, y, x, m_iMVy[y / MBSIZE][x / MBSIZE], m_iMVx[y / MBSIZE][x / MBSIZE], 'y');
        }
    }

    for(y = 0; y < _h2; y += BLOCKSIZE) {
        for(x = 0; x < _w2; x += BLOCKSIZE) {
            GetBlockRecon(_udata, _udataRef, y, x, m_iMVy[y / MBSIZEUV][x / MBSIZEUV] / 2, m_iMVx[y / MBSIZEUV][x / MBSIZEUV] / 2, 'u');
        }
    }

    for(y = 0; y < _h2; y += BLOCKSIZE) {
        for(x = 0; x < _w2; x += BLOCKSIZE) {
            GetBlockRecon(_vdata, _vdataRef, y, x, m_iMVy[y / MBSIZEUV][x / MBSIZEUV] / 2, m_iMVx[y / MBSIZEUV][x / MBSIZEUV] / 2, 'v');
        }
    }

}


void ICodec::GetBlockRecon(
    float **pfCurrFrame,    
    float **pfRefFrame,     
    int y0,     //y0, x0 give the upper-left corner of the block
    int x0, 
    int iMVy,   //iMVy, iMVx: MV.
    int iMVx,
	char yuv)
{
    int i, j;

    for (i = 0; i < BLOCKSIZE; i++) {
        for (j = 0; j < BLOCKSIZE; j++) {

			if( MBSIZE/2 == 8) {
				if(yuv == 'y'){
					if(y0 + i + iMVy >143 || y0 + i + iMVy < 0) {
						continue;
					}
				}else if(yuv == 'v' || yuv == 'u'){
					if(y0 + i + iMVy >71 || y0 + i + iMVy < 0) {
						continue;
					}
				}
			}

			//cout << "top left corner of x and y is (" << x0 << " , " << y0 << endl;
			//cout << "pfCurrFrame[" << y0 + i << "][" << x0 + j << "] is " << pfCurrFrame[y0 + i][x0 + j] << endl;
			//cout << "ppfRefFrame[" << y0 + i + iMVy<< "][" << x0 + j + iMVx << "] is " << pfRefFrame[y0 + i][x0 + j] << endl;

            pfCurrFrame[y0 + i][x0 + j] = pfCurrFrame[y0 + i][x0 + j] + pfRefFrame[y0 + i + iMVy][x0 + j + iMVx];
        }
    }
}


// Map signed integer to unsigned index for Golomb Rice encoding
// Positive numbers are mapped to odd indices.
// Negative numbers are mapped to even indices.
unsigned ICodec::ForwardGolombRiceIndexMapping(
    int iNumber)
{
    //Convert signed int to unsigned index:
    //iNumber:   0     1   -1  2   -2  3   -3 ...
    //uiIndex:   0     1    2  3    4  5    6 ...
    unsigned uiIndex;
    if (iNumber > 0) {
        uiIndex = 2 * iNumber - 1;
    } else {
        uiIndex = (unsigned) (-2 * iNumber);
    }
    return uiIndex;
}

// Map unsigned index back to signed integer for Golomb Rice decoding
// Odd indices are mapped to positive numbers.
// Even indices are mapped to negative numbers.
int ICodec::InverseGolombRiceIndexMapping(
    unsigned uiIndex)
{
    //Convert signed int to unsigned index:
    //uiIndex: 0     1    2  3    4  5    6 ...
    //iNumber: 0     1   -1  2   -2  3   -3 ...

    int iNumber;
    if (uiIndex & 1) {
        iNumber = (uiIndex + 1) >> 1;
    } else {
        iNumber = -1 * ((int) (uiIndex >> 1));
    }
    return iNumber;
}


//Copy image to the internal buffer
//and subtract 128.
void ICodec::SetImage(
	unsigned char *pcBuf)
{
  unsigned len = _w * _h, i;
  for(i = 0; i < len; i++) {
    _ydata[0][i] = (float) (pcBuf[i]);
  }

  pcBuf += len;
  for(i = 0; i < len / 4; i++) {
    _udata[0][i] = (float) (pcBuf[i]);
  }

  pcBuf += len / 4;
  for(i = 0; i < len / 4; i++) {
    _vdata[0][i] = (float) (pcBuf[i]);
  }
}

//Read internal image to external buffer
//After SetFrameBufPointer in CodeImage, the recently decoded image is in Ref buffer.
void ICodec::GetImage(
	unsigned char *pcBuf)
{
  int tmp;
  unsigned len = _w * _h, i;

  for(i = 0; i < len; i++) {
      //round to integer
      tmp = (int) (_ydataRef[0][i] + 0.5);

      //clipping to [0, 255]
      pcBuf[i] = (unsigned char) (tmp < 0 ? 0 : (tmp > 255 ? 255 : tmp));
  }

  pcBuf += len;
  for(i = 0; i < len / 4; i++) {
      //round to integer
      tmp = (int) (_udataRef[0][i] + 0.5);

      //clipping to [0, 255]
      pcBuf[i] = (unsigned char) (tmp < 0 ? 0 : (tmp > 255 ? 255 : tmp));
  }

  pcBuf += len / 4;
  for(i = 0; i < len / 4; i++) {
      //round to integer
      tmp = (int) (_vdataRef[0][i] + 0.5);

      //clipping to [0, 255]
      pcBuf[i] = (unsigned char) (tmp < 0 ? 0 : (tmp > 255 ? 255 : tmp));
  }

}


//-----------------------------------
// IEncoder class
//-----------------------------------

IEncoder::IEncoder(unsigned w,unsigned h)
  :ICodec(w, h)
{
  _ace = new ACEncoder();
  _out = new OFlow();
}

IEncoder::~IEncoder(void)
{
  delete _ace;
  delete _out;
}

// Perform unary code and binary arithmetic code:
// The first bit use Conext 0, and the rest bits use Context 1.
// Similar approach to H.264, since the first bit has more prob of 0.
void IEncoder::EncodeUnary(
    unsigned uiIndex)   //index to be encoded
{    
    if (uiIndex == 0) {
        //Codeword is 0: code with Context 0.
        _ace->codeSymbol(0, _context[0], _out);
    } else {
        //Other codeword has the fomat of 111...10
        //Code the first bit with Context 0.
        _ace->codeSymbol(1, _context[0], _out);
        uiIndex--;

        //Code other bits with Context 1.
        while (uiIndex > 0) {
            _ace->codeSymbol(1, _context[1], _out);
            uiIndex--;
        }

        //the last bit is 0, still use context 1
        _ace->codeSymbol(0, _context[1], _out);
    }    
}


// Perform Golomb-Rice code and binary arithmetic code:
// The first bit use Conext 0, and the rest bits use Context 1.
// Similar approach to H.264, since the first bit has more prob of 0.
void IEncoder::EncodeGolombRice(
    unsigned uiIndex,           //index to be encoded
    int iGRPara)                //Number of Golomb-Rice remainder bits
{    
    bool bBit;

    //Encode Group ID (with 1 or 2 contexts)
    unsigned uiGroup = uiIndex >> iGRPara;
    EncodeUnary(uiGroup);

    //encode remainder bits with Context 1
    for (int i = iGRPara - 1; i >= 0; i--) {
        bBit = (uiIndex >> i) & 1;
        _ace->codeSymbol(bBit, _context[1], _out);
    }
}


// Perform Exp-Golomb code and binary arithmetic code:
// The first bit use Conext 0, and the rest bits use Context 1.
// Similar approach to H.264, since the first bit has more prob of 0.
void IEncoder::EncodeExpGolomb(
    unsigned uiIndex)           //index to be encoded
{    
    bool bBit;

    if (uiIndex == 0) {
        EncodeUnary(0);
    } else {
        //Get Group ID
        unsigned uiGroup = 2;
        while ((int) uiIndex > (1 << uiGroup) - 2) {
            uiGroup++;
        }
        uiGroup = uiGroup - 1;

        //Encode Group ID (with 1 or 2 contexts)
        EncodeUnary(uiGroup);

        //encode offset within each group bits with Context 1
        unsigned uiOffset = uiIndex - (1 << uiGroup) + 1;
        for (int i = uiGroup - 1; i >= 0; i--) {
            bBit = (uiOffset >> i) & 1;
            _ace->codeSymbol(bBit, _context[1], _out);
        }
    }
}

//jennys section*************************************************************************
int getAVG(long int *Arr, int halfSIZE)
{
	long int average=0;
	long int sum = 0;
	for(int i = 0; i < halfSIZE*halfSIZE; i++)
	{
		//cout << Arr[i] << " ";
		sum += Arr[i];
	}
	average = sum/(halfSIZE*halfSIZE);

	return average;
}

bool CBsizeFLAG(long int *avgVL)
{
	long int sum = 0;
	for(int i = 1; i < 4; i++)
	{
		sum += avgVL[i] - avgVL[i - 1];
	}
	sum += avgVL[3] - avgVL[0];

	//cout << "\nsum = " << abs(sum) << endl;
	//cout << " Diff = " << abs(sum) - 20000000 << endl;

	if( abs(sum) < 20000000)
	{
		//cout << "flag = 0" << endl;
		return false;
	}
	else
	{
		//cout << "flag = 1" << endl;
		return true;
	}
}

bool SplitBlock(int *CTU, int fullSIZE)
{
	int halfSIZE = fullSIZE/2;
	long int avgVL[4]={0};
	long int *Arr1, *Arr2, *Arr3, *Arr4;
	Arr1 = new long int[halfSIZE*halfSIZE];
	Arr2 = new long int[halfSIZE*halfSIZE];
	Arr3 = new long int[halfSIZE*halfSIZE];
	Arr4 = new long int[halfSIZE*halfSIZE];

	//makes block 1
	for(int j = 0; j < halfSIZE; j++)
	{
		for(int i = 0; i < halfSIZE; i++)
		{
			Arr1[i + j*halfSIZE] = CTU[i + j*2*halfSIZE];
		}
	}
	
	avgVL[0] = getAVG(Arr1, halfSIZE);

	//makes block 2
	for(int j = 0; j < halfSIZE; j++)
	{
		for(int i = 0; i < halfSIZE; i++)
		{
			Arr2[i + j*halfSIZE] = CTU[i + (1 + j*2)*halfSIZE];
		}
	}

	avgVL[1] = getAVG(Arr2, halfSIZE);

	//makes block 3
	for(int j = 0; j < halfSIZE; j++)
	{
		for(int i = 0; i < halfSIZE; i++)
		{
			Arr3[i + j*halfSIZE] = CTU[i + (halfSIZE + j)*2*halfSIZE];
		}
	}

	avgVL[2] = getAVG(Arr3, halfSIZE);

	//makes block 4
	for(int j = 0; j < halfSIZE; j++)
	{
		for(int i = 0; i < halfSIZE; i++)
		{
			Arr4[i + j*halfSIZE] = CTU[i + (2*halfSIZE + 1 + j*2)*halfSIZE];
		}
	}

	avgVL[3] = getAVG(Arr4, halfSIZE);
	
	bool temp = false;
	temp = CBsizeFLAG(avgVL);
	//cout << "boolean after CBsizeFLAG function= " << temp << endl;
	return temp;
}
//***************************************************************************************

// Main function to encode an image
int IEncoder::codeImage(
    bool bIsIFrame,
	unsigned char *pcBitstreamBuf,   //output buffer
    float fQstep)                    //quantization step size,
{
  unsigned x, y;
  int offsetY = 0, offsetX = 0, counter = 0;
  bool *CBsizeFlag;
  int CTU[MAX_CONTEXT_NUM] = {0};
  CBsizeFlag = new bool[(_h / MBSIZE)*(_w / MBSIZE)];

  //Initialize output buffer
  _out->reset(pcBitstreamBuf);

  //Initialize arithmetic coder
  _ace->start();

 //Choose the valid CBsize, calculate the avg value of each pixel and compara
  //Split frame to CTU
  /*have 9*11 CTU(when MB is 16*16)  in total, each CTU goes through Jenny's function
   Jenny's fuction will return an array with flags of CBsize
   where MBsize^2 is 0, halfsize^2 is 1. 
   slipblock(CTUblk);*/
  
  for(y = 0; y < _h; y += MBSIZE){
    for(x = 0; x < _w; x += MBSIZE){

	  //cout << "\nBlock " << counter << "\n";

      //Making CTU
      for(offsetY = 0; offsetY < MBSIZE; offsetY ++){
	    for(offsetX = 0; offsetX < MBSIZE; offsetX ++){
		//  cout << "_ydata[" << (x+offsetX) + (y+offsetY)*_w << "] = CTU[" << offsetX + offsetY*MBSIZE << "]" << endl;
          CTU[offsetX + offsetY*MBSIZE] = (int)_ydata[ (x+offsetX) + (y+offsetY)*_w];
        }
      }
      //Get CBsizeflag
      CBsizeFlag[counter] = SplitBlock(CTU,MBSIZE);
      counter ++;

    }//for x
  }//for y
  
  
  if (!bIsIFrame) {
      //P frames: motion est, find prediction error
      MotionEst(CBsizeFlag);

      EncodeMV();

      GetPredError();
  }
  
  //encode Y component: _d contains reconstructed result after codeBlock().
  for(y = 0; y < _h; y += BLOCKSIZE) {
	  for(x = 0; x < _w; x += BLOCKSIZE) {
			codeBlock(_ydata, y, x, fQstep);
	  }
   }
  //encode U component
  for(y = 0; y < _h2; y += BLOCKSIZE) {
      for(x = 0; x < _w2; x += BLOCKSIZE) {
          codeBlock(_udata, y, x, fQstep);
      }
  }

  //encode V component
  for(y = 0; y < _h2; y += BLOCKSIZE) {
      for(x = 0; x < _w2; x += BLOCKSIZE) {
          codeBlock(_vdata, y, x, fQstep);
      }
  }

  if (!bIsIFrame) {
      //Get reconstructed frame, will be used as reference for the next frame
      GetReconstructedFrame();
  }

  //swap pointers:
  //the next frame will be written into the ref frame of the current frame,
  //and the current frame will become the ref frame of the next frame.
  SetFrameBufPointer(_ydataRef[0], _ydata[0]);
  
  _ace->stop(_out);
  return _out->bytesUsed();
}

//Encode an image block
void IEncoder::codeBlock(
    float **buf,    //Y/U/V buffer header,
    unsigned y,     //y index of upper-left corner of the block
    unsigned x,     //x index of upper-left corner of the block
    float fQstep)   //quantization step size,
{
    unsigned i, j, uiIndex;
    bool skipblk = false;

    //update m_fBlkBuf so that it points to the current block
    for (i = 0; i < BLOCKSIZE; i++) {
        m_fBlkBuf[i] = buf[y + i] + x;
    }
  
    //Forward DCT
    if (BLOCKSIZE == 4) {
        Transform::FDCT4(m_fBlkBuf);
    } else {
        Transform::FDCT8(m_fBlkBuf);
    }

    //Quantization
    Quant::QuantMidtread(m_fBlkBuf, BLOCKSIZE, fQstep);

    //check whether there is any non-zero coeff, 
    //send 1 bit flag to signal this: 1: skip, 0: not skip, encode all coeffs after this.
    //this is a simplified version of CBP, but we send it at the beginning of each block.
    skipblk = true;
    for (i = 0; i < BLOCKSIZE; i++) {
        for (j = 0; j < BLOCKSIZE; j++) {
            if (m_fBlkBuf[i][j] != 0) {
                skipblk = false;
                break;
            }
        }
        if (!skipblk) {
            break;
        }
    }
    //encode the skip bit using context 2.
    _ace->codeSymbol(skipblk, _context[2], _out);

    //encode all coeffs of m_fBlkBuf, not efficient since we encode all zeros at the end.
    if (!skipblk) {
        for (i = 0; i < BLOCKSIZE; i++) {
            for (j = 0; j < BLOCKSIZE; j++) {
                //convert signed int to unsigned index for G-R code
                uiIndex = ForwardGolombRiceIndexMapping((int) m_fBlkBuf[i][j]);

                //encode with unary code and arithmetic code
                EncodeUnary(uiIndex);
            }
        }
    }

    //Get recon for prediction
    //Dequantization
    Quant::DequantMidtread(m_fBlkBuf, BLOCKSIZE, fQstep);

    //IDCT
    if (BLOCKSIZE == 4) {
        Transform::IDCT4(m_fBlkBuf);
    } else {
        Transform::IDCT8(m_fBlkBuf);
    }

}

//find MV of each macro-block
void IEncoder::MotionEst(bool *CBsizeFlag)
{
    unsigned y, x;
	int counter=0;

    for(y = 0; y < _h; y += MBSIZE) {
        for(x = 0; x < _w; x += MBSIZE) {
			if(CBsizeFlag[counter] == 0){
				MBMotionEst(_ydata, _ydataRef, y, x, MBSIZE);
			}else{
				MBMotionEst(_ydata, _ydataRef, y, x, MBSIZE/2);
			}
        
		counter++;
		}
    }
	
	//cout << "motions estimation is done" << endl << endl;
}

void IEncoder::MBMotionEst(
	float **pfCurrFrame,
	float **pfRefFrame,
	int y,
	int x,
	//********************************Change by jerry
	int CBsize
	
	//********************************Change by jerry
	)
{
	int i, j, iSAD0, iSAD;
	int iMVy = 0, iMVx = 0, range = 0,offa = 0, offb = 0;
	int PBsize = 0;
	int block8[4] = {};
	int initialX = 0, initialY = 0;
	//************************************************************
	if (y >= 136 ) {
		//cout << "In MBEst, y is " << y << " and X is " << x << endl;
	}
	if (CBsize == 8)
	{
		
		PBsize = 4;
	}
	else
	{
		
		PBsize = CBsize;
	}
		//************************************************************
	switch (PBsize)
	{
	case 4:
		range = 6;
		break;
	case 16:
		range = 18;
		break;
	default:
		cout << "No proper CBsize is chosen";
	}
	
	//************************************************************
		
	iSAD0 = GetSAD(pfCurrFrame, pfRefFrame, y, x, 0, 0, 65535, PBsize);
	
	for (i = -range; i <= range; i++) {
       for (j = -range; j <= range; j++) {
		  
           //prevent MV from pointing to outside
           if (y + i  < 0 || y + i + PBsize - 1 >= (int) _h || x +  j < 0 || x + j +  PBsize - 1 >= (int) _w) 
		   {
               continue;
           }

           if (i == 0 && j == 0) 
		   {
               continue;
           }

           //get SAD
		   if (PBsize == 16)
		   {
			   iSAD = GetSAD(pfCurrFrame, pfRefFrame, y, x, i, j, iSAD0, PBsize);

			   //update min SAD
			   if (iSAD < iSAD0 * 0.925f)
			   {
				   iSAD0 = iSAD;
				   iMVy = i;
				   iMVx = j;
			   }

		   }
		   else
		   {
			  // cout << "at this moment, b4 entering BESTsad...... , offset y is " << y << " and b4 entering BESTsad...... offest X is " << x << endl;
		       BestSadForSmallestPB(pfCurrFrame, pfRefFrame, y, x, i, j, PBsize, iSAD0, block8);
			  // cout << "after best sad...., offset y is " << y << " after best sad.... and offest X is " << x << endl;

			   //cout << "iMVY is " << block8[0] << " and iMVX is " << block8[1] << endl << endl;
			  iMVy = block8[0];
			  iMVx = block8[1];
		   }
		   
		 

       }
    }
	//cout << "iMVY is " << block8[0] << " and iMVX is " << block8[1] << endl << endl;
	
   m_iMVy[y / MBSIZE][x / MBSIZE] = iMVy;
   m_iMVx[y / MBSIZE][x / MBSIZE] = iMVx;
  // cout << "m_iMVy[" << y / MBSIZE << "][" << x / MBSIZE << "] is " << m_iMVy[y / MBSIZE][x / MBSIZE] << endl;
 //  cout << "m_iMVx[" << y / MBSIZE << "][" << x / MBSIZE << "] is " << m_iMVx[y / MBSIZE][x / MBSIZE] << endl << endl;
   
	
}
// ************************************** Jerry****************************************************
void IEncoder::BestSadForSmallestPB(
	float **pfCurrFrame, 
	float **pfRefFrame, 
	int y, 
	int x, 
	int i, 
	int j, 
	int PBsize, 
	int iSAD0,
	int* block8)
{
	int smallblock[16] = {};
	int count = 0, mim = 0;
	int bestX = 0, bestY = 0;
	int initial_y = y;
	int initial_x = x;
	int MAX3a=0; 
	int MAX3b=0;

	for (y = initial_y; y < (int)_h; y += 4)
	{
		if(MAX3a > 3) break;
		
		for(x = initial_x; x<(int)_w; x += 4)
		{
			if(MAX3b > 3) break;

			if (y + i < 0 || y + i + PBsize - 1 >= (int)_h || x + j < 0 || x + j + PBsize - 1 >= (int)_w)
			{
				//cout << "y = " << y << ",      offa = " << offa << " ,    PBsize = " << PBsize << ",      i = " << i << ",      y + offa + i + PBsize - 1 = " << y + offa + i + PBsize - 1 << endl;
				continue;
			}

			if (i == 0 && j == 0)
			{
				continue;
			}
			smallblock[count] = GetSAD(pfCurrFrame, pfRefFrame, y, x, i, j, iSAD0, PBsize);
			if (smallblock[count] < iSAD0)
			{
				iSAD0 = smallblock[count];
			
				bestY = i;
				bestX = j;
			}
			count++;
			MAX3b++;
		}//for offb
		MAX3a++;
	}//for offa

	block8[0] = bestY;
	block8[1] = bestX;

}
// ********************************Jerry*********************************
//return the sum of absolute difference of two blocks with the given motion vectors
int IEncoder::GetSAD(
    float **pfCurrFrame,    //pointer to the current frame,
    float **pfRefFrame,     //pointer to the reference frame,
    int y,                  //(y, x) is the upper-left corner of the current block
    int x,                  
    int iMVy,               //(iMVy, iMVx) is the given MV
    int iMVx,
    int iSAD0,
	int PBsize
	)              //Minimum SAD so far, used for early termination.
				  
{
    int iDiff = 0, iSAD = 0;
	
    //HW: compute the SAD
    for (int g = 0; g < PBsize; g++)// g = 0 to 1 to 2 to 3
	{
        for (int h = 0; h < PBsize; h++) // h = 0 to 1 to 2 to 3
		{
			/*if ( y >= 136) {
				cout << "pfCurrFrame[y + g][x + h] " << "("<< y + g <<" , "<<x + h <<" )"<< endl << "pfRefFrame[y + g + iMVy][x + h + iMVx] " << "(" << y + g + iMVy << " , " << x + h + iMVx << " )"<< endl;
				cout << "at this moment , x is " << x << " and y is "<< y <<endl;
				cout << "and Value of y: " << y << " and value of g is " << g << endl;
				cout << "Value of pfCurrFrame: " << pfCurrFrame[y + g][x + h] << endl;
				cout << "and Value of pfRefFrame: " << pfRefFrame[y + g + iMVy][x + h + iMVx] << endl;
				
			}*/
			
			iDiff = (int) (pfCurrFrame[y + g][x + h] - pfRefFrame[y + g + iMVy][x + h + iMVx]);

            iSAD += iDiff > 0 ? iDiff : -iDiff;
			
            if (iSAD > iSAD0) 
			{
                return iSAD;
            }
        }
    }
	
    return iSAD;
}

//Save the MV information in a separate file for Matlab to plot.
void IEncoder::DumpMV(ofstream& DumpFile)
{
    DumpFile.write((const char *) m_iMVy[0], _w * _h / MBSIZE / MBSIZE * sizeof(int));
    DumpFile.write((const char *) m_iMVx[0], _w * _h / MBSIZE / MBSIZE * sizeof(int));
}

//*******************************************************************************

void IEncoder::EncodeMV()
{
    unsigned y, x;
    unsigned uiIndex;

    for(y = 0; y < m_iMBNumH; y++) {
        for(x = 0; x < m_iMBNumW; x++) {
			//cout << "in encodeMV"<<"m_iMVy[" << y << "][" << x <<"] and value of m_iMVy[y][x] is "<< m_iMVy[y][x] << endl;
            uiIndex = ForwardGolombRiceIndexMapping(m_iMVy[y][x]);
            EncodeUnary(uiIndex);
			//cout << "EncodeMV : EncodeUnary for y direction is done"<< endl;
			//cout << "in encodeMV" << "m_iMVx[" << y << "][" << x << "] and value of m_iMVx[y][x] is " << m_iMVy[y][x] << endl;
			uiIndex = ForwardGolombRiceIndexMapping(m_iMVx[y][x]);
            EncodeUnary(uiIndex);
			//cout << "EncodeMV : EncodeUnary for x direction is done" << endl;
         }
    }
	
}

//Called by encoder after motion est to get the prediction error.
//The error will then be encoded.
void IEncoder::GetPredError()
{
    unsigned y, x;

    for(y = 0; y < _h; y += BLOCKSIZE) {
        for(x = 0; x < _w; x += BLOCKSIZE) {
			//cout << "GetBlockPredError for Y has started " << endl;
            GetBlockPredError(_ydata, _ydataRef, y, x, m_iMVy[y / MBSIZE][x / MBSIZE], m_iMVx[y / MBSIZE][x / MBSIZE],'y');
        }
    }
	//cout << "GetBlockPredError for Y is done " << endl;


    for(y = 0; y < _h2; y += BLOCKSIZE) {
        for(x = 0; x < _w2; x += BLOCKSIZE) {
			//cout << "GetBlockPredError for U has started " << endl;
            GetBlockPredError(_udata, _udataRef, y, x, m_iMVy[y / MBSIZEUV][x / MBSIZEUV] / 2, m_iMVx[y / MBSIZEUV][x / MBSIZEUV] / 2,'u');
        }
    }
	//cout << "GetBlockPredError for U is done " << endl;


    for(y = 0; y < _h2; y += BLOCKSIZE) {
        for(x = 0; x < _w2; x += BLOCKSIZE) {
            GetBlockPredError(_vdata, _vdataRef, y, x, m_iMVy[y / MBSIZEUV][x / MBSIZEUV] / 2, m_iMVx[y / MBSIZEUV][x / MBSIZEUV] / 2,'v');
        }
    }
	//cout << "GetBlockPredError for V is done " << endl;

}


void IEncoder::GetBlockPredError(
    float **pfCurrFrame,    
    float **pfRefFrame,     
    int y0,     //y0, x0 give the upper-left corner of the block
    int x0, 
    int iMVy,   //iMVy, iMVx: MV.
    int iMVx,
	char yuv)
{
    //HW4: Get prediction error of the current frame
    int i, j;
	
	//cout << "in GetBlockPredError, iMVy and iMVx=" << iMVy << " and " << iMVx << endl;
	
    for (i = 0; i < BLOCKSIZE; i++) {
        for (j = 0; j < BLOCKSIZE; j++) {
			if( MBSIZE/2 == 8) {
				if(yuv == 'y'){
					if(y0 + i + iMVy >143) {
						continue;
					}
				}else if(yuv == 'v' || yuv == 'u'){
					if(y0 + i + iMVy >71) {
						continue;
					}
				}
			}
				//cout << "top left corner of x and y is (" << x0 << " , " << y0 << endl;
				//cout << "pfCurrFrame[" << y0 + i << "][" << x0 + j << "] is " << pfCurrFrame[y0 + i][x0 + j] << endl;
				//cout << "ppfRefFrame[" << y0 + i + iMVy<< "][" << x0 + j + iMVx << "] is " << pfRefFrame[y0 + i][x0 + j] << endl;
		
			pfCurrFrame[y0 + i][x0 + j] = pfCurrFrame[y0 + i][x0 + j] - pfRefFrame[y0 + i + iMVy][x0 + j + iMVx];
			if (y0 >= 136) {
			//cout << "BlockPredError is "<< pfCurrFrame[y0 + i][x0 + j] << endl;
			
			
			}
        }

    }
	
}

//-----------------------------------
// IDecoder class
//-----------------------------------
IDecoder::IDecoder(unsigned w,unsigned h)
  :ICodec(w,h)
{
  _acd=new ACDecoder();
  _in=new IFlow();
}

IDecoder::~IDecoder(void)
{
  delete _acd;
  delete _in;
}

// Decode unary code and binary arithmetic code:
// The first bit use Conext 0, and the rest bits use Context 1.
// Similar approach to H.264, since the first bit has more prob of 0.
unsigned IDecoder::DecodeUnary()
{
  unsigned uiIndex = 0;
  bool nextbit = _acd->decodeSymbol(_context[0], _in);

  if (nextbit == 1) {
    uiIndex = 0;
    do {
      nextbit = _acd->decodeSymbol(_context[1], _in);
      uiIndex++;
    } while (nextbit != 0);

  }

  return uiIndex;
}


// Decode Golomb-Rice code and binary arithmetic code:
// The first bit use Conext 0, and the rest bits use Context 1.
// Similar approach to H.264, since the first bit has more prob of 0.
unsigned IDecoder::DecodeGolombRice(
    int iGRPara)      //Number of Golomb-Rice remainder bits
{

    bool nextbit;

    //Decode Group ID
    unsigned uiIndex = DecodeUnary();

    //Decode the remainder bits
    for (unsigned char i = 0; i < iGRPara; i++) {
        nextbit = _acd->decodeSymbol(_context[1], _in);
        uiIndex = (uiIndex << 1) + (unsigned) nextbit;
    }

    return uiIndex;
}

// Decode Exp-Golomb code and binary arithmetic code:
// The first bit use Conext 0, and the rest bits use Context 1.
// Similar approach to H.264, since the first bit has more prob of 0.
unsigned IDecoder::DecodeExpGolomb()
{

    bool nextbit;

    //Decode Group ID
    unsigned uiGroup = DecodeUnary();

    //Decode the offset within each group
    unsigned uiIndex = 0;
    for (unsigned char i = 0; i < uiGroup; i++) {
        nextbit = _acd->decodeSymbol(_context[1], _in);
        uiIndex = (uiIndex << 1) + (unsigned) nextbit;
    }
    uiIndex += ((1 << uiGroup) - 1);

    return uiIndex;
}


// Main function to decode an image
int IDecoder::decodeImage(
    bool bIsIFrame,
	unsigned char *b,    //input buffer
    float fQstep)        //quantization step size,
{
  unsigned x, y;

  //Initialize _in to input buffer
  _in->reset(b);

  //Initialize arithmetic decoder
  _acd->start(_in);

  //HW4:
  if (!bIsIFrame) {
      //decode MV for P frames
      DecodeMV();
  }

  for(y = 0; y < _h; y += BLOCKSIZE) {
      for(x = 0; x < _w; x += BLOCKSIZE) {
          decodeBlock(_ydata, y, x, fQstep);
      }
  }

  for(y = 0; y < _h2; y += BLOCKSIZE) {
      for(x = 0; x < _w2; x += BLOCKSIZE) {
          decodeBlock(_udata, y, x, fQstep);
      }
  }

  for(y = 0; y < _h2; y += BLOCKSIZE) {
      for(x = 0; x < _w2; x += BLOCKSIZE) {
          decodeBlock(_vdata, y, x, fQstep);
      }
  }

  if (!bIsIFrame) {
      //add back prediction for P frames
      GetReconstructedFrame();
  }

  //swap pointers:
  //the next frame will be written into the ref frame of the current frame,
  //and the current frame will become the ref frame of the next frame.
  SetFrameBufPointer(_ydataRef[0], _ydata[0]);

  return _in->bytesUsed();
}


//Decode intra block
void IDecoder::decodeBlock(
    float **buf, 
    unsigned y,     //y index of upper-left corner of the block
    unsigned x,     //x index of upper-left corner of the block
    float fQstep)   //quantization step size,
{
    unsigned i, j, uiIndex;
    bool skipblk = false;

    //update m_fBlkBuf so that it points to the current block
    for (i = 0; i < BLOCKSIZE; i++) {
        m_fBlkBuf[i] = buf[y + i] + x;
    }

    //decode the skip bit using context 2.
    skipblk = _acd->decodeSymbol(_context[2], _in);

    if (skipblk) {
        //if skip bit = 1, reset block buffer, skip inverse quant and inverse DCT.
        for (i = 0; i < BLOCKSIZE; i++) {
            for (j = 0; j < BLOCKSIZE; j++) {
                m_fBlkBuf[i][j] = 0;
            }
        }
    } else {
        //decode all coeffs of the current block into m_fBlkBuf:
        for (i = 0; i < BLOCKSIZE; i++) {
            for (j = 0; j < BLOCKSIZE; j++) {
                uiIndex = DecodeUnary();

                //convert to signed int
                m_fBlkBuf[i][j] = (float) InverseGolombRiceIndexMapping(uiIndex);
            }
        }

        //Dequantization
        Quant::DequantMidtread(m_fBlkBuf, BLOCKSIZE, fQstep);

        //IDCT
        if (BLOCKSIZE == 4) {
            Transform::IDCT4(m_fBlkBuf);
        } else {
            Transform::IDCT8(m_fBlkBuf);
        }
    }
}


void IDecoder::DecodeMV()
{
    unsigned y, x;
    unsigned uiIndex;

    for(y = 0; y < m_iMBNumH; y++) {
        for(x = 0; x < m_iMBNumW; x++) {
            uiIndex = DecodeUnary();
            m_iMVy[y][x] = InverseGolombRiceIndexMapping(uiIndex);

            uiIndex = DecodeUnary();
            m_iMVx[y][x] = InverseGolombRiceIndexMapping(uiIndex);
         }
    }
}



