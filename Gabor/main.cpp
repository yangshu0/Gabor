#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

// my_gabor.cpp : 定义控制台应用程序的入口点。  
//  

#include<iostream>  
#include "math.h"  
#define  PI 3.1415926  
#define  N 4  

using namespace std;
using namespace cv;

void m_filer(unsigned char *src, int height, int width, double *mask_rel, double *mask_img, int mW, int mH, int k)
{
	Mat tmp = Mat::zeros(Size(width, height), CV_8UC1);
	double a, b, c;
	char res[20];       //保存的图像名称  
	unsigned char* tmpData = tmp.data;

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			a = 0.0;
			b = 0.0;
			c = 0.0;
			//去掉靠近边界的行列，无法滤波，超出范围  
			if (i > int(mH / 2) && i < height - int(mH / 2) && j > int(mW) && j < width - int(mW / 2))
			{
				for (int m = 0; m < mH; m++)
				{
					for (int n = 0; n < mW; n++)
					{
						//printf("%f\n",src[(i+m-int(mH/2))*width+(j+n-int(mW))]);  
						a += src[(i + m - int(mH / 2))*width + (j + n - int(mW))] * mask_rel[m*mW + n];
						b += src[(i + m - int(mH / 2))*width + (j + n - int(mW))] * mask_img[m*mW + n];
						//printf("%f,%f\n",a,b);  
					}
				}
			}
			c = sqrt(a*a + b*b);
			c /= mW * mH;
			tmpData[i*width + j] = (unsigned char)c;
		}
	}
	sprintf(res, "result%d.jpg", k);
	imshow("tmp", tmp);
	waitKey(0);
	//cvSaveImage(res, tmp);
}

int main(int argc, char* argv[])
{
	Mat src;
	double *rel, *img, *src_data, xtmp, ytmp, tmp1, tmp2, tmp3, re, im;
	double Theta, sigma, Gamma, Lambda, Phi;        //公式中的5个参数  
	int gabor_height, gabor_width, x, y;

	src = imread("test.jpg", CV_LOAD_IMAGE_GRAYSCALE);

	gabor_height = 10;
	gabor_width = 10;
	Gamma = 1.0;
	Lambda = 10.0;
	sigma = 100;
	Phi = 0;

	rel = new double[sizeof(double)*gabor_width*gabor_height];//实数部分  
	img = new double[sizeof(double)*gabor_width*gabor_height];//虚数部分  
	
	unsigned char* data = src.data;

	//构造gabor函数  
	for (int k = 0; k < N; k++)                    //定义N方向  
	{
		Theta = PI*((double)k / N);
		for (int i = 0; i < gabor_height; i++) //定义模版大小  
		{
			for (int j = 0; j < gabor_width; j++)
			{
				x = j - gabor_width / 2;
				y = i - gabor_height / 2;

				xtmp = (double)x*cos(Theta) + (double)y*sin(Theta);
				ytmp = (double)y*cos(Theta) - (double)x*sin(Theta);

				tmp1 = exp(-(pow(xtmp, 2) + pow(ytmp*Gamma, 2)) / (2 * pow(sigma, 2)));
				tmp2 = cos(2 * PI*xtmp / Lambda + Phi);
				tmp3 = sin(2 * PI*xtmp / Lambda + Phi);

				re = tmp1*tmp2;
				im = tmp1*tmp3;

				rel[i*gabor_width + j] = re;
				img[i*gabor_width + j] = im;
				//printf("%f,%f\n",re,im);  
			}
		}
		//用不同方向的GABOR函数对图像滤波并保存图片  
		m_filer(data, src.rows, src.cols, rel, img, 10, 10, k);
	}
	imshow("lena", src);
	waitKey(0);
	free(rel); 
	free(img);
	return 0;
}