#pragma once

#include "fft.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

class FormantSeparator//��������
{
private:
	constexpr static int FFTSize = 1024;
	int windowSize = FFTSize;
	int hopSize = FFTSize / 4;

	float inbuf[FFTSize];
	float outbuf1[FFTSize * 2];
	float fftdatre1[FFTSize];
	float fftdatim1[FFTSize];
	int posout = 0;
	int pos = 0;

	float normv = 1.0;

	float window[FFTSize];
	void UpdataWindow()
	{
		float energy = 0;
		for (int j = 0; j < windowSize; ++j)
		{
			float v = (1.0 - cosf(2.0 * M_PI * j / windowSize)) * 0.5;
			window[j] = v;
			energy += v;
		}
		float n = (FFTSize / 4) / hopSize;
		normv = 1.0 / (FFTSize * n);
	}

	//formant separation
	//����������
	constexpr static int SliderWindowSize = 1024 + FFTSize;
	float sliderwindow[SliderWindowSize];//������ ƽ��ֵ�˲���
	float sliderbuf[SliderWindowSize];
	//������ض�����
	float ebufre[FFTSize];
	float ebufim[FFTSize];

	int sliderlen = 5;//�ضϳ��Ȼ��߻���������
	int isFormantUpdata = 0;
protected:
	void ProcessSTFT(float* rev1, float* imv1, int numBins, int hopSize)
	{  //������ȡ���׵Ĺ���������Ч������ܺ�
		isFormantUpdata = 1;
		memset(sliderwindow, 0, sizeof(sliderwindow));
		float sum = 0;
		for (int i = 0; i < numBins; ++i)
		{
			float mag = sqrtf(rev1[i] * rev1[i] + imv1[i] * imv1[i]);
			float lnv = logf(mag + 1e-8);//����
			sliderwindow[i + sliderlen] = lnv;
			sum += lnv;
			sum -= sliderwindow[i];
			sliderbuf[i] = sum;
		}
		for (int i = numBins; i < numBins + sliderlen; ++i)
		{
			sliderwindow[i + sliderlen] = 0;
			sum -= sliderwindow[i];
			sliderbuf[i] = sum;
		}
		float lastv = 0;
		for (int i = 0; i < numBins; ++i)
		{
			//float mag = sqrtf(rev1[i] * rev1[i] + imv1[i] * imv1[i]);
			//float lnv = logf(mag + 1e-8);//����

			float v = sliderbuf[i + sliderlen / 2] / sliderlen;
			//float dv = v - lastv;//б��
			//lastv = v;

			//float mix = fabs(dv) * 10.0;//б��ԽС��Խ�Ƿ���߹�
			//if (mix > 1.0)mix = 1.0;
			//v = lnv * (1.0 - mix) + v * mix;//���

			float expv = expf(v) + 0.01;//�����Ƶ��
			sliderbuf[i] = expv;//�洢�����

			//rev1��������
			rev1[i] = rev1[i] / expv;				//���빲�����ź�
			imv1[i] = imv1[i] / expv;
		}

		//������ضϷ���Ҳ������
		/*
		for (int i = 0; i < numBins; ++i)
		{
			float mag = sqrtf(rev1[i] * rev1[i] + imv1[i] * imv1[i]);
			float lnv = logf(mag + 1e-24);//����
			ebufre[i] = lnv;
			ebufim[i] = 0;
		}
		for (int i = 0; i < FFTSize / 2; ++i)
		{
			ebufre[FFTSize - i - 1] = ebufre[i];
			ebufim[FFTSize - i - 1] = -ebufim[i];
		}
		fft_f32(ebufre, ebufim, FFTSize, -1);//ifft
		for (int i = sliderlen; i < FFTSize; ++i)//�ض�
		{
			ebufre[i] = 0;
			ebufim[i] = 0;
		}
		fft_f32(ebufre, ebufim, FFTSize, 1);//fft

		for (int i = 0; i < numBins; ++i)
		{
			//float mage = sqrtf(ebufre[i] * ebufre[i] + ebufim[i] * ebufim[i]) / FFTSize;
			float mage = fabs(ebufre[i]) / FFTSize;
			float expv = expf(mage + 1e-24);//�����Ƶ��
			sliderbuf[i] = expv;//�洢�����

			//rev1��������
			rev2[i] = rev1[i] / expv;				//���빲�����ź�
			imv2[i] = imv1[i] / expv;
		}
		for (int i = 0; i < numBins; ++i)
		{
			int index = format * i;
			if (index < numBins)
			{
				rev1[i] = rev2[i] * sliderbuf[index];//�Ը���Ӧ��format
				imv1[i] = imv2[i] * sliderbuf[index];
			}
			else
			{
				rev1[i] = 0;
				imv1[i] = 0;
			}
		}*/
	}

public:
	FormantSeparator() {
		memset(inbuf, 0, sizeof(inbuf));
		memset(outbuf1, 0, sizeof(outbuf1));
		memset(fftdatre1, 0, sizeof(fftdatre1));
		memset(fftdatim1, 0, sizeof(fftdatim1));
		UpdataWindow();
	}
	void SetWindowAndHop(float window, float hop)//window = 1.0,hop = 0.25 is good
	{
		windowSize = window * FFTSize;
		hopSize = hop * FFTSize;
		UpdataWindow();
	}
	int GetWindowSize() { return windowSize; }
	int GetHopSize() { return hopSize; }
	int GetFFTSize() { return FFTSize; }

	void SetSliderLen(float sliderLen)//����������
	{
		this->sliderlen = sliderLen * 25 + 2;
	}

	void ProcessBlock(const float* in, float* out1, float numSamples)
	{
		for (int i = 0; i < numSamples; ++i)
		{
			inbuf[pos] = in[i];
			pos++;
			if (pos >= FFTSize)pos = 0;
			if (pos % hopSize == 0)//������ˣ�
			{
				for (int j = 0; j < windowSize; ++j)//��
				{
					fftdatre1[j] = inbuf[(pos + j) % FFTSize] * window[j];
					fftdatim1[j] = 0;
				}
				for (int j = windowSize; j < FFTSize; ++j)//��
				{
					fftdatre1[j] = 0;
					fftdatim1[j] = 0;
				}

				fft_f32(fftdatre1, fftdatim1, FFTSize, 1);//fft

				ProcessSTFT(fftdatre1, fftdatim1, FFTSize / 2, hopSize);

				for (int j = 0; j < FFTSize / 2; ++j)
				{
					fftdatre1[FFTSize - j - 1] = -fftdatre1[j];
					fftdatim1[FFTSize - j - 1] = fftdatim1[j];
				}
				fft_f32(fftdatre1, fftdatim1, FFTSize, -1);//ifft

				for (int j = 0; j < windowSize; ++j)//�Ӵ�
				{
					fftdatre1[j] *= window[j];
				}
				int pos2 = posout;
				for (int j = 0; j < windowSize; ++j)
				{
					pos2++;
					outbuf1[pos2 % (FFTSize * 2)] += fftdatre1[j];
				}
				for (int j = 0; j < hopSize; ++j)
				{
					pos2++;
					outbuf1[pos2 % (FFTSize * 2)] = 0;
				}
			}

			out1[i] = outbuf1[posout] * normv;
			posout++;
			if (posout >= FFTSize * 2)posout = 0;
		}
	}

	int IsFormantUpdata()//���������ظ����µ�
	{
		return isFormantUpdata;
	}
	int GetFormant(float*& formant)
	{
		formant = sliderbuf;
		isFormantUpdata = 0;
		return FFTSize / 2;//numBins
	}
};