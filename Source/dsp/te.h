#pragma once

#include "fft.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

class TransientExtractor
{

private:
	constexpr static int FFTSize = 1024;
	int windowSize = FFTSize;
	int hopSize = FFTSize / 8;

	float inbuf[FFTSize];
	float outbuf1[FFTSize * 2];//瞬态
	float outbuf2[FFTSize * 2];//稳态
	float fftdatre1[FFTSize];
	float fftdatim1[FFTSize];
	float fftdatre2[FFTSize];
	float fftdatim2[FFTSize];
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

	//phase vocoder & transient extraction
	float pitch = 1.0, matchf = 1.0;

	float lastphase[FFTSize / 2] = { 0 };
	float lastfreq[FFTSize / 2] = { 0 };
	float dfreqs[FFTSize / 2] = { 0 };//瞬时频率变化值

	float lastmag[FFTSize / 2] = { 0 };//上次的幅度值
	float dmags[FFTSize / 2] = { 0 };//瞬时幅度变化值
	float lastmagsum = 0;

	float dfreqThreshold = 0.5;//频率变化阈值
	float diffThreshold = 0.5;//周围差异阈值
	float dmagThreshold = 0.5;//幅度变化阈值
	float dmagoffset = 0.0;//幅度偏移量
protected:
	void ProcessSTFT(float* rev1, float* imv1, float* rev2, float* imv2, int numBins, int hopSize)
	{
		float osamp = (float)FFTSize / hopSize;
		float omega = 2.0 * M_PI * osamp * hopSize / FFTSize;
		int binfreq = 48000.0 / FFTSize;

		float magsum = 0.0;//总能量
		for (int j = 0; j < numBins; ++j)
		{
			float mag = sqrtf(rev1[j] * rev1[j] + imv1[j] * imv1[j]);
			magsum += mag;
		}
		for (int j = 0; j < numBins; ++j)
		{
			float mag = sqrtf(rev1[j] * rev1[j] + imv1[j] * imv1[j]);
			float phase = atan2f(imv1[j], rev1[j]);
			float freq = phase - lastphase[j];//对时间的导数是频率
			lastphase[j] = phase;

			freq -= j * omega;
			freq = fmod(freq + M_PI, -2.0 * M_PI) + M_PI;
			freq = osamp * freq / (2.0 * M_PI);
			freq = (long)j * binfreq + freq * binfreq * matchf;

			//瞬时频率
			dfreqs[j] = freq - lastfreq[j];
			lastfreq[j] = freq;

			//瞬时幅度
			//dmags[j] = mag / lastmagsum - lastmag[j] / magsum;
			//dmags[j] = mag - lastmag[j] / magsum;
			dmags[j] = mag;
			lastmag[j] = mag;
		}
		lastmagsum = magsum;

		for (int j = 1; j < numBins - 1; ++j)
		{
			float diff = (fabs(dfreqs[j] - dfreqs[j - 1]) + fabs(dfreqs[j + 1] - dfreqs[j])) / 4.0;//如果各bin频率变化相近，则认为是稳态
			float dfreq = fabs(dfreqs[j]);
			float dmag = fabs(dmags[j]);
			/*
			if (fabs(dfreqs[j]) < dfreqThreshold)//稳态
			{
				if (diff < diffThreshold)//如果周围也是稳态
				{
					rev2[j] = rev1[j];
					imv2[j] = imv1[j];
					rev1[j] = 0;
					imv1[j] = 0;
				}
			}*/
			//float mix = (diff * diffThreshold + dfreq * dfreqThreshold) / binfreq / 2.0;
			//if (mix < 1.0)//稳态
			
			if (dmag * dmagThreshold > magsum / numBins)//稳态
			{
				rev2[j] = rev1[j];
				imv2[j] = imv1[j];
				rev1[j] = 0;
				imv1[j] = 0;
			}
			/*
			float mix = ((dmag + dmagoffset) * dmagThreshold) / (magsum / numBins);//越大越稳态
			if (mix > 1.0)mix = 1.0;
			if (mix < 0.0)mix = 0.0;
			rev2[j] = rev1[j] * mix;
			imv2[j] = imv1[j] * mix;
			rev1[j] = rev1[j] * (1.0 - mix);
			imv1[j] = imv1[j] * (1.0 - mix);*/
		}
	}
public:
	TransientExtractor() {
		memset(inbuf, 0, sizeof(inbuf));
		memset(outbuf1, 0, sizeof(outbuf1));
		memset(outbuf2, 0, sizeof(outbuf2));
		memset(fftdatre1, 0, sizeof(fftdatre1));
		memset(fftdatim1, 0, sizeof(fftdatim1));
		memset(fftdatre2, 0, sizeof(fftdatre2));
		memset(fftdatim2, 0, sizeof(fftdatim2));
		memset(lastphase, 0, sizeof(lastphase));
		memset(lastfreq, 0, sizeof(lastfreq));
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

	void SetDFreqThreshold(float dfreqthreshold, float diffthreshold, float dmagthreshold, float dmagoffset)//瞬态检测阈值(0->1)
	{
		dfreqThreshold = dfreqthreshold * 100.0;
		diffThreshold = diffthreshold * 100.0;
		dmagThreshold = expf((dmagthreshold - 0.5) * 12.0);
		this->dmagoffset = (dmagoffset - 0.5) * 10.0;
	}

	void ProcessBlock(const float* in, float* out1, float* out2, float numSamples)
	{
		for (int i = 0; i < numSamples; ++i)
		{
			inbuf[pos] = in[i];
			pos++;
			if (pos >= FFTSize)pos = 0;
			if (pos % hopSize == 0)//填充完了！
			{
				for (int j = 0; j < windowSize; ++j)//灌
				{
					fftdatre1[j] = inbuf[(pos + j) % FFTSize] * window[j];
					fftdatim1[j] = 0;
					fftdatre2[j] = 0;
					fftdatim2[j] = 0;
				}
				for (int j = windowSize; j < FFTSize; ++j)//灌
				{
					fftdatre1[j] = 0;
					fftdatim1[j] = 0;
					fftdatre2[j] = 0;
					fftdatim2[j] = 0;
				}

				fft_f32(fftdatre1, fftdatim1, FFTSize, 1);//fft

				ProcessSTFT(fftdatre1, fftdatim1, fftdatre2, fftdatim2, FFTSize / 2, hopSize);

				for (int j = 0; j < FFTSize / 2; ++j)
				{
					fftdatre1[FFTSize - j - 1] = -fftdatre1[j];
					fftdatim1[FFTSize - j - 1] = fftdatim1[j];
					fftdatre2[FFTSize - j - 1] = -fftdatre2[j];
					fftdatim2[FFTSize - j - 1] = fftdatim2[j];
				}
				fft_f32(fftdatre1, fftdatim1, FFTSize, -1);//ifft
				fft_f32(fftdatre2, fftdatim2, FFTSize, -1);//ifft

				for (int j = 0; j < windowSize; ++j)//加窗
				{
					fftdatre1[j] *= window[j];
					fftdatre2[j] *= window[j];
				}
				int pos2 = posout;
				for (int j = 0; j < windowSize; ++j)
				{
					pos2++;
					outbuf1[pos2 % (FFTSize * 2)] += fftdatre1[j];
					outbuf2[pos2 % (FFTSize * 2)] += fftdatre2[j];
				}
				for (int j = 0; j < hopSize; ++j)
				{
					pos2++;
					outbuf1[pos2 % (FFTSize * 2)] = 0;
					outbuf2[pos2 % (FFTSize * 2)] = 0;
				}
			}

			out1[i] = outbuf1[posout] * normv;
			out2[i] = outbuf2[posout] * normv;
			posout++;
			if (posout >= FFTSize * 2)posout = 0;
		}
	}
};