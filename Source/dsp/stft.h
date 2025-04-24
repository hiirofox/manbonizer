#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include "fft.h"

class STFT
{
private:
	constexpr static int FFTSize = 1024;
	int windowSize = FFTSize;
	int hopSize = FFTSize / 32;

	float inbuf[FFTSize];
	float outbuf[FFTSize * 2];
	float fftdatre[FFTSize];
	float fftdatim[FFTSize];
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
protected:
	virtual void ProcessSTFT(float* re, float* im, int numBins, int hopSize) {};//继承它。处理一帧信号的频率成分
public:
	STFT() {
		memset(inbuf, 0, sizeof(inbuf));
		memset(outbuf, 0, sizeof(outbuf));
		memset(fftdatre, 0, sizeof(fftdatre));
		memset(fftdatim, 0, sizeof(fftdatim));
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
	void ProcessBlock(const float* in, float* out, float numSamples)
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
					fftdatre[j] = inbuf[(pos + j) % FFTSize] * window[j];
					fftdatim[j] = 0;
				}
				for (int j = windowSize; j < FFTSize; ++j)//灌
				{
					fftdatre[j] = 0;
					fftdatim[j] = 0;
				}

				fft_f32(fftdatre, fftdatim, FFTSize, 1);//fft

				ProcessSTFT(fftdatre, fftdatim, FFTSize / 2, hopSize);

				for (int j = 0; j < FFTSize / 2; ++j)
				{
					fftdatre[FFTSize - j - 1] = -fftdatre[j];
					fftdatim[FFTSize - j - 1] = fftdatim[j];
					//fftdatre[FFTSize - j - 1] = 0;
					//fftdatim[FFTSize - j - 1] = 0;
				}
				fft_f32(fftdatre, fftdatim, FFTSize, -1);//ifft

				for (int j = 0; j < windowSize; ++j)//加窗
				{
					fftdatre[j] *= window[j];
				}
				int pos2 = posout;
				for (int j = 0; j < windowSize; ++j)
				{
					pos2++;
					outbuf[pos2 % (FFTSize * 2)] += fftdatre[j];
				}
				for (int j = 0; j < hopSize; ++j)
				{
					pos2++;
					outbuf[pos2 % (FFTSize * 2)] = 0;
				}
			}

			out[i] = outbuf[posout] * normv;
			posout++;
			if (posout >= FFTSize * 2)posout = 0;
		}
	}
};
