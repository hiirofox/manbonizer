#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "fft.h"

class ApplyFormant//共振峰分离
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

	//apply formant
	float* formant = NULL;
	int numFormantBins = 0;

	float format = 1.0;

protected:
	void ProcessSTFT(float* rev1, float* imv1, int numBins, int hopSize)
	{
		if (formant)
		{
			for (int i = 0; i < numBins; ++i)
			{
				int index = format * i;
				if (index < numFormantBins)
				{
					rev1[i] *= formant[index];
					imv1[i] *= formant[index];
				}
				else
				{
					rev1[i] = 0;
					imv1[i] = 0;
				}
			}
		}
	}

public:
	ApplyFormant() {
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

	void SetFormatValue(float formatPitch)
	{
		this->format = formatPitch;
	}
	void SetFormantData(float* formantDataPtr)
	{
		this->formant = formantDataPtr;
	}

	void ProcessBlock(const float* in, float* out1, float numSamples)
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
				}
				for (int j = windowSize; j < FFTSize; ++j)//灌
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

				for (int j = 0; j < windowSize; ++j)//加窗
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
};