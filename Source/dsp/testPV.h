#pragma once

#include "stft.h"
#include <string.h>

class PhaseVocoder :public STFT
{
private:
	constexpr static int MaxDelay = 800;//sample = MaxDelay*hopSize
	constexpr static int FFTSize = 1024;


	float pitch = 1.0, matchf = 1.0;

	float lastphase[FFTSize / 2] = { 0 };
	float resynphase[FFTSize / 2] = { 0 };
	float resynmagn[FFTSize / 2] = { 0 };
	float resynfreq[FFTSize / 2] = { 0 };

	inline float warp_pi(float x)
	{/*
		x /= 2.0 * M_PI;
		x -= roundf(x);
		return x * 2.0 * M_PI;*/
		//return fmodf(x + M_PI, 2.0 * M_PI) - M_PI;
		return fmodf(x, 2.0 * M_PI);
	}

	void ProcessSTFT(float* rev, float* imv, int numBins, int hopSize) override
	{
		float osamp = (float)FFTSize / hopSize;
		float omega = 2.0 * M_PI * osamp * hopSize / FFTSize;
		int binfreq = 48000 / FFTSize;

		for (int j = 0; j < numBins; ++j)
		{
			resynmagn[j] = 0;
			resynfreq[j] = 0;
		}
		for (int j = 0; j < numBins; ++j)
		{
			float mag = sqrtf(rev[j] * rev[j] + imv[j] * imv[j]);
			float phase = atan2f(imv[j], rev[j]);
			float freq = phase - lastphase[j];//对时间的导数是频率
			lastphase[j] = phase;

			freq -= j * omega;
			freq = fmod(freq + M_PI, -2.0 * M_PI) + M_PI;
			freq = osamp * freq / (2.0 * M_PI);
			freq = (long)j * binfreq + freq * binfreq * matchf;

			int index = j * pitch;
			if (index < numBins)
			{
				resynfreq[index] = freq * pitch;
				resynmagn[index] += mag;
			}
		}
		for (int j = 0; j < numBins; ++j)
		{
			rev[j] = resynmagn[j] * cosf(resynphase[j]);
			imv[j] = resynmagn[j] * sinf(resynphase[j]);
			float tmp = resynfreq[j] - (float)j * binfreq;
			tmp /= binfreq;
			tmp = 2.0 * M_PI * tmp / (float)osamp;
			tmp += j * omega;
			resynphase[j] += tmp;

		}
	}
public:
	PhaseVocoder()
	{
		memset(lastphase, 0, sizeof(lastphase));
		memset(resynmagn, 0, sizeof(resynmagn));
		memset(resynfreq, 0, sizeof(resynfreq));
		memset(resynphase, 0, sizeof(resynphase));
	}
	void SetPitch(float pitch, float matchf)
	{
		this->pitch = pitch;
		this->matchf = matchf;
	}
};