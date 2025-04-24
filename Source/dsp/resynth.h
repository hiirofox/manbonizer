#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

struct QuadOsc//嫌慢可以打表cos
{
	float k1, k2;
	float u, v, a;
	inline void UpdataFreq(float freq)
	{
		k1 = tanf(freq * M_PI);
		k2 = 2.0 * k1 / (1.0 + k1 * k1);
		//选做
		//float r = sqrtf(p->u * p->u + p->v * p->v);
		//p->u /= r;
		//p->v /= r;
	}
	inline void UpdataAmp(float amp)
	{
		a = amp;
	}
	inline float Process()
	{
		float tmp = u - k1 * v;
		v += k2 * tmp;
		u = tmp - k1 * v;
		return u * a;
	}
	inline void ResetPhase()
	{
		u = 1.0;//cos
		v = 0.0;//sin
		a = 0.0;
	}
};



class Resynth
{

private:
	constexpr static int FFTSize = 1024;
	int windowSize = FFTSize;
	int hopSize = FFTSize / 8;

	float inbuf[FFTSize];
	float fftdatre1[FFTSize];
	float fftdatim1[FFTSize];
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

	//phase vocoder(加法合成)
	float pitch = 1.0, matchf = 1.0;
	float lastphase[FFTSize / 2] = { 0 };
	QuadOsc oscs[FFTSize / 2];
protected:
	void ProcessSTFT(float* rev1, float* imv1, int numBins, int hopSize)
	{
		float osamp = (float)FFTSize / hopSize;
		float omega = 2.0 * M_PI * osamp * hopSize / FFTSize;
		int binfreq = 48000.0 / FFTSize;

		for (int j = 0; j < numBins; ++j)
		{
			float mag = sqrtf(rev1[j] * rev1[j] + imv1[j] * imv1[j]);
			float phase = atan2f(imv1[j], rev1[j]);
			float freq = phase - lastphase[j];//对时间的导数是频率
			lastphase[j] = phase;

			freq -= j * omega;
			freq = fmod(freq + M_PI, -2.0 * M_PI) + M_PI;
			freq = osamp * freq / (2.0 * M_PI);
			freq = (long)j * binfreq + freq * binfreq * matchf;//瞬时频率
			if (j * pitch < numBins)
			{
				oscs[j].UpdataFreq(freq / 48000.0 * pitch);
				oscs[j].UpdataAmp(mag);
			}
		}
	}
public:
	Resynth() {
		memset(inbuf, 0, sizeof(inbuf));
		memset(fftdatre1, 0, sizeof(fftdatre1));
		memset(fftdatim1, 0, sizeof(fftdatim1));
		memset(lastphase, 0, sizeof(lastphase));
		for (int i = 0; i < FFTSize / 2; ++i)
		{
			oscs[i].ResetPhase();
		}
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

	void SetMatchF(float matchF)//是否需要相位信息
	{
		this->matchf = matchF;
	}
	void SetPitch(float pitch)//建议先设置好matchF再设置pitch
	{
		float unitFreq = 48000 / hopSize;
		float globalFreq = 440.0;//1.0=A4(440hz)
		float realPitch = pitch * matchf + (pitch / unitFreq * globalFreq) * (1.0 - matchf);
		this->pitch = realPitch;
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
				ProcessSTFT(fftdatre1, fftdatim1, FFTSize / 2, hopSize);//analize
			}

			float sum = 0;
			for (int j = 0; j < FFTSize / 2; ++j)
			{
				sum += oscs[j].Process();
			}
			sum /= FFTSize / 2;
			out1[i] = sum;
		}
	}
};