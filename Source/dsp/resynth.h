#pragma once

#define _USE_MATH_DEFINES
#include <math.h>


struct QuadOsc//嫌慢可以打表cos
{
	float k1 = 0, k2 = 0;
	float u = 0, v = 0, a = 0;
	inline void SetTabel(float* tabel)
	{
	}
	inline void UpdataFreq(float freq)
	{
		k1 = tanf(freq * M_PI);
		k2 = 2.0 * k1 / (1.0 + k1 * k1);
		//选做
		float r = sqrtf(u * u + v * v);
		u /= r;
		v /= r;
	}
	inline void UpdataAmp(float amp)
	{
		a = amp;
	}
	inline float Process()
	{
		float tmp = u - k1 * v;
		v = v + k2 * tmp;
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

/*
struct QuadOsc
{
	float f, t, a;
	inline void SetTabel(float* tabel)
	{
	}
	inline void UpdataFreq(float freq)
	{
		f = freq;
	}
	inline void UpdataAmp(float amp)
	{
		a = amp;
	}
	inline float Process()
	{
		t += f;
		t -= (int)t;
		return cosf(t * 2.0 * M_PI) * a;
	}
	inline void ResetPhase()
	{
		t = 0;
	}

};*/
/*
struct QuadOsc
{

	unsigned int dt = 0, phase = 0;
	float a = 0;
	float* QuadOscCosTabel = NULL;
	inline void SetTabel(float* tabel)
	{
		QuadOscCosTabel = tabel;
	}
	inline void UpdataFreq(float freq)
	{
		dt = freq * 4294967296;
	}
	inline void UpdataAmp(float amp)
	{
		a = amp;
	}
	inline float Process()
	{
		phase += dt;
		return QuadOscCosTabel[phase >> 16] * a;
	}
	inline void ResetPhase()
	{
		phase = 0;
	}
};*/
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

	float QuadOscCosTabel[65536];

	float amp[FFTSize / 2] = { 0 };//幅度值
	float lastAmp[FFTSize / 2] = { 0 };//上次的幅度值
	float mixv = 0.0, mixdt = 0.0;
protected:
	float warp_pi(float x)
	{
		x /= 2.0 * M_PI;
		x -= roundf(x);
		return x * 2.0 * M_PI;
	}
	void ProcessSTFT(float* rev1, float* imv1, int numBins, int hopSize)
	{
		float osamp = (float)FFTSize / hopSize;
		float omega = 2.0 * M_PI * osamp * hopSize / FFTSize;
		int binfreq = 48000.0 / FFTSize;

		for (int j = 0; j < numBins; ++j)
		{
			lastAmp[j] = amp[j];
		}

		for (int j = 0; j < numBins; ++j)
		{
			float thetal = atan2f(rev1[j], imv1[j]);
			float mag = sqrtf(rev1[j] * rev1[j] + imv1[j] * imv1[j]);
			float dtl = thetal - lastphase[j];
			lastphase[j] = thetal;

			float wl = warp_pi(dtl + (float)j * M_PI * hopSize / FFTSize * 2) / hopSize;
			float fl = -wl * matchf / (2.0 * M_PI) * FFTSize + j;
			fl = fl / numBins * pitch / 2.0;

			amp[j] = 0;
			if (fl > 0.0 && fl < 1.0)
			{
				oscs[j].UpdataFreq(fl);
				//oscs[j].UpdataAmp(mag);
				amp[j] = mag;
			}
		}
		mixv = 0;
		mixdt = 1.0 / hopSize;
	}
public:
	Resynth() {
		memset(inbuf, 0, sizeof(inbuf));
		memset(fftdatre1, 0, sizeof(fftdatre1));
		memset(fftdatim1, 0, sizeof(fftdatim1));
		memset(lastphase, 0, sizeof(lastphase));
		UpdataWindow();
		for (int i = 0; i < 65536; ++i)
		{
			QuadOscCosTabel[i] = cosf((float)i * 2.0 * M_PI / 65536.0);
		}
		for (int i = 0; i < FFTSize / 2; ++i)
		{
			oscs[i].SetTabel(QuadOscCosTabel);
		}
		for (int i = 0; i < FFTSize / 2; ++i)
		{
			oscs[i].ResetPhase();
		}
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
		float unitFreq = 48000.0 / hopSize;
		float globalFreq = 440.0 * powf(2.0, 3.0 / 12.0);//1.0=C5
		float realPitch = pitch * matchf + (pitch / unitFreq * globalFreq) * (1.0 - matchf);
		this->pitch = realPitch;
		//this->pitch = pitch;
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

			//updata amp

			for (int j = 0; j < FFTSize / 2; ++j)
			{
				float ampv = lastAmp[j] * (1.0 - mixv) + amp[j] * mixv;
				oscs[j].UpdataAmp(ampv);
			}
			mixv += mixdt;

			//resynth
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