#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

/*
struct QuadOsc//�������Դ��cos
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
		//ѡ��
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
*/
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
struct QuadOsc//���Ǵ���
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

struct QuadOscs//�������Դ��cos
{
	constexpr static int NumOscs = 1024;
	float k1[NumOscs] = { 0 }, k2[NumOscs] = { 0 };
	float u[NumOscs] = { 0 }, v[NumOscs] = { 0 }, a[NumOscs] = { 0 };
	inline void SetTabel(float* tabel)
	{
	}
	inline void UpdataFreq(int bin, float freq)
	{
		/*
		k1 = tanf(freq * M_PI);
		k2 = 2.0 * k1 / (1.0 + k1 * k1);
		//ѡ��
		float r = sqrtf(u * u + v * v);
		u /= r;
		v /= r;*/
		k1[bin] = tanf(freq * M_PI);
		k2[bin] = 2.0 * k1[bin] / (1.0 + k1[bin] * k1[bin]);
		//ѡ��
		float r = sqrtf(u[bin] * u[bin] + v[bin] * v[bin]);
		u[bin] /= r;
		v[bin] /= r;
	}
	inline void UpdataAmp(int bin, float amp)
	{
		a[bin] = amp;
	}
	inline float Process(int numBins)
	{/*
		float tmp = u - k1 * v;
		v = v + k2 * tmp;
		u = tmp - k1 * v;
		return u * a;*/
		float sum = 0;
		for (int j = 0; j < numBins; ++j)
		{
			float tmp = u[j] - k1[j] * v[j];
			v[j] = v[j] + k2[j] * tmp;
			u[j] = tmp - k1[j] * v[j];
			sum += u[j] * a[j];
		}
		return sum;
	}
	inline void ResetPhase()
	{/*
		u = 1.0;//cos
		v = 0.0;//sin
		a = 0.0;*/
		for (int j = 0; j < NumOscs; ++j)
		{
			u[j] = 1.0;//cos
			v[j] = 0.0;//sin
			a[j] = 0.0;
		}
	}
};

struct TableQuadOscs//���Ǵ���
{
	constexpr static int NumOscs = 1024;
	unsigned int dt[NumOscs] = { 0 }, phase[NumOscs] = { 0 };
	float a[NumOscs] = { 0 };
	float QuadOscCosTabel[65536] = { 0 };
	inline void InitTabel()
	{
		for (int i = 0; i < 65536; ++i)
		{
			QuadOscCosTabel[i] = cosf((float)i * 2.0 * M_PI / 65536.0);
		}
	}
	inline void UpdataFreq(int bin, float freq)
	{
		dt[bin] = freq * 4294967296;
	}
	inline void UpdataAmp(int bin, float amp)
	{
		a[bin] = amp;
	}
	inline float Process(int numBins)
	{
		float sum = 0;
		for (int i = 0; i < numBins; ++i)
		{
			phase[i] += dt[i];
			sum += QuadOscCosTabel[phase[i] >> 16] * a[i];
		}
		return sum;
	}
	inline void ResetPhase()
	{
		for (int i = 0; i < NumOscs; ++i)
		{
			phase[i] = rand() % 1000;
		}
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
	int pos = 0, hopPos = 0;

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

	//phase vocoder(�ӷ��ϳ�)
	float pitch = 1.0, matchf = 1.0;
	float lastphase[FFTSize / 2] = { 0 };
	//QuadOsc oscs[FFTSize / 2];
	TableQuadOscs oscs;

	//float QuadOscCosTabel[65536];//��

	float amp[FFTSize / 2] = { 0 };//����ֵ
	float lastAmp[FFTSize / 2] = { 0 };//�ϴεķ���ֵ
	float mixv = 0.0, mixdt = 0.0;
protected:
	inline float warp_pi(float x)
	{
		x /= 2.0 * M_PI;
		x -= roundf(x);
		return x * 2.0 * M_PI;
	}
	void ProcessSTFT(float* rev1, float* imv1, int numBins, int hopSize)
	{
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
				//oscs[j].UpdataFreq(fl);
				oscs.UpdataFreq(j, fl);
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
		/*for (int i = 0; i < 65536; ++i)//��
		{
			QuadOscCosTabel[i] = cosf((float)i * 2.0 * M_PI / 65536.0);
		}*/
		/*
		for (int i = 0; i < FFTSize / 2; ++i)//��
		{
			oscs[i].SetTabel(QuadOscCosTabel);
		}*/
		/*
		for (int i = 0; i < FFTSize / 2; ++i)
		{
			oscs[i].ResetPhase();
		}*/
		oscs.InitTabel();
		oscs.ResetPhase();
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

	void SetMatchF(float matchF)//ƥ����λ��Ϣ
	{
		this->matchf = matchF;
	}
	void SetPitch(float pitch)//���������ú�matchF������pitch
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
			hopPos++;
			if (pos >= FFTSize)pos = 0;
			if (hopPos >= hopSize)//������ˣ�
			{
				hopPos = 0;
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
				ProcessSTFT(fftdatre1, fftdatim1, FFTSize / 2, hopSize);//analize
			}

			//updata amp

			for (int j = 0; j < FFTSize / 2; ++j)
			{
				float ampv = lastAmp[j] * (1.0 - mixv) + amp[j] * mixv;
				//oscs[j].UpdataAmp(ampv);
				oscs.UpdataAmp(j, ampv);
			}
			mixv += mixdt;

			//resynth

			float sum = 0;
			/*
			for (int j = 0; j < FFTSize / 2; ++j)
			{
				sum += oscs[j].Process();
			}*/
			sum = oscs.Process(FFTSize / 2);
			sum /= FFTSize / 2;
			out1[i] = sum;
		}
	}
};