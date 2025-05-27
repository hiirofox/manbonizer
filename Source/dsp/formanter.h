#pragma once

#include "fft.h"
#include "elliptic-blep.h"

class BlockResample
{
private:
	signalsmith::blep::EllipticBlep<float> blep;
public:
	void Process(const float* in, float* out, int numSamplesIn, int numSamplesOut, float stretch, float when)
	{
		blep.reset();
		float t = when - (int)when;

		int numOut = numSamplesIn * stretch;
		if (numOut > numSamplesOut)
		{
			numOut = numSamplesOut;
		}
		float k = 1.0 / stretch;
		int posIn = 0;
		for (int i = 0; i < numOut; ++i)
		{
			blep.step();
			t += k;
			while ((int)t)
			{
				t -= 1.0;
				float samplesInPast = t * stretch;
				float v = in[posIn++];
				blep.add(-v, 0, samplesInPast);
			}
			out[i] = blep.get();
		}
		for (int i = numOut; i < numSamplesOut; ++i)
		{
			out[i] = 0.0f;
		}
	}
};

class Formanter
{
private:
	constexpr static int MaxBufferInSize = 1024;//fft长度也是用这个
	constexpr static int MaxBufferOutSize = MaxBufferInSize * 2;//fft长度也是用这个
	constexpr static int MaxBlockSize = MaxBufferInSize / 2;
	float block[MaxBlockSize];
	float bufin[MaxBufferInSize];
	int posin = 0;
	float bufout[MaxBufferOutSize];
	int posout = 0;

	float pitch = 440.0 / 48000.0, formant = 1.0;


	float re1[MaxBufferInSize];
	float im1[MaxBufferInSize];
	float re2[MaxBufferInSize];
	float im2[MaxBufferInSize];
	int GetMaxCorrIndex()//找最大相关(fft)
	{
		int len1 = MaxBufferInSize;
		for (int i = 0; i < len1; ++i) {
			//float w = window2[i];
			//re1[i] = copybuf[i] * w;
			re1[i] = bufin[i];
			im1[i] = 0.0f;
		}
		for (int i = 0; i < MaxBlockSize; ++i) {
			re2[i] = block[i];
			im2[i] = 0.0f;
		}
		for (int i = MaxBlockSize; i < MaxBufferInSize; ++i) {
			re2[i] = 0.0f;
			im2[i] = 0.0f;
		}
		fft_f32(re1, im1, MaxBufferInSize, 1);
		fft_f32(re2, im2, MaxBufferInSize, 1);
		for (int i = 0; i < MaxBufferInSize; ++i) {
			float re1v = re1[i];
			float im1v = im1[i];
			float re2v = re2[i];
			float im2v = -im2[i];//这个得共轭
			re1[i] = re1v * re2v - im1v * im2v;
			im1[i] = re1v * im2v + im1v * re2v;
		}
		fft_f32(re1, im1, MaxBufferInSize, -1);
		int index = 0;
		float max = -9999999999;
		for (int i = 0; i < MaxBlockSize; ++i) {
			float r = re1[i];
			if (r > max) {
				max = r;
				index = i;
			}
		}
		return index;
	}

	signalsmith::blep::EllipticBlep<float> blit;//用来产生无混叠冲激串
	BlockResample resampler;//用来进行重采样
	float step = 0.0;//step += pitch; (per sample)
	float re3[MaxBufferInSize * 2];
	float im3[MaxBufferInSize * 2];
	float re4[MaxBufferInSize * 2];
	float im4[MaxBufferInSize * 2];
	void Process()
	{
		int index = GetMaxCorrIndex();
		for (int i = 0; i < MaxBlockSize; ++i)
		{
			block[i] = bufin[i + index];
		}
		for (int i = 0; i < MaxBufferInSize; ++i)
		{
			blit.step();
			step += pitch;
			float sum = 0;
			if (step >= 1.0)
			{
				step -= (int)step;
				blit.add(-1.0, 0, step / pitch);
				sum += 1.0;
			}
			//float w = 0.5 - 0.5 * cosf(2.0 * M_PI * i / MaxBufferInSize);//加窗
			float w = 1.0;
			re3[i] = blit.get() * w;
			//re3[i] = sum * w;
			im3[i] = 0;
		}
		for (int i = MaxBufferInSize; i < MaxBufferInSize * 2; ++i)
		{
			re3[i] = 0;
			im3[i] = 0;
		}

		for (int i = 0; i < MaxBlockSize; ++i)
		{
			float w = 0.5 - 0.5 * cosf(2.0 * M_PI * i / MaxBlockSize);//加窗
			//float w = 1.0;
			im4[i] = block[i] * w;
		}
		resampler.Process(im4, re4, MaxBlockSize, MaxBufferInSize, 1.0 / formant, 0);//重采样
		for (int i = 0; i < MaxBlockSize; ++i)
		{
			im4[i] = 0;
		}
		for (int i = MaxBlockSize; i < MaxBufferInSize * 2; ++i)
		{
			re4[i] = 0;
			im4[i] = 0;
		}

		fft_f32(re3, im3, MaxBufferInSize * 2, 1);
		fft_f32(re4, im4, MaxBufferInSize * 2, 1);
		for (int i = 1; i < MaxBufferInSize * 2; ++i)
		{
			float re3v = re3[i];
			float im3v = im3[i];
			float re4v = re4[i];
			float im4v = im4[i];
			re4[i] = re3v * re4v - im3v * im4v;
			im4[i] = re3v * im4v + im3v * re4v;
		}
		re3[0] = 0;
		im3[0] = 0;
		re4[0] = 0;
		im4[0] = 0;
		fft_f32(re4, im4, MaxBufferInSize * 2, -1);

		for (int i = 0; i < MaxBufferInSize * 2; ++i)
		{
			re4[i] /= (float)(MaxBufferInSize * 2.0);
			//re4[i] *= 0.5 - 0.5 * cosf(2.0 * M_PI * i / (MaxBufferInSize * 2.0));//加窗
		}
		int pos = posout;
		for (int i = 0; i < MaxBufferInSize * 2; ++i)
		{
			bufout[pos++] += re4[i];
			if (pos >= MaxBufferOutSize)
			{
				pos = 0;
			}
		}
	}
public:
	Formanter()
	{
		memset(bufin, 0, sizeof(bufin));
		memset(bufout, 0, sizeof(bufout));
		memset(block, 0, sizeof(block));
		memset(re1, 0, sizeof(re1));
		memset(im1, 0, sizeof(im1));
		memset(re2, 0, sizeof(re2));
		memset(im2, 0, sizeof(im2));
		memset(re3, 0, sizeof(re3));
		memset(im3, 0, sizeof(im3));
		memset(re4, 0, sizeof(re4));
		memset(im4, 0, sizeof(im4));
	}
	void SetPitch(float freq)
	{
		pitch = freq / 48000.0f;
	}
	void SetFormant(float formant)
	{
		this->formant = formant;
	}
	void ProcessBlock(const float* in, float* out, int numSamples)
	{
		for (int i = 0; i < numSamples; ++i)
		{
			bufin[posin] = in[i];
			out[i] = bufout[posout];
			bufout[posout] = 0;

			posin++;
			if (posin >= MaxBufferInSize)
			{
				posin = 0;
				Process();
			}
			posout++;
			if (posout >= MaxBufferOutSize)
			{
				posout = 0;
			}
		}
	}
};

class Formanter_FFT
{
private:
	constexpr static int MaxBufferInSize = 1024;//fft长度也是用这个
	constexpr static int MaxBufferOutSize = MaxBufferInSize * 2;//fft长度也是用这个
	int hopSize = MaxBufferInSize / 2;
	int posHop = 0;

	float bufin[MaxBufferInSize];//re
	float bufout[MaxBufferOutSize];

	float re1[MaxBufferInSize];
	float im1[MaxBufferInSize];
	float re2[MaxBufferInSize];
	float im2[MaxBufferInSize];
	float buf_3[MaxBufferInSize];//共振峰拉伸用的

	int posin = 0;
	int posout = 0;
	float pitch = 1.0, formant = 1.0;

	signalsmith::blep::EllipticBlep<float> blit;//用来产生无混叠冲激串
	float step = 0.0;
public:
	Formanter_FFT()
	{
		memset(bufin, 0, sizeof(bufin));
		memset(bufout, 0, sizeof(bufout));
		memset(re1, 0, sizeof(re1));
		memset(im1, 0, sizeof(im1));
		memset(re2, 0, sizeof(re2));
		memset(im2, 0, sizeof(im2));
		memset(buf_3, 0, sizeof(buf_3));
	}
	void SetPitch(float freq)
	{
		pitch = freq / 48000.0;
	}
	void SetFormant(float formant)
	{
		this->formant = formant;
	}
	void ProcessBlock(const float* in, float* out, int numSamples)
	{
		for (int i = 0; i < numSamples; ++i)
		{
			bufin[posin] = in[i];
			out[i] = bufout[posout];
			bufout[posout] = 0;
			posin++, posHop++;
			posout++;
			if (posin >= MaxBufferInSize) posin = 0;
			if (posHop >= hopSize)
			{
				posHop = 0;
				for (int i = 0; i < MaxBufferInSize; ++i)
				{
					float w = 0.5f - 0.5f * cosf(2.0f * M_PI * i / MaxBufferInSize);//加窗
					re1[i] = bufin[(posin + i) % MaxBufferInSize] * w;
					im1[i] = 0.0f;
				}
				fft_f32(re1, im1, MaxBufferInSize, 1);

				for (int i = 0; i < MaxBufferInSize / 2; ++i)//丢掉相位信息
				{
					float r = sqrtf(re1[i] * re1[i] + im1[i] * im1[i]);
					re1[i] = r;
					im1[i] = 0;
				}

				int len = (float)MaxBufferInSize / 2 * formant;
				if (len > MaxBufferInSize / 2)len = MaxBufferInSize / 2;
				for (int i = 0; i < len; ++i)//共振峰拉伸
				{
					float x = (float)i / formant;
					int index = x;
					float mix = x - index;
					buf_3[i] = re1[index] * (1.0 - mix) + re1[index + 1];
				}
				for (int i = len; i < MaxBufferInSize / 2; ++i)
				{
					buf_3[i] = 0.0f;
				}
				////
				blit.reset();
				float t = step;
				for (int i = 0; i < MaxBufferInSize; ++i)
				{
					blit.step();
					t += pitch;
					if (t >= 1.0f)
					{
						t -= (int)t;
						blit.add(-1.0f, 0, t / pitch);
					}
					float w = 0.5f - 0.5f * cosf(2.0f * M_PI * i / MaxBufferInSize);//加窗
					re2[i] = blit.get() * w;
					im2[i] = 0.0f;
				}
				step += pitch * hopSize;
				step -= (int)step;
				fft_f32(re2, im2, MaxBufferInSize, 1);

				for (int i = 0; i < MaxBufferInSize / 2; ++i)
				{
					re1[i] = re2[i] * buf_3[i];
					im1[i] = im2[i] * buf_3[i];
				}
				for (int i = 0; i < MaxBufferInSize / 2; ++i)
				{
					re1[MaxBufferInSize - i - 1] = re1[i];
					im1[MaxBufferInSize - i - 1] = -im1[i];//共轭
				}
				fft_f32(re1, im1, MaxBufferInSize, -1);
				for (int i = 0; i < MaxBufferInSize; ++i)
				{
					float w = 0.5f - 0.5f * cosf(2.0f * M_PI * i / MaxBufferInSize);//加窗
					re1[i] *= w / MaxBufferInSize;
				}

				posout = 0;
				for (int i = 0; i < MaxBufferInSize; ++i)
				{
					bufout[i] = bufout[i + hopSize];
				}
				for (int i = MaxBufferInSize - hopSize; i < MaxBufferOutSize; ++i)
				{
					bufout[i] = 0.0f;
				}
				for (int i = 0; i < MaxBufferInSize; ++i)
				{
					bufout[i] += re1[i];
				}
			}
		}
	}
};