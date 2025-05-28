#pragma once

#include "fft.h"
#include "elliptic-blep.h"

class BlockResample
{
private:
	signalsmith::blep::EllipticBlep<float> blep;
	float buf[8192];
public:
	BlockResample()
	{
		memset(buf, 0, sizeof(buf));
	}

	void Process1(const float* in, float* out, int numSamplesIn, int numSamplesOut, float stretch, float when)
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

	void Process2(const float* in, float* out, int numSamplesIn, int numSamplesOut, float stretch, float when)
	{
		float y0 = 0, y1 = 0, y2 = 0, y3 = 0, y4 = 0, y5 = 0;

		for (int j = 0; j < numSamplesIn; ++j)
		{
			buf[j + 3] = in[j];
		}
		buf[0] = buf[1] = buf[2] = 0;
		buf[numSamplesIn + 3] = buf[numSamplesIn + 4] = buf[numSamplesIn + 5] = 0;
		int len = (float)numSamplesIn * stretch;
		if (len > numSamplesOut)len = numSamplesOut;
		for (int j = 0; j < len; ++j)
		{
			float indexf = (float)j / stretch;
			int index = indexf;


			// 正常情况，取六个点
			y0 = buf[index + 0];
			y1 = buf[index + 1];
			y2 = buf[index + 2];
			y3 = buf[index + 3];
			y4 = buf[index + 4];
			y5 = buf[index + 5];


			// 计算五次Lagrange基函数
			float t = indexf - index;
			float t_plus_2 = t + 2;
			float t_plus_1 = t + 1;
			float t_minus_1 = t - 1;
			float t_minus_2 = t - 2;
			float t_minus_3 = t - 3;

			float L0 = (t_plus_1 * t * t_minus_1 * t_minus_2 * t_minus_3) / (-120.0f);
			float L1 = (t_plus_2 * t * t_minus_1 * t_minus_2 * t_minus_3) / 24.0f;
			float L2 = (t_plus_2 * t_plus_1 * t_minus_1 * t_minus_2 * t_minus_3) / (-12.0f);
			float L3 = (t_plus_2 * t_plus_1 * t * t_minus_2 * t_minus_3) / 12.0f;
			float L4 = (t_plus_2 * t_plus_1 * t * t_minus_1 * t_minus_3) / (-24.0f);
			float L5 = (t_plus_2 * t_plus_1 * t * t_minus_1 * t_minus_2) / 120.0f;

			out[j] = y0 * L0 + y1 * L1 + y2 * L2 + y3 * L3 + y4 * L4 + y5 * L5;
		}
		for (int j = len; j < numSamplesOut; ++j)
		{
			out[j] = 0;
		}
	}
	void Process(const float* in, float* out, int numSamplesIn, int numSamplesOut, float stretch, float when)
	{
		if (stretch < 1)Process1(in, out, numSamplesIn, numSamplesOut, stretch, when);
		else Process2(in, out, numSamplesIn, numSamplesOut, stretch, when);
		//Process2(in, out, numSamplesIn, numSamplesOut, stretch, when);
	}
};

class Formanter
{
private:
	constexpr static int MaxBufferInSize = 1024;//fft长度也是用这个
	constexpr static int MaxBufferOutSize = MaxBufferInSize * 2;//fft长度也是用这个
	constexpr static int MaxBlockSize = MaxBufferInSize / 3;
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

		resampler.Process(block, re4, MaxBlockSize, MaxBufferInSize, 1.0 / formant, 0);//重采样
		int len = (float)MaxBlockSize / formant;
		if (len > MaxBlockSize) len = MaxBlockSize;
		for (int i = 0; i < len; ++i)
		{
			float w = 0.5 - 0.5 * cosf(2.0 * M_PI * i / len);//加窗
			re4[i] *= w;
			im4[i] = 0;
		}
		for (int i = len; i < MaxBufferInSize * 2; ++i)
		{
			re4[i] = 0;
			im4[i] = 0;
		}

		fft_f32(re3, im3, MaxBufferInSize * 2, 1);
		fft_f32(re4, im4, MaxBufferInSize * 2, 1);
		for (int i = 0; i < MaxBufferInSize * 2; ++i)
		{
			float re3v = re3[i];
			float im3v = im3[i];
			float re4v = re4[i];
			float im4v = im4[i];
			re4[i] = re3v * re4v - im3v * im4v;
			im4[i] = re3v * im4v + im3v * re4v;
		}
		int fv = (1.0 - formant) * MaxBufferInSize;
		if (fv < 0) fv = 0;
		if (fv > MaxBufferInSize) fv = MaxBufferInSize;
		for (int i = MaxBufferInSize - fv; i < MaxBufferInSize + fv; ++i)
		{
			re4[i] = 0;
			im4[i] = 0;//滤除高频
		}
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

				for (int i = 0; i < MaxBufferInSize; ++i)
				{
					//re1[i] = re2[i] * buf_3[i];
					//im1[i] = im2[i] * buf_3[i];
					re1[i] = buf_3[i];
					im1[i] = buf_3[i];
				}
				//for (int i = 0; i < MaxBufferInSize / 2; ++i)
				//{
				//	re1[MaxBufferInSize - i - 1] = re1[i];
				//	im1[MaxBufferInSize - i - 1] = -im1[i];//共轭
				//}
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

class Formanter3 //psola
{
private:
	constexpr static int MaxBufferSize = 1024;
	constexpr static int MaxBlockSize = MaxBufferSize / 2;
	float bufin[MaxBufferSize];
	float bufout[MaxBufferSize];
	int pos = 0;
	float posHop = 0, step = 440.0 / 48000.0;

	float searchbuf[MaxBufferSize];
	float block[MaxBlockSize];

	int blockSize = MaxBlockSize, range = MaxBlockSize;
	int hopSize = MaxBufferSize / 2;//不固定

	int maxIndexLen = 0;//两次最大相关位置的差，用于预测下一次最大相关位置的区间。

	BlockResample resampler;
	float block2[MaxBlockSize];
	float formant = 1.0;

	constexpr static int JumpPoint = 2;
	int GetMaxCorrIndex()
	{
		float maxv = 99999999999;
		int maxIdx = 0;
		for (int i = 0; i < range; i += JumpPoint)
		{
			float sum = 0;
			for (int j = 0; j < blockSize; j += 4)
			{
				float a = searchbuf[i + j];
				float b = block[j];
				float c = a - b;
				sum += c * c;
			}
			if (sum < maxv)
			{
				maxv = sum;
				maxIdx = i;
			}
		}
		int start = maxIdx - JumpPoint + 1;
		int end = maxIdx + JumpPoint - 1;
		if (start < 0)start = 0;
		if (end > range)end = range;
		for (int i = start; i < end; i += 2)
		{
			float sum = 0;
			for (int j = 0; j < blockSize; j += 4)
			{
				float a = searchbuf[i + j];
				float b = block[j];
				float c = a - b;
				sum += c * c;
			}
			if (sum < maxv)
			{
				maxv = sum;
				maxIdx = i;
			}
		}

		return maxIdx;
	}
public:
	Formanter3()
	{
		memset(bufin, 0, sizeof(bufin));
		memset(bufout, 0, sizeof(bufout));
		memset(block, 0, sizeof(block));
	}
	void SetPitch(float freq)
	{
		step = freq / 48000.0;
		blockSize = 1.0 / step * 1.5 + 64;
		range = 1.0 / step * 32.0;
		if (blockSize > MaxBlockSize)blockSize = MaxBlockSize;
		if (range > MaxBlockSize)range = MaxBlockSize;
	}
	void SetFormant(float formant)
	{
		this->formant = formant;
	}
	void ProcessBlock(const float* in, float* out, int numSamples)
	{
		for (int i = 0; i < numSamples; ++i)
		{
			bufin[pos] = in[i];
			out[i] = bufout[pos];
			bufout[pos] = 0;
			pos = (pos + 1) % MaxBufferSize;
			posHop += step;
			if (posHop >= 1.0)
			{
				posHop -= (int)posHop;
				int hopSize = 1.0 / step;
				int start = (maxIndexLen + hopSize) % MaxBufferSize - range / 2;
				if (start < 0)start = 0;//0 -> MaxBufferSize-range-blockSize
				if (start > MaxBufferSize - range - blockSize)start = MaxBufferSize - range - blockSize;

				for (int j = 0; j < range + blockSize; ++j)
				{
					searchbuf[j] = bufin[(pos + j + start) % MaxBufferSize];
				}

				int index = GetMaxCorrIndex();
				maxIndexLen = (pos + start + index) % MaxBufferSize;

				for (int j = 0; j < blockSize; ++j)
				{
					block[j] = searchbuf[j + index];
				}

				resampler.Process2(block, block2, blockSize, blockSize, 1.0 / formant, 0);
				float len = blockSize / formant;
				if (len > blockSize)len = blockSize;
				for (int j = 0; j < len; ++j)
				{
					float w = 0.5 - 0.5 * cosf(2.0 * M_PI * j / len);
					block2[j] *= w;
				}

				for (int j = 0; j < blockSize; ++j)
				{
					bufout[(pos + j) % MaxBufferSize] += block2[j];
				}
			}
		}
	}
};

class Formanter4
{
private:
	constexpr static int MaxBufferSize = 2048;
	constexpr static int MaxBlockSize = MaxBufferSize / 2;
	float bufin[MaxBufferSize];
	float bufout[MaxBufferSize];
	int pos = 0;
	float posHop = 0, step = 440.0 / 48000.0;

	float re1[MaxBlockSize];
	float im1[MaxBlockSize];

	int blockSize = MaxBlockSize, range = MaxBlockSize;
	int hopSize = MaxBufferSize / 2;//不固定

public:
	Formanter4()
	{
		memset(bufin, 0, sizeof(bufin));
		memset(bufout, 0, sizeof(bufout));
		memset(re1, 0, sizeof(re1));
		memset(im1, 0, sizeof(im1));
	}
	void SetPitch(float freq)
	{
		step = freq / 48000.0;
		blockSize = 1.0 / step * 2.0;
		if (blockSize > MaxBlockSize)blockSize = MaxBlockSize;
	}
	void SetFormant(float formant)
	{

	}
	void ProcessBlock(const float* in, float* out, int numSamples)
	{
		for (int i = 0; i < numSamples; ++i)
		{
			bufin[pos] = in[i];
			out[i] = bufout[pos];
			bufout[pos] = 0;
			pos = (pos + 1) % MaxBufferSize;
			posHop += step;
			if (posHop >= 1.0)
			{
				posHop -= (int)posHop;
				//int FFTN = 2 * (1 << (int)(ceilf(log2f(blockSize))));
				//if (FFTN > MaxBlockSize)FFTN = MaxBlockSize;
				int FFTN = MaxBlockSize;
				for (int j = 0; j < FFTN; ++j)
				{
					float w = 0.5 - 0.5 * cosf(2.0 * M_PI * j / FFTN);
					re1[j] = bufin[(pos + j) % MaxBufferSize] * w;
					im1[j] = 0;
				}
				fft_f32(re1, im1, FFTN, 1);
				for (int j = 0; j < FFTN; ++j)
				{
					float r = sqrtf(re1[j] * re1[j] + im1[j] * im1[j]);
					re1[j] = r;
					im1[j] = 0;
				}
				fft_f32(re1, im1, FFTN, -1);
				for (int j = 0; j < blockSize; ++j)
				{
					float w = 0.5 - 0.5 * cosf(2.0 * M_PI * j / blockSize);
					bufout[(pos + j) % MaxBufferSize] += re1[j] * w / FFTN;
				}
			}
		}
	}
};