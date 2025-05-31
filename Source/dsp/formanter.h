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

	float window[MaxBufferSize];

	constexpr static int JumpPoint = 2;
	int GetMaxCorrIndex()
	{
		float maxv = 99999999999;
		int maxIdx = 0;
		for (int i = 0; i < range; i += JumpPoint)
		{
			float sum = 0;
			for (int j = 0; j < blockSize; j += 3)
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
			for (int j = 0; j < blockSize; j += 3)
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
	int lastlen;
	void UpdataWindow(int len)
	{
		if (lastlen == len)return;
		lastlen = len;
		for (int i = 0; i < len; ++i)
		{
			window[i] = 0.5 - 0.5 * cosf(2.0 * M_PI * i / len);
		}
	}
public:
	Formanter3()
	{
		memset(bufin, 0, sizeof(bufin));
		memset(bufout, 0, sizeof(bufout));
		memset(block, 0, sizeof(block));
	}
	void Reset()
	{
		//memset(bufin, 0, sizeof(bufin));
		memset(bufout, 0, sizeof(bufout));
		memset(block, 0, sizeof(block));
		//pos = 0;
		posHop = 0;
	}
	void FillBuffer(const float* in, int numSamples)
	{
		for (int i = 0; i < numSamples; ++i)
		{
			bufin[pos] = in[i];
			bufout[pos] = 0;
			pos = (pos + 1) % MaxBufferSize;
		}
	}
	void SetPitch(float freq, float blocksize = 1.0, float samplerate = 48000)
	{
		step = freq / samplerate;
		blockSize = (1.0 / step * 2.0 + 96) * blocksize;
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
		float len = blockSize / formant;//每个周期，块的实际长度
		if (len > blockSize)len = blockSize;
		UpdataWindow(len);
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

				for (int j = 0; j < len; ++j)
				{
					float w = window[j];
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
