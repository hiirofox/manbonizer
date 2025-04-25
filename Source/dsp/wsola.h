#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

class Wsola
{
private:
	constexpr static int BufSize = 16384;//应该留有足够大的空间来寻找最大自相关
	constexpr static int WindowSize = 32;//窗长
	float buf[BufSize] = { 0 };
	float bufOut[BufSize] = { 0 };//输出缓冲区
	int pos = 0;

	float maxcorr = 0;//最大自相关值
	int maxpos = 0;//最大自相关位置
	int lastmaxpos = 0;//上次最大自相关位置
	int findpos = 0;//寻找最大自相关位置
	//两次自相关位置可以大致推断出音高，进而实现pitch shift

	float smallblock[WindowSize] = { 0 };//小块

	float phase = 0;

	float pitch = 440.0 / 48000.0;
	float matchf = 1.0;//匹配频率(类似phase vocoder的保留/丢弃相位信息)
	float format = 1.0;//音色(拉伸块长度实现拉伸共振峰)
public:
	Wsola()
	{
		memset(buf, 0, sizeof(buf));
		memset(smallblock, 0, sizeof(smallblock));
	}
	float ProcessSample(float in)
	{
		buf[pos] = in;//填充缓冲区
		pos++;
		if (pos >= BufSize)pos = 0;

		int start = (BufSize + pos - 1 - WindowSize) % BufSize;//求一次自相关
		float corr = 0;
		for (int i = 0; i < WindowSize; ++i)
		{
			corr += smallblock[i] * buf[(start + i) % BufSize];
		}

		if (corr > maxcorr)//更新最大自相关
		{
			maxcorr = corr;
			maxpos = start;
		}

		phase += pitch;
		if (phase >= 1.0)//此时更新smallblock
		{
			phase -= (int)phase;
			for (int i = 0; i < WindowSize; ++i)
			{
				float window = (1.0 - cosf(2.0 * M_PI * i / WindowSize)) * 0.5;
				smallblock[i] = buf[(maxpos + i) % BufSize] * window;
			}
			maxcorr = 0;
			lastmaxpos = maxpos;
			maxpos = start;

			//将smallblock叠加进输出缓冲区
			for (int i = WindowSize / 2; i < WindowSize; i++)
			{
				bufOut[(lastmaxpos + i) % BufSize] = 0;
			}
			for (int i = 0; i < WindowSize; i++)
			{
				bufOut[(lastmaxpos + i) % BufSize] += smallblock[i];
			}
		}

		return bufOut[pos];
	}
	float ProcessBlock(const float* in, float* out, int numSamples)
	{
		for (int i = 0; i < numSamples; ++i)
		{
			out[i] = ProcessSample(in[i]);
		}
		return 0;
	}

	void SetPitch(float pitch)
	{
		this->pitch = pitch;
	}
	void SetMatchF(float matchf)
	{
		this->matchf = matchf;
	}
	void SetFormat(float format)
	{
		this->format = format;
	}
};