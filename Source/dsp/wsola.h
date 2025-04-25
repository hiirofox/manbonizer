#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

class Wsola
{
private:
	constexpr static int BufSize = 16384;//Ӧ�������㹻��Ŀռ���Ѱ����������
	constexpr static int WindowSize = 32;//����
	float buf[BufSize] = { 0 };
	float bufOut[BufSize] = { 0 };//���������
	int pos = 0;

	float maxcorr = 0;//��������ֵ
	int maxpos = 0;//��������λ��
	int lastmaxpos = 0;//�ϴ���������λ��
	int findpos = 0;//Ѱ����������λ��
	//���������λ�ÿ��Դ����ƶϳ����ߣ�����ʵ��pitch shift

	float smallblock[WindowSize] = { 0 };//С��

	float phase = 0;

	float pitch = 440.0 / 48000.0;
	float matchf = 1.0;//ƥ��Ƶ��(����phase vocoder�ı���/������λ��Ϣ)
	float format = 1.0;//��ɫ(����鳤��ʵ�����칲���)
public:
	Wsola()
	{
		memset(buf, 0, sizeof(buf));
		memset(smallblock, 0, sizeof(smallblock));
	}
	float ProcessSample(float in)
	{
		buf[pos] = in;//��仺����
		pos++;
		if (pos >= BufSize)pos = 0;

		int start = (BufSize + pos - 1 - WindowSize) % BufSize;//��һ�������
		float corr = 0;
		for (int i = 0; i < WindowSize; ++i)
		{
			corr += smallblock[i] * buf[(start + i) % BufSize];
		}

		if (corr > maxcorr)//������������
		{
			maxcorr = corr;
			maxpos = start;
		}

		phase += pitch;
		if (phase >= 1.0)//��ʱ����smallblock
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

			//��smallblock���ӽ����������
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