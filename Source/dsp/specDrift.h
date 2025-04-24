#pragma once

#include "stft.h"

class SpecDrift : public STFT
{
private:
	constexpr static int MaxDelay = 800;//sample = MaxDelay*hopSize
	constexpr static int FFTSize = 1024;
	float bufre[MaxDelay][FFTSize / 2];
	float bufim[MaxDelay][FFTSize / 2];
	int pos = 0;
	float ldelay = 0, rdelay = 0, t0 = 0, fb = 0;
	void ProcessSTFT(float* re, float* im, int numBins, int hopSize) override
	{
		for (int i = 0; i < numBins; ++i)
		{
			// 计算当前频点的总延迟样本数
			float t = ldelay + (rdelay - ldelay) * ((float)i / numBins) + t0;

			// 分解为块延迟和相位延迟（基于hopSize）
			int t1 = (int)(t / hopSize);
			float t2 = t - t1 * hopSize; // 剩余样本数

			// 计算相位旋转量，分母使用FFTSize
			float phase = 2.0f * M_PI * t2 * i / FFTSize;
			float mulre = cosf(phase);
			float mulim = -sinf(phase);

			// 应用旋转并写入缓冲区
			int writePos = (pos + t1) % MaxDelay;
			float rev = re[i] + bufre[pos][i] * fb;
			float imv = im[i] + bufre[pos][i] * fb;
			bufre[writePos][i] = rev * mulre - imv * mulim;
			bufim[writePos][i] = rev * mulim + imv * mulre;

			// 读取当前pos的数据
			re[i] = bufre[pos][i];
			im[i] = bufim[pos][i];
		}
		pos = (pos + 1) % MaxDelay; // 更新环形缓冲区位置
	}
public:
	SpecDrift()
	{
		memset(bufre, 0, sizeof(bufre));
		memset(bufim, 0, sizeof(bufim));
	}
	void SetDelay(float ldelay, float rdelay, float t0, float fb)
	{
		int hopSize = GetHopSize();
		this->ldelay = ldelay * hopSize * MaxDelay;
		this->rdelay = rdelay * hopSize * MaxDelay;
		this->t0 = t0 * hopSize * MaxDelay;
		this->fb = fb;
	}
};
class SpecDriftPhaseVocoder : public STFT
{
private:
	constexpr static int MaxDelay = 800;//sample = MaxDelay*hopSize
	constexpr static int FFTSize = 1024;
	float ldelay = 0, rdelay = 0, t0 = 0, fb = 0;

	float lastphase[FFTSize / 2] = { 0 };
	float resynphase[FFTSize / 2] = { 0 };
	float resynmagn[MaxDelay][FFTSize / 2] = { 0 };
	float resynfreq[MaxDelay][FFTSize / 2] = { 0 };
	int pos = 0;

	inline float warp_pi(float x)
	{/*
		x /= 2.0 * M_PI;
		x -= roundf(x);
		return x * 2.0 * M_PI;*/
		return fmodf(x, 2.0 * M_PI);
	}

	void ProcessSTFT(float* re, float* im, int numBins, int hopSize) override
	{
		int hopSize_div = FFTSize / hopSize;
		float expct = 2.0 * M_PI * hopSize / FFTSize;
		float binfreq = 48000.0 / FFTSize;

		for (int j = 0; j < numBins; ++j)
		{
			resynmagn[pos][j] = 0;
			resynfreq[pos][j] = 0;
		}
		for (int j = 0; j < numBins; ++j) {
			float rev = re[j];
			float imv = im[j];
			float magn = 2.0 * sqrtf(rev * rev + imv * imv);
			float phase = atan2f(imv, rev);
			float delta = phase - lastphase[j];
			lastphase[j] = phase;
			delta -= (float)j * expct;
			delta = warp_pi(delta);
			delta = hopSize_div * delta / (2.0 * M_PI);
			delta = (float)j * binfreq + delta * binfreq;

			// 计算当前频点的总延迟样本数
			float t = ldelay + (rdelay - ldelay) * ((float)j / numBins) + t0;
			// 分解为块延迟和相位延迟（基于hopSize）
			int t1 = (int)(t / hopSize);
			float t2 = t - t1 * hopSize; // 剩余样本数

			int index = (pos + t1 + 1) % MaxDelay;
			resynmagn[index][j] += magn;
			resynfreq[index][j] = delta;


			//magn = resynmagn[pos][j];//resynth
			//float deltaphase = (resynfreq[pos][j] - binfreq * j) / binfreq;
			//magn = resynmagn[pos][j] * (t2)+resynmagn[(pos + 1) % MaxDelay][j] * (1.0 - t2);//线性插值
			magn = resynmagn[(pos + 1) % MaxDelay][j];//线性插值
			float deltaphase = (resynfreq[(pos + 1) % MaxDelay][j] - binfreq * j) / binfreq;
			deltaphase = 2.0 * M_PI * deltaphase / hopSize_div;

			re[j] = cosf(resynphase[j]) * magn;
			im[j] = sinf(resynphase[j]) * magn;

			deltaphase += expct * j;
			resynphase[j] = warp_pi(resynphase[j] + deltaphase);

		}
		pos = (pos + 1) % MaxDelay;
	}
public:
	SpecDriftPhaseVocoder()
	{
		memset(lastphase, 0, sizeof(lastphase));
		memset(resynmagn, 0, sizeof(resynmagn));
		memset(resynfreq, 0, sizeof(resynfreq));
		memset(resynphase, 0, sizeof(resynphase));
	}
	void SetDelay(float ldelay, float rdelay, float t0, float fb)
	{
		int hopSize = GetHopSize();
		this->ldelay = ldelay * hopSize * MaxDelay;
		this->rdelay = rdelay * hopSize * MaxDelay;
		this->t0 = t0 * hopSize * MaxDelay;
		this->fb = fb;
	}
};