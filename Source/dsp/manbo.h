#pragma once

#include "formanter.h"

class Manbonizer
{
private:
	constexpr static int MaxPolyNum = 16;
	const float C4Freq = 440.0 * powf(2.0, (3.0 - 24.0) / 12.0);

	Formanter3 fml[MaxPolyNum];
	Formanter3 fmr[MaxPolyNum];
	int lastNote[MaxPolyNum];
	float noteV[MaxPolyNum];
	int noteState[MaxPolyNum];
	float noteVelo[MaxPolyNum];
	float notePan[MaxPolyNum];

	float corr, formant, glide, mix;

	float tmp1l[8192];
	float tmp1r[8192];
	float tmp2l[8192];
	float tmp2r[8192];
public:
	Manbonizer()
	{
		for (int i = 0; i < MaxPolyNum; ++i)
		{
			lastNote[i] = 0;
			noteState[i] = -1;
			noteVelo[i] = 0;
			notePan[i] = 0;
		}
		memset(tmp1l, 0, sizeof(tmp1l));
		memset(tmp1r, 0, sizeof(tmp1r));
		memset(tmp1l, 0, sizeof(tmp1l));
		memset(tmp1r, 0, sizeof(tmp1r));
	}
	void NoteOn(int note, float velo, float pan)
	{
		for (int i = 0; i < MaxPolyNum; ++i)
		{
			if (noteState[i] == -1)
			{
				noteState[i] = note;
				noteVelo[i] = velo;
				notePan[i] = pan;
				fml[i].Reset();
				fmr[i].Reset();
				return;
			}
		}
		noteState[0] = note;
		noteVelo[0] = velo;
		notePan[0] = pan;
	}
	void NoteOff(int note)
	{
		for (int i = 0; i < MaxPolyNum; ++i)
		{
			if (noteState[i] == note)
			{
				lastNote[i] = noteState[i];
				noteState[i] = -1;
				noteVelo[i] = 0;
				notePan[i] = 0;
				return;
			}
		}
	}
	void SetParam(float corr, float formant, float glide, float mix)
	{
		this->corr = corr;
		this->formant = formant;
		this->glide = glide;
		this->mix = mix;
	}
	void ProcessBlock(const float* inl, const float* inr, float* outl, float* outr, int numSamples)
	{
		for (int j = 0; j < numSamples; ++j)
		{
			tmp2l[j] = 0;
			tmp2r[j] = 0;
		}
		for (int i = 0; i < MaxPolyNum; ++i)
		{
			if (noteState[i] != -1)
			{
				float sign = 1;
				if (lastNote[i] <= noteState[i] && noteV[i] >= noteState[i])noteV[i] = noteState[i], sign = 0;
				else if (lastNote[i] <= noteState[i])sign = 1.0;
				if (lastNote[i] >= noteState[i] && noteV[i] <= noteState[i])noteV[i] = noteState[i], sign = 0;
				else if (lastNote[i] >= noteState[i])sign = -1.0;
				noteV[i] += glide * sign;

				float note = noteV[i];
				float freq = C4Freq * powf(2.0, (float)(note - 48) / 12.0);
				float fm = 1.0 * corr + (freq / C4Freq) * (1.0 - corr);
				fm *= formant;
				fml[i].SetPitch(freq);
				fmr[i].SetPitch(freq);
				fml[i].SetFormant(fm);
				fmr[i].SetFormant(fm);
				fml[i].ProcessBlock(inl, tmp1l, numSamples);
				fmr[i].ProcessBlock(inr, tmp1r, numSamples);

				float velo = noteVelo[i];
				float panl = cosf(notePan[i] * M_PI / 2.0);
				float panr = sinf(notePan[i] * M_PI / 2.0);
				if (panl < 0)panl = 0;
				if (panr < 0)panr = 0;
				panl *= velo;
				panr *= velo;
				for (int j = 0; j < numSamples; ++j)
				{
					tmp2l[j] += tmp1l[j] * panl;
					tmp2r[j] += tmp1r[j] * panr;
				}
			}
			else
			{
				fml[i].FillBuffer(inl, numSamples);//仅填充缓冲，没啥其他作用。
				fmr[i].FillBuffer(inr, numSamples);
			}
		}
		float dry = cosf(mix * M_PI / 2.0);
		float wet = sinf(mix * M_PI / 2.0);
		for (int j = 0; j < numSamples; ++j)
		{
			outl[j] = inl[j] * dry + tmp2l[j] * wet;
			outr[j] = inr[j] * dry + tmp2r[j] * wet;
		}
	}
};