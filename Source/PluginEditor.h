/*
  ==============================================================================

	This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"
#include "ui/LM_slider.h"
#include "ui/specturmUI.h"
//==============================================================================
/**
*/
class LModelAudioProcessorEditor : public juce::AudioProcessorEditor, juce::Timer
{
public:
	LModelAudioProcessorEditor(LModelAudioProcessor&);
	~LModelAudioProcessorEditor() override;

	//==============================================================================
	void paint(juce::Graphics&) override;
	void resized() override;
	void timerCallback() override;//注意，这里面不定期会给fft填东西。

private:
	// This reference is provided as a quick way for your editor to
	// access the processor object that created it.
	LModelAudioProcessor& audioProcessor;

	LMKnob K_Corr;
	LMKnob K_Formant;
	LMKnob K_Glide;
	LMKnob K_Mix;
	LMKnob K_Randpan;
	LMKnob K_Keep;

	SpectrumUI sp;

	juce::ComponentBoundsConstrainer constrainer;  // 用于设置宽高比例
	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(LModelAudioProcessorEditor)
};

