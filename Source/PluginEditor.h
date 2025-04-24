/*
  ==============================================================================

	This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"
#include "ui/LM_slider.h"
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
	void timerCallback() override;//ע�⣬�����治���ڻ��fft�����

private:
	// This reference is provided as a quick way for your editor to
	// access the processor object that created it.
	LModelAudioProcessor& audioProcessor;

	LMKnob K_FQT, K_DFT, K_DMT, K_SV, K_TV;
	juce::ComponentBoundsConstrainer constrainer;  // �������ÿ�߱���
	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(LModelAudioProcessorEditor)
};

