/*
  ==============================================================================

	This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
LModelAudioProcessorEditor::LModelAudioProcessorEditor(LModelAudioProcessor& p)
	: AudioProcessorEditor(&p), audioProcessor(p)
{
	// Make sure that before the constructor has finished, you've set the
	// editor's size to whatever you need it to be.
	setResizable(true, true); // �����ڵ�����С

	setOpaque(false);  // �����ڱ߿��������

	//setResizeLimits(64 * 11, 64 * 5, 10000, 10000); // ������С���Ϊ300x200�������Ϊ800x600
	setSize(64 * 13, 64 * 4);
	setResizeLimits(64 * 4, 64 * 4, 64 * 13, 64 * 4);

	//constrainer.setFixedAspectRatio(11.0 / 4.0);  // ����Ϊ16:9����
	//setConstrainer(&constrainer);  // �󶨴��ڵĿ������

	K_ALEN.setText("AnalyzeLen");
	K_ALEN.ParamLink(audioProcessor.GetParams(), "alen");
	addAndMakeVisible(K_ALEN);
	K_PITCH.setText("Pitch");
	K_PITCH.ParamLink(audioProcessor.GetParams(), "pitch");
	addAndMakeVisible(K_PITCH);
	K_MATCHF.setText("MatchF");
	K_MATCHF.ParamLink(audioProcessor.GetParams(), "matchf");
	addAndMakeVisible(K_MATCHF);
	K_FORMAT.setText("Format");
	K_FORMAT.ParamLink(audioProcessor.GetParams(), "format");
	addAndMakeVisible(K_FORMAT);

	startTimerHz(30);

}

LModelAudioProcessorEditor::~LModelAudioProcessorEditor()
{
}

//==============================================================================
void LModelAudioProcessorEditor::paint(juce::Graphics& g)
{
	g.fillAll(juce::Colour(0x00, 0x00, 0x00));

	g.fillAll();
	g.setFont(juce::Font("FIXEDSYS", 17.0, 1));
	g.setColour(juce::Colour(0xff00ff00));;

	int w = getBounds().getWidth(), h = getBounds().getHeight();

	//g.drawText("L-MODEL Vibron", juce::Rectangle<float>(32, 16, w, 16), 1);
}

void LModelAudioProcessorEditor::resized()
{
	juce::Rectangle<int> bound = getBounds();
	int x = bound.getX(), y = bound.getY(), w = bound.getWidth(), h = bound.getHeight();
	auto convXY = juce::Rectangle<int>::leftTopRightBottom;

	K_ALEN.setBounds(32 + 64 * 0, 32 + 64 * 0, 64, 64);
	K_PITCH.setBounds(32 + 64 * 1, 32 + 64 * 0, 64, 64);
	K_MATCHF.setBounds(32 + 64 * 2, 32 + 64 * 0, 64, 64);
	K_FORMAT.setBounds(32 + 64 * 3, 32 + 64 * 0, 64, 64);
}

void LModelAudioProcessorEditor::timerCallback()
{
	repaint();
}
