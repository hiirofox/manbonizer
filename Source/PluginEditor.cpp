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
	setResizable(true, true); // 允许窗口调整大小

	setOpaque(false);  // 允许在边框外面绘制

	//setResizeLimits(64 * 11, 64 * 5, 10000, 10000); // 设置最小宽高为300x200，最大宽高为800x600
	setSize(64 * 9, 64 * 3);
	setResizeLimits(64 * 9, 64 * 3, 64 * 13, 64 * 3);

	//constrainer.setFixedAspectRatio(11.0 / 4.0);  // 设置为16:9比例
	//setConstrainer(&constrainer);  // 绑定窗口的宽高限制

	K_Corr.setText("corr");
	K_Corr.ParamLink(audioProcessor.GetParams(), "corr");
	addAndMakeVisible(K_Corr);
	K_Formant.setText("formant");
	K_Formant.ParamLink(audioProcessor.GetParams(), "formant");
	addAndMakeVisible(K_Formant);
	K_Glide.setText("glide");
	K_Glide.ParamLink(audioProcessor.GetParams(), "glide");
	addAndMakeVisible(K_Glide);
	K_Mix.setText("mix");
	K_Mix.ParamLink(audioProcessor.GetParams(), "mix");
	addAndMakeVisible(K_Mix);
	K_Randpan.setText("randpan");
	K_Randpan.ParamLink(audioProcessor.GetParams(), "randpan");
	addAndMakeVisible(K_Randpan);
	K_Keep.setText("keepForm");
	K_Keep.ParamLink(audioProcessor.GetParams(), "keep");
	addAndMakeVisible(K_Keep);

	sp.SetBuffer(p.outbuf, &p.pos, p.OutBufferLen);
	addAndMakeVisible(sp);
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

	g.drawText("L-MODEL Manbonizer", juce::Rectangle<float>(32, 16, w, 16), 1);
}

void LModelAudioProcessorEditor::resized()
{
	juce::Rectangle<int> bound = getBounds();
	int x = bound.getX(), y = bound.getY(), w = bound.getWidth(), h = bound.getHeight();
	auto convXY = juce::Rectangle<int>::leftTopRightBottom;

	K_Corr.setBounds(w - 32 - 64 * 3, 32 + 64 * 0, 64, 64);
	K_Formant.setBounds(w - 32 - 64 * 2, 32 + 64 * 0, 64, 64);
	K_Keep.setBounds(w - 32 - 64 * 1, 32 + 64 * 0, 64, 64);
	K_Randpan.setBounds(w - 32 - 64 * 3, 32 + 64 * 1, 64, 64);
	K_Glide.setBounds(w - 32 - 64 * 2, 32 + 64 * 1, 64, 64);
	K_Mix.setBounds(w - 32 - 64 * 1, 32 + 64 * 1, 64, 64);

	sp.setBounds(32, 32, w - 64 * 4 - 32, 64 * 2);
}

void LModelAudioProcessorEditor::timerCallback()
{
	repaint();
}
