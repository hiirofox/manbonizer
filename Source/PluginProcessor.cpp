/*
  ==============================================================================

	This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"


//==============================================================================
LModelAudioProcessor::LModelAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
	: AudioProcessor(BusesProperties()
#if ! JucePlugin_IsMidiEffect
#if ! JucePlugin_IsSynth
		.withInput("Input", juce::AudioChannelSet::stereo(), true)
#endif
		.withOutput("Output", juce::AudioChannelSet::stereo(), true)
#endif
	)
#endif
{

}


juce::AudioProcessorValueTreeState::ParameterLayout LModelAudioProcessor::createParameterLayout()
{
	juce::AudioProcessorValueTreeState::ParameterLayout layout;

	layout.add(std::make_unique<juce::AudioParameterFloat>("alen", "alen", 0, 1, 0.5));
	layout.add(std::make_unique<juce::AudioParameterFloat>("pitch", "pitch", 0, 1, 0.5));
	layout.add(std::make_unique<juce::AudioParameterFloat>("matchf", "matchf", 0, 1, 1.0));
	layout.add(std::make_unique<juce::AudioParameterFloat>("format", "format", 0, 1, 0.5));

	return layout;
}

LModelAudioProcessor::~LModelAudioProcessor()
{
}

//==============================================================================
const juce::String LModelAudioProcessor::getName() const
{
	return JucePlugin_Name;
}

bool LModelAudioProcessor::acceptsMidi() const
{
#if JucePlugin_WantsMidiInput
	return true;
#else
	return false;
#endif
}

bool LModelAudioProcessor::producesMidi() const
{
#if JucePlugin_ProducesMidiOutput
	return true;
#else
	return false;
#endif
}

bool LModelAudioProcessor::isMidiEffect() const
{
#if JucePlugin_IsMidiEffect
	return true;
#else
	return false;
#endif
}

double LModelAudioProcessor::getTailLengthSeconds() const
{
	return 0.0;
}

int LModelAudioProcessor::getNumPrograms()
{
	return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
	// so this should be at least 1, even if you're not really implementing programs.
}

int LModelAudioProcessor::getCurrentProgram()
{
	return 0;
}

void LModelAudioProcessor::setCurrentProgram(int index)
{
}

const juce::String LModelAudioProcessor::getProgramName(int index)
{
	return "Virbon";
}

void LModelAudioProcessor::changeProgramName(int index, const juce::String& newName)
{

}

//==============================================================================
void LModelAudioProcessor::prepareToPlay(double sampleRate, int samplesPerBlock)
{
}

void LModelAudioProcessor::releaseResources()
{
	// When playback stops, you can use this as an opportunity to free up any
	// spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool LModelAudioProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
#if JucePlugin_IsMidiEffect
	juce::ignoreUnused(layouts);
	return true;
#else
	// This is the place where you check if the layout is supported.
	// In this template code we only support mono or stereo.
	// Some plugin hosts, such as certain GarageBand versions, will only
	// load plugins that support stereo bus layouts.
	if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
		&& layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
		return false;

	// This checks if the input layout matches the output layout
#if ! JucePlugin_IsSynth
	if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
		return false;
#endif

	return true;
#endif
}
#endif


void LModelAudioProcessor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midiMessages)
{
	int isMidiUpdata = 0;
	juce::MidiMessage MidiMsg;//先处理midi事件
	int MidiTime;
	juce::MidiBuffer::Iterator MidiBuf(midiMessages);
	while (MidiBuf.getNextEvent(MidiMsg, MidiTime))
	{
		if (MidiMsg.isNoteOn())
		{
			int note = MidiMsg.getNoteNumber() - 24;
		}
		if (MidiMsg.isNoteOff())
		{
			int note = MidiMsg.getNoteNumber() - 24;
		}
	}
	midiMessages.clear();

	const int numSamples = buffer.getNumSamples();
	float* wavbufl = buffer.getWritePointer(0);
	float* wavbufr = buffer.getWritePointer(1);
	const float* recbufl = buffer.getReadPointer(0);
	const float* recbufr = buffer.getReadPointer(1);

	float SampleRate = getSampleRate();

	float alen = *Params.getRawParameterValue("alen");
	float pitch = *Params.getRawParameterValue("pitch");
	float matchf = *Params.getRawParameterValue("matchf");
	float format = *Params.getRawParameterValue("format");

	//pvl.SetSliderLen(dft);
	//pvr.SetSliderLen(dft);
	//pvl.SetFormat(dmt);
	//pvr.SetFormat(dmt);
	//pvl.ProcessBlock(recbufl, buf1l, buf2l, numSamples);
	//pvr.ProcessBlock(recbufr, buf1r, buf2r, numSamples);

	//rsl.SetMatchF(tv);
	//rsr.SetMatchF(tv);
	//rsl.SetPitch(powf(2.0, (sv - 0.5) * 2.0 * 2.0));//4个8度
	//rsr.SetPitch(powf(2.0, (sv - 0.5) * 2.0 * 2.0));//4个8度
	//rsl.ProcessBlock(recbufl, buf1l, numSamples);
	//rsr.ProcessBlock(recbufr, buf1r, numSamples);

	//formsep->resynth->applyformant

	pvl.SetSliderLen(alen);
	pvr.SetSliderLen(alen);
	pvl.ProcessBlock(recbufl, buf1l, numSamples);//buf1:补偿了共振峰的信号
	pvr.ProcessBlock(recbufr, buf1r, numSamples);

	rsl.SetMatchF(matchf);
	rsr.SetMatchF(matchf);
	rsl.SetPitch(powf(2.0, (pitch - 0.5) * 2.0 * 2.0));//4个8度
	rsr.SetPitch(powf(2.0, (pitch - 0.5) * 2.0 * 2.0));//4个8度
	rsl.ProcessBlock(buf1l, buf2l, numSamples);//buf2:shift过后的信号
	rsr.ProcessBlock(buf1r, buf2r, numSamples);

	float* fdatl, * fdatr;
	afl.SetFormantLen(pvl.GetFormant(fdatl));
	afr.SetFormantLen(pvr.GetFormant(fdatr));
	afl.SetFormantData(fdatl);
	afr.SetFormantData(fdatr);

	afl.SetFormatValue(powf(2.0, (format - 0.5) * 2.0 * 2.0));
	afr.SetFormatValue(powf(2.0, (format - 0.5) * 2.0 * 2.0));
	afl.ProcessBlock(buf2l, buf1l, numSamples);//buf1:应用共振峰的信号
	afr.ProcessBlock(buf2r, buf1r, numSamples);

	for (int i = 0; i < numSamples; ++i)
	{
		wavbufl[i] = buf1l[i];
		wavbufr[i] = buf1r[i];
	}
}

//==============================================================================
bool LModelAudioProcessor::hasEditor() const
{
	return true; // (change this to false if you choose to not supply an editor)
}

juce::AudioProcessorEditor* LModelAudioProcessor::createEditor()
{
	return new LModelAudioProcessorEditor(*this);

	//return new juce::GenericAudioProcessorEditor(*this);//自动绘制(调试)
}

//==============================================================================

void LModelAudioProcessor::getStateInformation(juce::MemoryBlock& destData)
{
	// You should use this method to store your parameters in the memory block.
	// You could do that either as raw data, or use the XML or ValueTree classes
	// as intermediaries to make it easy to save and load complex data.

	// 创建一个 XML 节点
	juce::XmlElement xml("LMEQ_Settings");

	/*juce::MemoryBlock eqDataBlock;
	eqDataBlock.append(&manager, sizeof(ResonatorManager));
	juce::String base64Data = eqDataBlock.toBase64Encoding();
	xml.setAttribute("VIB_MANAGER", base64Data);//添加resonance数据
	*/
	auto state = Params.copyState();
	xml.setAttribute("Knob_Data", state.toXmlString());//添加旋钮数据

	juce::String xmlString = xml.toString();
	destData.append(xmlString.toRawUTF8(), xmlString.getNumBytesAsUTF8());
}

void LModelAudioProcessor::setStateInformation(const void* data, int sizeInBytes)
{
	// You should use this method to restore your parameters from this memory block,
	// whose contents will have been created by the getStateInformation() call.
	  // 将 data 转换为字符串以解析 XML
	juce::String xmlString(static_cast<const char*>(data), sizeInBytes);

	// 解析 XML
	std::unique_ptr<juce::XmlElement> xml(juce::XmlDocument::parse(xmlString));
	if (xml == nullptr || !xml->hasTagName("LMEQ_Settings"))
	{
		DBG("Error: Unable to load XML settings");
		return;
	}

	/*
	juce::String base64Data = xml->getStringAttribute("VIB_MANAGER");
	juce::MemoryBlock eqDataBlock;
	eqDataBlock.fromBase64Encoding(base64Data);
	if (eqDataBlock.getData() != NULL)
	{
		std::memcpy(&manager, eqDataBlock.getData(), sizeof(ResonatorManager));
	}
	*/
	auto KnobDataXML = xml->getStringAttribute("Knob_Data");
	Params.replaceState(juce::ValueTree::fromXml(KnobDataXML));
}

//==============================================================================
// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
	return new LModelAudioProcessor();
}
