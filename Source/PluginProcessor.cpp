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

	layout.add(std::make_unique<juce::AudioParameterFloat>("corr", "corr", 0, 1, 1));
	layout.add(std::make_unique<juce::AudioParameterFloat>("formant", "formant", -1, 1, 0));
	layout.add(std::make_unique<juce::AudioParameterFloat>("glide", "glide", 0, 1, 0));
	layout.add(std::make_unique<juce::AudioParameterFloat>("mix", "mix", 0, 1, 1));

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
	juce::MidiMessage MidiMsg;//�ȴ���midi�¼�
	int MidiTime;
	juce::MidiBuffer::Iterator MidiBuf(midiMessages);
	while (MidiBuf.getNextEvent(MidiMsg, MidiTime))
	{
		if (MidiMsg.isController() && MidiMsg.getControllerNumber() == 10) { // 10 = Pan������
			const int channel = MidiMsg.getChannel(); // ͨ��ת������0��
			if (channel >= 0 && channel < 16) {
				channelPan[channel] = (float)(MidiMsg.getControllerValue()) / 127.0;
			}
		}
		if (MidiMsg.isNoteOn())
		{
			int note = MidiMsg.getNoteNumber();
			float velo = (float)MidiMsg.getVelocity() / 127.0;
			float pan = channelPan[MidiMsg.getChannel()];
			manbo.NoteOn(note, velo, pan);
		}
		if (MidiMsg.isNoteOff())
		{
			int note = MidiMsg.getNoteNumber();
			manbo.NoteOff(note);
		}
	}
	midiMessages.clear();

	const int numSamples = buffer.getNumSamples();
	float* wavbufl = buffer.getWritePointer(0);
	float* wavbufr = buffer.getWritePointer(1);
	const float* recbufl = buffer.getReadPointer(0);
	const float* recbufr = buffer.getReadPointer(1);

	float SampleRate = getSampleRate();

	float corr = *Params.getRawParameterValue("corr");
	float formant = *Params.getRawParameterValue("formant");
	float glide = *Params.getRawParameterValue("glide");
	float mix = *Params.getRawParameterValue("mix");

	auto toexp = [](float x, float n) {
		float sign = x < 0 ? -1.0f : 1.0f;
		x = fabsf(x);
		return (expf((x - 1.0) * n) - expf(-n)) / (1.0 - expf(-n)) * sign;
		};

	formant = powf(2.0, formant * 2.0);
	glide = toexp((1.0 - glide) * 2, 3);

	manbo.SetParam(corr, formant, glide, mix);
	manbo.ProcessBlock(recbufl, recbufr, wavbufl, wavbufr, numSamples);
}

//==============================================================================
bool LModelAudioProcessor::hasEditor() const
{
	return true; // (change this to false if you choose to not supply an editor)
}

juce::AudioProcessorEditor* LModelAudioProcessor::createEditor()
{
	return new LModelAudioProcessorEditor(*this);

	//return new juce::GenericAudioProcessorEditor(*this);//�Զ�����(����)
}

//==============================================================================

void LModelAudioProcessor::getStateInformation(juce::MemoryBlock& destData)
{
	// You should use this method to store your parameters in the memory block.
	// You could do that either as raw data, or use the XML or ValueTree classes
	// as intermediaries to make it easy to save and load complex data.

	// ����һ�� XML �ڵ�
	juce::XmlElement xml("LMEQ_Settings");

	/*juce::MemoryBlock eqDataBlock;
	eqDataBlock.append(&manager, sizeof(ResonatorManager));
	juce::String base64Data = eqDataBlock.toBase64Encoding();
	xml.setAttribute("VIB_MANAGER", base64Data);//���resonance����
	*/
	auto state = Params.copyState();
	xml.setAttribute("Knob_Data", state.toXmlString());//�����ť����

	juce::String xmlString = xml.toString();
	destData.append(xmlString.toRawUTF8(), xmlString.getNumBytesAsUTF8());
}

void LModelAudioProcessor::setStateInformation(const void* data, int sizeInBytes)
{
	// You should use this method to restore your parameters from this memory block,
	// whose contents will have been created by the getStateInformation() call.
	  // �� data ת��Ϊ�ַ����Խ��� XML
	juce::String xmlString(static_cast<const char*>(data), sizeInBytes);

	// ���� XML
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
