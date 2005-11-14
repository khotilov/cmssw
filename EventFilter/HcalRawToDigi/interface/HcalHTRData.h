/* -*- C++ -*- */
#ifndef HcalHTRData_H
#define HcalHTRData_H

/**  \class HcalHTRData
 *
 *  Interpretive class for HcalHTRData
 *  Since this class requires external specification of the length of the data, it is implemented 
 *  as an interpreter, rather than a cast-able header class.
 *
 *  $Date: 2005/07/26 15:10:51 $
 *  $Revision: 1.1 $
 *  \author J. Mans - UMD
 */

class HcalHTRData {
 public:
  static const int CHANNELS_PER_SPIGOT        ;// = 24;
  static const int MAXIMUM_SAMPLES_PER_CHANNEL;// = 20;

  HcalHTRData();
  ~HcalHTRData() { if (m_ownData!=0) delete [] m_ownData; }
  HcalHTRData(int version_to_create);
  HcalHTRData(const unsigned short* data, int length);
  HcalHTRData(const HcalHTRData&);

  HcalHTRData& operator=(const HcalHTRData&);
  void adoptData(const unsigned short* data, int length);

  /** \brief Get the version number of this event */
  inline int getFormatVersion() const { return m_formatVersion; }

  /** \brief Get a pointer to the raw data */
  inline const unsigned short* getRawData() const { return m_rawConst; }

  /** \brief Get the length of the raw data */
  inline const int getRawLength() const { return m_rawLength; }

  /** \brief Check for a good event 
   Requires a minimum length, matching wordcount and length, not an empty event */
  bool check() const;

  /** \brief Obtain the starting and ending pointers for external unpacking of the data 
      \param daq_first Pointer to a pointer to the start of the DAQ data
      \param daq_last Pointer to a pointer to the end of the DAQ data
      \param tp_first Pointer to a pointer to the start of the TP data
      \param tp_last Pointer to a pointer to the end of the TP data
  */
  void dataPointers(const unsigned short** daq_first, 
		    const unsigned short** daq_last, 
		    const unsigned short** tp_first, 
		    const unsigned short** tp_last); 


  /** \brief Unpack the HTR data into TP and DAQ data sorted by channel 
      \param daq_lengths unsigned char[24] of lengths.  High bit set indicates error with this channel
      \param daq_samples unsigned short [24*20] of data
      \param tp_lengths  unsigned char[24] of lengths
      \param tp_samples  unsigned short [24*20] of data
  */
  void unpack(unsigned char* daq_lengths, unsigned short* daq_samples,
	      unsigned char* tp_lengths, unsigned short* tp_samples) const;

  /** \brief Unpack special histogramming mode data
      \param fiber 
      \param fiberchan
      \param capid Capacitor id for which to extract a histogram 
      \param histogram unsigned int[32] into which the data should be deposited
  */
  bool unpackHistogram(int fiber, int fiberchan, int capid, unsigned short* histogram) const;

  /** \brief Unpack the HTR data into TP and DAQ data sorted by channel 
      \param daq_lengths unsigned char[24] of lengths
      \param daq_samples unsigned short [24*20] of data
      \param tp_lengths  unsigned char[24] of lengths
      \param tp_samples  unsigned short [24*20] of data
  */
  void pack(unsigned char* daq_lengths, unsigned short* daq_samples,
            unsigned char* tp_lengths, unsigned short* tp_samples, bool do_capid=false);
  /** \brief pack header and trailer (call _after_ pack)*/
  void packHeaderTrailer(int L1Anumber, int bcn, int submodule, int orbitn, int pipeline, int firmwareRev=0);

  /** \brief Get the HTR event number */
  inline unsigned int getL1ANumber() const { return (m_rawConst[0]&0xFF)+(m_rawConst[1]<<8); }
  /** \brief Get the HTR bunch number */
  inline unsigned int getBunchNumber() const { return (m_rawConst[4]&0xFFF); }
  /** \brief Get the HTR orbit number */
  unsigned int getOrbitNumber() const;
  /** \brief Get the HTR submodule number */
  unsigned int getSubmodule() const;
  /** \brief Is this event a calibration-stream event? */
  bool isCalibrationStream() const;
  /** \brief Is this event a pattern-ram event? */
  bool isPatternRAMEvent() const;
  /** \brief Is this event a histogram event? (do not call standard unpack in this case!!!!!) */
  bool isHistogramEvent() const;
  /** \brief Get the fiber numbers for the data present in this event (only in histogram mode!) */
  void getHistogramFibers(int& a, int& b) const;
  /** \brief Get the pipeline length used for this event */
  unsigned int getPipelineLength() const;
  /** \brief Get the HTR firmware version */
  unsigned int getFirmwareRevision() const;
  /** \brief Get the errors word */
  inline unsigned int getErrorsWord() const { return m_rawConst[2]&0xFF; }
  /** \brief Get the number of daq data samples when not zero-suppressed */
  int getNDD() const;
  /** \brief Get the number of trigger data samples when not zero-suppressed */
  int getNTP() const;
  /** \brief Get the number of presamples in daq data */
  int getNPS() const;
  /** \brief Was there an error on the given fiber for this event (only in histogram mode!) */
  bool wasHistogramError(int ifiber) const; 

private:
  void determineSectionLengths(int& tpWords, int& daqWords, int& headerWords, int& trailerWords) const;
  void determineStaticLengths(int& headerWords, int& trailerWords) const;
  int m_formatVersion;
  int m_rawLength;
  const unsigned short* m_rawConst;
  unsigned short* m_ownData;
};

#endif
