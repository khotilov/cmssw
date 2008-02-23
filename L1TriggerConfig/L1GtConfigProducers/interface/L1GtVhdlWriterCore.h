#ifndef L1GtConfigProducers_L1GtVhdlWriterCore_h
#define L1GtConfigProducers_L1GtVhdlWriterCore_h

/**
 * \class L1GtVhdlWriterCore
 *
 *
 * \Description writes the actual VHDL code and stores global information like the common header.
 *  This class is using parameter maps containing the substitution parameters and their content.
 *
 * Implementation:
 *    <TODO: enter implementation details>
 *
 * \author Philipp Wagner
 *
 * $Date: 2008/02/21 21:59:50 $
 * $Revision: 1.5 $
 *
 */

// system include files

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

// CMSSW headers
#include "L1GtVhdlTemplateFile.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"
#include "CondFormats/L1TObjects/interface/L1GtFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"

#include "CondFormats/L1TObjects/interface/L1GtMuonTemplate.h"
#include "CondFormats/L1TObjects/interface/L1GtCaloTemplate.h"
#include "CondFormats/L1TObjects/interface/L1GtEnergySumTemplate.h"

#include "L1TriggerConfig/L1GtConfigProducers/interface/L1GtVhdlWriterBitManager.h"
#include "L1TriggerConfig/L1GtConfigProducers/interface/L1GtVhdlWriterCore.h"

class L1GtVhdlWriterCore
{

    public:

        /// constructor
        L1GtVhdlWriterCore(const std::string &templatesDirectory, const std::string &outputDirectory, const bool &debug);

        /// destructor
        ~L1GtVhdlWriterCore();

        /// returns all condition of the same class. Conditions belong to one class if they are matching
        /// in type (Type1s, Type2s..) category (CondMuon, CondCalo,)and are defined for the same object (Mu, fwdJet..)
        /// \param map is the ConditionMap that shall be filtered
        /// \param conditionToIntegerMap is the condition is added to this conversion map
        /// \param outputMap here the conditions matching the requirements are stored
        bool returnConditionsOfOneClass(const L1GtConditionType &type,  const L1GtConditionCategory &category, const L1GtObject &object, const ConditionMap &map,  ConditionMap &outputMap);

        /// builds a parameter map that can be used by the L1GtVhdlWriterCore class
        /// the output is written into the parameter &muonParameters;
        /// \param conditionToIntegerMap is the condition is added to this conversion map
        /// \param muonParameters here the routine stores the content of the subsitution parameters
        /// \param conditionMap containing the input
        void buildMuonParameterMap(const unsigned short int &condChip, std::map<std::string,int> &conditionToIntegerMap,
            std::map<std::string,std::string> &muonParameters, const std::vector<ConditionMap> &conditionMap);

        bool buildCaloParameterMap(const unsigned short int &condChip, std::map<std::string,int> &conditionToIntegerMap,
            std::map<std::string,std::string> &muonParameters,const L1GtObject &caloObject, const std::vector<ConditionMap> &conditionMap);

        bool buildEnergySumParameter(const unsigned short int &condChip,  const L1GtObject &object, std::map<std::string,int> &conditionToIntegerMap,
            std::string &energySumParameter, const std::vector<ConditionMap> &conditionMap);

        /// builds the substitution parameters for the cond_chip.vhd
        /// \condChip the condition chip that will be processed
        /// \param conditionToIntegerMap this has to be a already FILLED conversion map.
        ///  therefore this routine has to be called after all those which are adding information
        ///  to the conditionToIntegerMap map (buildMuonParameterMap, buildCaloParameterMap...)
        /// \param templates in this map the final content for the subsitution parameters is stored in
        /// VHDL template file format
        bool buildCondChipParameters(const unsigned short int &condChip, std::map<std::string,int> &conditionToIntegerMap,
            const std::vector<ConditionMap> &conditionMap, std::map<std::string, L1GtVhdlTemplateFile> &templates, std::map<std::string, std::string> &commonParams);

        std::string retNumberOfConditionsString(const std::string &typeStr, const int &number);
        
        bool findObjectType(const L1GtObject &object, ConditionMap &map);

        /// builds the parameters particle_common for the cond_chip.vhd's
        bool buildCommonParameter(L1GtVhdlTemplateFile &particle,const L1GtObject &object, const L1GtConditionCategory &category, std::string &parameterStr, const ConditionMap &conditionMap);

        /// calculates the integer value for jet counts conditions and furthermore counts
        /// how many jet counts of one type are in trigger menu
        void addJetCountsToCond2IntMap(const int chip, const std::vector<ConditionMap> &conditionMap, std::map<std::string,int> &conditionToIntegerMap);

        /// processes algorithm map delivered by parser, replaces condition names by types and serial numbers,
        /// splits the map in two seperate ones for the two condition chips
        bool processAlgorithmMap(const AlgorithmMap &algorithmMap,  std::vector<ConditionMap> &conditionMap,
            std::map<std::string,int> &conditionToIntegerMap, std::map<int, std::string> &algorithmsChip1, std::map<int, std::string> &algorithmsChip2);

        /// builds the common header for all files
        void buildCommonHeader(std::map<std::string,std::string> &headerParameters , const std::vector<std::string> &connectedChannels);

        /// opens a new template file and inserts the common header
        L1GtVhdlTemplateFile openVhdlFileWithCommonHeader(const std::string &filename, const std::string &outputFilename);

        /// returns the common header
        L1GtVhdlTemplateFile retrunCommonHeader();

        /// prints the common header
        void printCommonHeader();

        /// produces the firmware code
        bool makeFirmware(std::vector<ConditionMap> &conditionMap,const AlgorithmMap &algorithmMap);

        /// builds muon setup files
        void writeMuonSetupVhdl(std::map<std::string,std::string> &muonParameters, const std::string &particle, unsigned short int &condChip);

        /// builds cond_chip.vhds
        void writeConditionChipSetup(std::map<std::string, L1GtVhdlTemplateFile> templates, const std::map<std::string, std::string> &common, const unsigned short int &chip);

        /// writes def_val_pkg.vhd
        void writeDefValPkg();
        
        /// builds etm setup files
        void writeEtmSetup(std::string &etmString, const int &condChip);

        /// builds the prealgo_and_or setup
        void writeAlgoSetup(std::map<int, std::string> &algorithmsChip1, std::map<int, std::string> &algorithmsChip2);

        /// builds the two quartus setup files. This routine is called in buildCommonHeader!
        void writeQsfSetupFiles(const std::string &version);

        /// converts a integer into a string
        static std::string int2str(const int &integerValue);

        /// checks weather value searchValue exists in a <string,int> map, saves it in &intVal if it exists and returns false if not
        bool getIntVal(const std::map<std::string,int> &map,  const std::string &searchValue, int &intVal);

        /// returns the templates path
        std::string gtTemplatesPath();

        /// adds a string to intern message buffer
        void msg(const std::string &message);

        /// returns intern message buffer
        std::vector<std::string> getMsgBuf();

        /// builds the index for the cond_chip vhds. (They have to start at 1 so one has to be added
        /// to the value in conversion map)
        std::string index4CondChipVhd(int intval);

        /// returns conditions of a class in output map and counts the occurance of all conditions of a
        /// certain type at the same time with the result being stored in numberOfConditions_ in the
        /// following format: CONSTANT nr_muon_3 : integer := 0; Tis is the format for the cond_chip_pkg files.
        void countCondsAndAdd2NumberVec(const L1GtConditionType &type,  const L1GtConditionCategory &category, const L1GtObject &object, const ConditionMap &map, ConditionMap &outputMap,const int &condChip);

        void initializeDeltaConditions();

        /// produces a control output file for condition to integer conversion
        void writeCond2intMap2File(const std::map<std::string,int> &conditionToIntegerMap);

        /// for debuggin
        void printConditionsOfCategory(const L1GtConditionCategory &category, const ConditionMap &map);

        void writeCondChipPkg(const int &chip);
        
        std::map<std::string,int> getCond2IntMap();

    private:

        /// stores to condition name to integer conversion table
        std::map<std::string,int> conditionToIntegerMap_;
        
        /// converts L1GtConditionType to firmware string
        std::map<L1GtConditionType,std::string> condType2Str_;

        /// converts L1GtObject to calo_nr
        std::map<L1GtObject,std::string> caloType2Int_;

        /// converts L1GtObject to string
        std::map<L1GtObject,std::string> objType2Str_;

        /// list of all possible calo objects
        std::vector<L1GtObject> caloObjects_;

        /// list of all possible esums objects
        std::vector<L1GtObject> esumObjects_;

        /// templates directory
        std::string vhdlDir_;

        /// output directory
        std::string outputDir_;

        /// bit manager for bit operations
        L1GtVhdlWriterBitManager bm_;

        /// common header for all files
        L1GtVhdlTemplateFile commonHeader_;

        std::vector<std::string> internMessageBuf_;

        /// vector containing the initialization of all conditions
        std::vector<std::vector<std::string> > numberOfConditions_;

        /// class will produce some additional debugging output if set
        bool debugMode_;

};
#endif                                            /*L1GtConfigProducers_L1GtVhdlWriterCore_h*/
