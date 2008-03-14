#ifndef PhysicsTools_PatAlgos_interface_IsoDepositIsolator_h
#define PhysicsTools_PatAlgos_interface_IsoDepositIsolator_h

#include "PhysicsTools/PatAlgos/interface/BaseIsolator.h"
#include "DataFormats/MuonReco/interface/MuIsoDeposit.h"


namespace pat { namespace helper {
class IsoDepositIsolator : public BaseIsolator {
    public:
        typedef edm::ValueMap<reco::MuIsoDeposit> Isolation;
 
        IsoDepositIsolator() {}
        IsoDepositIsolator(const edm::ParameterSet &conf) ;
        virtual ~IsoDepositIsolator() {}
        virtual void beginEvent(const edm::Event &event) ;
        virtual void endEvent() ;

        virtual std::string description() const ;
    protected:
        edm::Handle<Isolation> handle_;
        float deltaR_;
        virtual float getValue(const edm::ProductID &id, size_t index) const {
            return handle_->get(id, index).depositWithin(deltaR_);
        }
}; // class IsoDepositIsolator
} } // namespaces

#endif
