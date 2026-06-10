/********* Begin Register Information **********
{
    "name": "free_peaks_multiparty",
    "identifier": "FreePeaksMultiParty",
    "tags": [ "continuous", "multi_objective", "multiparty", "multimodal" ]
}
*********** End Register Information **********/

#ifndef OFEC_FREE_PEAKS_MULTIPARTY_H
#define OFEC_FREE_PEAKS_MULTIPARTY_H

#include "free_peaks.h"
#include <array>

namespace ofec {
    class FreePeaksMultiParty : public FreePeaks {
        OFEC_CONCRETE_INSTANCE(FreePeaksMultiParty)

    public:
        struct PartyPeakSpec {
            std::string name;
            std::vector<Real> center;
            Real height = 90;
            std::string shape = "s1";
            Real condition = 1;
            Real bias = 0;
            Real rotation = 0;
            Real linkage = 0;
            Real asymmetry = 0;
            Real deceptive = 0;
            Real radius = 32;
        };

        struct PartySpec {
            std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>> tree;
            std::vector<PartyPeakSpec> peaks;
        };

        struct SuiteSpec {
            std::string name;
            std::string feature;
            std::array<PartySpec, 2> parties;
            std::vector<std::vector<Real>> shared_optima;
        };

        static const std::vector<SuiteSpec>& suite();
        const SuiteSpec& currentSpec() const;
        int suiteId() const { return m_suite_id; }
        std::string partyRegionName(size_t party_id, const std::vector<Real>& x) const;

    protected:
        void addInputParameters();
        void initialize_(Environment* env) override;
        void evaluate(const VariableBase& vars, std::vector<Real>& objs, std::vector<Real>& cons) const override;
        void calculateOptima();
        void updateOptima(Environment* env) override;
		free_peaks::Subproblem* subproblem(const std::string& subspace_name) const override;
		const std::vector<std::pair<Real, Real>>& subspaceBox(const std::string& subspace_name) const override;

    private:
        void buildAssignedSuite(Environment* env);
        SuiteSpec makeRandomSuiteSpec();
        void generateRandomPeaks(SuiteSpec& spec);
        void setPartyKDTree(size_t party_id, const std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>>& tree_data);
        void setPartySubproblem(size_t party_id, const PartyPeakSpec& peak_spec);
        static std::array<Real, 2> point(Real x, Real y);

    private:
        int m_suite_id = 1;
        size_t m_problem_dimension = 2;
        SuiteSpec m_current_spec;
        std::vector<free_peaks::SubspaceKDTree> m_party_trees;
    };
}

#endif
