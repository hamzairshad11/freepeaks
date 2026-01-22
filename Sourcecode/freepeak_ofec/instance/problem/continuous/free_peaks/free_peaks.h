/********* Begin Register Information **********
{
	"name": "free_peaks",
	"identifier": "FreePeaks",
	"tags": [ "continuous"]
}
*********** End Register Information **********/

#ifndef OFEC_FREE_PEAKS_H
#define OFEC_FREE_PEAKS_H

#include "../../../../core/problem/continuous/continuous.h"
#include "subspace_kdtree.h"

namespace ofec {
#define CAST_FPs(pro) dynamic_cast<FreePeaks*>(pro)

	class FreePeaks : virtual public Continuous {
		OFEC_CONCRETE_INSTANCE(FreePeaks)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void evaluate(const VariableBase &vars, std::vector<Real> &objs, std::vector<Real> &cons) const override;
		void updateOptima(Environment *env) override;
		void readFromFile(const std::string &file_path);

		std::string m_generation_type, m_file_name;
		free_peaks::SubspaceKDTree m_subspace_tree;
		Real m_objective_accuracy;
		static bool ms_registered;

	public:
		void updateCandidates(const SolutionBase &sol, std::vector<std::unique_ptr<SolutionBase>> &candidates) const;
		bool isSolved(const std::vector<std::unique_ptr<SolutionBase>> &candidates) const;
		std::vector<Real> errorToPeaks(const std::vector<std::unique_ptr<SolutionBase>> &candidates) const;

		void writeToFile(const std::string &file_path) const;
		const free_peaks::SubspaceKDTree& subspaceTree() const { return m_subspace_tree; }
		free_peaks::Subproblem* subproblem(const std::string &subspace_name) const;
		std::string generationType() const { return m_generation_type; }
		std::string fileName() const { return m_file_name; }
	};
}

#endif // !OFEC_FREE_PEAKS_H
